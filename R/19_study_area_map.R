# =============================================================================
# 19_study_area_map.R  --  the study area (Figure 2 in the manuscript; step 18 numbers it)
# =============================================================================
# WHY THIS STEP EXISTS
#   The manuscript needs a study-area figure. This step draws it -- CONUS
#   Albers, state outlines, region hulls, locator inset, scale bar, north
#   arrow -- from the site coordinate file (data/site_coordinates.csv in the
#   repository, 90 sites in the three regions; outside the repository layout a
#   120-site project file, filtered to the three regions), drawn with this
#   pipeline's PAL_REGION so the map and every other figure share one palette.
#
# WHY IT CARRIES PANEL (b)
#   A map alone shows where the sites are and nothing else. Panel (b) puts the
#   condition gradient each region spans beside it, so the range of B-IBI
#   scores behind every within-region test is visible rather than buried in
#   prose. Together the two panels carry the study system, so Table 1 is no
#   longer needed in the main text.
#
# INPUTS   the site coordinate file above; IN$bibi (site B-IBI and condition class).
# OUTPUT   output/figures/fig1_study_area.{png,pdf,svg}; step 18 copies it as Figure 2.
# =============================================================================

say("STEP 19 -- study area map (Figure 2)")
suppressPackageStartupMessages({
  library(sf); library(maps); library(ggspatial); library(cowplot); library(patchwork)
})

SITES_CSV <- if (REPO_LAYOUT) file.path(DATA, "site_coordinates.csv") else
             file.path(dirname(CH2_ROOT), "Chapter 3 - Land Use + Managment Assessmnet", "R_Poster", "sites_120.csv")
stopifnot("site coordinate file missing" = file.exists(SITES_CSV))

# The three regions of this study. The project coordinate file labels this study's Piedmont "Piedmont — mixed";
# its urban Piedmont sites belong to a companion study and are dropped by the region filter.
sites <- read_csv(SITES_CSV, show_col_types = FALSE) %>%
  mutate(Region = dplyr::recode(Region, "Piedmont — mixed" = "Piedmont")) %>%
  filter(Region %in% REGIONS) %>%
  mutate(Region = factor(Region, levels = REGIONS))
say(sprintf("  %d sites across %d regions", nrow(sites), dplyr::n_distinct(sites$Region)))

bibi <- read_csv(IN$bibi, show_col_types = FALSE) %>%
  select(SiteID, BIBI) %>% distinct(SiteID, .keep_all = TRUE)
sites <- left_join(sites, bibi, by = "SiteID")
stopifnot(sum(is.na(sites$BIBI)) == 0)

ALBERS <- 5070
pts <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326,
                remove = FALSE) %>% st_transform(ALBERS)
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  st_set_crs(4326) %>% st_make_valid() %>% st_transform(ALBERS)

hulls <- pts %>% group_by(Region) %>%
  summarise(do_union = TRUE, .groups = "drop") %>%
  st_convex_hull() %>% st_buffer(9000)

# 25 km of air, not 60. The three regions sit in an L, so a generous buffer
# bought a large empty quadrant in the south-west and nothing else.
bb <- st_bbox(st_buffer(pts, 25000))

p_map <- ggplot() +
  geom_sf(data = states, fill = "white", colour = "grey60", linewidth = 0.35) +
  geom_sf(data = hulls, aes(fill = Region), colour = NA, alpha = 0.15,
          show.legend = FALSE) +
  geom_sf(data = pts, aes(colour = Region), size = 1.7, shape = 21,
          fill = "white", stroke = 0.7) +
  geom_sf(data = pts, aes(colour = Region), size = 0.95, show.legend = FALSE) +
  scale_colour_manual(values = PAL_REGION, name = NULL) +
  scale_fill_manual(values = PAL_REGION, guide = "none") +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE, crs = st_crs(ALBERS)) +
  # ggspatial draws the scale bar and the north arrow with its own font setting,
  # not the ggplot theme's, so these two need FIG_FAMILY passed explicitly or
  # they alone come out sans while the rest of the figure is Times.
  annotation_scale(location = "bl", width_hint = 0.26, height = unit(0.12, "cm"),
                   text_cex = 0.55, text_family = FIG_FAMILY,
                   pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.85, "cm"), width = unit(0.7, "cm"),
                         pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                         style = north_arrow_orienteering(text_size = 6,
                                                          text_family = FIG_FAMILY)) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, stroke = 0.9))) +
  theme_void(base_size = FIG_BASE_PT, base_family = FIG_FAMILY) +
  theme(panel.background = element_rect(fill = "#F4F7F9", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 0.6),
        legend.position = "bottom", legend.text = element_text(size = FIG_BASE_PT - 0.5),
        legend.margin = margin(t = -2), legend.key.height = unit(0.30, "cm"))

inset <- ggplot() +
  geom_sf(data = states, fill = "white", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = st_as_sfc(st_bbox(c(xmin = bb[["xmin"]], xmax = bb[["xmax"]],
                                     ymin = bb[["ymin"]], ymax = bb[["ymax"]]),
                                   crs = st_crs(ALBERS))),
          fill = NA, colour = "#B03A00", linewidth = 0.7) +
  coord_sf(expand = FALSE, crs = st_crs(ALBERS)) +
  theme_void(base_family = FIG_FAMILY) +
  theme(panel.background = element_rect(fill = "white", colour = "grey70",
                                        linewidth = 0.5))

panel_a <- ggdraw(p_map) + draw_plot(inset, x = 0.005, y = 0.70, width = 0.26, height = 0.26)

# ---- panel (b): one panel, every site -----------------------------------------
# WHY b AND c ARE NOW ONE PANEL.
#   Separately they showed the same 90 sites twice: a range bar that hid how the
#   sites are distributed inside it, and a stacked count that hid where along the
#   gradient those counts fall. Stacking each site as a dot at its own B-IBI
#   value turns the two into one mini-histogram per region: the width of the
#   stack is the span, the height is how many sites sit there, and the shape is
#   the condition class. Shapes match the response figures, so a reader learns
#   them once, and n is equal by design (30 per region), so it belongs in the
#   caption rather than repeated three times on the panel.
site <- read_csv(IN$bibi, show_col_types = FALSE) %>%
  filter(Region %in% REGIONS) %>%
  mutate(Region     = factor(Region, levels = rev(REGIONS)),
         BIBI_Class = factor(BIBI_Class, levels = CLASS_LEVELS))

# Stack sites that share a bin, so density along the gradient is visible as
# height. Bin width is set in B-IBI units, not by n, so the three regions stay
# directly comparable.
BINW <- 3.5
site <- site %>%
  mutate(bin = floor(BIBI / BINW)) %>%
  group_by(Region, bin) %>%
  arrange(BIBI, .by_group = TRUE) %>%
  mutate(stack = dplyr::row_number()) %>%
  ungroup()
DY <- min(0.105, 0.72 / max(site$stack))
site <- site %>%
  mutate(ypos = as.numeric(Region) + 0.10 + (stack - 1) * DY)

gr <- site %>% group_by(Region) %>%
  summarise(lo = min(BIBI), hi = max(BIBI), med = median(BIBI), .groups = "drop") %>%
  mutate(ybar = as.numeric(Region) + 0.02)

panel_b <- ggplot(site) +
  geom_segment(data = gr, aes(x = lo, xend = hi, y = ybar, yend = ybar),
               colour = "grey82", linewidth = 2.2, lineend = "round") +
  geom_segment(data = gr, aes(x = med, xend = med, y = ybar - 0.055,
                              yend = ybar + 0.055),
               colour = "grey30", linewidth = 0.55) +
  geom_point(aes(BIBI, ypos, colour = Region, shape = BIBI_Class), size = 1.6) +
  scale_colour_manual(values = PAL_REGION, guide = "none") +
  scale_shape_manual(values = CLASS_SHAPES, name = "Condition class") +
  scale_y_continuous(breaks = seq_along(levels(site$Region)),
                     labels = levels(site$Region),
                     limits = c(0.55, max(site$ypos) + 0.14),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = 0.02)) +
  labs(x = BIBI_AXIS_LABEL, y = NULL) +
  guides(shape = guide_legend(nrow = 2, title.position = "top",
                              override.aes = list(size = 2))) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y  = element_text(size = FIG_BASE_PT, face = "bold",
                                    vjust = 0, margin = margin(r = 3)),
        axis.ticks.y = element_blank(),
        axis.title.x      = element_text(margin = margin(t = 3, b = 4)),
        legend.position   = "bottom",
        legend.title      = element_text(size = FIG_BASE_PT - 1, hjust = 0.5,
                                         margin = margin(b = 2)),
        legend.text       = element_text(size = FIG_BASE_PT - 1),
        legend.key.size   = unit(0.28, "cm"),
        legend.margin     = margin(t = -1))

fig1 <- panel_a + panel_b + plot_layout(widths = c(1.45, 1))
fig1 <- tag_panels(fig1)
save_fig(fig1, "fig1_study_area", FIG_W_DOUBLE, 4.1)
say("STEP 19 complete")
