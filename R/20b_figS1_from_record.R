# =============================================================================
# 20b_figS1_from_record.R   (added 12 October 2026)
# Figure S1, the chironomid check, drawn from step 20's own tables and the site matrices; nothing is recomputed.
#   (a) the share of the trait-set density and of the biomass made up by the four chironomid taxa, against biological
#       condition, within region (the mechanism behind Discussion 4.2 and 4.3; coefficients are Table S7 part A);
#   (b) functional dispersion and overlap redundancy, weighted by density, with every taxon and with the four removed
#       (Table S7 part B). Rao's quadratic entropy nearly duplicates functional dispersion (Appendix S1.6) and is not drawn.
# The former panel (b), trait states under the fourth-corner analysis, went to the trait-state manuscript
# (R/supplement_asset_map.csv). Step 20 also writes the old figure name; this step overwrites it with this version.
# INPUTS   den_mat, bio_mat and dat (step 01); 20_chironomid_share.csv and 20_chironomid_functional_metrics.csv (step 20).
# OUTPUT   output/tables/20b_chironomid_share_by_site.csv; output/figures/figS5_chironomid_removal.{png,pdf,svg}, also
#          copied as Figure_S1 into output/manuscript/figures/ and output/manuscript/supplement_final/ when those folders exist
#          (on a first run they do not yet, the copies fail with a warning, and steps 18 and 18b write the figure: 18 as Figure_S3,
#          18b as supplement_final/Figure_S1).
# =============================================================================
say("STEP 20b -- Figure S1 from the step 20 record")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(mgcv) })
CHIRONOMIDAE <- c("Chironomini", "Orthocladiinae", "Tanypodinae", "Tanytarsini")
stopifnot("a chironomid group is missing from the trait-set matrix" = all(CHIRONOMIDAE %in% colnames(den_mat)), all(CHIRONOMIDAE %in% colnames(bio_mat)))
RL <- c("Valley & Ridge" = "Valley and Ridge", "Piedmont" = "Piedmont", "Coastal Plain" = "Coastal Plain")

# ---- (a) share by site, read from the same matrices step 20 used; asserted against step 20's table
st <- dat %>% select(SiteID, Region, BIBI); st <- st[match(rownames(den_mat), st$SiteID), ]
stopifnot(setequal(st$SiteID, rownames(den_mat)), setequal(st$SiteID, rownames(bio_mat)))
bio_mat2 <- bio_mat[rownames(den_mat), colnames(den_mat), drop = FALSE]
share_site <- bind_rows(
  tibble(SiteID = st$SiteID, Region = st$Region, BIBI = st$BIBI, Weighting = "Density", share = 100 * rowSums(den_mat[, CHIRONOMIDAE]) / rowSums(den_mat)),
  tibble(SiteID = st$SiteID, Region = st$Region, BIBI = st$BIBI, Weighting = "Biomass", share = 100 * rowSums(bio_mat2[, CHIRONOMIDAE]) / rowSums(bio_mat2)))
write_csv(share_site, file.path(OUT_TAB, "20b_chironomid_share_by_site.csv"))
sh_rec <- read_csv(file.path(OUT_TAB, "20_chironomid_share.csv"), show_col_types = FALSE)
chk <- share_site %>% group_by(Weighting, Region) %>% summarise(rho = suppressWarnings(cor(share, BIBI, method = "spearman")), .groups = "drop") %>%
  left_join(sh_rec %>% select(Weighting, Region, rec = `rho with B-IBI`), by = c("Weighting", "Region"))
stopifnot("the site-level shares do not reproduce step 20's coefficients" = all(abs(chk$rho - chk$rec) < 1e-9))
say("  site-level shares reproduce all 6 coefficients of 20_chironomid_share.csv")
fa <- share_site %>% mutate(Region = factor(Region, REGIONS), Weighting = factor(Weighting, c("Density", "Biomass")))
lev <- c(Density = "weighted by density", Biomass = "weighted by biomass")
pa <- ggplot(fa, aes(BIBI, share, colour = Region, fill = Region)) +
  geom_smooth(aes(linetype = Weighting), method = "gam", formula = y ~ s(x, k = 3), method.args = list(method = "REML"), linewidth = .6, alpha = .18) +
  geom_point(aes(shape = Weighting), size = 1.6, stroke = .6) +
  facet_wrap(~Region, nrow = 1, labeller = labeller(Region = RL)) +
  scale_colour_manual(values = REGION_COLS, limits = REGIONS, guide = "none") + scale_fill_manual(values = REGION_COLS, limits = REGIONS, guide = "none") +
  scale_shape_manual(values = c(Density = 16, Biomass = 1), labels = lev) + scale_linetype_manual(values = c(Density = "solid", Biomass = "dashed"), labels = lev) +
  scale_x_bibi(BIBI_AXIS_LABEL) + scale_y_continuous(breaks = seq(0, 100, 25)) + coord_cartesian(ylim = c(0, 100)) +  # coord_cartesian, not limits =, which clipped the smooths
  labs(y = "Four chironomid taxa, share of the trait-set total (%)") +
  theme_pub() + theme(legend.position = "bottom", legend.title = element_blank(), panel.spacing.x = unit(1.2, "lines"))

# ---- (b) step 20's table, Rao's quadratic entropy dropped
met_tab <- read_csv(file.path(OUT_TAB, "20_chironomid_functional_metrics.csv"), show_col_types = FALSE)
COL2 <- c("Every taxon" = "#1A1A1A", "Four chironomid taxa removed" = "#8C8C8C")
fm <- met_tab %>% filter(Weighting == "Density", Metric %in% c("FDis_den", "FR_RD_den")) %>%
  mutate(Metric = factor(c(FDis_den = "Functional dispersion", FR_RD_den = "Overlap redundancy")[Metric], c("Functional dispersion", "Overlap redundancy")),
         Region = factor(RL[Region], rev(RL[REGIONS]))) %>%
  tidyr::pivot_longer(c(rho_all, rho_removed, met_all, met_removed), names_to = c(".value", "set"), names_pattern = "(rho|met)_(all|removed)") %>%
  mutate(set = factor(ifelse(set == "all", "Every taxon", "Four chironomid taxa removed"), names(COL2)))
stopifnot(nrow(fm) == 12)
pb <- ggplot(fm, aes(rho, Region)) + geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") +
  geom_line(aes(group = interaction(Region, Metric)), linewidth = .5, colour = "grey55") +
  geom_point(aes(colour = set, shape = met), size = 2.4, stroke = .8, position = position_nudge(y = ifelse(fm$set == "Every taxon", 0.08, -0.08))) + facet_wrap(~Metric, nrow = 1) +
  scale_colour_manual(values = COL2) + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "Meets the criterion", `FALSE` = "Does not")) +
  scale_x_continuous(limits = c(-.8, .8), breaks = seq(-.6, .6, .3)) + labs(x = "Spearman ρ with B-IBI, weighted by density", y = NULL) +
  theme_pub() + theme(legend.position = "bottom", legend.title = element_blank()) + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))

fig <- (pa / pb) + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(plot.tag = element_text(face = "bold", size = FIG_BASE_PT + 1))
save_fig(fig, "figS5_chironomid_removal", FIG_W_DOUBLE, 6.6)
file.copy(file.path(OUT_FIG, paste0("figS5_chironomid_removal", c(".png", ".pdf", ".svg"))),
          file.path(OUT, "manuscript", "figures", paste0("Figure_S1", c(".png", ".pdf", ".svg"))), overwrite = TRUE)
file.copy(file.path(OUT_FIG, paste0("figS5_chironomid_removal", c(".png", ".pdf", ".svg"))),   # the supplement_final copy the manifest names
          file.path(OUT, "manuscript", "supplement_final", paste0("Figure_S1", c(".png", ".pdf", ".svg"))), overwrite = TRUE)
say("STEP 20b complete")
