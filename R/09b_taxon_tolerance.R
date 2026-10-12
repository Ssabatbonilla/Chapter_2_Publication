# =============================================================================
# 09b_taxon_tolerance.R
# =============================================================================
# THE QUESTION THIS STEP ANSWERS
#   Step 09 screens every taxon against the B-IBI gradient within region, and
#   finds the response one-sided: almost every taxon that meets FDR <= 0.10
#   INCREASES in density (and biomass) as condition improves. Are those the
#   pollution-sensitive taxa, as hypothesis 1's mechanism says (stressors reduce
#   sensitive taxa), or are they generalists that simply become more numerous?
#   The answer decides whether H1 can be stated through density ("sensitive taxa
#   are more numerous in better streams") without restating how HBI is built.
#
# HOW
#   Tolerance values are the family-level Chessie values the site HBI was
#   computed from (IN$hbi_tv, one per taxon), classed as the index classes them:
#   sensitive TV < 3, tolerant TV > 7, intermediate otherwise (Buchanan et al.
#   2011, Table 6). Each tested taxon carries its class and its body-size,
#   flight and respiration affinities (IN$traits), so the tolerance of the
#   responders can be read beside their traits. Increasers are compared with the
#   taxa tested in the same region and weighting that did not respond (Wilcoxon
#   on TV, run when each side has at least three values), because the responders
#   can only be drawn from the tested pool.
#
# OUTPUT
#   output/tables/09b_taxon_tolerance_screen.csv   one row per tested taxon; its increasers and decreasers
#                                                  are supplement Table S8 part B (legacy Table S24, step 18)
#   output/tables/09b_tolerance_summary.csv        one row per region x weighting (log record; no later step reads it)
#
# Added 21 September 2026.
# =============================================================================

say("STEP 09b -- tolerance values of the taxa that respond to condition")

tol_screen <- read_csv(file.path(OUT_TAB, "09_taxon_response_screen.csv"), show_col_types = FALSE)

tv <- read_csv(IN$hbi_tv, show_col_types = FALSE) %>%
  distinct(TAXANAME, Family, TV) %>%
  group_by(TAXANAME) %>%
  summarise(Family = first(Family), n_tv = n_distinct(TV[!is.na(TV)]),
            TV = suppressWarnings(as.numeric(first(TV[!is.na(TV)]))), .groups = "drop")
if (any(tv$n_tv > 1)) stop("Step 09b: a taxon carries more than one tolerance value in ", basename(IN$hbi_tv))
say("  tolerance values: %d taxa in %s, %d without a value",
    nrow(tv), basename(IN$hbi_tv), sum(is.na(tv$TV)))

tv_class <- function(x) case_when(
  is.na(x)                    ~ "no value",
  x <  TV_SENSITIVE_BELOW     ~ "sensitive",
  x >  TV_TOLERANT_ABOVE      ~ "tolerant",
  TRUE                        ~ "intermediate")

trait_cols <- c("Size_Small", "Size_Medium", "Size_Large",
                "Flying_Weak", "Flying_Strong",       # Twardochleb names since 22 Sep 2026; Respiration_Plastron is
                "Respiration_Gills", "Respiration_Tegument", "Respiration_Plastron")   # Twardochleb's plastron or spiracle
traits <- readRDS(IN$traits) %>% select(TAXANAME, any_of(trait_cols))

out <- tol_screen %>%
  left_join(tv %>% select(TAXANAME, Family, TV), by = "TAXANAME") %>%
  mutate(TV_class = tv_class(TV)) %>%
  left_join(traits, by = "TAXANAME") %>%
  mutate(Region = factor(Region, levels = REGIONS)) %>%
  arrange(Region, Currency, Response, desc(rho)) %>%
  mutate(Region = as.character(Region))
write_csv(out, file.path(OUT_TAB, "09b_taxon_tolerance_screen.csv"))

summ <- map_dfr(REGIONS, function(rg) map_dfr(c("Density", "Biomass"), function(cu) {
  s    <- out %>% filter(Region == rg, Currency == cu)
  inc  <- s %>% filter(Response == "increaser")
  dec  <- s %>% filter(Response == "decreaser")
  back <- s %>% filter(Response == "no response")
  cnt  <- function(d, k) sum(d$TV_class == k)
  w <- if (sum(!is.na(inc$TV)) >= 3 && sum(!is.na(back$TV)) >= 3)
         suppressWarnings(wilcox.test(inc$TV, back$TV)$p.value) else NA_real_
  tibble(Region = rg, Currency = cu,
         n_tested = nrow(s), n_increaser = nrow(inc), n_decreaser = nrow(dec),
         inc_sensitive = cnt(inc, "sensitive"), inc_intermediate = cnt(inc, "intermediate"),
         inc_tolerant = cnt(inc, "tolerant"), inc_no_value = cnt(inc, "no value"),
         dec_sensitive = cnt(dec, "sensitive"), dec_intermediate = cnt(dec, "intermediate"),
         dec_tolerant = cnt(dec, "tolerant"), dec_no_value = cnt(dec, "no value"),
         back_sensitive = cnt(back, "sensitive"), back_intermediate = cnt(back, "intermediate"),
         back_tolerant = cnt(back, "tolerant"), back_no_value = cnt(back, "no value"),
         median_TV_increasers = median(inc$TV, na.rm = TRUE),
         median_TV_background = median(back$TV, na.rm = TRUE),
         p_wilcox_TV = w)
}))
write_csv(summ, file.path(OUT_TAB, "09b_tolerance_summary.csv"))

say("  responders by tolerance class (sensitive TV < %g, tolerant TV > %g):", TV_SENSITIVE_BELOW, TV_TOLERANT_ABOVE)
for (i in seq_len(nrow(summ))) {
  z <- summ[i, ]
  say("    %-15s %-8s increasers %2d (sens %d, int %d, tol %d, none %d) | decreasers %d | median TV %.1f vs background %.1f (Wilcoxon p %s)",
      z$Region, z$Currency, z$n_increaser, z$inc_sensitive, z$inc_intermediate, z$inc_tolerant,
      z$inc_no_value, z$n_decreaser, z$median_TV_increasers, z$median_TV_background,
      ifelse(is.na(z$p_wilcox_TV), "n/a", sprintf("%.3f", z$p_wilcox_TV)))
}
say("STEP 09b complete")
