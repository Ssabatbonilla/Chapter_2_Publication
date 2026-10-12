# =============================================================================
# 02_concentration_metrics.R — the concentration family (NEW in v2)
# =============================================================================
# WHY THIS STEP EXISTS
#   An earlier draft of the manuscript claimed that biomass is concentrated in a
#   few large-bodied taxa and that this is a response to biological condition.
#
#   Nothing in the 19 screened metrics of the four original families measures concentration. There is no
#   Simpson, no Berger-Parker, no percent-dominant-taxon in the response set.
#   The only place dominance appears as a magnitude anywhere in the paper is the
#   circle sizing on the trait-space figure, which has no test behind it.
#
#   This step supplies the missing test.
#
# WHY NOT SIMPSON DOMINANCE
#   Simpson_den and Simpson_bio are already in the site table (computed by 00c)
#   and deliberately excluded from screening. Keep them excluded -- but for a better reason than "QC field":
#   Simpson dominance is 1 - Simpson diversity, and Simpson diversity D is the
#   DENOMINATOR of FR_RD = 1 - Q/D. Screening it would re-import exactly the
#   arithmetic-dependency problem that already makes FR_DB_RD non-independent
#   of FR_RD_den and FR_RD_bio.
#
#   Gini and N50 are independent of the redundancy family. Berger-Parker was
#   tested and is too noisy at n = 30 -- nothing passed in any region.
#
# WHY A SEPARATE BH FAMILY
#   Benjamini-Hochberg is applied within family. Putting these five metrics in
#   their own family means adding them CANNOT move any published FDR value in
#   the four existing families. That is a deliberate design choice, not an
#   oversight.
#
# OUTPUT
#   conc : 90 rows, one per site, five concentration metrics (screened by step 03:
#          the inequality rows of Table 2 and Table S4 part A; also read by step 17)
#   output/tables/02_concentration_site_level.csv
# =============================================================================

say("STEP 02 -- concentration metrics")

# ---- Metric definitions ------------------------------------------------------

#' Gini coefficient of a community's taxon share vector.
#' 0 = every taxon holds an equal share; 1 = one taxon holds everything.
#' Higher = more concentrated.
gini_shares <- function(x) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n < 2) return(NA_real_)
  x <- sort(x / sum(x))
  sum((2 * seq_len(n) - n - 1) * x) / n
}

#' Number of taxa needed to reach 50% of community density or biomass.
#' Higher = LESS concentrated (it takes more taxa to make up half the total).
n50_taxa <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(NA_real_)
  p <- sort(x / sum(x), decreasing = TRUE)
  as.numeric(which(cumsum(p) >= 0.5)[1])
}

# ---- Compute on the SAME matrices FD::dbFD uses ------------------------------
# den_mat / bio_mat are the 90-site trait-eligible matrices from step 01.
# Using them, rather than the raw site x taxon file, is what makes these metrics
# directly comparable to FDiv, FR_RD and the rest.

conc <- tibble(SiteID = rownames(den_mat)) %>%
  mutate(
    Gini_den = apply(den_mat, 1, gini_shares),
    Gini_bio = apply(bio_mat[SiteID, , drop = FALSE], 1, gini_shares),
    N50_den  = apply(den_mat, 1, n50_taxa),
    N50_bio  = apply(bio_mat[SiteID, , drop = FALSE], 1, n50_taxa)
  ) %>%
  # The weighting analogue of FR_DB_RD: how much more (or less) concentrated the
  # community is by biomass than by density.
  mutate(Gini_DB_DOM = Gini_bio - Gini_den)

stopifnot(nrow(conc) == 90)
stopifnot(!any(is.na(conc$Gini_den)), !any(is.na(conc$Gini_bio)))

write_csv(conc, file.path(OUT_TAB, "02_concentration_site_level.csv"))

say("  computed %d concentration metrics for %d sites",
    length(CONCENTRATION_FAMILY), nrow(conc))
say("  Gini_den range %.3f-%.3f | Gini_bio range %.3f-%.3f",
    min(conc$Gini_den), max(conc$Gini_den), min(conc$Gini_bio), max(conc$Gini_bio))
say("STEP 02 complete")
