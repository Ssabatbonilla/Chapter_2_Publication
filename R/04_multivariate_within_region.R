# =============================================================================
# 04_multivariate_within_region.R — within-region and cross-region multivariate tests
# =============================================================================
# WHY THIS STEP EXISTS
#   Section 3.1 reports whether community structure differs among regions and,
#   within each region, whether it varies with biological condition, weighted by
#   density and by biomass.
#
# METHOD
#   Bray-Curtis dissimilarity on square-root-transformed site-by-taxon matrices
#   (TRANSFORM_PRIMARY, as Methods 2.5 describes), with the untransformed matrices
#   run beside them (TRANSFORM_COMPARE) so that the effect of the transformation
#   is a table rather than an assumption. Within each region and weighting:
#   ANOSIM by condition class (N_PERM_ANOSIM permutations), PERMANOVA by class
#   and by the continuous B-IBI (N_PERM), BETADISPER by class (N_PERM). Across
#   regions: PERMANOVA with Region, B-IBI and their interaction, and BETADISPER
#   by region. The earlier analysis (v1) ran the untransformed matrices and
#   described them as transformed; under the transformation the Piedmont
#   separates by class under both weightings, and section 3.1 says so.
#
# CAUTION
#   Condition classes are the ICPRB ratings. On the rarefied-median basis two
#   classes hold a single site (Valley & Ridge Excellent; Coastal Plain Good), so
#   class-based tests are weakly supported there and BETADISPER within region is
#   not defined for those regions. The continuous-B-IBI PERMANOVA is run alongside
#   because it has no class-size constraint; it is the inference Results report.
#
# OUTPUT
#   output/tables/04_within_region_multivariate.csv
#   output/tables/04_cross_region_multivariate.csv    <- section 3.1 paragraph 1
#   output/tables/04_transformation_sensitivity.csv   <- transformed against untransformed, and the v1 reference
#   Step 07 combines the first two into 07_S4_multivariate_tests.csv, which step 18 writes as supplement
#   Table S5 (legacy Table S9; the class-based within-region tests are left out there).
# =============================================================================

say("STEP 04 -- within-region and cross-region multivariate tests")
set.seed(SEED)

meta <- dat %>% select(SiteID, Region, BIBI, BIBI_Class)

run_region <- function(rg, mat, weighting, transform) {
  m <- meta %>% filter(Region == rg)
  X <- mat[m$SiteID, , drop = FALSE]
  X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]
  Z <- if (transform == "sqrt") sqrt(X) else X
  D <- vegdist(Z, method = "bray")

  cls   <- factor(m$BIBI_Class)
  n_cls <- length(levels(cls))
  small <- any(table(cls) < 2)
  out   <- list()

  nm <- suppressWarnings(metaMDS(D, k = 2, trymax = 50, trace = 0))
  out[[length(out)+1]] <- tibble(Test = "NMDS", Term = "stress",
                                 Statistic = nm$stress, R2 = NA_real_, p = NA_real_)

  if (n_cls >= 2) {
    an <- anosim(D, cls, permutations = N_PERM_ANOSIM)
    out[[length(out)+1]] <- tibble(Test = "ANOSIM", Term = "B-IBI class",
                                   Statistic = an$statistic, R2 = NA_real_, p = an$signif)
    pc <- adonis2(D ~ cls, permutations = N_PERM, by = "terms")
    out[[length(out)+1]] <- tibble(Test = "PERMANOVA", Term = "B-IBI class",
                                   Statistic = pc$F[1], R2 = pc$R2[1], p = pc$`Pr(>F)`[1])
  }

  pb <- adonis2(D ~ m$BIBI, permutations = N_PERM, by = "terms")
  out[[length(out)+1]] <- tibble(Test = "PERMANOVA", Term = "B-IBI continuous",
                                 Statistic = pb$F[1], R2 = pb$R2[1], p = pb$`Pr(>F)`[1])

  if (n_cls >= 2 && !small) {
    bd <- betadisper(D, cls); pt <- permutest(bd, permutations = N_PERM)
    out[[length(out)+1]] <- tibble(Test = "BETADISPER", Term = "B-IBI class",
                                   Statistic = pt$tab$F[1], R2 = NA_real_,
                                   p = pt$tab$`Pr(>F)`[1])
  } else {
    out[[length(out)+1]] <- tibble(Test = "BETADISPER", Term = "B-IBI class",
                                   Statistic = NA_real_, R2 = NA_real_, p = NA_real_)
  }

  bind_rows(out) %>%
    mutate(Region = rg, Weighting = weighting, Transform = transform,
           n_sites = nrow(m), n_taxa = ncol(X), n_classes = n_cls,
           singleton_class = small, .before = 1)
}

grid <- expand_grid(Region = REGIONS,
                    Weighting = c("Density", "Biomass"),
                    Transform = unique(c(TRANSFORM_PRIMARY, TRANSFORM_COMPARE)))

multi <- pmap_dfr(grid, function(Region, Weighting, Transform) {
  run_region(Region, if (Weighting == "Density") den_mat else bio_mat, Weighting, Transform)
}) %>% arrange(Transform, Weighting, Region, Test, Term)

write_csv(multi, file.path(OUT_TAB, "04_within_region_multivariate.csv"))

# ---- Cross-region model -----------------------------------------------------
# WHY THIS IS HERE
#   Section 3.1's first paragraph reports the cross-region model: Region as
#   the dominant axis of variation in community structure, and the Region x
#   B-IBI interaction. That model is what justifies analyzing within region at
#   all. It is run here under the same transformation as the within-region
#   tests, so both paragraphs of section 3.1 rest on one basis.
#
#   BETADISPER by region is run alongside because Methods 2.5 says dispersion
#   heterogeneity is checked before a PERMANOVA centroid effect is interpreted.
#   With three regions of very different trait-pool breadth that check is not
#   a formality.

run_global <- function(mat, weighting, transform) {
  X <- mat[meta$SiteID, , drop = FALSE]
  X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]
  Z <- if (transform == "sqrt") sqrt(X) else X
  D <- vegdist(Z, method = "bray")

  ad  <- adonis2(D ~ Region * BIBI, data = meta, permutations = N_PERM, by = "terms")
  trm <- setdiff(rownames(ad), c("Residual", "Total"))
  i   <- match(trm, rownames(ad))
  out <- tibble(Test = "PERMANOVA", Term = trm,
                Statistic = ad$F[i], R2 = ad$R2[i], p = ad$`Pr(>F)`[i])

  bd <- betadisper(D, factor(meta$Region))
  pt <- permutest(bd, permutations = N_PERM)
  out <- bind_rows(out, tibble(Test = "BETADISPER", Term = "Region",
                               Statistic = pt$tab$F[1], R2 = NA_real_,
                               p = pt$tab$`Pr(>F)`[1]))

  out %>% mutate(Weighting = weighting, Transform = transform,
                 n_sites = nrow(meta), n_taxa = ncol(X), .before = 1)
}

global_grid <- expand_grid(Weighting = c("Density", "Biomass"),
                           Transform = unique(c(TRANSFORM_PRIMARY, TRANSFORM_COMPARE)))

cross <- pmap_dfr(global_grid, function(Weighting, Transform) {
  run_global(if (Weighting == "Density") den_mat else bio_mat, Weighting, Transform)
}) %>% arrange(Transform, Weighting, Test, Term)

write_csv(cross, file.path(OUT_TAB, "04_cross_region_multivariate.csv"))

say("  cross-region model (D ~ Region * B-IBI), primary transform '%s':", TRANSFORM_PRIMARY)
for (w in c("Density", "Biomass")) {
  r <- cross %>% filter(Weighting == w, Transform == TRANSFORM_PRIMARY, Test == "PERMANOVA")
  for (i in seq_len(nrow(r)))
    say("    %-8s %-14s F = %6.2f  R2 = %.3f  p = %.3f", w, r$Term[i], r$Statistic[i], r$R2[i], r$p[i])
  b <- cross %>% filter(Weighting == w, Transform == TRANSFORM_PRIMARY, Test == "BETADISPER")
  say("    %-8s BETADISPER     F = %6.2f              p = %.3f%s", w, b$Statistic[1], b$p[1],
      ifelse(b$p[1] <= 0.05,
             "   <<< dispersion heterogeneous: read Region as spread AND centroid",
             ""))
}

# ---- Reference: the v1 ANOSIM values, computed untransformed on the observed index and its classes ------
# They can reproduce only on that basis (CH2_BIBI_MEDIAN unset); on the published median basis the classes
# differ, so the comparison is not applicable and the column is NA.
published_31 <- tribble(
  ~Region,          ~Weighting, ~R_published, ~R2_published,
  "Coastal Plain",  "Density",   0.235,        0.083,
  "Coastal Plain",  "Biomass",   0.263,        0.051,
  "Piedmont",       "Density",   0.137,        0.101,
  "Piedmont",       "Biomass",   0.301,        0.090,
  "Valley & Ridge", "Density",   0.068,        0.024,
  "Valley & Ridge", "Biomass",   0.101,        0.020
)

sens <- multi %>%
  filter(Test %in% c("ANOSIM", "PERMANOVA"), Term %in% c("B-IBI class", "B-IBI continuous")) %>%
  select(Region, Weighting, Transform, Test, Term, Statistic, R2, p) %>%
  pivot_wider(names_from = Transform, values_from = c(Statistic, R2, p))

anosim_chk <- multi %>% filter(Test == "ANOSIM", Transform == TRANSFORM_COMPARE) %>%
  select(Region, Weighting, R_none = Statistic, p_none = p) %>%
  left_join(multi %>% filter(Test == "ANOSIM", Transform == TRANSFORM_PRIMARY) %>%
              select(Region, Weighting, R_sqrt = Statistic, p_sqrt = p),
            by = c("Region","Weighting")) %>%
  left_join(published_31 %>% select(Region, Weighting, R_published),
            by = c("Region","Weighting")) %>%
  mutate(reproduces_published = if (nzchar(Sys.getenv("CH2_BIBI_MEDIAN"))) NA else abs(R_none - R_published) < 0.005,
         sig_flips_under_sqrt = (p_none > 0.05) != (p_sqrt > 0.05))

write_csv(anosim_chk, file.path(OUT_TAB, "04_transformation_sensitivity.csv"))

say("  ANOSIM | primary = '%s' (as Methods 2.5 describes) | untransformed = '%s' | v1 reference",
    TRANSFORM_PRIMARY, TRANSFORM_COMPARE)
for (i in seq_len(nrow(anosim_chk))) {
  a <- anosim_chk[i, ]
  say("    %-15s %-8s sqrt R=%.3f (p=%.3f) | none R=%.3f (p=%.3f) | v1 %.3f  %s%s",
      a$Region, a$Weighting, a$R_sqrt, a$p_sqrt, a$R_none, a$p_none, a$R_published,
      ifelse(is.na(a$reproduces_published), "", ifelse(a$reproduces_published, "v1 reproduces", "*** v1 DOES NOT REPRODUCE ***")),
      ifelse(a$sig_flips_under_sqrt, "   <<< significance differs between transformations", ""))
}

if (all(is.na(anosim_chk$reproduces_published))) {
  say("  v1 reference values were computed on the observed index and its classes; not applicable on the published basis")
} else if (!all(anosim_chk$reproduces_published)) {
  say("  WARNING: the v1 ANOSIM values do not reproduce untransformed on the observed basis. Investigate.")
} else {
  say("  v1 reference reproduces: all six ANOSIM values reproduce untransformed on the observed basis")
}
if (any(anosim_chk$sig_flips_under_sqrt)) {
  chg <- anosim_chk %>% filter(sig_flips_under_sqrt)
  say("  significance at 0.05 differs between the transformations for: %s",
      paste(sprintf("%s %s", chg$Region, chg$Weighting), collapse = "; "))
}
say("STEP 04 complete")
