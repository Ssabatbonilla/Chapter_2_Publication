# =============================================================================
# 17_null_model_ses.R
# =============================================================================
# THE ALTERNATIVE EXPLANATION THIS STEP TESTS
#
#   The paper's central functional finding is that trait-space geometry
#   responds under DENSITY weighting and not under biomass. Step 02 then
#   found that CONCENTRATION -- Gini and N50 -- behaves identically:
#   density-only, flat in biomass, in the same two regions.
#
#   Every geometry metric is abundance-weighted. So there is a live rival
#   explanation for that finding: the geometry response may BE the evenness change that
#   Gini and N50 already measure, re-expressed in trait units. If so, the paper
#   reports one phenomenon twice and calls the second one functional.
#
#   A trait-label-shuffle null separates them. Hold the community matrix
#   exactly as observed -- same richness, same abundance distribution, same
#   evenness -- and permute which taxon carries which trait profile. Anything
#   the observed metric retains over that null is attributable to WHICH TRAITS
#   the abundant taxa actually have, not to how abundance is distributed.
#
#     SES = (observed - mean(null)) / sd(null)
#
#   If rho(SES, B-IBI) holds up, the geometry response is a trait phenomenon.
#   If it collapses toward zero while raw rho stays strong, the geometry result
#   restates the concentration result.
#
# WHY THIS COULD NOT BE RUN BEFORE 2026-08-27
#   A null model recomputes the metric thousands of times, which requires
#   knowing exactly how the metric was computed. That was not reproducible
#   until step 13. This step uses the step 13 geometry and verifies against it
#   before permuting anything.
#
# SCOPE
#   RaoQ and FDis only. Both are computed on the uncorrected Gower distance (see
#   step 13), so both have exact analytic forms and N_PERM_SES permutations cost
#   seconds rather than the ~40 minutes dbFD would take. FDiv and FRic are
#   convex-hull quantities and are left out; if they are wanted later the same
#   machinery works, at roughly 3.5 hours as a background job.
#
# INPUTS
#   den_mat and dat (step 01); IN$traits_fuz, IN$traits_bin and IN$trait_blocks
#   (the step 13 trait space); 02_concentration_site_level.csv (step 02).
#   Density weighting only.
#
# OUTPUT
#   output/tables/17_ses_vs_raw.csv          (Table S6; legacy Table S15 in step 18)
#   output/tables/17_ses_site_level.csv
# =============================================================================

say("STEP 17 -- null-model SES: is the geometry response a trait signal or evenness?")
suppressPackageStartupMessages({ library(FD); library(ade4) })
set.seed(SEED)

# ---- 1. Rebuild the step 13 trait space -------------------------------------
.std <- function(df) {
  if (!"TAXANAME" %in% names(df)) df$TAXANAME <- df$Genus
  if (!"Genus"    %in% names(df)) df$Genus    <- df$TAXANAME
  df %>% mutate(Genus = stringr::str_squish(as.character(Genus)),
                TAXANAME = stringr::str_squish(as.character(TAXANAME))) %>%
    filter(!is.na(TAXANAME), TAXANAME != "") %>% distinct(TAXANAME, .keep_all = TRUE)
}
tf <- .std(readRDS(IN$traits_fuz)); tbn <- .std(readRDS(IN$traits_bin))
tbk <- readRDS(IN$trait_blocks)
idc <- c("Genus", "TAXANAME")
cl  <- intersect(setdiff(names(tbn), idc), setdiff(names(tf), idc))
cl  <- unique(cl[vapply(tbn[cl], is.numeric, logical(1)) & vapply(tf[cl], is.numeric, logical(1))])
tx  <- intersect(colnames(den_mat), intersect(tbn$TAXANAME, tf$TAXANAME))
.bd <- function(df) { d <- df %>% filter(TAXANAME %in% tx) %>% distinct(TAXANAME, .keep_all = TRUE)
  d <- d[match(tx, d$TAXANAME), , drop = FALSE]; m <- as.matrix(d[, cl, drop = FALSE])
  rownames(m) <- d$TAXANAME; storage.mode(m) <- "double"; m }   # NA kept: a trait group with no record stays missing (22 Sep 2026)
MF <- .bd(tf); MB <- .bd(tbn)
kp <- (colSums(MB, na.rm = TRUE) > 0) | (colSums(MF, na.rm = TRUE) > 0); MF <- MF[, kp, drop = FALSE]
gw <- tbk$gowdis_w[colnames(MF)]

# CRITICAL, and not what the dbFD documentation implies: RaoQ and FDis are
# computed from the UNCORRECTED Gower distance. Read FD::dbFD source --
#   line 577  RaoQ <- divc(data.frame(t(a)), x.dist, ...)
#   line 580  disp <- fdisp(x.dist, a)
# both take x.dist, NOT the corrected x.dist2. The sqrt/cailliez correction
# only ever reaches FDiv, FEve and FRic, which go through the PCoA at line 303.
# That is also why step 13 found RaoQ and FDis identical under both
# corrections -- the correction is irrelevant to them.
#
# divc() squares the distance internally, so RaoQ = sum_ij p_i p_j d_ij^2 / 2.
gd <- suppressWarnings(gowdis(as.data.frame(MF), w = gw))
G2 <- as.matrix(gd)^2

# ---- 2. Analytic RaoQ and FDis, verified against dbFD before use ------------
A <- den_mat[dat$SiteID, rownames(MF), drop = FALSE]   # MF is taxa x modalities
P <- A / rowSums(A)

raoq_of <- function(M) rowSums((P %*% M) * P) / 2
fdis_of <- function(d)  FD::fdisp(d, as.matrix(A))$FDis
obs_rao <- raoq_of(G2); obs_fdis <- fdis_of(gd)

ref <- suppressWarnings(FD::dbFD(x = MF, a = A,
                                 w = gw, stand.x = TRUE, corr = "cailliez", w.abun = TRUE,   # RaoQ and FDis do not depend on the correction; sqrt fails with NA groups
                                 calc.FRic = FALSE, calc.FDiv = FALSE, calc.CWM = FALSE,
                                 m = 10, messages = FALSE))
chk <- c(RaoQ = max(abs(obs_rao  - ref$RaoQ[dat$SiteID])),
         FDis = max(abs(obs_fdis - ref$FDis[dat$SiteID])))
say("  analytic vs dbFD: RaoQ max diff = %.2e | FDis max diff = %.2e", chk[1], chk[2])
if (any(chk > 1e-8)) {
  say("  ANALYTIC FORM DOES NOT MATCH dbFD -- the null would be computed in the")
  say("  wrong geometry, so it is not run. Nothing written. Fix this first.")
  say("STEP 17 ABORTED")
} else {

# ---- 3. Trait-label shuffle null --------------------------------------------
NP <- N_PERM_SES                           # declared in 00_config.R (27 Sep 2026)
nR <- matrix(NA_real_, nrow(P), NP); nF <- nR
gdm <- as.matrix(gd)
for (k in seq_len(NP)) {
  i <- sample(ncol(P))                      # permute taxon labels on the traits
  nR[, k] <- raoq_of(G2[i, i])
  Dp <- gdm[i, i]; dimnames(Dp) <- dimnames(gdm)
  nF[, k] <- fdis_of(as.dist(Dp))
  if (k %% 250 == 0) say("    ... %d of %d permutations", k, NP)
}
ses <- function(o, n) (o - rowMeans(n)) / apply(n, 1, sd)
s_rao <- ses(obs_rao, nR); s_fdis <- ses(obs_fdis, nF)

site_tab <- tibble(SiteID = dat$SiteID, Region = dat$Region, BIBI = dat$BIBI,
                   RaoQ_den = obs_rao, SES_RaoQ = s_rao,
                   FDis_den = obs_fdis, SES_FDis = s_fdis)
write_csv(site_tab, file.path(OUT_TAB, "17_ses_site_level.csv"))

cmp <- purrr::map_dfr(REGIONS, function(rg) {
  s <- dat$Region == rg
  purrr::map_dfr(list(list("RaoQ_den", obs_rao, s_rao), list("FDis_den", obs_fdis, s_fdis)),
    function(z) {
      a <- suppressWarnings(cor.test(dat$BIBI[s], z[[2]][s], method = "spearman", exact = FALSE))
      b <- suppressWarnings(cor.test(dat$BIBI[s], z[[3]][s], method = "spearman", exact = FALSE))
      tibble(Region = rg, Metric = z[[1]],
             rho_raw = unname(a$estimate), p_raw = a$p.value,
             rho_SES = unname(b$estimate), p_SES = b$p.value,
             # meaningless when the raw coefficient is ~0, so do not report it
             pct_retained = ifelse(abs(unname(a$estimate)) < RHO_MIN, NA_real_,
                                   100 * abs(unname(b$estimate)) / abs(unname(a$estimate))),
             verdict = case_when(
               abs(unname(a$estimate)) < RHO_MIN            ~ "raw not responsive",
               abs(unname(b$estimate)) >= SES_RETAIN_SURVIVES * abs(unname(a$estimate)) ~ "trait signal survives",
               abs(unname(b$estimate)) >= SES_RETAIN_PARTLY * abs(unname(a$estimate)) ~ "PARTLY evenness-driven",
               TRUE                                          ~ "COLLAPSES -- evenness, not traits"))
    })
})
write_csv(cmp, file.path(OUT_TAB, "17_ses_vs_raw.csv"))

say("")
say("  raw rho vs SES rho (trait-label shuffle, %d permutations):", NP)
for (i in seq_len(nrow(cmp)))
  say("    %-15s %-9s raw %+.3f (p %.3f) -> SES %+.3f (p %.3f) | %s retained | %s",
      cmp$Region[i], cmp$Metric[i], cmp$rho_raw[i], cmp$p_raw[i],
      cmp$rho_SES[i], cmp$p_SES[i],
      ifelse(is.na(cmp$pct_retained[i]), "  n/a", sprintf("%3.0f%%", cmp$pct_retained[i])),
      cmp$verdict[i])

conc <- read_csv(file.path(OUT_TAB, "02_concentration_site_level.csv"), show_col_types = FALSE)
say("")
say("  cross-check -- how strongly does raw geometry track numerical evenness?")
for (rg in REGIONS) {
  s  <- dat$Region == rg
  n50 <- conc$N50_den[match(dat$SiteID[s], conc$SiteID)]
  say("    %-15s rho(RaoQ, N50_den) = %+.3f | rho(FDis, N50_den) = %+.3f", rg,
      cor(obs_rao[s], n50, method = "spearman", use = "complete.obs"),
      cor(obs_fdis[s], n50, method = "spearman", use = "complete.obs"))
}
say("STEP 17 complete")
}
