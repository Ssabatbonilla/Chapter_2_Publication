# =============================================================================
# 15_rlq_fourthcorner.R
# =============================================================================
# NOT REPORTED IN THIS PAPER
#   The RLQ and fourth-corner results left this paper for a companion trait-state manuscript (legacy Table 4,
#   Tables S16, S17 and S19, and Figure 6; step 18b copies them to output/manuscript/trait_state_manuscript/).
#   The step still runs because later steps read its tables: step 20 asserts its fourth-corner rerun against
#   15_fourthcorner_per_block*.csv and stops without them, step 15b reproduces the same statistics, and steps 06
#   and 18 build the companion manuscript's figure and tables from them.
#
# WHAT IT ASKS
#   RLQ with the fourth-corner test relates trait states to the condition gradient directly, a constrained
#   analysis (Doledec et al. 1996; Dray & Legendre 2008; Dray et al. 2014; Chessman 2015).
#
# INPUTS
#   den_mat or bio_mat and dat (step 01); IN$traits (the fuzzy trait file); IN$site_metrics (the CWM columns of
#   section 4); N_PERM_FC, N_SEEDS_FC, FDR_MAX and SEED from 00_config.R.
#
# 24 SEPTEMBER 2026. Land use removed from the environment table; see the note above run_region().
#
# 22 SEPTEMBER 2026. The trait matrix is the Twardochleb fuzzy file; permutations are N_PERM_FC (9999) with the
# seed written to 15_permutation_settings.csv; the RLQ trait table is a fuzzy correspondence analysis (dudi.fca on
# prep.fuzzy.var), the ordination for fuzzy-coded traits (Chevenet et al. 1994; Doledec et al. 1996).
#
# WHAT IT DOES (as run for the paper; the log prints every count)
#
#   PRIMARY: PER-BLOCK. One RLQ and fourth-corner analysis per region-by-trait-group
#   block, on every taxon that has a record for that group (a trait group with no
#   record is NA for a taxon, so a complete-case analysis would discard the taxon
#   for every group). The fourth-corner runs ade4's model 6: model 2 (sites
#   permuted) and model 4 (taxa permuted), each Benjamini-Hochberg adjusted within
#   the block, and the larger adjusted p reported for each state (Dray et al.
#   2014). N_PERM_FC permutations under N_SEEDS_FC seeds; a state passes only when
#   its adjusted p meets FDR_MAX under every seed, and is borderline when it passes
#   under some. Every state carries its carrier count (present taxa with a nonzero
#   affinity); every passing state gets a leave-one-taxon-out.
#
#   SENSITIVITIES, reported beside the primary: (a) the correction across all
#   states of a region (Benjamini-Hochberg on the model-6 p, per seed); (b) one
#   multivariate test per trait group (randtest.rlq under model 6), corrected
#   across the region's groups, with the states in a fuzzy correspondence analysis
#   and, separately, standardized; (c) the complete-case run (taxa with every group
#   resolved), corrected across the region's states; (d) the coverage confound test
#   on every passing state; (e) the community-weighted mean of each named state
#   against the fourth-corner statistic, in this run's weighting. The difference
#   between the weightings is tested in 15b.
#
#   A state that passes model 2 and fails model 4 is carried by abundance shifts
#   among taxa whose traits are unremarkable relative to the pool; that is not
#   trait-based environmental filtering. Trait associations are reported as
#   standing-stock composition, not as filtering.
#
# OUTPUT (per weighting; the biomass files carry the suffix _biomass; legacy numbers, none in this paper)
#   output/tables/15_fourthcorner_per_block.csv           (legacy Table 4, Figure 6 and Table S16; read by steps 15b and 20)
#   output/tables/15_fourthcorner_per_block_by_seed.csv   (every per-block test under each seed)
#   output/tables/15_fourthcorner_regionwide.csv          (sensitivity a; legacy Table S16)
#   output/tables/15_fourthcorner_pergroup.csv, 15_fourthcorner_pergroup_by_seed.csv   (sensitivity b; legacy Table S17)
#   output/tables/15_leave_one_taxon_out.csv              (legacy Table S19)
#   output/tables/15_fourthcorner_coverage_confound.csv   (sensitivity d)
#   output/tables/15_rlq_global_tests.csv, 15_fourthcorner_by_model.csv   (sensitivity c; diagnostic)
#   output/tables/15_cwm_vs_fourthcorner.csv              (sensitivity e; diagnostic)
#   output/tables/15_permutation_settings.csv             (permutations, seeds, rules and versions)
# =============================================================================

say("STEP 15 -- RLQ and fourth-corner: constrained trait-environment analysis")
suppressPackageStartupMessages(library(ade4))
set.seed(SEED)
say("  permutations: %d per test (N_PERM_FC), seed %d", N_PERM_FC, SEED)

# ---- 0. Currency switch -----------------------------------------------------
# WHY THIS EXISTS. This step was written on the density matrix with no comment
# explaining the choice, in a paper whose H3 is that currency changes what is
# detected. Nothing in RLQ requires counts: the L table enters through a
# correspondence analysis, which needs a non-negative matrix, and bio_mat has
# identical row and column names to den_mat, no NAs and no zero margins. Density stays PRIMARY on substantive grounds -- correspondence
# analysis gives large-bodied taxa disproportionate leverage on the axes, and
# biomass carries length-mass regression error that counts do not -- so the
# biomass run is a secondary analysis, reported as such.
#
#   Rscript run_all.R 15                          BOTH weightings, biomass then density
#                                                 (the runner loops; since 21 Sept 2026)
#   CH2_FC_CURRENCY=biomass Rscript run_all.R 15   biomass only (writes *_biomass.csv)
#   CH2_FC_CURRENCY=density Rscript run_all.R 15   density only
#
# Before 21 September 2026 the default was density only, so the documented run never rewrote the
# biomass files. The runner prints a banner naming the weightings. The density run writes the
# unsuffixed filenames; the biomass run writes the same files with the suffix _biomass.
CURRENCY <- Sys.getenv("CH2_FC_CURRENCY", "density")
if (!CURRENCY %in% c("density", "biomass"))
  stop("CH2_FC_CURRENCY must be 'density' or 'biomass', got: ", CURRENCY)
W   <- if (CURRENCY == "density") den_mat else bio_mat
SUF <- if (CURRENCY == "density") "" else "_biomass"
stopifnot(identical(dimnames(W), dimnames(den_mat)))
say("  currency: %s (%d sites x %d taxa)", CURRENCY, nrow(W), ncol(W))

# ---- 1. Build R, L, Q -------------------------------------------------------
tr_raw <- as.data.frame(readRDS(IN$traits))
rownames(tr_raw) <- tr_raw$TAXANAME
# Twardochleb fuzzy traits, 22 Sep 2026: Q is the Twardochleb fuzzy file; its identifier columns are Genus and TAXANAME.
Q_raw  <- as.matrix(tr_raw[colnames(W), setdiff(names(tr_raw), c("Genus", "TAXANAME")), drop = FALSE])
storage.mode(Q_raw) <- "double"
blk    <- sub("_.*$", "", colnames(Q_raw))

# Row-normalize within each trait block. A block with no record is NA and stays NA; is.na() catches
# both NA and NaN, so a missing block is never read as a zero affinity.
Q_norm <- Q_raw
for (b in unique(blk)) {
  i <- blk == b
  Q_norm[, i] <- Q_raw[, i, drop = FALSE] / rowSums(Q_raw[, i, drop = FALSE])
}
incomplete <- apply(Q_norm, 1, function(r) any(is.na(r)))
say("  taxa with >= 1 unresolved trait block: %d of %d (%.1f%% of %s)",
    sum(incomplete), nrow(Q_norm),
    100 * sum(colSums(W)[incomplete]) / sum(W), CURRENCY)

# 24 SEPTEMBER 2026: the environment table is the index alone. Until this date the RLQ and the
# complete-case fourth-corner carried land-use PC1 beside the index, so the global RLQ test and the complete-case
# family (82 tests per region) included land use. No analysis in this paper uses land use; it is the gradient of a companion study.
# With one environmental variable the RLQ has one axis, so rlq_axis1_pct is 100 by construction.
env <- dat

# ---- 2. Per-region RLQ + fourth-corner, complete-case taxa ------------------
run_region <- function(rg) {
  m <- env %>% filter(Region == rg)
  L <- W[m$SiteID, , drop = FALSE]
  L <- L[, colSums(L) > 0, drop = FALSE]
  L <- L[, !incomplete[colnames(L)], drop = FALSE]
  L <- L[, colSums(L) > 0, drop = FALSE]
  Q <- as.data.frame(Q_norm[colnames(L), , drop = FALSE])
  R <- data.frame(BIBI = m$BIBI)                     # the index alone (24 Sep 2026)

  # randtest.rlq() re-resolves BOTH its dudi arguments AND the data frames
  # those dudis were built from, by name, via eval.parent(). All of them must
  # be visible in the global environment or the call fails with
  # "object 'cR' not found" / "object 'R' not found". Assign, test, remove.
  assign("L", L, envir = globalenv())
  assign("R", R, envir = globalenv())
  assign("Q", Q, envir = globalenv())
  assign("cL", dudi.coa(as.data.frame(L), scannf = FALSE, nf = 2), envir = globalenv())
  assign("cR", dudi.pca(R, row.w = cL$lw, scannf = FALSE, nf = 2), envir = globalenv())
  # Fuzzy-coded traits enter RLQ through a fuzzy correspondence analysis (22 Sep 2026; was dudi.pca
  # on the affinities). A state no complete-case taxon of the region carries has zero weight in FCA and is
  # dropped from the ordination only; the fourth-corner below tests the full Q.
  qk <- colSums(Q) > 0
  qb <- rle(blk[qk]); stopifnot(!any(duplicated(qb$values)))           # blocks contiguous
  assign("Qf", prep.fuzzy.var(Q[, qk, drop = FALSE], col.blocks = setNames(qb$lengths, qb$values),
                              row.w = cL$cw), envir = globalenv())
  assign("cQ", dudi.fca(Qf, scannf = FALSE, nf = 2), envir = globalenv())
  rl <- rlq(cR, cL, cQ, scannf = FALSE, nf = 2)
  rt <- randtest(rl, nrepet = N_PERM_FC, modeltype = 6)
  p_rlq <- if (inherits(rt, "krandtest")) max(rt$pvalue) else rt$pvalue
  rm(list = c("L", "R", "Q", "Qf", "cL", "cR", "cQ"), envir = globalenv())

  glob <- tibble(Region = rg, n_sites = nrow(L), n_taxa = ncol(L), q_ordination = "dudi.fca",
                 n_states_in_fca = sum(qk),
                 pct_zeros = 100 * mean(L == 0),
                 rlq_axis1_pct = 100 * rl$eig[1] / sum(rl$eig),
                 rlq_p_model6  = p_rlq)

  # Fourth-corner under each model separately, then combined.
  fcs <- lapply(c(2, 4, 6), function(mt)
    fourthcorner(R, as.data.frame(L), Q, modeltype = mt, nrepet = N_PERM_FC,
                 p.adjust.method.G = "fdr", p.adjust.method.D = "fdr"))
  names(fcs) <- paste0("model", c(2, 4, 6))

  tabs <- purrr::imap_dfr(fcs, function(fc, nm) {
    tibble(Region = rg, model = nm, term = fc$tabD$names,
           stat = fc$tabD$obs, p_raw = fc$tabD$pvalue,
           p_fdr = fc$tabD$adj.pvalue)
  })
  list(glob = glob, tabs = tabs)
}

res <- lapply(REGIONS, run_region); names(res) <- REGIONS
rlq_glob <- bind_rows(lapply(res, `[[`, "glob"))
fc_all   <- bind_rows(lapply(res, `[[`, "tabs"))

write_csv(rlq_glob, file.path(OUT_TAB, sprintf("15_rlq_global_tests%s.csv", SUF)))
write_csv(tibble(setting = c("permutations per test", "seed", "per-block seeds", "per-block pass rule", "fourth-corner model (per block)", "FDR adjustment (per block)",
                             "RLQ trait ordination", "ade4 version", "R version", "weighting", "primary family (L37)", "sensitivity tiers (L37)", "leave-one-taxon-out (FL1)"),
                 value   = c(N_PERM_FC, SEED, sprintf("%d to %d", SEED, SEED + N_SEEDS_FC - 1),
                             sprintf("FDR <= %.2f under all %d seeds", FDR_MAX, N_SEEDS_FC), "6 (models 2 and 4 combined)",
                             "ade4 model 6: models 2 (sites permuted) and 4 (taxa permuted) each Benjamini-Hochberg adjusted within the region-by-trait-group block, the larger adjusted p reported (Dray et al. 2014)",
                             "dudi.fca on prep.fuzzy.var", as.character(packageVersion("ade4")),
                             paste(R.version$major, R.version$minor, sep = "."), CURRENCY,
                             "within each region-by-trait-group block (the tier of hypothesis 3's state-by-state predictions)",
                             "across all states of a region (Benjamini-Hochberg on the model-6 p, per seed, pass under every seed); one multivariate test per trait group (fuzzy correspondence analysis, and states standardized), BH across the region's groups",
                             "every carrier of each passing state removed in turn; the most influential removal rerun under model 6 at N_PERM_FC (pipeline seed)")),
          file.path(OUT_TAB, sprintf("15_permutation_settings%s.csv", SUF)))
write_csv(fc_all,   file.path(OUT_TAB, sprintf("15_fourthcorner_by_model%s.csv", SUF)))

say("")
say("  RLQ global test (models 2 and 4 combined):")
for (i in seq_len(nrow(rlq_glob)))
  say("    %-15s L = %2d x %3d (%.0f%% zeros) | axis 1 = %.1f%% | p = %.3f %s",
      rlq_glob$Region[i], rlq_glob$n_sites[i], rlq_glob$n_taxa[i],
      rlq_glob$pct_zeros[i], rlq_glob$rlq_axis1_pct[i], rlq_glob$rlq_p_model6[i],
      ifelse(rlq_glob$rlq_p_model6[i] <= 0.05, "", "n.s."))

say("")
say("  fourth-corner, tests at raw p <= 0.05 (chance expectation in brackets):")
smry <- fc_all %>% group_by(Region, model) %>%
  summarise(n_tests = n(), n_raw = sum(p_raw <= 0.05),
            n_fdr = sum(p_fdr <= FDR_MAX), .groups = "drop")
for (i in seq_len(nrow(smry)))
  say("    %-15s %-7s %2d of %2d [%.1f] | surviving FDR: %d",
      smry$Region[i], smry$model[i], smry$n_raw[i], smry$n_tests[i],
      0.05 * smry$n_tests[i], smry$n_fdr[i])

# ---- 3. Per-block analysis, which keeps the taxa complete-case discards -----
# Five seeds (22 Sep 2026). At N_PERM_FC permutations some states still passed FDR_MAX
# under one seed and not another, each with its adjusted p within about 0.02 of the line; more permutations cannot
# settle a p that sits on it. So the per-block analysis runs under seeds SEED to SEED + N_SEEDS_FC - 1, each set
# explicitly (identical results run in parallel or not), and a state PASSES only when it meets FDR_MAX under every
# seed. The statistic does not depend on the seed (asserted); only the p-values do. The primary columns (p_raw,
# p_fdr) are the SEED run; n_seeds_pass, p_fdr_min, p_fdr_max and verdict carry the rule.
per_block_one_seed <- function(k) {
  set.seed(SEED + k)
  purrr::map_dfr(REGIONS, function(rg) {
    m <- env %>% filter(Region == rg)
    purrr::map_dfr(unique(blk), function(b) {
      i  <- blk == b
      ok <- !apply(Q_norm[, i, drop = FALSE], 1, function(r) any(is.na(r)))
      L  <- W[m$SiteID, , drop = FALSE]
      L  <- L[, colSums(L) > 0 & ok[colnames(L)], drop = FALSE]
      if (ncol(L) < 5) return(NULL)
      Q  <- as.data.frame(Q_norm[colnames(L), i, drop = FALSE])
      R  <- data.frame(BIBI = m$BIBI)
      fc <- fourthcorner(R, as.data.frame(L), Q, modeltype = 6, nrepet = N_PERM_FC,
                         p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
      ncar <- colSums(Q > 0)[sub("^BIBI / ", "", fc$tabD$names)]                 # carrier taxa: present taxa with a nonzero affinity
      tibble(seed = SEED + k, Region = rg, block = b, n_taxa = ncol(L),
             pct_cov = 100 * sum(W[m$SiteID, colnames(L)]) /
                       sum(W[m$SiteID, ]),
             term = fc$tabD$names, stat = fc$tabD$obs,
             p_raw = fc$tabD$pvalue, p_fdr = fc$tabD$adj.pvalue, n_carriers = as.integer(ncar))
    })
  })
}
seeds_k <- seq_len(N_SEEDS_FC) - 1
pb_runs <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(seeds_k, per_block_one_seed, mc.cores = min(N_SEEDS_FC, max(1, parallel::detectCores() - 1)))
} else lapply(seeds_k, per_block_one_seed)
if (!all(vapply(pb_runs, is.data.frame, logical(1)))) stop("Step 15: a per-block seed run failed: ", paste(pb_runs[!vapply(pb_runs, is.data.frame, logical(1))], collapse = " | "))
pb_all <- bind_rows(pb_runs)
stopifnot(identical(sort(unique(pb_all$seed)), SEED + seeds_k))
stat_spread <- pb_all %>% group_by(Region, block, term) %>% summarise(s = max(stat) - min(stat), n = n(), .groups = "drop")
stopifnot(max(stat_spread$s) < 1e-12, all(stat_spread$n == N_SEEDS_FC))       # statistic never depends on the seed
write_csv(pb_all, file.path(OUT_TAB, sprintf("15_fourthcorner_per_block_by_seed%s.csv", SUF)))
per_block <- pb_all %>% filter(seed == SEED) %>% select(-seed) %>%
  left_join(pb_all %>% group_by(Region, block, term) %>%
              summarise(n_seeds_pass = sum(p_fdr <= FDR_MAX), p_fdr_min = min(p_fdr), p_fdr_max = max(p_fdr), .groups = "drop"),
            by = c("Region", "block", "term")) %>%
  mutate(seeds = sprintf("%d to %d", SEED, SEED + N_SEEDS_FC - 1),
         verdict = case_when(n_seeds_pass == N_SEEDS_FC ~ "passes under every seed",
                             n_seeds_pass > 0            ~ sprintf("borderline, passes under %d of %d seeds", n_seeds_pass, N_SEEDS_FC),
                             TRUE                        ~ "does not pass"))
say("  per-block fourth-corner under %d seeds: %d states pass under every seed, %d borderline", N_SEEDS_FC,
    sum(per_block$n_seeds_pass == N_SEEDS_FC), sum(per_block$n_seeds_pass > 0 & per_block$n_seeds_pass < N_SEEDS_FC))
pb_out <- per_block
names(pb_out)[names(pb_out) == "pct_cov"] <- paste0("pct_", CURRENCY, "_covered")
write_csv(pb_out, file.path(OUT_TAB, sprintf("15_fourthcorner_per_block%s.csv", SUF)))

# ---- 3b. The two stricter multiplicity tiers, reported beside the block family (24 Sep 2026) ---------------
# Within each region-by-trait-group block is the primary family, the tier at which hypothesis 3's state-by-state
# predictions were made. Two stricter tiers are written as sensitivities: (a) the same bivariate tests
# corrected across all states of a region, per seed, a state passing only under every seed; (b) one multivariate
# fourth-corner test per trait group (the RLQ total inertia, Dray & Legendre 2008; randtest.rlq under model 6, the
# larger of the site- and taxon-permutation p, ter Braak et al. 2012), corrected across the region's groups, per seed.
# (b) is run with the group's states in a fuzzy correspondence analysis (dudi.fca, the rule fixed first) and with
# the states standardized (dudi.pca with abundance row weights, the weighting ade4::fourthcorner2 uses), because the
# fuzzy-correspondence weighting gives a rare state most of the inertia and loses power under taxon permutation;
# the standardized run is the declared sensitivity statistic.
rw <- pb_all %>% group_by(seed, Region) %>% mutate(p_fdr_regionwide = p.adjust(p_raw, "BH"), n_tests_regionwide = n()) %>% ungroup() %>%
  group_by(Region, block, term, stat) %>%
  summarise(n_tests_regionwide = first(n_tests_regionwide), n_seeds_pass_regionwide = sum(p_fdr_regionwide <= FDR_MAX),
            p_fdr_regionwide_min = min(p_fdr_regionwide), p_fdr_regionwide_max = max(p_fdr_regionwide), .groups = "drop") %>%
  mutate(verdict_regionwide = case_when(n_seeds_pass_regionwide == N_SEEDS_FC ~ "passes under every seed",
                                        n_seeds_pass_regionwide > 0 ~ sprintf("borderline, passes under %d of %d seeds", n_seeds_pass_regionwide, N_SEEDS_FC),
                                        TRUE ~ "does not pass"))
write_csv(rw, file.path(OUT_TAB, sprintf("15_fourthcorner_regionwide%s.csv", SUF)))
say("  region-wide correction (all states of a region, %d seeds): %d states pass under every seed", N_SEEDS_FC, sum(rw$n_seeds_pass_regionwide == N_SEEDS_FC))

pergroup_one_seed <- function(k, q_ord) {
  set.seed(SEED + k); out <- list()
  for (rg in REGIONS) {
    m <- env %>% filter(Region == rg)
    for (b in unique(blk)) {
      i <- blk == b; ok <- !apply(Q_norm[, i, drop = FALSE], 1, function(r) any(is.na(r)))
      L0 <- W[m$SiteID, , drop = FALSE]; L0 <- L0[, colSums(L0) > 0 & ok[colnames(L0)], drop = FALSE]
      if (ncol(L0) < 5) next
      Qb <- as.data.frame(Q_norm[colnames(L0), i, drop = FALSE]); qk <- colSums(Qb) > 0
      if (sum(qk) < 2) { out[[length(out) + 1]] <- tibble(seed = SEED + k, q_ordination = q_ord, Region = rg, block = b, n_taxa = ncol(L0), n_states = sum(qk),
                                                          trace = NA_real_, p_sites = NA_real_, p_taxa = NA_real_, p_raw = NA_real_, note = "fewer than two states with a carrier"); next }
      # randtest.rlq re-resolves its dudi arguments and their data frames by name in the global environment (section 2 above)
      assign("L", as.data.frame(L0), envir = globalenv()); assign("R", data.frame(BIBI = m$BIBI), envir = globalenv()); assign("Q", Qb, envir = globalenv())
      assign("cL", dudi.coa(L, scannf = FALSE, nf = 2), envir = globalenv())
      assign("cR", dudi.pca(R, row.w = cL$lw, scannf = FALSE, nf = 2), envir = globalenv())
      if (q_ord == "fca") {
        assign("Qf", prep.fuzzy.var(Q[, qk, drop = FALSE], col.blocks = setNames(sum(qk), b), row.w = cL$cw), envir = globalenv())
        assign("cQ", dudi.fca(Qf, scannf = FALSE, nf = 2), envir = globalenv())
      } else {
        assign("Qf", Q[, qk, drop = FALSE], envir = globalenv())
        assign("cQ", dudi.pca(Qf, row.w = cL$cw, scannf = FALSE, nf = 2), envir = globalenv())
      }
      rl <- rlq(cR, cL, cQ, scannf = FALSE, nf = 2); rt <- randtest(rl, nrepet = N_PERM_FC, modeltype = 6)
      pv <- rt$pvalue
      out[[length(out) + 1]] <- tibble(seed = SEED + k, q_ordination = q_ord, Region = rg, block = b, n_taxa = ncol(L0), n_states = sum(qk), trace = sum(rl$eig),
                                       p_sites = if (length(pv) >= 2) pv[1] else NA_real_, p_taxa = if (length(pv) >= 2) pv[2] else NA_real_, p_raw = max(pv), note = "")
      rm(list = c("L", "R", "Q", "Qf", "cL", "cR", "cQ"), envir = globalenv())
    }
  }
  bind_rows(out)
}
pg_runs <- lapply(c("fca", "pca"), function(q_ord) {
  r <- if (.Platform$OS.type == "unix") parallel::mclapply(seeds_k, pergroup_one_seed, q_ord = q_ord, mc.cores = min(N_SEEDS_FC, max(1, parallel::detectCores() - 1))) else lapply(seeds_k, pergroup_one_seed, q_ord = q_ord)
  if (!all(vapply(r, is.data.frame, logical(1)))) stop("Step 15: a per-group seed run failed (", q_ord, ")")
  bind_rows(r)
})
pg_all <- bind_rows(pg_runs) %>% group_by(seed, q_ordination, Region) %>% mutate(p_fdr = p.adjust(p_raw, "BH"), n_tests = sum(!is.na(p_raw))) %>% ungroup()
write_csv(pg_all, file.path(OUT_TAB, sprintf("15_fourthcorner_pergroup_by_seed%s.csv", SUF)))
pg_sum <- pg_all %>% group_by(q_ordination, Region, block, n_taxa, n_states, trace) %>%
  summarise(n_tests = first(n_tests), n_seeds_pass = sum(p_fdr <= FDR_MAX, na.rm = TRUE), p_raw_min = min(p_raw), p_raw_max = max(p_raw),
            p_fdr_min = min(p_fdr), p_fdr_max = max(p_fdr), .groups = "drop") %>%
  mutate(q_ordination = ifelse(q_ordination == "fca", "fuzzy correspondence analysis (the rule fixed first)", "states standardized (the declared sensitivity)"),
         verdict = case_when(n_seeds_pass == N_SEEDS_FC ~ "passes under every seed", n_seeds_pass > 0 ~ sprintf("borderline, passes under %d of %d seeds", n_seeds_pass, N_SEEDS_FC), TRUE ~ "does not pass")) %>%
  arrange(q_ordination, Region, block)
write_csv(pg_sum, file.path(OUT_TAB, sprintf("15_fourthcorner_pergroup%s.csv", SUF)))
for (q in unique(pg_sum$q_ordination)) say("  one test per trait group, %s: %d groups pass under every seed (%s)", q, sum(pg_sum$n_seeds_pass[pg_sum$q_ordination == q] == N_SEEDS_FC),
    paste(with(pg_sum[pg_sum$q_ordination == q & pg_sum$n_seeds_pass == N_SEEDS_FC, ], paste(Region, block)), collapse = "; "))

# ---- 3c. Leave-one-taxon-out for every state passing under every seed ------------------------------
# The bivariate statistic is recomputed here exactly as ade4::fourthcorner computes it (asserted against tabD$obs), so
# every carrier taxon of a passing state can be removed in turn without a permutation run; the most influential removal
# is then rerun through fourthcorner under model 6 at N_PERM_FC permutations (the pipeline seed) to see whether the state
# still meets FDR_MAX within its block. A one-carrier state has no statistic without its taxon: it is lost with it.
fc_block <- function(L, x, Qb) {
  w <- L / sum(L); wi <- rowSums(w); wj <- colSums(w)
  zx <- (x - sum(wi * x)) / sqrt(sum(wi * (x - sum(wi * x))^2))
  zq <- apply(Qb, 2, function(q) { v <- sum(wj * (q - sum(wj * q))^2); if (v < 1e-15) rep(0, length(q)) else (q - sum(wj * q)) / sqrt(v) })
  as.vector(crossprod(colSums(w * zx), zq))
}
loo <- purrr::pmap_dfr(per_block %>% filter(n_seeds_pass == N_SEEDS_FC) %>% select(Region, block, term, stat, p_fdr), function(Region, block, term, stat, p_fdr) {
  m <- env %>% filter(Region == !!Region); i <- blk == block; ok <- !apply(Q_norm[, i, drop = FALSE], 1, function(r) any(is.na(r)))
  L0 <- W[m$SiteID, , drop = FALSE]; L0 <- L0[, colSums(L0) > 0 & ok[colnames(L0)], drop = FALSE]
  Qb <- Q_norm[colnames(L0), i, drop = FALSE]; state <- sub("^BIBI / ", "", term)
  full <- fc_block(L0, m$BIBI, Qb); names(full) <- colnames(Qb); stopifnot(abs(full[[state]] - stat) < 1e-10)
  carriers <- colnames(L0)[Qb[, state] > 0]
  if (length(carriers) == 1) return(tibble(Region, block, term, stat, n_carriers = 1L, stat_min_loo = 0, stat_max_loo = 0, taxon_most_influential = carriers,
                                           stat_without_it = 0, delta = -stat, p_fdr_without_it = NA_real_, verdict = "one carrier: the state is lost with its taxon (L36)"))
  loo_stats <- vapply(carriers, function(cx) { keep <- setdiff(colnames(L0), cx); v <- fc_block(L0[, keep, drop = FALSE], m$BIBI, Qb[keep, , drop = FALSE]); v[match(state, colnames(Qb))] }, numeric(1))
  worst <- names(which.max(abs(loo_stats - stat)))
  keep <- setdiff(colnames(L0), worst); set.seed(SEED)
  fcw <- fourthcorner(data.frame(BIBI = m$BIBI), as.data.frame(L0[, keep, drop = FALSE]), as.data.frame(Qb[keep, , drop = FALSE]), modeltype = 6, nrepet = N_PERM_FC,
                      p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
  pw <- fcw$tabD$adj.pvalue[fcw$tabD$names == term]
  tibble(Region, block, term, stat, n_carriers = length(carriers), stat_min_loo = min(loo_stats), stat_max_loo = max(loo_stats), taxon_most_influential = worst,
         stat_without_it = loo_stats[[worst]], delta = loo_stats[[worst]] - stat, p_fdr_without_it = pw,
         verdict = ifelse(pw <= FDR_MAX, "passes without its most influential taxon (pipeline seed)", "fails without its most influential taxon (pipeline seed)"))
})
write_csv(loo, file.path(OUT_TAB, sprintf("15_leave_one_taxon_out%s.csv", SUF)))
say("  leave-one-taxon-out on %d passing states: %d pass without their most influential taxon, %d fail, %d have one carrier",
    nrow(loo), sum(grepl("^passes", loo$verdict)), sum(grepl("^fails", loo$verdict)), sum(loo$n_carriers == 1))

say("")
say("  per-block analysis keeps more taxa. %s coverage by block:", CURRENCY)
cov_b <- per_block %>% distinct(Region, block, n_taxa, pct_cov) %>%
  group_by(block) %>%
  summarise(min_cov = min(pct_cov), .groups = "drop") %>%
  arrange(min_cov)
for (i in seq_len(min(5, nrow(cov_b))))
  say("    %-12s worst-region %s coverage %.1f%%", cov_b$block[i], CURRENCY, cov_b$min_cov[i])

# ---- 4. Do the five named trait-state CWMs agree with the fourth-corner? ---
# WHICH FIVE. Five (modality, region) pairs, chosen for the process each
# regional signature implicates: Piedmont, litter breakdown and particle capture
# (shredder, collector-filterer); Coastal Plain, dispersal range and timing of the
# emergent adult subsidy (strong flying, weak flying, winter emergence). This is a
# named list, not a screening rule: no |rho| / FDR threshold reproduces it without a
# hand-picked trait-block map (tested 21 September 2026), so do not describe it as
# "the strongest CWMs". The list was fixed on 21 September 2026.
C8_ROWS <- tribble(
  ~modality,                     ~Region,
  "Feed_SH",                     "Piedmont",        # Twardochleb names since 22 Sep 2026 (were Feed_SH_Shredder,
  "Feed_CF",                     "Piedmont",        # Feed_CF_Collector-filterer, Flying_Strong.flying and
  "Flying_Strong",               "Coastal Plain",   # Flying_Weak.flying in the legacy binary file)
  "Flying_Weak",                 "Coastal Plain",
  "Emerge_Winter",               "Coastal Plain"
)

# WHICH NUMBER. cwm_rho is computed here, not typed: the within-region Spearman of the
# CWM column in THIS run's weighting (_CWM_den or _CWM_bio) against the same BIBI the
# fourth-corner above used, so the two statistics share a weighting and a B-IBI basis.
# The values this replaces were density coefficients on an earlier B-IBI file,
# compared against both the density and the biomass fourth-corner.
cwm_site <- read_csv(IN$site_metrics, show_col_types = FALSE) %>%
  mutate(SiteID = std_site(SiteID))
cwm_tag  <- if (CURRENCY == "density") "_CWM_den" else "_CWM_bio"
cwm_rho_one <- function(modality, rg) {
  col <- paste0(modality, cwm_tag)
  if (!col %in% names(cwm_site))
    stop("Step 15: no column ", col, " in ", basename(IN$site_metrics))
  m  <- env %>% filter(Region == rg) %>% select(SiteID, BIBI) %>%
    left_join(cwm_site %>% select(SiteID, cwm = all_of(col)), by = "SiteID")
  ok <- is.finite(m$BIBI) & is.finite(m$cwm)
  if (sum(ok) < 4) stop("Step 15: fewer than 4 sites with ", col, " in ", rg)
  cor(m$BIBI[ok], m$cwm[ok], method = "spearman")
}
say("  cwm_rho (%s): %s vs BIBI from %s", CURRENCY, basename(IN$site_metrics), basename(IN$bibi))
cmp <- C8_ROWS %>%
  mutate(cwm_rho = purrr::map2_dbl(modality, Region, cwm_rho_one)) %>%
  left_join(
    fc_all %>% filter(model == "model2", grepl("^BIBI", term)) %>%
      mutate(modality = sub("^BIBI / ", "", term)) %>%
      select(Region, modality, fc_stat = stat, fc_p_model2 = p_raw),
    by = c("Region", "modality")) %>%
  left_join(
    fc_all %>% filter(model == "model4", grepl("^BIBI", term)) %>%
      mutate(modality = sub("^BIBI / ", "", term)) %>%
      select(Region, modality, fc_p_model4 = p_raw),
    by = c("Region", "modality")) %>%
  mutate(agreement = case_when(
    is.na(fc_stat)                            ~ "not testable",
    sign(fc_stat) != sign(cwm_rho)            ~ "CONTRADICTS",
    fc_p_model2 <= 0.05 & fc_p_model4 <= 0.05 ~ "corroborated, both models",
    fc_p_model2 <= 0.05                       ~ "corroborated under site permutation only",
    fc_p_model2 <= 0.10                       ~ "directional, marginal",
    TRUE                                      ~ "not corroborated"))
write_csv(cmp, file.path(OUT_TAB, sprintf("15_cwm_vs_fourthcorner%s.csv", SUF)))

say("")
say("  the five named trait-state CWMs, against the fourth-corner:")
for (i in seq_len(nrow(cmp)))
  say("    %-28s %-14s CWM %+.2f | 4C %s p2=%s p4=%s -> %s",
      cmp$modality[i], cmp$Region[i], cmp$cwm_rho[i],
      ifelse(is.na(cmp$fc_stat[i]), "  --  ", sprintf("%+.3f", cmp$fc_stat[i])),
      ifelse(is.na(cmp$fc_p_model2[i]), " -- ", sprintf("%.3f", cmp$fc_p_model2[i])),
      ifelse(is.na(cmp$fc_p_model4[i]), " -- ", sprintf("%.3f", cmp$fc_p_model4[i])),
      cmp$agreement[i])
# ---- 5. Coverage confound, applied to every per-block FDR survivor ----------
# Step 12 retired the REL composites because they correlated with trait
# coverage at r = 0.93-0.99 and collapsed under partial correlation. The same
# question has to be asked of any fourth-corner result whose trait group is not
# recorded for every taxon present (the coverage lines printed above give each
# group's worst-region coverage).

# FDR_MAX (0.10) is the paper's declared admission standard, set in 00_config.R
# and used by the metric screen. This step hardcoded 0.05 until 2 September 2026.
surv <- per_block %>% filter(n_seeds_pass == N_SEEDS_FC)       # passes under every seed
conf <- purrr::pmap_dfr(list(surv$Region, surv$block, surv$term), function(rg, b, tm) {
  m  <- env %>% filter(Region == rg)
  i  <- blk == b
  ok <- !apply(Q_norm[, i, drop = FALSE], 1, function(r) any(is.na(r)))
  X  <- W[m$SiteID, , drop = FALSE]
  keep_tx <- colnames(X)[colSums(X) > 0 & ok[colnames(X)]]
  # per-site coverage = share of this run's weighting (density or biomass) held by taxa resolved for THIS block
  cov <- rowSums(X[, keep_tx, drop = FALSE]) / rowSums(X)
  # Since the trait rebuild (22 Sep 2026) most blocks have a record for every taxon present, so coverage is 1 at
  # every site: nothing is missing, no coverage confound is possible, and the correlations below are undefined.
  complete <- all(abs(cov - 1) < 1e-12)
  modality <- sub("^BIBI / ", "", tm)
  P  <- X[, keep_tx, drop = FALSE]; P <- P / rowSums(P)
  cwm <- as.numeric(P %*% Q_norm[keep_tx, modality])
  r_cov <- if (complete) NA_real_ else cor(cwm, cov, method = "spearman")
  rho   <- cor(m$BIBI, cwm, method = "spearman")
  # partial rho of cwm ~ BIBI given coverage, on ranks
  rb <- if (complete) NA_real_ else cor(m$BIBI, cov, method = "spearman")
  part <- (rho - r_cov * rb) / sqrt((1 - r_cov^2) * (1 - rb^2))
  tibble(Region = rg, block = b, modality = modality,
         cov_pct = 100 * mean(cov), r_cwm_coverage = r_cov,
         rho = rho, rho_partial = part,
         pct_lost = 100 * (1 - abs(part) / abs(rho)),
         verdict = case_when(complete                   ~ "complete coverage -- no confound possible",
                             abs(part) >= 0.7 * abs(rho) ~ "survives -- real trait signal",
                             TRUE                        ~ "COLLAPSES -- coverage artifact"))
})
cf_out <- conf
names(cf_out)[names(cf_out) == "cov_pct"] <- paste0(CURRENCY, "_coverage")
write_csv(cf_out, file.path(OUT_TAB, sprintf("15_fourthcorner_coverage_confound%s.csv", SUF)))

say("")
say("  coverage confound test on every per-block FDR survivor:")
for (i in seq_len(nrow(conf)))
  say("    %-14s %-26s cov %.0f%% | r(cwm,cov)=%+.2f | rho %+.3f -> %+.3f partial (%.0f%% lost)  %s",
      conf$Region[i], conf$modality[i], conf$cov_pct[i], conf$r_cwm_coverage[i],
      conf$rho[i], conf$rho_partial[i], conf$pct_lost[i], conf$verdict[i])

# The verdict is COMPUTED, not hardcoded. An earlier version stated the density
# conclusions as fixed prose ("flight CORROBORATED in both models"), which
# contradicted the CWM table printed directly above it -- that table reports
# flight as directional/marginal under the complete-case models, while flight
# survives FDR only in the per-block analysis. Hardcoded verdicts drift.
say("")
say("  VERDICT -- currency = %s", CURRENCY)
say("    complete-case: %d of %d state tests survive FDR",
    sum(fc_all$p_fdr <= FDR_MAX), nrow(fc_all))
say("    per-block    : %d modalities pass FDR under every seed, %d coverage-robust (complete coverage or surviving the partial test); %d borderline",
    nrow(surv), sum(grepl("^(survives|complete coverage)", conf$verdict)), sum(per_block$n_seeds_pass > 0 & per_block$n_seeds_pass < N_SEEDS_FC))
if (nrow(surv)) {
  say("    surviving modalities:")
  for (i in seq_len(nrow(surv)))
    say("      %-14s %-28s stat %+.3f  FDR %.4f", surv$Region[i],
        sub("^BIBI / ", "", surv$term[i]), surv$stat[i], surv$p_fdr[i])
}
say("    CWM agreement (complete-case models 2 and 4):")
for (i in seq_len(nrow(cmp)))
  say("      %-28s %-14s %s", cmp$modality[i], cmp$Region[i], cmp$agreement[i])
say("    Report trait associations as standing-stock composition, not trait filtering.")
say("STEP 15 complete")
