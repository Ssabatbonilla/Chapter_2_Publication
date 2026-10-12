# =============================================================================
# 15b_weighting_difference.R -- does a trait state's fourth-corner statistic differ between the weightings? (24 Sep 2026)
# =============================================================================
# NOT REPORTED IN THIS PAPER
#   This test belongs to the fourth-corner analysis, which left this paper for a companion trait-state manuscript
#   (legacy Table S20 and Figure S2; step 18b copies them to output/manuscript/trait_state_manuscript/). The step
#   still runs because step 06 draws legacy Figure S2 from 15b_weighting_difference.csv, and step 18 reads that file
#   and 15b_weighting_difference_settings.csv for legacy Table S20 and the Figure S2 caption.
#
# INPUTS
#   den_mat, bio_mat and dat (step 01); IN$traits; step 15's 15_fourthcorner_per_block.csv and
#   15_fourthcorner_per_block_biomass.csv, against which the observed statistics are asserted.
#
# WHY THIS STEP EXISTS
#   Hypothesis 3 says density- and biomass-weighted responses differ. For a trait state, "passes weighted by density
#   and not weighted by biomass" compares two significance verdicts, and a difference in significance is not itself a
#   significant difference (Gelman & Stern, 2006). This step tests the difference directly, at N_PERM_FC resamples
#   under N_SEEDS_FC seeds.
#
# THE TEST
#   For each region, trait group and state, the fourth-corner statistic (the inflated-table correlation between the
#   index and the state's affinity, exactly step 15's) is computed weighted by density and weighted by biomass, and the
#   difference (biomass minus density) is bootstrapped by resampling streams and taxa INDEPENDENTLY, WITH REPLACEMENT,
#   ON EACH DRAW, both statistics computed on the same resample (the pigeonhole bootstrap; Owen, 2007). Resampling
#   streams alone understates the uncertainty of a trait-environment statistic because taxa are random units too
#   (ter Braak, Peres-Neto & Dray, 2017); resampling both together in one permutation is a different procedure that
#   does not control type I error, which is why the wording is "independently, on each draw".
#   Reported per state: the observed difference, the bootstrap 95% percentile interval (range over the seeds), the
#   share of draws in which the difference had the opposite sign to the observed one, a two-sided percentile p
#   (p = 2 x (1 + k)/(B + 1), k the smaller count on either side of zero with exact zeros counted half; 26 Sep 2026), the
#   Benjamini-Hochberg adjustment of that p across the region's contrasts per seed, and the pass rule (FDR <= FDR_MAX,
#   0.10, under every seed).
#   In a two-state group the two states' differences are equal and opposite, so the family counts ONE state per such
#   group (its first state); the mirror state is listed with its statistics and no adjusted p.
#   A one-carrier state: the share of draws in which no carrier taxon was drawn is reported beside its p,
#   because the approximation is weakest where one taxon dominates a state (Owen & Eckles, 2012).
#
# OUTPUT
#   output/tables/15b_weighting_difference.csv           one row per region x state
#   output/tables/15b_weighting_difference_by_seed.csv   the per-seed interval and p
#   output/tables/15b_weighting_difference_settings.csv
# =============================================================================
say("STEP 15b -- density-against-biomass difference of each trait state's fourth-corner statistic (pigeonhole bootstrap)")
tr_raw <- as.data.frame(readRDS(IN$traits)); rownames(tr_raw) <- tr_raw$TAXANAME
Q_raw  <- as.matrix(tr_raw[colnames(den_mat), setdiff(names(tr_raw), c("Genus", "TAXANAME")), drop = FALSE]); storage.mode(Q_raw) <- "double"
blk    <- sub("_.*$", "", colnames(Q_raw)); Q_norm <- Q_raw
for (b in unique(blk)) { i <- blk == b; Q_norm[, i] <- Q_raw[, i, drop = FALSE] / rowSums(Q_raw[, i, drop = FALSE]) }
stopifnot(identical(dimnames(den_mat), dimnames(bio_mat)))

# the fourth-corner statistic for one quantitative variable and a block's states, as step 15 (ade4::fourthcorner) computes it:
# w_ij = L_ij / L++, x and each state standardized with the marginal weights; a state no drawn taxon has is 0, as ade4 reports it
fc_block <- function(L, x, Qb) {
  w <- L / sum(L); wi <- rowSums(w); wj <- colSums(w)
  zx <- (x - sum(wi * x)) / sqrt(sum(wi * (x - sum(wi * x))^2))
  zq <- apply(Qb, 2, function(q) { v <- sum(wj * (q - sum(wj * q))^2); if (v < 1e-15) rep(0, length(q)) else (q - sum(wj * q)) / sqrt(v) })
  as.vector(crossprod(colSums(w * zx), zq))
}
block_taxa <- function(s, b) {
  i <- blk == b; ok <- !apply(Q_norm[, i, drop = FALSE], 1, function(r) any(is.na(r)))
  L <- den_mat[s, , drop = FALSE]; keep <- colSums(L) > 0 & ok[colnames(L)]
  stopifnot(identical(unname(keep), unname(colSums(bio_mat[s, , drop = FALSE]) > 0 & ok[colnames(L)])))   # the same taxa are present by density and by biomass
  list(cols = colnames(L)[keep], Qb = Q_norm[colnames(L)[keep], i, drop = FALSE])
}
# observed statistics, asserted against step 15's per-block tables (the statistic never depends on the seed)
pb15 <- bind_rows(read_csv(file.path(OUT_TAB, "15_fourthcorner_per_block.csv"), show_col_types = FALSE) %>% mutate(weighting = "density"),
                  read_csv(file.path(OUT_TAB, "15_fourthcorner_per_block_biomass.csv"), show_col_types = FALSE) %>% mutate(weighting = "biomass"))
obs <- list()
for (rg in REGIONS) {
  m <- dat %>% filter(Region == rg); s <- m$SiteID; x <- m$BIBI
  for (b in unique(blk)) {
    st <- block_taxa(s, b); if (length(st$cols) < 5) next
    vd <- fc_block(den_mat[s, st$cols, drop = FALSE], x, st$Qb); vb <- fc_block(bio_mat[s, st$cols, drop = FALSE], x, st$Qb)
    obs[[length(obs) + 1]] <- tibble(Region = rg, block = b, term = paste0("BIBI / ", colnames(st$Qb)), n_taxa = length(st$cols),
                                     n_carriers = colSums(st$Qb > 0), stat_den = vd, stat_bio = vb)
  }
}
obs <- bind_rows(obs) %>% mutate(diff = stat_bio - stat_den)
chk <- obs %>% inner_join(pb15 %>% select(weighting, Region, term, stat) %>% pivot_wider(names_from = weighting, values_from = stat), by = c("Region", "term"))
stopifnot(nrow(chk) == nrow(obs), max(abs(chk$stat_den - chk$density)) < 1e-10, max(abs(chk$stat_bio - chk$biomass)) < 1e-10)
say("  observed statistics reproduce step 15's per-block tables for %d region-by-state rows", nrow(obs))
# one contrast per two-state group: the first state of the block
two_state <- obs %>% count(Region, block) %>% filter(n == 2) %>% pull(block) %>% unique()
first_of <- obs %>% group_by(Region, block) %>% slice(1) %>% ungroup() %>% mutate(keep = TRUE) %>% select(Region, block, term, keep)
obs <- obs %>% left_join(first_of, by = c("Region", "block", "term")) %>% mutate(keep = ifelse(block %in% two_state, coalesce(keep, FALSE), TRUE))

one_seed <- function(k) {
  set.seed(SEED + k); B <- N_PERM_FC; rows <- list()
  for (rg in REGIONS) {
    m <- dat %>% filter(Region == rg); s <- m$SiteID; x0 <- m$BIBI
    for (b in unique(blk)) {
      st <- block_taxa(s, b); tx <- st$cols; if (length(tx) < 5) next
      Qb <- st$Qb
      D <- matrix(NA_real_, B, ncol(Qb)); nocar <- matrix(FALSE, B, ncol(Qb))
      for (r in seq_len(B)) {
        ix <- sample.int(length(s), replace = TRUE); j <- sample(tx, replace = TRUE)
        Ld <- den_mat[s[ix], j, drop = FALSE]; Lb <- bio_mat[s[ix], j, drop = FALSE]
        ok <- colSums(Ld) > 0 & colSums(Lb) > 0
        if (sum(ok) < 2) next
        Qj <- Qb[j[ok], , drop = FALSE]
        D[r, ] <- fc_block(Lb[, ok, drop = FALSE], x0[ix], Qj) - fc_block(Ld[, ok, drop = FALSE], x0[ix], Qj)
        nocar[r, ] <- colSums(Qj > 0) == 0
      }
      o <- obs$diff[match(paste(rg, paste0("BIBI / ", colnames(Qb))), paste(obs$Region, obs$term))]
      rows[[length(rows) + 1]] <- tibble(seed = SEED + k, Region = rg, block = b, term = paste0("BIBI / ", colnames(Qb)),
        boot_mean = colMeans(D, na.rm = TRUE), lo95 = apply(D, 2, quantile, 0.025, na.rm = TRUE), hi95 = apply(D, 2, quantile, 0.975, na.rm = TRUE),
        share_opposite_sign = vapply(seq_len(ncol(D)), function(cc) { d <- D[, cc]; d <- d[is.finite(d)]; if (o[cc] >= 0) mean(d <= 0) else mean(d >= 0) }, numeric(1)),
        # 26 Sep 2026: the two-sided percentile-bootstrap p is twice the smaller
        # count of resampled differences on either side of zero, a difference of exactly zero counted half on each side (Rousselet
        # et al., 2021, p. 4), with one added to the count and to the number of resamples so no p is zero (Phipson & Smyth, 2010).
        # Was min(1, 2 * min(mean(d <= 0), mean(d >= 0)) + 1/(B + 1)): the +1 after doubling has no source, and a one-carrier
        # state's exact zeros (draws without its carrier) were counted in both tails.
        p_boot = vapply(seq_len(ncol(D)), function(cc) { d <- D[, cc]; d <- d[is.finite(d)]; t0 <- sum(d == 0)
                        k <- min(sum(d < 0) + t0 / 2, sum(d > 0) + t0 / 2); min(1, 2 * (1 + k) / (length(d) + 1)) }, numeric(1)),
        # 27 Sep 2026: the count the p is built from, so it can be recomputed from this table
        k_smaller_side = vapply(seq_len(ncol(D)), function(cc) { d <- D[, cc]; d <- d[is.finite(d)]; t0 <- sum(d == 0); min(sum(d < 0) + t0 / 2, sum(d > 0) + t0 / 2) }, numeric(1)),
        n_zero = vapply(seq_len(ncol(D)), function(cc) { d <- D[, cc]; sum(d[is.finite(d)] == 0) }, numeric(1)),
        n_draws = colSums(is.finite(D)), share_draws_no_carrier = colMeans(nocar))
    }
  }
  bind_rows(rows)
}
t0 <- Sys.time()
runs <- if (.Platform$OS.type == "unix") parallel::mclapply(seq_len(N_SEEDS_FC) - 1, one_seed, mc.cores = min(N_SEEDS_FC, max(1, parallel::detectCores() - 1))) else lapply(seq_len(N_SEEDS_FC) - 1, one_seed)
if (!all(vapply(runs, is.data.frame, logical(1)))) stop("Step 15b: a seed run failed: ", paste(runs[!vapply(runs, is.data.frame, logical(1))], collapse = " | "))
by_seed <- bind_rows(runs) %>% left_join(obs %>% select(Region, term, keep), by = c("Region", "term")) %>%
  group_by(seed, Region) %>% mutate(p_fdr = ifelse(keep, p.adjust(ifelse(keep, p_boot, NA), "BH"), NA_real_)) %>% ungroup()
say("  %d resamples x %d seeds in %.1f min", N_PERM_FC, N_SEEDS_FC, as.numeric(difftime(Sys.time(), t0, units = "mins")))
summ <- by_seed %>% group_by(Region, term) %>%
  summarise(boot_mean = mean(boot_mean), lo95_min = min(lo95), lo95_max = max(lo95), hi95_min = min(hi95), hi95_max = max(hi95),
            share_opposite_sign_min = min(share_opposite_sign), share_opposite_sign_max = max(share_opposite_sign),
            p_boot_min = min(p_boot), p_boot_max = max(p_boot), p_fdr_min = suppressWarnings(min(p_fdr)), p_fdr_max = suppressWarnings(max(p_fdr)),
            n_seeds_pass = sum(p_fdr <= FDR_MAX, na.rm = TRUE), share_draws_no_carrier = mean(share_draws_no_carrier), n_draws_min = min(n_draws), .groups = "drop") %>%
  mutate(across(c(p_fdr_min, p_fdr_max), ~ ifelse(is.finite(.x), .x, NA_real_)))
out <- obs %>% left_join(summ, by = c("Region", "term")) %>%
  mutate(interval_excludes_zero_all_seeds = (lo95_max > 0) | (hi95_min < 0),
         verdict = case_when(!keep ~ "mirror of the first state of its group; not counted",
                             n_seeds_pass == N_SEEDS_FC ~ "differs under every seed",
                             n_seeds_pass > 0 ~ sprintf("borderline, passes under %d of %d seeds", n_seeds_pass, N_SEEDS_FC),
                             TRUE ~ "not shown to differ")) %>%
  arrange(Region, block, term)
write_csv(out, file.path(OUT_TAB, "15b_weighting_difference.csv"))
write_csv(by_seed, file.path(OUT_TAB, "15b_weighting_difference_by_seed.csv"))
write_csv(tibble(setting = c("resamples per seed", "seeds", "resampling", "statistic", "family", "pass rule", "p-value rule", "R version"),
                 value = c(N_PERM_FC, sprintf("%d to %d", SEED, SEED + N_SEEDS_FC - 1),
                           "streams and taxa resampled independently, with replacement, on each draw; both statistics on the same resample (Owen, 2007; ter Braak et al., 2017)",
                           "fourth-corner statistic weighted by biomass minus weighted by density, per state",
                           "Benjamini-Hochberg across a region's contrasts, one state per two-state group",
                           sprintf("FDR <= %.2f under all %d seeds", FDR_MAX, N_SEEDS_FC),
                           "two-sided percentile bootstrap p = 2 (1 + k) / (B + 1), k the smaller count of resampled differences on either side of zero with exact zeros counted half (Rousselet et al. 2021; Phipson and Smyth 2010)",
                           R.version.string)),
          file.path(OUT_TAB, "15b_weighting_difference_settings.csv"))
i_max <- which.max(abs(out$diff) * out$keep)
say("  contrasts counted: %d; differ under every seed: %d; borderline: %d; largest |difference| %s %s (%+.3f, interval %.3f to %.3f)",
    sum(out$keep), sum(out$keep & out$verdict == "differs under every seed"), sum(out$keep & grepl("borderline", out$verdict)),
    out$Region[i_max], out$term[i_max], out$diff[i_max], out$lo95_min[i_max], out$hi95_max[i_max])
say("STEP 15b complete")
