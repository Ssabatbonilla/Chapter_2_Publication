# =============================================================================
# 16_beta_turnover_nestedness.R
# =============================================================================
# NOT REPORTED IN THIS PAPER
#   Beta-diversity partitioning moved to its own paper, and none of its results
#   are reported here. The step still runs because step 18 reads
#   16_cross_region_pool_beta.csv for the number of taxa recorded in each region
#   and shared between regions (the repository table regional_pools_shared_taxa).
#
# WHAT IT DOES
#   Partitions the dissimilarity between sites into turnover and
#   nestedness-resultant components (the Baselga framework), on taxa and in trait
#   space (Villeger et al. 2013), within each region, and asks whether each
#   component tracks the condition gradient (Mantel test against the difference
#   in B-IBI). Functional beta partitioning is PRESENCE-BASED, unlike the
#   abundance-weighted metrics of the paper.
#
# INPUTS
#   den_mat and dat (step 01); IN$traits_fuz, IN$traits_bin and IN$trait_blocks
#   (the step 13 trait space).
#
# OUTPUT
#   output/tables/16_beta_partition_by_region.csv
#   output/tables/16_beta_gradient_mantel.csv
#   output/tables/16_beta_multisite_by_condition.csv  (mean PAIRWISE beta within each condition third)
#   output/tables/16_cross_region_pool_beta.csv       (read by step 18)
# =============================================================================

say("STEP 16 -- beta diversity partitioned into turnover and nestedness")
suppressPackageStartupMessages({ library(betapart); library(FD); library(ade4) })
set.seed(SEED)

N_AXES <- 3   # functional hulls need n_taxa > n_axes at every site; 3 is the
              # usual compromise and keeps every site testable here.

# PERFORMANCE NOTE -- do not reintroduce functional.beta.multi().
#   It computes an N-WAY convex hull intersection, which scales roughly 2.7x
#   per added site: 4 sites 7.2 s, 5 sites 19.4 s, 6 sites over 60 s. On the
#   10-site condition groups used below it does not finish in hours.
#   functional.beta.pair() on all 30 sites takes 1.1 s, so the condition-group
#   summary is computed by averaging the PAIRWISE components within each group
#   instead. Same quantity, and it terminates.
#   The taxonomic components are averaged the same way, from beta.pair().

# ---- 1. Trait space, built with the step 13 recipe --------------------------
.std <- function(df) {
  if (!"TAXANAME" %in% names(df)) df$TAXANAME <- df$Genus
  if (!"Genus"    %in% names(df)) df$Genus    <- df$TAXANAME
  df %>% mutate(Genus = stringr::str_squish(as.character(Genus)),
                TAXANAME = stringr::str_squish(as.character(TAXANAME))) %>%
    filter(!is.na(TAXANAME), TAXANAME != "") %>% distinct(TAXANAME, .keep_all = TRUE)
}
# The binary file is read only to intersect column and taxon names with the fuzzy file (the either-coding
# retention rule of step 13); no binary value enters the distance.
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
gd <- suppressWarnings(gowdis(as.data.frame(MF), w = gw))
X  <- as.matrix(dudi.pco(suppressWarnings(cailliez(gd)), scannf = FALSE, nf = N_AXES)$li)
say("  trait space: %d taxa x %d axes (cailliez, weighted Gower -- step 13 recipe)",
    nrow(X), ncol(X))

# ---- 2. Per-region partitioning ---------------------------------------------
one_region <- function(rg) {
  m  <- dat %>% filter(Region == rg) %>% arrange(SiteID)
  PA <- (den_mat[m$SiteID, , drop = FALSE] > 0) * 1
  PA <- PA[, colSums(PA) > 0, drop = FALSE]
  Xr <- X[colnames(PA), , drop = FALSE]

  tb_pair <- beta.pair(PA, index.family = "sorensen")
  fb_pair <- functional.beta.pair(PA, Xr, index.family = "sorensen")

  bd <- as.numeric(dist(m$BIBI))   # not used below: mt() passes dist(m$BIBI) directly
  mt <- function(d) {
    r <- suppressWarnings(vegan::mantel(as.dist(d), dist(m$BIBI),
                                        method = "spearman", permutations = N_PERM))
    c(r$statistic, r$signif)
  }
  parts <- list(
    tax_turnover   = tb_pair$beta.sim,  tax_nestedness = tb_pair$beta.sne,
    tax_total      = tb_pair$beta.sor,
    fun_turnover   = fb_pair$funct.beta.sim, fun_nestedness = fb_pair$funct.beta.sne,
    fun_total      = fb_pair$funct.beta.sor)

  summ <- purrr::imap_dfr(parts, function(d, nm) {
    mm <- mt(d)
    tibble(Region = rg, component = nm, mean_beta = mean(as.numeric(d)),
           mantel_r = mm[1], mantel_p = mm[2])
  })
  nest_share <- tibble(
    Region = rg,
    tax_pct_nestedness = 100 * mean(as.numeric(tb_pair$beta.sne)) /
                               mean(as.numeric(tb_pair$beta.sor)),
    fun_pct_nestedness = 100 * mean(as.numeric(fb_pair$funct.beta.sne)) /
                               mean(as.numeric(fb_pair$funct.beta.sor)))

  # Condition-group summary. Mean of the PAIRWISE components among sites within
  # each group -- see the performance note above for why not functional.beta.multi().
  q  <- quantile(m$BIBI, c(1/3, 2/3))
  sub_mean <- function(d, idx) mean(as.matrix(d)[idx, idx][lower.tri(diag(length(idx)))])
  multi <- purrr::map_dfr(c("bottom third", "top third"), function(g) {
    idx <- if (g == "bottom third") which(m$BIBI <= q[1]) else which(m$BIBI >= q[2])
    P   <- PA[idx, , drop = FALSE]; P <- P[, colSums(P) > 0, drop = FALSE]
    tibble(Region = rg, condition_group = g, n_sites = length(idx), n_taxa = ncol(P),
           tax_turnover   = sub_mean(tb_pair$beta.sim, idx),
           tax_nestedness = sub_mean(tb_pair$beta.sne, idx),
           tax_total      = sub_mean(tb_pair$beta.sor, idx),
           fun_turnover   = sub_mean(fb_pair$funct.beta.sim, idx),
           fun_nestedness = sub_mean(fb_pair$funct.beta.sne, idx),
           fun_total      = sub_mean(fb_pair$funct.beta.sor, idx))
  })
  list(summ = summ, nest = nest_share, multi = multi)
}

out <- lapply(REGIONS, one_region); names(out) <- REGIONS
mantel_tab <- bind_rows(lapply(out, `[[`, "summ"))
nest_tab   <- bind_rows(lapply(out, `[[`, "nest"))
multi_tab  <- bind_rows(lapply(out, `[[`, "multi"))

write_csv(mantel_tab, file.path(OUT_TAB, "16_beta_gradient_mantel.csv"))
write_csv(nest_tab,   file.path(OUT_TAB, "16_beta_partition_by_region.csv"))
write_csv(multi_tab,  file.path(OUT_TAB, "16_beta_multisite_by_condition.csv"))

# ---- 4. Cross-region pool comparison ----------------------------------------
# WHY THIS IS HERE
#   The number of taxa recorded in each region, and shared between regions,
#   comes from this table (step 18 reads it). The log below flags any pair of
#   regions whose pools differ more by nestedness than by turnover.
#
#   This is a POOL-level comparison -- a region is present or absent for a
#   taxon, not a site -- so it answers a different question from the
#   within-region partitioning above and does not replace it.

pool <- t(sapply(REGIONS, function(rg) {
  sid <- dat %>% filter(Region == rg) %>% pull(SiteID)
  (colSums(den_mat[sid, , drop = FALSE] > 0) > 0) * 1
}))
rownames(pool) <- REGIONS
pool <- pool[, colSums(pool) > 0, drop = FALSE]

pb  <- beta.pair(pool, index.family = "sorensen")
prs <- combn(REGIONS, 2)

pool_tab <- tibble(
  Region_A = prs[1, ], Region_B = prs[2, ],
  taxa_A   = as.integer(rowSums(pool)[prs[1, ]]),
  taxa_B   = as.integer(rowSums(pool)[prs[2, ]]),
  shared   = as.integer(mapply(function(a, b) sum(pool[a, ] > 0 & pool[b, ] > 0),
                               prs[1, ], prs[2, ])),
  turnover   = as.numeric(pb$beta.sim),
  nestedness = as.numeric(pb$beta.sne),
  total      = as.numeric(pb$beta.sor)) %>%
  mutate(pct_nestedness = 100 * nestedness / total,
         turnover_exceeds_nestedness = turnover > nestedness)

write_csv(pool_tab, file.path(OUT_TAB, "16_cross_region_pool_beta.csv"))

n_all3 <- sum(colSums(pool) == length(REGIONS))
say("")
say("  cross-region POOL comparison (regional presence-absence, %d taxa):", ncol(pool))
say("    taxa recorded in all three regions: %d of %d", n_all3, ncol(pool))
for (i in seq_len(nrow(pool_tab)))
  say("    %-15s vs %-15s turnover %.3f | nestedness %.3f | %s",
      pool_tab$Region_A[i], pool_tab$Region_B[i],
      pool_tab$turnover[i], pool_tab$nestedness[i],
      ifelse(pool_tab$turnover_exceeds_nestedness[i],
             "turnover exceeds nestedness",
             "*** NESTEDNESS EXCEEDS TURNOVER ***"))

say("")
say("  share of total dissimilarity that is NESTEDNESS (role loss rather than")
say("  role replacement), averaged over all site pairs:")
for (i in seq_len(nrow(nest_tab)))
  say("    %-15s taxonomic %5.1f%% | functional %5.1f%%", nest_tab$Region[i],
      nest_tab$tax_pct_nestedness[i], nest_tab$fun_pct_nestedness[i])

say("")
say("  does each component track the condition gradient? (Mantel vs |delta B-IBI|)")
for (rg in REGIONS) {
  r <- mantel_tab %>% filter(Region == rg, grepl("^fun", component))
  say("    %-15s %s", rg,
      paste(sprintf("%s r=%+.3f p=%.3f", sub("fun_", "", r$component),
                    r$mantel_r, r$mantel_p), collapse = " | "))
}

say("")
say("  mean pairwise beta, bottom vs top third of the within-region gradient:")
for (i in seq_len(nrow(multi_tab)))
  say("    %-15s %-13s functional: turnover %.3f | nestedness %.3f | total %.3f",
      multi_tab$Region[i], multi_tab$condition_group[i],
      multi_tab$fun_turnover[i], multi_tab$fun_nestedness[i], multi_tab$fun_total[i])
say("STEP 16 complete")
