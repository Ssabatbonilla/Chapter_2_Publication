# =============================================================================
# 20_chironomid_sensitivity.R   (added 8 October 2026)
# =============================================================================
# WHY THIS STEP EXISTS
#
#   Four chironomid groups (Chironomini, Orthocladiinae, Tanypodinae, Tanytarsini), identified to subfamily or tribe, hold
#   a third to two thirds of the trait-set density, and each carries one trait profile (respiration is coded identically for all four; dispersal, flying
#   strength and synchrony are coded at group level, so every genus in a group carries its group's code). A density-weighted trait result can therefore be a
#   statement about how much of the community is chironomid. This step asks, from the run of record's own inputs:
#
#     1. how large the chironomid share is, and how it tracks condition (Table S7 part A; legacy Table S26);
#     2. whether the density- and biomass-weighted functional metrics keep their association with condition when the four
#        groups are removed, against removing the four most abundant other taxa and 1,000 random sets of four (Table S7 part B; legacy Table S27);
#     3. whether each trait state passing under all five seeds keeps passing (same call, seeds and permutations as step 15;
#        legacy Table S28, not in this paper: it left with the fourth-corner analysis for the trait-state manuscript);
#     4. how far a genus-level split of the four groups, under fixed, varying and condition-tracking mixes and with database
#        profiles or the profile of a random other taxon, moves Rao's Q and overlap redundancy (Table S7 part C; legacy Table S29).
#
#   Nothing here changes a number in the main pipeline outputs; it only reads them. Every recomputed quantity is first asserted against
#   its published value: the fourth-corner statistic and the pass count of every state, the p value of every state under the
#   pipeline seeds, and the site-level functional metrics.
#
# INPUTS
#   den_mat, bio_mat and dat (step 01); IN$traits, IN$trait_blocks, IN$site_metrics; 03_screen_all_metrics.csv (step 03);
#   15_fourthcorner_per_block.csv and 15_fourthcorner_per_block_biomass.csv (step 15), against which section 3 is asserted.
#
# OUTPUT
#   output/tables/20_*.csv; output/figures/figS5_chironomid_removal.{png,pdf,svg} (Figure S1, legacy Figure S3), which
#   step 20b overwrites with the final two-panel version.
# =============================================================================

say("STEP 20 -- sensitivity of the density-weighted trait results to the Chironomidae")
suppressPackageStartupMessages({ library(FD); library(ade4); library(parallel); library(patchwork) })

CHIRONOMIDAE <- c("Chironomini", "Orthocladiinae", "Tanypodinae", "Tanytarsini")
stopifnot("a chironomid group is missing from the trait-set matrix" = all(CHIRONOMIDAE %in% colnames(den_mat)))
N_RANDOM_SETS <- 1000L; N_GENUS_DRAWS <- 1000L; SIM_SEED <- SEED + 100L
# Genera that carry a profile in the Twardochleb fuzzy file, assigned to the four groups by subfamily or tribe. The Podonominae
# and Diamesinae genera of the file (Boreochlus, Parochlus, Paraboreochlus, Diamesa, Lasiodiamesa) belong to none of the four groups.
GENERA_OF_GROUP <- list(
  Tanypodinae = c("Labrundinia", "Procladius", "Ablabesmyia", "Paramerina", "Larsia", "Nilotanypus", "Pentaneura"),
  Tanytarsini = c("Micropsectra", "Tanytarsus", "Rheotanytarsus"),
  Chironomini = c("Phaenopsectra", "Dicrotendipes", "Chironomus", "Polypedilum", "Stenochironomus", "Demicryptochironomus", "Apedilum",
                  "Paratendipes", "Pseudochironomus", "Cryptochironomus", "Stictochironomus", "Microtendipes", "Paracladopelma", "Harnischia"),
  Orthocladiinae = c("Orthocladius", "Eukiefferiella", "Rheocricotopus", "Brillia", "Psectrocladius", "Smittia", "Heterotrissocladius",
                     "Corynoneura", "Parametriocnemus", "Krenosmittia", "Thienemanniella", "Nanocladius", "Tvetenia", "Chaetocladius",
                     "Hydrobaenus", "Heleniella"))

sites <- dat %>% select(SiteID, Region, BIBI)
stopifnot(setequal(sites$SiteID, rownames(den_mat)))
sites <- sites[match(rownames(den_mat), sites$SiteID), ]
sp_cor <- function(a, b) { r <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)); c(rho = unname(r$estimate), p = r$p.value) }
CURRENCIES <- list(Density = den_mat, Biomass = bio_mat)

# ---- 1. The chironomid share and condition (Table S7 part A; legacy Table S26) -------------------------------------------------------------------------
share_tab <- purrr::map_dfr(names(CURRENCIES), function(cur) { W <- CURRENCIES[[cur]]
  purrr::map_dfr(REGIONS, function(rg) { i <- sites$Region == rg; tot <- rowSums(W[i, ]); ch <- rowSums(W[i, CHIRONOMIDAE]); r <- sp_cor(ch / tot, sites$BIBI[i])
    tibble(Weighting = cur, Region = rg, `Chironomidae, share of the trait-set total (%)` = 100 * sum(ch) / sum(tot),
           `Taxa to reach half the total` = { w <- sort(colSums(W[i, ]), decreasing = TRUE); which(cumsum(w) / sum(w) >= .5)[1] },
           `Of those, chironomid taxa` = { w <- sort(colSums(W[i, ]), decreasing = TRUE); k <- which(cumsum(w) / sum(w) >= .5)[1]; sum(names(w)[1:k] %in% CHIRONOMIDAE) },
           rho = r[["rho"]], p = r[["p"]]) }) }) %>%
  group_by(Weighting) %>% mutate(`FDR p` = p.adjust(p, "BH")) %>% ungroup() %>% rename(`rho with B-IBI` = rho, `raw p` = p) %>%
  arrange(factor(Weighting, c("Density", "Biomass")), factor(Region, REGIONS))
write_csv(share_tab, file.path(OUT_TAB, "20_chironomid_share.csv"))
say("  chironomid share of trait-set density: %s", paste(sprintf("%s %.1f%%", share_tab$Region[share_tab$Weighting == "Density"], share_tab[[3]][share_tab$Weighting == "Density"]), collapse = "; "))

# ---- 2. Functional metrics with the four groups kept and removed (Table S7 part B; legacy Table S27) ---------------------------------------------------
tf <- as.data.frame(readRDS(IN$traits)); if (!"TAXANAME" %in% names(tf)) tf$TAXANAME <- tf$Genus
tf <- tf[!duplicated(tf$TAXANAME), ]; rownames(tf) <- tf$TAXANAME
tblk <- readRDS(IN$trait_blocks)
cl <- setdiff(names(tf), c("Genus", "TAXANAME")); cl <- cl[vapply(tf[cl], is.numeric, logical(1))]
gw <- tblk$gowdis_w[cl]; stopifnot(all(is.finite(gw)))
tx <- colnames(den_mat)
MF <- as.matrix(tf[tx, cl]); storage.mode(MF) <- "double"
G_full <- as.matrix(suppressWarnings(FD::gowdis(as.data.frame(MF), w = gw)))
site_tab <- read_csv(IN$site_metrics, show_col_types = FALSE) %>% mutate(SiteID = std_site(SiteID))
site_tab <- site_tab[match(rownames(den_mat), site_tab$SiteID), ]
rng_full <- apply(MF, 2, function(v) diff(range(v, na.rm = TRUE))); rng_wo <- apply(MF[setdiff(tx, CHIRONOMIDAE), ], 2, function(v) diff(range(v, na.rm = TRUE)))
stopifnot("removing the four groups changes a Gower range; the removed geometry would differ from the published one" = all(abs(rng_full - rng_wo) < 1e-12))

fd_metrics <- function(W, drop = character()) {
  keep <- setdiff(tx, drop); X <- as.data.frame(MF[keep, ]); G <- G_full[keep, keep]
  a <- W[, keep]; P <- a / rowSums(a)
  Q <- rowSums((P %*% G) * P); Simp <- 1 - rowSums(P^2)
  fd <- suppressWarnings(FD::dbFD(x = X, a = a, w = gw, stand.x = TRUE, corr = "none", calc.FRic = FALSE, calc.FDiv = FALSE, calc.CWM = FALSE, messages = FALSE))
  tibble(SiteID = rownames(a), RaoQ = rowSums((P %*% G^2) * P) / 2, FDis = unname(fd$FDis[rownames(a)]), FR_RD = 1 - Q / Simp)
}
full <- lapply(CURRENCIES, fd_metrics); without <- lapply(CURRENCIES, fd_metrics, drop = CHIRONOMIDAE)
for (cur in names(CURRENCIES)) { sfx <- c(Density = "_den", Biomass = "_bio")[[cur]]
  for (m in c("RaoQ", "FDis", "FR_RD"))
    if (!(max(abs(full[[cur]][[m]] - site_tab[[paste0(m, sfx)]])) < 1e-9)) stop(sprintf("recomputed %s%s does not reproduce the published site value", m, sfx)) }
metric_sets <- list(Density = full$Density, Biomass = full$Biomass); removed_sets <- list(Density = without$Density, Biomass = without$Biomass)
# the redundancy divergence is the biomass-weighted minus the density-weighted overlap redundancy
div_kept <- full$Biomass$FR_RD - full$Density$FR_RD; div_removed <- without$Biomass$FR_RD - without$Density$FR_RD
stopifnot(max(abs(div_kept - site_tab$FR_DB_RD)) < 1e-9)

screen <- read_csv(file.path(OUT_TAB, "03_screen_all_metrics.csv"), show_col_types = FALSE)
# FDR of a replaced p within its own family and region, as step 03 corrects (BH within family, within region)
refdr <- function(rg, metric, p_new) { fam <- screen$Family[screen$Metric == metric][1]
  s <- screen %>% filter(Family == fam, Region == rg); stopifnot(isTRUE(all.equal(p.adjust(s$p_raw, "BH"), s$p_FDR)))
  s$p_raw[s$Metric == metric] <- p_new; p.adjust(s$p_raw, "BH")[s$Metric == metric] }
topk <- function(W, rg, k = 4) { d <- colSums(W[sites$Region == rg, ]); names(sort(d[setdiff(names(d), CHIRONOMIDAE)], decreasing = TRUE))[seq_len(k)] }
set.seed(SIM_SEED)
rows <- list()
for (cur in names(CURRENCIES)) for (m in c("RaoQ", "FDis", "FR_RD")) for (rg in REGIONS) {
  i <- sites$Region == rg; x <- sites$BIBI[i]; mm <- paste0(m, c(Density = "_den", Biomass = "_bio")[[cur]])
  k <- sp_cor(metric_sets[[cur]][[m]][i], x); r <- sp_cor(removed_sets[[cur]][[m]][i], x)
  pub <- screen %>% filter(Metric == mm, Region == rg)
  stopifnot(abs(pub$rho - k[["rho"]]) < 1e-9)
  t4 <- fd_metrics(CURRENCIES[[cur]], topk(CURRENCIES[[cur]], rg))[[m]][i]
  rows[[length(rows) + 1]] <- tibble(Weighting = cur, Metric = mm, Region = rg, rho_all = k[["rho"]], fdr_all = pub$p_FDR, rho_removed = r[["rho"]], p_removed = r[["p"]],
    fdr_removed = refdr(rg, mm, r[["p"]]), rho_top4_other = sp_cor(t4, x)[["rho"]], top4_other = paste(topk(CURRENCIES[[cur]], rg), collapse = ", ")) }
# divergence
for (rg in REGIONS) { i <- sites$Region == rg; x <- sites$BIBI[i]; k <- sp_cor(div_kept[i], x); r <- sp_cor(div_removed[i], x)
  pub <- screen %>% filter(Metric == "FR_DB_RD", Region == rg); stopifnot(abs(pub$rho - k[["rho"]]) < 1e-9)
  rows[[length(rows) + 1]] <- tibble(Weighting = "Divergence (biomass minus density)", Metric = "FR_DB_RD", Region = rg, rho_all = k[["rho"]], fdr_all = pub$p_FDR,
    rho_removed = r[["rho"]], p_removed = r[["p"]], fdr_removed = refdr(rg, "FR_DB_RD", r[["p"]]), rho_top4_other = NA_real_, top4_other = NA_character_) }
met_tab <- bind_rows(rows)
# random sets of four other taxa (Rao's Q and overlap redundancy, which are exact for any subset; functional dispersion is not rerun for random sets)
rand_tab <- purrr::map_dfr(c("Density", "Biomass"), function(cur) { W <- CURRENCIES[[cur]]
  purrr::map_dfr(REGIONS, function(rg) { i <- sites$Region == rg; x <- sites$BIBI[i]; occ <- colSums(W[i, ]); pool <- setdiff(names(occ)[occ > 0], CHIRONOMIDAE); pw <- occ[pool] / sum(occ[pool])
    ref <- metric_sets[[cur]]; refr <- removed_sets[[cur]]
    draws <- t(replicate(N_RANDOM_SETS, { d <- sample(pool, 4, prob = pw); keep <- setdiff(tx, d); a <- W[i, keep]; P <- a / rowSums(a); G <- G_full[keep, keep]
      Q <- rowSums((P %*% G) * P); c(RaoQ = sp_cor(rowSums((P %*% G^2) * P) / 2, x)[["rho"]], FR_RD = sp_cor(1 - Q / (1 - rowSums(P^2)), x)[["rho"]]) }))
    purrr::map_dfr(c("RaoQ", "FR_RD"), function(m) { base <- sp_cor(ref[[m]][i], x)[["rho"]]; chi <- sp_cor(refr[[m]][i], x)[["rho"]]
      tibble(Weighting = cur, Metric = paste0(m, c(Density = "_den", Biomass = "_bio")[[cur]]), Region = rg, random4_p05 = unname(quantile(draws[, m], .05)), random4_p95 = unname(quantile(draws[, m], .95)),
             random4_as_extreme = mean(abs(draws[, m] - base) >= abs(chi - base))) }) }) })
met_tab <- met_tab %>% left_join(rand_tab, by = c("Weighting", "Metric", "Region")) %>%
  mutate(met_all = abs(rho_all) >= 0.30 & fdr_all <= FDR_MAX, met_removed = abs(rho_removed) >= 0.30 & fdr_removed <= FDR_MAX)
write_csv(met_tab, file.path(OUT_TAB, "20_chironomid_functional_metrics.csv"))
say("  density-weighted RaoQ, FDis, FR_RD meeting the criterion: %d of 9 with every taxon, %d of 9 with the four groups removed",
    sum(met_tab$met_all[met_tab$Weighting == "Density"]), sum(met_tab$met_removed[met_tab$Weighting == "Density"]))

# ---- 3. Fourth-corner states, kept and removed (legacy Table S28; not in this paper) -----------------------------------------------------------------------
tr_raw <- as.data.frame(readRDS(IN$traits)); rownames(tr_raw) <- tr_raw$TAXANAME
Q_raw <- as.matrix(tr_raw[colnames(den_mat), setdiff(names(tr_raw), c("Genus", "TAXANAME"))]); storage.mode(Q_raw) <- "double"
blk <- sub("_.*$", "", colnames(Q_raw)); Q_norm <- Q_raw
for (b in unique(blk)) { i <- blk == b; Q_norm[, i] <- Q_raw[, i, drop = FALSE] / rowSums(Q_raw[, i, drop = FALSE]) }
fc_one_seed <- function(k, W, drop) {
  set.seed(SEED + k)
  purrr::map_dfr(REGIONS, function(rg) { s <- sites$Region == rg
    purrr::map_dfr(unique(blk), function(b) {
      i <- blk == b; ok <- !apply(Q_norm[, i, drop = FALSE], 1, function(r) any(is.na(r)))
      L <- W[s, , drop = FALSE]; L <- L[, colSums(L) > 0 & ok[colnames(L)] & !(colnames(L) %in% drop), drop = FALSE]
      if (ncol(L) < 5) return(NULL)
      Qb <- as.data.frame(Q_norm[colnames(L), i, drop = FALSE]); R <- data.frame(BIBI = sites$BIBI[s])
      fc <- fourthcorner(R, as.data.frame(L), Qb, modeltype = 6, nrepet = N_PERM_FC, p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
      tibble(seed = SEED + k, Region = rg, block = b, term = fc$tabD$names, stat = fc$tabD$obs, p_raw = fc$tabD$pvalue, p_fdr = fc$tabD$adj.pvalue) }) }) }
fc_run <- function(W, drop) bind_rows(parallel::mclapply(seq_len(N_SEEDS_FC) - 1, fc_one_seed, W = W, drop = drop, mc.cores = min(N_SEEDS_FC, max(1, parallel::detectCores() - 1))))
fc_all <- purrr::map_dfr(names(CURRENCIES), function(cur) bind_rows(
  fc_run(CURRENCIES[[cur]], character()) %>% mutate(Weighting = cur, chironomids = "kept"),
  fc_run(CURRENCIES[[cur]], CHIRONOMIDAE) %>% mutate(Weighting = cur, chironomids = "removed")))
fc_sum <- fc_all %>% group_by(Weighting, chironomids, Region, block, term) %>%
  summarise(stat = first(stat), n_seeds_pass = sum(p_fdr <= FDR_MAX), p_fdr_max = max(p_fdr), .groups = "drop") %>% mutate(state = sub("^BIBI / ", "", term))
rec <- bind_rows(read_csv(file.path(OUT_TAB, "15_fourthcorner_per_block.csv"), show_col_types = FALSE) %>% mutate(Weighting = "Density"),
                 read_csv(file.path(OUT_TAB, "15_fourthcorner_per_block_biomass.csv"), show_col_types = FALSE) %>% mutate(Weighting = "Biomass"))
ctl <- fc_sum %>% filter(chironomids == "kept") %>% inner_join(rec %>% select(Weighting, Region, block, term, nsp = n_seeds_pass, pm = p_fdr_max, st = stat), by = c("Weighting", "Region", "block", "term"))
stopifnot("the kept run does not reproduce the pass counts of step 15" = nrow(ctl) == nrow(rec) && all(ctl$n_seeds_pass == ctl$nsp),
          "the kept run does not reproduce the p values of step 15" = max(abs(ctl$p_fdr_max - ctl$pm)) < 1e-9,
          "the kept run does not reproduce the statistics of step 15" = max(abs(ctl$stat - ctl$st)) < 1e-9)
say("  fourth-corner control: %d state-tests reproduce step 15 (pass counts, p values, statistics)", nrow(ctl))
st_wide <- fc_sum %>% select(Weighting, Region, block, state, chironomids, stat, n_seeds_pass, p_fdr_max) %>%
  pivot_wider(names_from = chironomids, values_from = c(stat, n_seeds_pass, p_fdr_max)) %>%
  mutate(pass_all = n_seeds_pass_kept == N_SEEDS_FC, pass_removed = n_seeds_pass_removed == N_SEEDS_FC)
# the share of the state's weighted mass held by the four groups, regional mean over streams
mass_share <- purrr::map_dfr(names(CURRENCIES), function(cur) { W <- CURRENCIES[[cur]]
  purrr::map_dfr(REGIONS, function(rg) { i <- sites$Region == rg; P <- W[i, ] / rowSums(W[i, ]); ch <- colnames(P) %in% CHIRONOMIDAE
    purrr::map_dfr(colnames(Q_norm), function(st) { q <- Q_norm[, st]; u <- !is.na(q)
      tibble(Weighting = cur, Region = rg, state = st, chiro_part_of_mass = mean(P[, u & ch, drop = FALSE] %*% q[u & ch]) / mean(P[, u, drop = FALSE] %*% q[u])) }) }) })
state_tab <- st_wide %>% filter(pass_all | pass_removed) %>% left_join(mass_share, by = c("Weighting", "Region", "state")) %>%
  mutate(verdict_all = ifelse(pass_all, "passes", ifelse(n_seeds_pass_kept > 0, "borderline", "does not pass")),
         verdict_removed = ifelse(pass_removed, "passes", ifelse(n_seeds_pass_removed > 0, "borderline", "does not pass"))) %>%
  arrange(factor(Weighting, c("Density", "Biomass")), factor(Region, REGIONS), block, state)
write_csv(state_tab, file.path(OUT_TAB, "20_chironomid_fourthcorner_states.csv"))
write_csv(fc_all, file.path(OUT_TAB, "20_chironomid_fourthcorner_by_seed.csv"))
cnt <- st_wide %>% group_by(Weighting, Region) %>% summarise(passing_all = sum(pass_all), passing_removed = sum(pass_removed), lost = sum(pass_all & !pass_removed), gained = sum(!pass_all & pass_removed), .groups = "drop")
write_csv(cnt, file.path(OUT_TAB, "20_chironomid_fourthcorner_counts.csv"))
for (cur in c("Density", "Biomass")) say("  fourth-corner states passing under all five seeds, %s: %d with every taxon; %d lost and %d gained with the four groups removed",
    tolower(cur), sum(cnt$passing_all[cnt$Weighting == cur]), sum(cnt$lost[cnt$Weighting == cur]), sum(cnt$gained[cnt$Weighting == cur]))

# ---- 4. Genus-level split of the four groups (Table S7 part C; legacy Table S29) ---------------------------------------------------------------------------
gens <- unlist(lapply(GENERA_OF_GROUP, function(v) intersect(v, rownames(tf)))); stopifnot(!anyDuplicated(gens), length(gens) >= 30)
grp_of <- rep(names(GENERA_OF_GROUP), vapply(GENERA_OF_GROUP, function(v) length(intersect(v, rownames(tf))), integer(1))); names(grp_of) <- gens
others <- setdiff(tx, CHIRONOMIDAE); allt <- c(others, CHIRONOMIDAE, gens)
MFa <- as.matrix(tf[allt, cl]); storage.mode(MFa) <- "double"
Ga <- as.matrix(suppressWarnings(FD::gowdis(as.data.frame(MFa), w = gw))); dimnames(Ga) <- list(allt, allt); Ga2 <- Ga^2
gmetrics <- function(P) { Q <- rowSums((P %*% Ga) * P); S <- 1 - rowSums(P^2); cbind(RaoQ = rowSums((P %*% Ga2) * P) / 2, FR_RD = 1 - Q / S) }
set.seed(SIM_SEED)
sim_rows <- list()
for (rg in REGIONS) {
  i <- sites$Region == rg; x <- sites$BIBI[i]; z <- as.numeric(scale(x)); A <- matrix(0, sum(i), length(allt), dimnames = list(rownames(den_mat)[i], allt)); A[, c(others, CHIRONOMIDAE)] <- den_mat[i, c(others, CHIRONOMIDAE)]
  M0 <- gmetrics(A / rowSums(A)); base <- c(RaoQ = sp_cor(M0[, "RaoQ"], x)[["rho"]], FR_RD = sp_cor(M0[, "FR_RD"], x)[["rho"]])
  oth <- others[colSums(A[, others]) > 0]
  one <- function(mode, beta = 0) {
    B <- A; B[, CHIRONOMIDAE] <- 0; t_g <- rnorm(length(gens)); names(t_g) <- gens
    fixed <- lapply(names(GENERA_OF_GROUP), function(g) { v <- rgamma(sum(grp_of == g), 1); v / sum(v) }); names(fixed) <- names(GENERA_OF_GROUP)
    Gm <- Ga; Gm2 <- Ga2; cols <- allt                                   # every column; the four group columns are zero after the split
    if (mode == "donor") { dd <- sample(oth, length(gens), replace = TRUE); names(dd) <- gens; idx <- c(others, dd); Gm <- Ga[idx, idx]; Gm2 <- Gm^2; cols <- c(others, gens) }
    for (g in names(GENERA_OF_GROUP)) { gg <- names(grp_of)[grp_of == g]; d <- A[, g]
      W <- if (mode %in% c("fixed", "donor")) matrix(fixed[[g]], nrow(A), length(gg), byrow = TRUE) else if (mode == "site") { m <- matrix(rgamma(nrow(A) * length(gg), 1), nrow(A)); m / rowSums(m) }
           else { m <- matrix(rep(fixed[[g]], each = nrow(A)), nrow(A)) * exp(beta * outer(z, t_g[gg])); m / rowSums(m) }
      B[, gg] <- d * W }
    P <- B[, cols] / rowSums(B[, cols]); Q <- rowSums((P %*% Gm) * P); S <- 1 - rowSums(P^2)
    c(RaoQ = sp_cor(rowSums((P %*% Gm2) * P) / 2, x)[["rho"]], FR_RD = sp_cor(1 - Q / S, x)[["rho"]]) }
  for (sc in list(list("Profiles of the four taxa (as analyzed)", NA, 0), list("Fixed genus mix at every stream", "fixed", 0), list("Genus mix varies randomly among streams", "site", 0),
                  list("Genus mix tracks condition, moderate", "track", 1), list("Genus mix tracks condition, strong", "track", 2), list("Genera take the profile of a random other taxon", "donor", 0))) {
    if (is.na(sc[[2]])) { for (m in c("RaoQ", "FR_RD")) sim_rows[[length(sim_rows) + 1]] <- tibble(Region = rg, Scenario = sc[[1]], Metric = m, rho_median = base[[m]], rho_p05 = NA_real_, rho_p95 = NA_real_, share_in_baseline_direction_abs_ge_030 = NA_real_); next }
    r <- t(replicate(N_GENUS_DRAWS, one(sc[[2]], sc[[3]])))
    for (m in c("RaoQ", "FR_RD")) sim_rows[[length(sim_rows) + 1]] <- tibble(Region = rg, Scenario = sc[[1]], Metric = m, rho_median = median(r[, m]), rho_p05 = unname(quantile(r[, m], .05)), rho_p95 = unname(quantile(r[, m], .95)),
      share_in_baseline_direction_abs_ge_030 = mean(sign(r[, m]) == sign(base[[m]]) & abs(r[, m]) >= .30)) } }
sim_tab <- bind_rows(sim_rows) %>% mutate(Metric = c(RaoQ = "Rao's quadratic entropy", FR_RD = "Overlap redundancy")[Metric], Weighting = "Density") %>%
  arrange(factor(Region, REGIONS), Metric, match(Scenario, unique(Scenario)))
write_csv(sim_tab, file.path(OUT_TAB, "20_chironomid_genus_mix_simulation.csv"))
say("  genus-level simulation: %d scenarios x %d draws x %d regions", length(unique(sim_tab$Scenario)) - 1, N_GENUS_DRAWS, length(REGIONS))

# ---- 5. Figure S1 (legacy Figure S3; written as figS5_chironomid_removal and overwritten by step 20b) ----------------------------------------------------------
RL <- c("Valley & Ridge" = "Valley and Ridge", "Piedmont" = "Piedmont", "Coastal Plain" = "Coastal Plain")
COL2 <- c("Every taxon" = "#1A1A1A", "Four chironomid taxa removed" = "#8C8C8C")   # one neutral pair: blue is the Piedmont's color elsewhere
fm <- met_tab %>% filter(Weighting == "Density") %>%
  mutate(Metric = factor(c(FDis_den = "Functional dispersion", RaoQ_den = "Rao's quadratic entropy", FR_RD_den = "Overlap redundancy")[Metric], c("Functional dispersion", "Rao's quadratic entropy", "Overlap redundancy")),
         Region = factor(RL[Region], rev(RL[REGIONS]))) %>%
  tidyr::pivot_longer(c(rho_all, rho_removed, met_all, met_removed), names_to = c(".value", "set"), names_pattern = "(rho|met)_(all|removed)") %>%
  mutate(set = factor(ifelse(set == "all", "Every taxon", "Four chironomid taxa removed"), names(COL2)))
pa <- ggplot(fm, aes(rho, Region)) + geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") +
  geom_line(aes(group = interaction(Region, Metric)), linewidth = .5, colour = "grey55") + geom_point(aes(colour = set, shape = met), size = 2.4, stroke = .8, position = position_nudge(y = ifelse(fm$set == "Every taxon", 0.08, -0.08))) + facet_wrap(~Metric, nrow = 1) +
  scale_colour_manual(values = COL2) + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "Meets the criterion", `FALSE` = "Does not")) +
  scale_x_continuous(limits = c(-.8, .8), breaks = seq(-.6, .6, .3)) + labs(x = "Spearman ρ with B-IBI, weighted by density", y = NULL, title = "(a) Functional metrics") +
  theme_pub() + theme(legend.position = "bottom", legend.title = element_blank()) + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
# mirror, fs and pb build the former panel (b), the trait states with and without the four groups. That panel left this
# paper with the fourth-corner analysis (12 Oct 2026): pb is built but never drawn or saved, and no output depends on it.
mirror <- c("Disperse_Low", "Flying_Weak", "Respiration_Gills", "Synchro_Poorly")
fs <- state_tab %>% filter(!state %in% mirror) %>% transmute(Weighting = factor(Weighting, c("Density", "Biomass")), Region = factor(RL[Region], RL[REGIONS]),
        lab = paste(TRAIT_GROUP_SHORT[block], TRAIT_STATE_LAB[state], sep = ": "), `Every taxon` = stat_kept, `Four chironomid groups removed` = stat_removed,
        p1 = pass_all, p2 = pass_removed) %>%
  tidyr::pivot_longer(c(`Every taxon`, `Four chironomid groups removed`), names_to = "set", values_to = "stat") %>%
  mutate(pass = ifelse(set == "Every taxon", p1, p2), set = factor(set, names(COL2)), lab = factor(lab, rev(sort(unique(lab)))))
pb <- ggplot(fs, aes(stat, lab)) + geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") + geom_line(aes(group = interaction(lab, Region, Weighting)), linewidth = .5, colour = "grey55") +
  geom_point(aes(colour = set, shape = pass), size = 2.1, stroke = .8, position = position_nudge(y = ifelse(fs$set == "Every taxon", 0.17, -0.17))) + facet_grid(Weighting ~ Region, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = COL2) + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "Passes under all five seeds", `FALSE` = "Does not pass")) +
  scale_x_continuous(limits = c(-.55, .55), breaks = c(-.4, -.2, 0, .2, .4)) + labs(x = "Fourth-corner statistic", y = NULL, title = "(b) Trait states that passed under either analysis") +
  theme_pub() + theme(legend.position = "bottom", legend.title = element_blank(), axis.text.y = element_text(size = 7)) + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
# 12 Oct 2026: panel (b) left with the fourth-corner analysis (R/supplement_asset_map.csv); the one remaining panel takes no "(a)" title.
# The figure is also redrawn from this step's table by R/20b_figS1_from_record.R without recomputing anything.
save_fig(pa + labs(title = NULL), "figS5_chironomid_removal", FIG_W_DOUBLE, 3.2)   # overwritten by step 20b (two panels: the chironomid share, then these metrics without Rao)
say("STEP 20 complete")
