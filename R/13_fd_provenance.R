# =============================================================================
# 13_fd_provenance.R
# =============================================================================
# WHY THIS STEP EXISTS
#
#   Step 01 proves the SCREEN is faithful to the metric table: it recomputes
#   Spearman coefficients from the stored values in the site table. What it
#   cannot do -- structurally -- is prove the metric table itself is right,
#   because it reads the metrics rather than deriving them.
#
#   Until 2026-08-27 nothing closed that gap, and three attempts to close it by
#   hand failed (a hand-written unweighted Gower; dbFD on the legacy binary
#   trait file; a saved dbFD object). This step closes it. The canonical metrics are reproduced to machine
#   precision from raw inputs, and the run halts if that ever stops being true.
#
# WHAT THE ANSWER TURNED OUT TO BE
#
#   1. The trait matrices are the PROXY-AWARE pair, not traits_sel.rds:
#        traits_fuzzy_proxy_main.rds  +  traits_binary_proxy_main.rds
#      The upstream build ran with ANALYSIS_VERSION = "main", and its trait-file
#      selector (pick_trait_file(), not part of this repository) resolves to
#      *_proxy_main. The binary matrix matters even for the fuzzy result,
#      because the column-retention rule is
#        keep column if nonzero in EITHER coding.
#
#   2. Gower weights come from trait_blocks.rds$gowdis_w (1 / block size),
#      NOT from parsing the column names. Both give 41 weights; they are not
#      the same 41.
#
#   3. Every canonical metric is read from the Cailliez-corrected run. RaoQ and
#      FDis do not depend on the correction: in FD::dbFD (version 1.0-12) both are
#      computed on the UNCORRECTED Gower distance (RaoQ <- divc(..., x.dist, ...);
#      disp <- fdisp(x.dist, a)), and only FDiv, FEve and FRic pass through the
#      corrected principal coordinates. The sqrt correction cannot be applied to
#      this trait matrix (a trait group with no record is NA), so the correction
#      sensitivity in section 4 is written as a one-row note.
#
#   4. FRic reproduces in the fuzzy-pool geometry (section 3b verifies it to
#      machine precision). The earlier statement that FRic did not reproduce
#      applied to the retired binary-pool run.
#
# OUTPUT
#   output/tables/13_fd_provenance_verification.csv   repository table fd_metric_provenance.csv (legacy Table S8)
#   output/tables/13_fd_metrics_recomputed.csv
#   output/tables/13_redundancy_provenance.csv
#   output/tables/13_correction_sensitivity.csv
#   (the last three are records; no later step reads them)
# =============================================================================

say("STEP 13 -- FD metric provenance: derive the canonical metrics from raw inputs")

suppressPackageStartupMessages({ library(FD); library(stringr) })

# ---- 1. Rebuild the trait matrices exactly as the upstream build did ----------
# Mirrors the trait-matrix standardization of the upstream Step 3 build (the
# script that produced the site table; not part of this repository), so the same columns are retained.

.std_trait_df <- function(df) {
  if (!"TAXANAME" %in% names(df)) df$TAXANAME <- df$Genus
  if (!"Genus"    %in% names(df)) df$Genus    <- df$TAXANAME
  df %>%
    mutate(Genus    = str_squish(as.character(Genus)),
           TAXANAME = str_squish(as.character(TAXANAME))) %>%
    filter(!is.na(TAXANAME), TAXANAME != "") %>%
    distinct(TAXANAME, .keep_all = TRUE)
}

tr_fuz <- .std_trait_df(readRDS(IN$traits_fuz))
tr_bin <- .std_trait_df(readRDS(IN$traits_bin))
tblk   <- readRDS(IN$trait_blocks)

id_cols   <- c("Genus", "TAXANAME")
cols_fuz  <- setdiff(names(tr_fuz), id_cols)
cols_bin  <- setdiff(names(tr_bin), id_cols)
trait_cols <- intersect(cols_bin, cols_fuz)
trait_cols <- unique(trait_cols[
  vapply(tr_bin[trait_cols], is.numeric, logical(1)) &
  vapply(tr_fuz[trait_cols], is.numeric, logical(1))])

fd_taxa <- intersect(union(colnames(den_mat), colnames(bio_mat)),
                     intersect(tr_bin$TAXANAME, tr_fuz$TAXANAME))

.build <- function(df) {
  d <- df %>% filter(TAXANAME %in% fd_taxa) %>% distinct(TAXANAME, .keep_all = TRUE)
  d <- d[match(fd_taxa, d$TAXANAME), , drop = FALSE]
  m <- as.matrix(d[, trait_cols, drop = FALSE])
  rownames(m) <- d$TAXANAME
  storage.mode(m) <- "double"
  m   # NA kept: an empty trait group is missing, not zero (Twardochleb fuzzy traits, 22 Sep 2026)
}
M_fuz <- .build(tr_fuz); M_bin <- .build(tr_bin)

# Retain a column if it is nonzero in EITHER coding. Dropping this rule changes
# the retained set and the metrics stop reproducing.
keep  <- (colSums(M_bin, na.rm = TRUE) > 0) | (colSums(M_fuz, na.rm = TRUE) > 0)
M_fuz <- M_fuz[, keep, drop = FALSE]; M_bin <- M_bin[, keep, drop = FALSE]

gw <- tblk$gowdis_w[colnames(M_bin)]
if (!(length(gw) == ncol(M_bin) && all(is.finite(gw)))) {
  sl <- tblk$state_lookup
  g  <- sl$Trait_group[match(colnames(M_bin), sl$ColName)]
  gw <- 1 / as.numeric(table(g)[g]); names(gw) <- colnames(M_bin)
  say("  NOTE: gowdis_w incomplete; used the state_lookup fallback")
}

say("  trait matrices rebuilt: %d taxa x %d modalities | weights from %s",
    nrow(M_fuz), ncol(M_fuz), basename(IN$trait_blocks))

# ---- 2. Run dbFD under both corrections, both currencies --------------------
run_fd <- function(a, corr) {
  a <- a[, intersect(colnames(a), rownames(M_fuz)), drop = FALSE]
  suppressWarnings(FD::dbFD(
    x = M_fuz, a = a, w = gw, stand.x = TRUE, corr = corr, w.abun = TRUE,
    calc.FRic = TRUE, calc.FDiv = TRUE, stand.FRic = TRUE, calc.CWM = FALSE,
    m = 10, messages = FALSE))
}

# The sqrt correction fails once empty groups are NA (the distance stays non-Euclidean); RaoQ and FDis do
# not depend on the correction (section 4), so every canonical metric is now read from the Cailliez run.
.try_fd <- function(a, corr) tryCatch(run_fd(a, corr), error = function(e) NULL)
fd <- list(
  den_sqrt = .try_fd(den_mat, "sqrt"),    den_cail = run_fd(den_mat, "cailliez"),
  bio_sqrt = .try_fd(bio_mat, "sqrt"),    bio_cail = run_fd(bio_mat, "cailliez")
)

# ---- 3. Verify each canonical metric against the site metric table ----------
# Every metric is compared with the Cailliez run (header, point 3).
CANON <- tribble(
  ~Metric,     ~Currency, ~fd_name, ~corr,
  "RaoQ_den",  "den",     "RaoQ",   "cailliez",
  "FDis_den",  "den",     "FDis",   "cailliez",
  "FDiv_den",  "den",     "FDiv",   "cailliez",
  "FEve_den",  "den",     "FEve",   "cailliez",
  "RaoQ_bio",  "bio",     "RaoQ",   "cailliez",
  "FDis_bio",  "bio",     "FDis",   "cailliez",
  "FDiv_bio",  "bio",     "FDiv",   "cailliez",
  "FEve_bio",  "bio",     "FEve",   "cailliez"
)

sites <- dat$SiteID
verify13 <- CANON %>%
  rowwise() %>%
  mutate(
    n_compared = {
      obj <- fd[[paste0(Currency, "_", substr(corr, 1, 4))]]
      v   <- if (is.null(obj[[fd_name]])) rep(NA_real_, length(sites))
             else obj[[fd_name]][match(sites, names(obj[[fd_name]]))]
      p   <- dat[[Metric]][match(sites, dat$SiteID)]
      sum(is.finite(v) & is.finite(p))
    },
    max_abs_diff = {
      obj <- fd[[paste0(Currency, "_", substr(corr, 1, 4))]]
      v   <- if (is.null(obj[[fd_name]])) rep(NA_real_, length(sites))
             else obj[[fd_name]][match(sites, names(obj[[fd_name]]))]
      p   <- dat[[Metric]][match(sites, dat$SiteID)]
      ok  <- is.finite(v) & is.finite(p)
      if (!any(ok)) NA_real_ else max(abs(v[ok] - p[ok]))
    },
    # A metric that was never computed must FAIL. max(numeric(0)) is -Inf and
    # -Inf < tol is TRUE, which would turn "not computed" into "REPRODUCED".
    status = dplyr::case_when(
      n_compared < length(sites) ~ "INCOMPLETE",
      !is.finite(max_abs_diff)   ~ "NOT COMPUTED",
      max_abs_diff < 1e-8        ~ "REPRODUCED",
      TRUE                       ~ "MISMATCH")
  ) %>% ungroup()

recomputed <- tibble(SiteID = sites)
for (i in seq_len(nrow(CANON))) {
  obj <- fd[[paste0(CANON$Currency[i], "_", substr(CANON$corr[i], 1, 4))]]
  recomputed[[CANON$Metric[i]]] <-
    obj[[CANON$fd_name[i]]][match(sites, names(obj[[CANON$fd_name[i]]]))]
}

write_csv(verify13,   file.path(OUT_TAB, "13_fd_provenance_verification.csv"))
write_csv(recomputed, file.path(OUT_TAB, "13_fd_metrics_recomputed.csv"))

for (i in seq_len(nrow(verify13))) {
  say("    %-9s corr=%-8s n=%2d  max abs diff = %-9s %s", verify13$Metric[i],
      verify13$corr[i], verify13$n_compared[i],
      ifelse(is.finite(verify13$max_abs_diff[i]),
             sprintf("%.2e", verify13$max_abs_diff[i]), "--"),
      verify13$status[i])
}

n_bad <- sum(verify13$status != "REPRODUCED")
if (n_bad > 0) {
  print(as.data.frame(verify13 %>% filter(status == "MISMATCH")))
  stop("STEP 13 FAILED: ", n_bad, " of ", nrow(verify13),
       " FD metrics no longer reproduce from raw inputs.\n",
       "  An input file or the trait matrix construction has drifted.\n",
       "  Report: output/tables/13_fd_provenance_verification.csv")
}

say("  VERIFIED: all %d canonical FD metrics reproduced from raw trait and", nrow(verify13))
say("  community matrices, to machine precision.")

# ---- 3b. The redundancy family (added 18 September 2026) ---------------------
# The checks above do not cover FR_RD, FR_SF or FRic. FR_RD = 1 - Q/D with Q the
# linear Rao's Q on the Gower distance of the geometry in force (fuzzy,
# group-weighted when CH2_RD_GEOMETRY=fuzzy; the original binary, unweighted
# otherwise) and D Simpson's index on the trait-eligible community matrix.
# FRic_den in the original geometry is the upstream binary-pool dbFD run:
# identical binary profiles collapsed into one representative (abundances
# summed), sqrt correction, group weights, m = min(taxa - 1, 10), stand.FRic;
# in the fuzzy geometry it is the fuzzy-pool run (Cailliez, no collapsing). FR_SF_den_log10 =
# log10(Richness_FD / FRic_den + 1), richness = trait-eligible taxa present.
geom_fuzzy <- identical(Sys.getenv("CH2_RD_GEOMETRY"), "fuzzy")
G_rd <- if (geom_fuzzy) {
  as.matrix(suppressWarnings(FD::gowdis(as.data.frame(M_fuz), w = gw)))
} else {
  as.matrix(FD::gowdis(as.data.frame(M_bin)))
}
.P <- function(a) { a <- a[sites, rownames(G_rd), drop = FALSE]; a / rowSums(a) }
.q <- function(P) rowSums((P %*% G_rd) * P)
P_den13 <- .P(den_mat); P_bio13 <- .P(bio_mat)
D_den13 <- vegan::diversity(den_mat[sites, rownames(G_rd), drop = FALSE], index = "simpson")
D_bio13 <- vegan::diversity(bio_mat[sites, rownames(G_rd), drop = FALSE], index = "simpson")
FR_RD_den13 <- 1 - .q(P_den13) / D_den13; FR_RD_bio13 <- 1 - .q(P_bio13) / D_bio13

# FRic_den: dbFD on the deduplicated binary matrix, as upstream
.keys <- apply(M_bin, 1, paste, collapse = "_"); .rep <- rownames(M_bin)[!duplicated(.keys)]
a_dd <- sapply(.rep, function(r) { mem <- rownames(M_bin)[.keys == .keys[match(r, rownames(M_bin))]]
                                   rowSums(den_mat[sites, intersect(mem, colnames(den_mat)), drop = FALSE]) })
fd_bin <- tryCatch(suppressWarnings(FD::dbFD(x = M_bin[.rep, , drop = FALSE], a = a_dd, w = gw, stand.x = TRUE,
                                             corr = "sqrt", w.abun = TRUE, calc.FRic = TRUE, calc.FDiv = FALSE,
                                             calc.CWM = FALSE, stand.FRic = TRUE, m = min(length(.rep) - 1, 10L), messages = FALSE)),
                   error = function(e) NULL)
FRic13_bin <- if (is.null(fd_bin)) rep(NA_real_, length(sites)) else fd_bin$FRic[match(sites, names(fd_bin$FRic))]
# In the fuzzy geometry FRic_den is the fuzzy-pool run (Cailliez, group weights, no collapsing): the same dbFD
# call that reproduces FDiv and FEve above, so its FRic is taken from there.
FRic13_fuz <- fd$den_cail$FRic[match(sites, names(fd$den_cail$FRic))]
FRic13 <- if (geom_fuzzy) FRic13_fuz else FRic13_bin
FRic13_label <- if (geom_fuzzy) "fuzzy pool dbFD, cailliez" else "binary pool dbFD, sqrt, deduplicated"
Rich13 <- rowSums(den_mat[sites, rownames(G_rd), drop = FALSE] > 0)
FR_SF13 <- log10(Rich13 / FRic13 + 1)

.cmp <- function(v, m) { p <- dat[[m]][match(sites, dat$SiteID)]; ok <- is.finite(v) & is.finite(p)
  tibble(Metric = m, Currency = ifelse(grepl("bio", m), "bio", "den"),
         fd_name = m, corr = if (geom_fuzzy) "fuzzy weighted Gower" else "binary unweighted Gower",
         n_compared = sum(ok), max_abs_diff = if (any(ok)) max(abs(v[ok] - p[ok])) else NA_real_) }
red13 <- bind_rows(.cmp(FR_RD_den13, "FR_RD_den"), .cmp(FR_RD_bio13, "FR_RD_bio"),
                   .cmp(FR_RD_bio13 - FR_RD_den13, "FR_DB_RD"),
                   .cmp(FRic13, "FRic_den") %>% mutate(corr = FRic13_label),
                   .cmp(FR_SF13, "FR_SF_den_log10") %>% mutate(corr = "Richness_FD / FRic_den, log10(x + 1)")) %>%
  mutate(status = dplyr::case_when(n_compared < length(sites) ~ "INCOMPLETE",
                                   !is.finite(max_abs_diff) ~ "NOT COMPUTED",
                                   max_abs_diff < 1e-8 ~ "REPRODUCED", TRUE ~ "MISMATCH"))
write_csv(red13, file.path(OUT_TAB, "13_redundancy_provenance.csv"))
verify13 <- bind_rows(verify13, red13)
write_csv(verify13, file.path(OUT_TAB, "13_fd_provenance_verification.csv"))
for (i in seq_len(nrow(red13)))
  say("    %-16s %-38s n=%2d  max abs diff = %-9s %s", red13$Metric[i], red13$corr[i], red13$n_compared[i],
      ifelse(is.finite(red13$max_abs_diff[i]), sprintf("%.2e", red13$max_abs_diff[i]), "--"), red13$status[i])
if (any(red13$status != "REPRODUCED"))
  stop("STEP 13 FAILED: the redundancy family does not reproduce in the geometry in force (CH2_RD_GEOMETRY=",
       Sys.getenv("CH2_RD_GEOMETRY"), ").")

# ---- 4. Correction sensitivity ----------------------------------------------
# How much does the choice of correction matter? Same logic as step 03's BH-scope sensitivity and step
# 04's transformation sensitivity -- compute both, report the difference, let
# the decision be made on numbers.

if (is.null(fd$den_sqrt) || is.null(fd$bio_sqrt)) {
  say("  correction sensitivity not run: the sqrt correction fails when empty trait groups are NA")
  write_csv(tibble(note = "sqrt correction not applicable: empty trait groups are NA (22 Sep 2026)"), file.path(OUT_TAB, "13_correction_sensitivity.csv"))
} else {
corr_sens <- purrr::map_dfr(c("den", "bio"), function(cy) {
  purrr::map_dfr(c("RaoQ", "FDis", "FDiv", "FEve"), function(nm) {
    purrr::map_dfr(REGIONS, function(rg) {
      s  <- dat$SiteID[dat$Region == rg]
      b  <- dat$BIBI[dat$Region == rg]
      gv <- function(cr) {
        o <- fd[[paste0(cy, "_", substr(cr, 1, 4))]][[nm]]
        o[match(s, names(o))]
      }
      r1 <- suppressWarnings(cor(b, gv("sqrt"),     method = "spearman", use = "complete.obs"))
      r2 <- suppressWarnings(cor(b, gv("cailliez"), method = "spearman", use = "complete.obs"))
      tibble(Metric = paste0(nm, "_", cy), Region = rg,
             rho_sqrt = r1, rho_cailliez = r2, delta = r2 - r1)
    })
  })
})
write_csv(corr_sens, file.path(OUT_TAB, "13_correction_sensitivity.csv"))

inv <- corr_sens %>% filter(Metric %in% grep("^(RaoQ|FDis)", unique(Metric), value = TRUE))
mx  <- max(abs(inv$delta), na.rm = TRUE)
mvd <- corr_sens %>% filter(abs(delta) >= 0.10)

say("")
say("  CORRECTION SENSITIVITY (sqrt vs cailliez, Spearman vs B-IBI):")
say("    RaoQ and FDis   : max |delta| = %.4f across %d region x weighting tests", mx, nrow(inv))
say("    FDiv and FEve   : %d of %d tests move by >= 0.10", nrow(mvd), nrow(corr_sens) - nrow(inv))
for (i in seq_len(nrow(mvd)))
  say("        %-9s %-15s %+.3f -> %+.3f", mvd$Metric[i], mvd$Region[i],
      mvd$rho_sqrt[i], mvd$rho_cailliez[i])
say("    tests where |rho| crosses %.2f between the corrections: %d", RHO_MIN,
    sum((abs(corr_sens$rho_sqrt) >= RHO_MIN) != (abs(corr_sens$rho_cailliez) >= RHO_MIN), na.rm = TRUE))
say("  RaoQ and FDis are computed on the uncorrected Gower distance (FD::dbFD), so they cannot")
say("  depend on the correction; only FDiv, FEve and FRic pass through the corrected PCoA.")
}
say("STEP 13 complete")
