# =============================================================================
# 01_load_and_verify.R — load canonical inputs and PROVE they are the right ones
# =============================================================================
# WHY THIS STEP EXISTS
#   The single biggest risk in re-running this analysis is picking up a stale or
#   wrong input file and not noticing (see the B-IBI note in 00_config.R). So this
#   step does not just load. It checks every site's B-IBI against its ICPRB release
#   value (the provenance gate; the run stops on any difference), checks the shape
#   of the inputs (90 sites, 30 per region, every screened metric present, every
#   site in both matrices), and recomputes a table of reference correlations. The
#   reference values come from the earlier analysis on the observed index (chessie_bibi); on the
#   rarefied-median basis the manuscript uses they are expected to differ and are
#   reported, not enforced. On the observed basis they halt the run if they drift.
#
# INPUTS
#   every file in IN (00_config.R) and the three ICPRB scoring files in the data folder
#
# OUTPUT (read by every later step; no paper table or figure comes from this step directly)
#   dat            : 90 rows, one per site: SiteID, Region, BIBI, every screened metric
#   den_mat/bio_mat: site x taxon matrices of the trait-eligible taxa
#   output/tables/01_verification_report.csv
# =============================================================================

say("STEP 01 -- load and verify canonical inputs")

# ---- 1. Load ----------------------------------------------------------------
stopifnot(all(file.exists(unlist(IN))))

metrics_raw <- read_csv(IN$site_metrics, show_col_types = FALSE) %>%
  mutate(SiteID = std_site(SiteID)) %>%
  # The metrics file carries its own BIBI_score / BIBI_Class columns. Drop them
  # and take B-IBI only from IN$bibi, so there is exactly one B-IBI source in
  # the pipeline and the provenance gate below checks every value.
  select(-any_of(c("BIBI_score", "BIBI_Class", "BIBI", "Region")))

bibi_raw <- read_csv(IN$bibi, show_col_types = FALSE) %>%
  mutate(SiteID = std_site(SiteID)) %>%
  select(SiteID, BIBI, BIBI_Class)

den_mat <- readRDS(IN$den_mat)
bio_mat <- readRDS(IN$bio_mat)

# ---- 2. Provenance gate (installed 25 Sep 2026): every site's B-IBI equals its ICPRB release value
# Replaces the magnitude guard (BIBI < 1), which halted on NORW's genuine 0.3487. A value is checked against its
# source, never against a size it is expected to have. Median basis: median_score; default basis: chessie_bibi.
rel_col <- if (nzchar(Sys.getenv("CH2_BIBI_MEDIAN"))) "median_score" else "chessie_bibi"
icprb <- bind_rows(
  read_csv(file.path(DATA, "Macros_Typ-1_ChessieB-IBI_Scores_Metrics.csv"), show_col_types = FALSE) %>%
    transmute(SiteID = std_site(site_ID), BIBI_release = as.numeric(.data[[rel_col]])),
  read_csv(file.path(DATA, "Macros_Typ-2_ChessieB-IBI_Scores_Metrics.csv"), show_col_types = FALSE) %>%
    transmute(SiteID = std_site(site_id), BIBI_release = as.numeric(.data[[rel_col]])),
  read_csv(file.path(DATA, "Typ3_ChessieB-IBI_Scores_Metrics.csv"), show_col_types = FALSE) %>%
    transmute(SiteID = std_site(SITE_ID), BIBI_release = as.numeric(.data[[rel_col]])))
gate <- bibi_raw %>% left_join(icprb, by = "SiteID") %>%
  mutate(ok = is.finite(BIBI_release) & abs(BIBI - BIBI_release) < 1e-6)
if (!all(gate$ok)) {
  bad <- gate %>% filter(!ok)
  stop("B-IBI PROVENANCE GATE: ", nrow(bad), " site(s) differ from the ICPRB ", rel_col, ": ",
       paste(sprintf("%s %.6g (release %.6g)", bad$SiteID, bad$BIBI, bad$BIBI_release), collapse = "; "))
}
say("  B-IBI provenance gate passed: all %d sites equal the ICPRB %s (tolerance 1e-6)", nrow(gate), rel_col)
say("  B-IBI range %.1f-%.1f, n=%d",
    min(bibi_raw$BIBI), max(bibi_raw$BIBI), nrow(bibi_raw))

# ---- 3. Assemble ------------------------------------------------------------
dat <- metrics_raw %>%
  left_join(bibi_raw, by = "SiteID") %>%
  mutate(Region = as_region(Typology)) %>%
  select(SiteID, Region, BIBI, BIBI_Class,
         any_of(c(unlist(METRIC_FAMILIES, use.names = FALSE))))

# Structural guards -- cheap, and they catch join failures immediately.
stopifnot(nrow(dat) == 90)
stopifnot(!any(is.na(dat$BIBI)))
stopifnot(all(table(dat$Region) == 30))
stopifnot(identical(sort(rownames(den_mat)), sort(rownames(bio_mat))))
stopifnot(all(dat$SiteID %in% rownames(den_mat)))

missing_metrics <- setdiff(unlist(METRIC_FAMILIES, use.names = FALSE), names(dat))
if (length(missing_metrics)) stop("Missing metric columns: ",
                                  paste(missing_metrics, collapse = ", "))

say("  loaded %d sites x %d metrics; matrices %d x %d",
    nrow(dat), length(unlist(METRIC_FAMILIES)), nrow(den_mat), ncol(den_mat))

# ---- 4. Regression test against reference values -----------------------------
# These 22 rho values are from the earlier analysis on the observed index, after
# Decapoda were removed. Reproducing them on that basis proves the inputs are correct.
PUBLISHED <- tribble(
  ~Metric,           ~Region,          ~rho_published,
  "Richness_alpha",  "Piedmont",        0.297,
  "Richness_alpha",  "Coastal Plain",  -0.065,
  "Richness_alpha",  "Valley & Ridge",  0.298,
  "Shannon",         "Piedmont",        0.606,
  "Shannon",         "Coastal Plain",   0.402,
  "HBI_bio",         "Piedmont",       -0.473,
  "Total_Biomass",   "Piedmont",        0.443,
  "Total_Density",   "Coastal Plain",  -0.335,
  "FDiv_den",        "Coastal Plain",  -0.694,
  "FDiv_den",        "Piedmont",       -0.494,
  "FDiv_bio",        "Coastal Plain",  -0.208,
  "FDis_den",        "Coastal Plain",   0.543,
  "FDis_den",        "Piedmont",        0.432,
  "RaoQ_den",        "Coastal Plain",   0.531,
  "RaoQ_den",        "Piedmont",        0.433,
  "FR_RD_den",       "Coastal Plain",  -0.596,
  "FR_RD_den",       "Piedmont",       -0.271,
  "FR_RD_bio",       "Piedmont",        0.426,
  "FR_RD_bio",       "Coastal Plain",  -0.022,
  "FR_DB_RD",        "Coastal Plain",   0.563,
  "FR_DB_RD",        "Piedmont",        0.345,
  "FRic_den",        "Piedmont",        0.307
)

verify <- PUBLISHED %>%
  rowwise() %>%
  mutate(
    rho_computed = {
      s <- dat[dat$Region == Region, ]
      suppressWarnings(cor(s$BIBI, s[[Metric]], method = "spearman",
                           use = "complete.obs"))
    }
  ) %>%
  ungroup() %>%
  mutate(diff = rho_computed - rho_published,
         status = ifelse(abs(diff) < 0.0015, "OK", "MISMATCH"))

write_csv(verify, file.path(OUT_TAB, "01_verification_report.csv"))

n_bad <- sum(verify$status == "MISMATCH")
# On the median-score run the published values are EXPECTED not to reproduce:
# the axis itself changed. Report the divergence, do not halt. The guard still
# halts on the default (chessie) run, where a mismatch means a real defect.
if (n_bad > 0 && nzchar(Sys.getenv("CH2_BIBI_MEDIAN"))) {
  say("  RUNNING ON MEDIAN B-IBI: %d of %d reference values differ, as expected.", n_bad, nrow(verify))
  say("  The reference values were computed on chessie_bibi. Guard not enforced.")
} else if (n_bad > 0) {
  print(as.data.frame(verify %>% filter(status == "MISMATCH")))
  stop("STEP 01 FAILED: ", n_bad, " of ", nrow(verify),
       " reference values did not reproduce.\n",
       "  Do NOT proceed. An input file is stale or wrong.\n",
       "  Report written to output/tables/01_verification_report.csv")
}

if (n_bad == 0) say("  VERIFIED: %d of %d reference values reproduced (|diff| < 0.0015)", nrow(verify), nrow(verify))
say("STEP 01 complete")
