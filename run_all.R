#!/usr/bin/env Rscript
# =============================================================================
# run_all.R — Chapter 2 analysis pipeline v2
# =============================================================================
# RUN IT (the published analysis; both flags are required, see the banners below). ./run.sh does all of this.
#   cd <repository root>
#   Rscript R/00c_resolve_parent_rows.R; Rscript R/00a_trait_columns.R; Rscript R/00b_redundancy_geometry.R; Rscript R/00d_site_replicate_medians.R
#   CH2_BIBI_MEDIAN=1 CH2_RD_GEOMETRY=fuzzy Rscript run_all.R              # everything
#   CH2_BIBI_MEDIAN=1 CH2_RD_GEOMETRY=fuzzy Rscript run_all.R 01 02 03     # only these steps (and what they need)
#
# It is safe to run repeatedly. It reads the data files, writes only to ./output/,
# and never modifies the data folder (the four 00 steps above do, by design).
#
# Paper assets: output/manuscript/tables_final/ (Table 2) and output/manuscript/supplement_final/ (Tables S1 to S8,
# Figure S1), written by steps 18 and 18b.
#
# If step 01 fails, STOP. It means an input file is stale or wrong, and every
# number produced after it would be untrustworthy. Fix the input, do not skip.
# =============================================================================

.args     <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args, value = TRUE)
PIPE_DIR  <- if (length(.file_arg)) {
  # Rscript encodes spaces in --file= as "~+~". This project lives in a folder
  # tree whose names carry spaces, so the path is full of them.
  # Without this substitution the pipeline only runs from its own directory.
  .p <- gsub("~\\+~", " ", sub("^--file=", "", .file_arg[1]), perl = TRUE)
  normalizePath(dirname(.p))
} else normalizePath(".")

STEPS <- c(
  "01" = "01_load_and_verify.R",
  "02" = "02_concentration_metrics.R",
  "03" = "03_screen_metrics.R",
  "03b" = "03b_leverage_audit.R",       # leave-few-out leverage audit of the responsive associations, Table S3 part B (legacy Table S18) (24 Sep 2026)
  "04" = "04_multivariate_within_region.R",
  # "05" (land-use covariate models) removed 24 Sep 2026: no analysis in this paper models land use (step 18
  # reports the land-use gradient only as context in Table S1 part A); the script is not in this repository.
  "07" = "07_tables.R",
  # "08" (manuscript number audit) retired 27 Sep 2026: it compared the pipeline with coefficients typed from a
  # draft that no longer exists, and nothing read its two tables; the script is not in this repository.
  "09" = "09_taxon_responses_and_hidden_biomass.R",
  "09b" = "09b_taxon_tolerance.R",      # tolerance values of the responding taxa, Table S8 part B (legacy Table S24) (21 Sep 2026)
  "10" = "10_taxon_set_accounting.R",
  "11" = "11_taxon_convention_and_traitspace_quality.R",
  "12" = "12_evidence_sufficiency.R",
  "13" = "13_fd_provenance.R",
  "14" = "14_power_and_detectable_effects.R",
  "15" = "15_rlq_fourthcorner.R",
  "15b" = "15b_weighting_difference.R",  # density-against-biomass difference per trait state, pigeonhole bootstrap (24 Sep 2026); companion trait-state manuscript, read by step 06
  "06" = "06_figures.R",           # after 15 and 15b: step 06 draws the trait-state figures (legacy Figure 6 and Figure S2, companion manuscript) from their tables of THIS run
  "16" = "16_beta_turnover_nestedness.R",
  "17" = "17_null_model_ses.R",
  "19" = "19_study_area_map.R",   # before 18: 18 copies this map out as Figure 2
  "20" = "20_chironomid_sensitivity.R",   # 8 Oct 2026: the four chironomid groups removed, controls, genus-level split; Table S7 parts A to C (legacy Tables S26, S27, S29; S28 is trait-state), Figure S1 (legacy Figure S3)
  "20b" = "20b_figS1_from_record.R",   # 12 Oct 2026: Figure S1 redrawn from step 20's tables (share of the chironomids by site; functional dispersion and overlap redundancy with and without them); recomputes nothing
  "18" = "18_manuscript_assets.R",
  "18b" = "18b_supplement_renumber.R",  # 10 Oct 2026: final supplement numbering from R/supplement_asset_map.csv; reads the files step 18 wrote
  "21" = "21_session_record.R"          # 11 Oct 2026: R version, package versions and citations of this run (Appendix S1.8; output/manuscript/software/)
)

requested <- commandArgs(trailingOnly = TRUE)
if (!length(requested)) requested <- names(STEPS)
unknown <- setdiff(requested, names(STEPS))
if (length(unknown)) stop("Unknown step(s): ", paste(unknown, collapse = ", "),
                          "\nValid: ", paste(names(STEPS), collapse = " "))

# Dependency graph. Requesting a step pulls in everything it needs, so
# `Rscript run_all.R 07` does the right thing rather than failing obscurely.
DEPS <- list(
  "01" = character(0),
  "02" = "01",
  "03" = c("01","02"),
  "03b" = c("01","02","03"),
  "04" = "01",
  "06" = c("01","02","03","11","15","15b"),   # 11: the trait-space ordination behind fig4_traitspace_binned (drawn, not in the paper)
  "07" = c("01","02","03","04"),
  "09" = c("01","02","03"),
  "09b" = c("01","02","03","09"),
  "10" = c("01","02","03"),
  "11" = c("01","02","03"),
  "12" = c("01","02","03"),
  "13" = "01",
  "14" = c("01","02","03"),
  "15" = c("01","02","03"),
  "15b" = c("01","02","03","15"),
  "16" = c("01","02","03"),
  "17" = c("01","02","03"),
  "18" = c("01","02","03","03b","06","07","09","09b","10","11","12","13","14","15","15b","16","17","19","20","20b"),
  "19" = character(0),
  "20" = c("01","02","03","15"),
  "20b" = c("01","20"),
  "18b" = character(0)                           # reads step 18's output; run 18 first when assets change
)
expand_deps <- function(x) {
  repeat {
    y <- unique(c(x, unlist(DEPS[x], use.names = FALSE)))
    if (setequal(y, x)) break
    x <- y
  }
  x
}
requested <- expand_deps(requested)
requested <- names(STEPS)[names(STEPS) %in% requested]   # restore run order

# vegan's permutest/betadisper open a base graphics device, which leaves a
# stray Rplots.pdf next to the script. Send those to a null device instead.
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

source(file.path(PIPE_DIR, "R", "00_config.R"))

# Repository layout (a data/ folder beside R/): the two flags below are required, because without them the
# runner would write two retired analyses under the published file names. In the project tree the flags
# stay optional and the banners say which basis ran.
if ((REPO_LAYOUT || nzchar(Sys.getenv("CH2_DATA_DIR"))) && (!nzchar(Sys.getenv("CH2_BIBI_MEDIAN")) || !identical(Sys.getenv("CH2_RD_GEOMETRY"), "fuzzy")))
  stop("This copy runs the published analysis only. Run:\n",
       "  CH2_BIBI_MEDIAN=1 CH2_RD_GEOMETRY=fuzzy Rscript run_all.R ", paste(requested, collapse = " "),
       "\nor ./run.sh (see README.md).")

# ---------------------------------------------------------------------------
# WHICH B-IBI IS THIS RUN ON?  Read this line before trusting any output.
#
# The published Chapter 2 values are computed on the RAREFIED MEDIAN score, which
# is selected only by the environment variable CH2_BIBI_MEDIAN=1. Running without
# it silently screens on the observed chessie_bibi score instead and overwrites
# output/tables/ with a DIFFERENT ANALYSIS -- Valley & Ridge flips from the most
# responsive region to a null, and nothing errors.
#
# That happened on 8 September 2026. Always run:
#     CH2_BIBI_MEDIAN=1 Rscript run_all.R <steps>
# and check the banner below before using anything the run produced.
# ---------------------------------------------------------------------------
BIBI_BASIS <- if (nzchar(Sys.getenv("CH2_BIBI_MEDIAN")))
  "RAREFIED MEDIAN  (matches the manuscript)" else
  "observed chessie_bibi  <-- NOT the published basis"
say("=== Chapter 2 pipeline v2 | steps: %s ===", paste(requested, collapse = " "))
say("=== B-IBI basis: %s ===", BIBI_BASIS)
if (!nzchar(Sys.getenv("CH2_BIBI_MEDIAN"))) {
  say("!!! Output will NOT match the manuscript.")
  say("!!! Re-run as: CH2_BIBI_MEDIAN=1 Rscript run_all.R %s", paste(requested, collapse = " "))
  say("!!! unless you deliberately want the chessie baseline.")
}
# Second basis, added 18 September 2026 after a full rerun with CH2_BIBI_MEDIAN=1 alone
# silently rebuilt the BINARY redundancy geometry over the fuzzy one the manuscript now
# carries (overlap redundancy, redundancy divergence, FRic, FR_SF; see R/00b_redundancy_geometry.R).
# The fuzzy geometry is selected only by CH2_RD_GEOMETRY=fuzzy, and 00b must have been run first.
RD_GEOMETRY <- if (identical(Sys.getenv("CH2_RD_GEOMETRY"), "fuzzy"))
  "FUZZY group-weighted Gower  (matches the manuscript since 18 September 2026)" else
  "BINARY unweighted Gower  <-- NOT the manuscript basis"
say("=== Redundancy geometry: %s ===", RD_GEOMETRY)
if (!identical(Sys.getenv("CH2_RD_GEOMETRY"), "fuzzy")) {
  say("!!! Redundancy and geometry outputs will NOT match the manuscript.")
  say("!!! Re-run as: Rscript R/00b_redundancy_geometry.R; CH2_BIBI_MEDIAN=1 CH2_RD_GEOMETRY=fuzzy Rscript run_all.R %s", paste(requested, collapse = " "))
}
# Third basis, added 21 September 2026. Step 15 has a weighting switch, CH2_FC_CURRENCY. Its
# default used to be density only, and the documented command never wrote the biomass files, so
# the biomass fourth-corner outputs (the biomass rows of legacy Table 4, now in the companion trait-state manuscript) sat on
# the observed chessie B-IBI from 2 to 21 September. The runner now runs step 15 once per
# weighting unless CH2_FC_CURRENCY is set, which restricts it to that one weighting.
FC_CURRENCY_SET <- nzchar(Sys.getenv("CH2_FC_CURRENCY"))
if ("15" %in% requested) {
  say("=== Step 15 weightings: %s ===",
      if (FC_CURRENCY_SET) sprintf("%s ONLY  <-- the other weighting's step-15 files are NOT rewritten", Sys.getenv("CH2_FC_CURRENCY"))
      else "density AND biomass  (both written)")
  if (FC_CURRENCY_SET)
    say("!!! Unset CH2_FC_CURRENCY to write both weightings in one run.")
}
PR <- file.path(RA, "tables", "H1", "H1_parent_rows_provenance.csv")
if (file.exists(PR)) {
  pr <- read.csv(PR, stringsAsFactors = FALSE); n_res <- pr$max_abs_diff[pr$column == "parent-site cases resolved"]
  say("=== Parent rows: resolved by DPAC before any metric (R/00c; %s parent-site cases)  (matches the manuscript since 24 September 2026) ===", n_res)
} else {
  say("=== Parent rows: NOT resolved  <-- R/00c_resolve_parent_rows.R has not run; richness and every metric keep parent rows beside their children ===")
}
t0 <- Sys.time()

for (s in requested) {
  f <- file.path(PIPE_DIR, "R", STEPS[[s]])
  if (!file.exists(f)) { say("  step %s not present yet (%s) -- skipped", s, STEPS[[s]]); next }
  if (s == "15" && !FC_CURRENCY_SET) {
    # Biomass first, density last, so the session ends in the state a density-only run always left.
    for (cur in c("biomass", "density")) {
      Sys.setenv(CH2_FC_CURRENCY = cur)
      source(f)
    }
    Sys.unsetenv("CH2_FC_CURRENCY")
  } else source(f)
}

say("=== done in %.1f s | outputs in %s ===",
    as.numeric(difftime(Sys.time(), t0, units = "secs")), if (REPO_LAYOUT) "output/" else OUT)
