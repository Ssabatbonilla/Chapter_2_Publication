# =============================================================================
# 00d_site_replicate_medians.R -- site totals and site HBI as the median over the processed replicates (added 5 Oct 2026)
# =============================================================================
# Rule adopted 5 Oct 2026: no occurrence weighting. Every site x taxon density and biomass is the mean over the
# processed replicates, zeros included ("M0"); it feeds the community table and the matrices (00c, 00a, 00b) and through
# them Shannon, Bray-Curtis, Gini, N50 and the fourth corner. The SITE TOTALS (Total_Density, Total_Biomass) and the
# site HBI (HBI_rel_abund, HBI_den, HBI_bio) are the median over the processed replicates of each replicate's value
# ("M3"), read from Step2_SiteReplicateMedians.csv in the data folder (written upstream).
# Feeds the total density, total biomass and HBI rows of Table 2 and Table S4 part A.
#
# This step runs after 00b. It puts the median values into the default columns of both site tables
# (H1_site_summary_macro_metrics.csv and _fuzzyRD.csv) and keeps the mean versions beside them
# (Total_Density_mean, Total_Biomass_mean, HBI_rel_abund_mean, HBI_den_mean, HBI_bio_mean). The mean totals are
# checked against the replicate means in the same file (00c's parent resolution conserves site sums). Idempotent: a
# second run leaves the tables as the first left them. Without the median file this step changes nothing.
#
# Order: 00c, 00a, 00b, 00d, run_all.R (./run.sh). The upstream steps that built the data files are not in this repository.
# =============================================================================
suppressPackageStartupMessages({ library(readr); library(dplyr) })
PIPE_DIR <- if (exists("PIPE_DIR")) PIPE_DIR else normalizePath(".")
source(file.path(PIPE_DIR, "R", "00_paths.R"))
MED <- file.path(DATA, "Step2_SiteReplicateMedians.csv")
if (!file.exists(MED)) { cat("[00d] no Step2_SiteReplicateMedians.csv: Step 2 did not run the default method; site tables unchanged\n"); quit(status = 0) }
med <- read_csv(MED, show_col_types = FALSE) %>% mutate(SiteID = toupper(trimws(SiteID)))
stopifnot(all(med$density_method == "M3M0"), !anyDuplicated(med$SiteID))
cat("=== Site totals and site HBI: MEDIAN over the processed replicates (M3); taxa and matrices: MEAN (M0); means kept beside the medians ===\n")
SWAP <- c(Total_Density = "Total_Density_median", Total_Biomass = "Total_Biomass_median", HBI_rel_abund = "HBI_rel_abund_median",
          HBI_den = "HBI_den_median", HBI_bio = "HBI_bio_median")
for (f in file.path(H1, c("H1_site_summary_macro_metrics.csv", "H1_site_summary_macro_metrics_fuzzyRD.csv"))) {
  st <- read_csv(f, show_col_types = FALSE); key <- toupper(trimws(st$SiteID))
  stopifnot(all(key %in% med$SiteID)); m <- med[match(key, med$SiteID), ]
  for (col in names(SWAP)) {
    mc <- paste0(col, "_mean")
    if (!mc %in% names(st)) {                                   # first pass: the M0 value 00c computed becomes the mean column
      st[[mc]] <- st[[col]]
      st <- st %>% relocate(all_of(mc), .after = all_of(col))
    }
    st[[col]] <- m[[SWAP[[col]]]]
  }
  rel <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))
  d_den <- rel(st$Total_Density_mean, m$Total_Density_mean); d_bio <- rel(st$Total_Biomass_mean, m$Total_Biomass_mean)
  if (d_den > 1e-8 || d_bio > 1e-8) stop("[00d] the site table's mean totals differ from the replicate means in Step2_SiteReplicateMedians.csv (", signif(d_den, 3), ", ", signif(d_bio, 3), "): stale table?")
  write_csv(st, f)
  cat(sprintf("[00d] %s: Total_Density, Total_Biomass and HBI_* now replicate medians; *_mean beside them (mean totals match Step 2 to %.1e)\n",
              basename(f), max(d_den, d_bio)))
}
