# =============================================================================
# 18b_supplement_renumber.R
# =============================================================================
# WHAT THIS STEP DOES
#   Step 18 writes every table and figure under its legacy number (Table_S5, Table_S10, ... Figure_S3). On
#   10 October 2026 the supplement was consolidated into Tables S1 to S8 and Figure S1, and the fourth-corner
#   (trait-state) assets left the paper for a companion manuscript. This step applies R/supplement_asset_map.csv to
#   the files step 18 wrote, so the pipeline and the supplement agree on the numbers and nothing is renamed by hand.
#
#   Step 18 is left as it was: it still writes the legacy names, and the map is the one place the old and new
#   numbers meet.
#
# INPUTS
#   R/supplement_asset_map.csv; output/manuscript/tables/, figures/ and ASSET_MANIFEST.csv (step 18).
#
# OUTPUT (all under output/manuscript/)
#   supplement_final/            Table_S1A.csv ... Table_S8B.csv, Figure_S1.{png,pdf,svg}, SUPPLEMENT_MANIFEST.csv
#   trait_state_manuscript/      the assets cut from this paper (Table_S16, S17, S19, S20, S28, Table_4, Figure_S2, Figure_6)
#   tables_final/Table_2.csv     main-text Table 2: Tables 3a and 3b merged, the four retired metrics removed
# A table with no pipeline source (Table S1 part B and Table S2 part A, authored in the supplement) is listed in the
# manifest and not written; Table S2 part B (terms) is authored too, and the pipeline file of its "B and C" row is
# part C.
# =============================================================================
say("STEP 18b -- final supplement numbering from supplement_asset_map.csv")
MS      <- file.path(OUT, "manuscript")
MS_TAB  <- file.path(MS, "tables"); MS_FIG <- file.path(MS, "figures")
FINAL   <- file.path(MS, "supplement_final"); SPIN <- file.path(MS, "trait_state_manuscript"); TFIN <- file.path(MS, "tables_final")
for (d in c(FINAL, SPIN, TFIN)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(FINAL, full.names = TRUE))   # a renumbering must not leave files under the previous numbers
MAP <- readr::read_csv(file.path(PIPE_DIR, "R", "supplement_asset_map.csv"), show_col_types = FALSE, col_types = readr::cols(.default = "c"))
# legacy (step 18's manifest) is read but not used below.
legacy <- readr::read_csv(file.path(MS, "ASSET_MANIFEST.csv"), show_col_types = FALSE)

keep <- MAP[MAP$action == "keep", ]
nums <- sort(unique(as.integer(sub("^S", "", keep$new[grepl("^Table", keep$old)]))))
stopifnot("table numbers are not consecutive from S1" = identical(nums, seq_len(max(nums))))
stopifnot("a part letter is used twice inside one table" = !anyDuplicated(keep[grepl("^Table", keep$old), c("new", "part")]))

rows <- list()
for (i in seq_len(nrow(MAP))) {
  m <- MAP[i, ]; stem <- gsub(" ", "_", m$old); is_tab <- grepl("^Table", m$old)
  if (m$action == "delete") {
    ext <- if (is_tab) "csv" else c("png", "pdf", "svg"); dir_ <- if (is_tab) MS_TAB else MS_FIG
    for (e in ext) { f <- file.path(dir_, paste0(stem, ".", e)); if (file.exists(f)) file.copy(f, file.path(SPIN, basename(f)), overwrite = TRUE) }
    rows[[i]] <- tibble::tibble(Final = NA_character_, Old = m$old, Part = NA_character_, Action = "to the trait-state manuscript",
                                File = file.path("trait_state_manuscript", paste0(stem, if (is_tab) ".csv" else ".*")), Rows = NA_integer_, Cols = NA_integer_, Note = m$note)
    next
  }
  part <- if (is.na(m$part) || !nzchar(m$part)) "" else sub(".* ", "", m$part)       # "B and C": the pipeline file is the last part
  final <- paste0(if (is_tab) "Table_" else "Figure_", m$new, part)
  if (is_tab) {
    src <- file.path(MS_TAB, paste0(stem, ".csv"))
    if (file.exists(src)) { df <- readr::read_csv(src, show_col_types = FALSE); readr::write_csv(df, file.path(FINAL, paste0(final, ".csv")))
      rows[[i]] <- tibble::tibble(Final = final, Old = m$old, Part = m$part, Action = "renumbered", File = file.path("supplement_final", paste0(final, ".csv")), Rows = nrow(df), Cols = ncol(df), Note = m$note) }
    else rows[[i]] <- tibble::tibble(Final = final, Old = m$old, Part = m$part, Action = "authored in the supplement; no pipeline source", File = NA_character_, Rows = NA_integer_, Cols = NA_integer_, Note = m$note)
  } else {
    for (e in c("png", "pdf", "svg")) file.copy(file.path(MS_FIG, paste0(stem, ".", e)), file.path(FINAL, paste0(final, ".", e)), overwrite = TRUE)
    rows[[i]] <- tibble::tibble(Final = final, Old = m$old, Part = m$part, Action = "renumbered", File = file.path("supplement_final", paste0(final, ".*")), Rows = NA_integer_, Cols = NA_integer_, Note = m$note)
  }
}
man <- dplyr::bind_rows(rows)
readr::write_csv(man, file.path(FINAL, "SUPPLEMENT_MANIFEST.csv"))

# Main-text Table 2: Tables 3a and 3b merged, the four retired metrics removed (functional evenness and divergence met the criterion
# nowhere; Rao's quadratic entropy and space-filling redundancy repeat functional dispersion and richness; Table S4 parts A and C).
t3 <- dplyr::bind_rows(readr::read_csv(file.path(MS_TAB, "Table_3a.csv"), show_col_types = FALSE),
                       readr::read_csv(file.path(MS_TAB, "Table_3b.csv"), show_col_types = FALSE))
t3 <- t3[!grepl("Rao|Space-filling|FR_SF|evenness|divergence \\(FDiv\\)|FEve|FDiv", t3$Metric) | grepl("Redundancy divergence|Concentration divergence", t3$Metric), ]
t3 <- t3[order(match(t3$Region, c("Valley & Ridge", "Piedmont", "Coastal Plain")), -abs(t3$`Spearman rho`)), ]
readr::write_csv(t3, file.path(TFIN, "Table_2.csv"))
stopifnot("merged Table 2 should have 29 rows (the manuscript's)" = nrow(t3) == 29)

say("  supplement: %d tables (S1 to S%d, parts A to C) and %d figure(s); %d assets to trait_state_manuscript/; Table 2 has %d rows",
    length(nums), max(nums), sum(grepl("^Figure", keep$old)), sum(MAP$action == "delete"), nrow(t3))
say("STEP 18b complete")
