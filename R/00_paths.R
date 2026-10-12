# =============================================================================
# 00_paths.R -- where the data live; sourced by 00_config.R and by 00a, 00b and 00c
# =============================================================================
# Two layouts, one rule. In the project tree the pipeline folder sits beside Data/Data
# and reads from there. In the public repository a data/ folder sits beside R/ and
# mirrors the layout of Data/Data, so every path below it is the same and the scripts
# run unchanged. Whichever exists is used. CH2_DATA_DIR overrides both (a scratch copy
# of the data in the same layout). The pre-pipeline steps 00c, 00a, 00b and 00d write their
# derived files into the same folder they read from, beside the record copies.
stopifnot(exists("PIPE_DIR"))
CH2_ROOT <- dirname(PIPE_DIR)
REPO_LAYOUT <- !nzchar(Sys.getenv("CH2_DATA_DIR")) && dir.exists(file.path(PIPE_DIR, "data"))
DATA <- if (nzchar(Sys.getenv("CH2_DATA_DIR"))) normalizePath(Sys.getenv("CH2_DATA_DIR"), mustWork = TRUE) else
        if (REPO_LAYOUT) file.path(PIPE_DIR, "data") else file.path(CH2_ROOT, "Data", "Data")
RA <- file.path(DATA, "Regional Analysis")
H1 <- file.path(RA, "tables", "H1")
