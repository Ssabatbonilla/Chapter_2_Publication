# =============================================================================
# 21_session_record.R  --  the R version and the packages of THIS run, with their citations
# =============================================================================
# WHY THIS STEP EXISTS (11 October 2026)
#   The manuscript says which R version and which packages computed each number (Appendix S1.8), and a
#   journal asks for them to be cited. This step runs last, in the same R session as every other step, so what it
#   writes is what produced the tables and figures that the same run wrote:
#     output/logs/session_record_<stamp>.txt       R version, platform, locale, flags, sessionInfo()
#     output/tables/R_packages_and_versions.csv    package, version, role, whether this run loaded it
#     output/manuscript/software/R_software_citations.bib       BibTeX from citation()
#     output/manuscript/software/R_software_references_APA.txt   the same, as text, for the reference list
# =============================================================================

say("STEP 21 -- session record (R version, packages, citations)")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(file.path(OUT, "manuscript", "software"), showWarnings = FALSE, recursive = TRUE)

# role of each package in this chapter's pipeline (what a reader needs to know to repeat it)
ROLES <- tibble::tribble(
  ~package,      ~role,
  "FD",          "functional richness, evenness, divergence, dispersion and Rao's quadratic entropy (dbFD); step 03 onward",
  "ade4",        "Gower dissimilarity, Cailliez correction, principal coordinates analysis; RLQ and fourth-corner in step 15 (the held trait-state manuscript)",
  "vegan",       "PERMANOVA (adonis2), multivariate dispersion (betadisper), Bray-Curtis dissimilarity; step 04",
  "betapart",    "beta diversity partitioning (turnover and nestedness); step 16",
  "dplyr",       "data manipulation",
  "tidyr",       "data reshaping",
  "tibble",      "data frames",
  "readr",       "reading and writing tables",
  "stringr",     "string handling",
  "purrr",       "iteration over regions, metrics and permutations",
  "rlang",       "tidy evaluation used by the tidyverse packages",
  "ggplot2",     "figures",
  "patchwork",   "multi-panel figure layout",
  "cowplot",     "figure assembly",
  "scales",      "axis scales and labels",
  "ggspatial",   "scale bar and north arrow on the study-area map (Figure 2)",
  "sf",          "spatial data for the study-area map (Figure 2)",
  "maps",        "state boundaries for the study-area map (Figure 2)",
  "svglite",     "SVG graphics device",
  "ragg",        "PNG graphics device",
  "systemfonts", "font handling for the graphics devices",
  "stats",       "Spearman correlations, Wilcoxon tests, Benjamini-Hochberg adjustment, ordination helpers",
  "parallel",    "permutation loops",
  "grid",        "figure layout",
  "grDevices",   "graphics devices and colors",
  "utils",       "reading, writing and citation helpers"
)
base_pkgs <- c("stats", "parallel", "grid", "grDevices", "utils", "graphics", "methods", "datasets", "tools")
ip <- installed.packages()
ROLES$installed <- ROLES$package %in% rownames(ip)
ROLES$version   <- vapply(ROLES$package, function(p) if (p %in% rownames(ip)) as.character(packageVersion(p)) else NA_character_, "")
ROLES$base_R    <- ROLES$package %in% base_pkgs
ROLES$loaded_in_this_run <- ROLES$package %in% loadedNamespaces()
miss <- ROLES$package[!ROLES$installed]
if (length(miss)) say("  !!! package(s) named in the roles table but not installed: %s", paste(miss, collapse = ", "))
# loaded namespaces this run used that the roles table does not name (a package added to a script and not to this record)
extra <- setdiff(loadedNamespaces(), c(ROLES$package, base_pkgs, "compiler"))
direct <- c(ROLES$package[ROLES$loaded_in_this_run])
say("  R %s | %d packages named; %d loaded this run; namespaces loaded but not named here (dependencies): %d",
    paste(R.version$major, R.version$minor, sep = "."), nrow(ROLES), sum(ROLES$loaded_in_this_run), length(extra))
ROLES$citation_key <- paste0("pkg_", ROLES$package)
write.csv(ROLES, file.path(OUT_TAB, "R_packages_and_versions.csv"), row.names = FALSE)

# the session record
flags <- Sys.getenv(c("CH2_BIBI_MEDIAN", "CH2_RD_GEOMETRY", "CH2_FC_CURRENCY", "CH2_DENSITY_METHOD", "CH2_MISSING_REPS", "CH2_00C_STAMP"), unset = "(unset)")
py <- tryCatch(system2("python3", c("-c", shQuote("import sys;print(sys.version.split()[0])")), stdout = TRUE, stderr = FALSE), error = function(e) "not available")
sink(file.path(OUT_LOG, sprintf("session_record_%s.txt", stamp)))
cat("Chapter 2 pipeline v2: session record\n")
cat("written:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R:", R.version.string, "\n")
cat("platform:", R.version$platform, "| running under:", sessionInfo()$running, "\n")
cat("locale:", Sys.getlocale(), "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
cat("Rscript:", file.path(R.home("bin"), "Rscript"), "\n")
cat("Python (Figure 1 build script only):", paste(py, collapse = " "), "\n")
cat("analysis flags:\n"); print(flags)
cat("\nsteps requested this run:", paste(requested, collapse = " "), "\n\n")
print(sessionInfo())
sink()

# citations: R itself, then each installed, non-base package the roles table names (whether or not this run loaded it)
cits <- list(R = citation())
for (p in ROLES$package[ROLES$installed & !ROLES$base_R]) cits[[p]] <- tryCatch(citation(p), error = function(e) NULL)
bib <- unlist(lapply(names(cits), function(n) if (is.null(cits[[n]])) NULL else {
  ver <- if (n == "R") paste(R.version$major, R.version$minor, sep = ".") else ROLES$version[ROLES$package == n]
  ce <- cits[[n]]
  unlist(lapply(seq_along(ce), function(i) {
    # a unique key per entry (a package can carry several: FD two, ade4 five, sf two), the version in a comment so a reader sees what was run
    key <- if (n == "R") "RCoreTeam" else sprintf("pkg_%s_%d", n, i)
    b <- as.character(toBibtex(ce[i]))
    b <- sub("^@([A-Za-z]+)\\{[^,]*,", sprintf("@\\1{%s,", key), b)
    c(b, sprintf("%% %s version used: %s (citation entry %d of %d)", n, ver, i, length(ce)), "")
  }))
}))
writeLines(bib, file.path(OUT, "manuscript", "software", "R_software_citations.bib"))
txt <- unlist(lapply(names(cits), function(n) if (is.null(cits[[n]])) NULL else {
  ver <- if (n == "R") paste(R.version$major, R.version$minor, sep = ".") else ROLES$version[ROLES$package == n]
  c(sprintf("[%s %s]", n, ver), format(cits[[n]], style = "text"), "")
}))
writeLines(txt, file.path(OUT, "manuscript", "software", "R_software_references_APA.txt"))
say("  wrote session_record_%s.txt, R_packages_and_versions.csv, R_software_citations.bib and R_software_references_APA.txt", stamp)
say("STEP 21 complete")
