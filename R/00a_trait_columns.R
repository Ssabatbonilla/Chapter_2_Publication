# =============================================================================
# 00a_trait_columns.R  (added 22 September 2026; run before 00b by ./run.sh, not part of run_all.R)
# =============================================================================
# Rebuilds the trait-derived columns of the site table from the Twardochleb fuzzy file, under the
# trait rules adopted on 22 September 2026: every trait is fuzzy-coded from
# Twardochleb et al. (2021) or melded from related Twardochleb genera; a trait group with no record is
# filled from the taxon's family in Twardochleb, and left NA, never zero, where the family has none.
#
#   Functional diversity. FD::dbFD exactly as step 13 reproduces the published values (the fuzzy and
#   binary proxy matrices, the either-coding column rule, trait_blocks$gowdis_w, stand.x, w.abun, m = 10),
#   with the Cailliez correction for every metric. dbFD computes RaoQ and FDis on the uncorrected Gower
#   distance, so they are identical under either correction (step 13, section 4); the sqrt correction
#   fails once empty groups are NA. Columns: RaoQ_, FDis_, FDiv_, FEve_ (den, bio), FRic_den_fuzPool
#   and FRic_den; RaoQ_gower_ in the binary geometry, which 00b checks and then replaces.
#   Community-weighted means. For each modality, sum_i w_i a_i / sum_i w_i over the taxa with a record
#   for that trait group, per site, weighted by density (_CWM_den) and by biomass (_CWM_bio), named by
#   the fuzzy columns (for example Feed_SH_CWM_den). The legacy-named CWM and REL columns, built
#   upstream from the legacy binary file (a taxon with several binary states counted in full in each, so
#   feeding, habit and respiration shares exceeded 1), are dropped.
#   Dropped as stale and unread by the pipeline: the binPool, binReg and fuzReg variants and every
#   fuzPool variant except FRic_den_fuzPool, which 00b reads.
#   Consistency check. With CH2_00A_VERIFY=1 every recomputed functional diversity column must
#   reproduce the site table this step reads (since 00c that table is 00c's output, computed by the same
#   dbFD call, so this is a self-consistency check; 1e-8, binary RaoQ_gower 1e-12), and nothing is
#   written but the provenance table.
#
# Reads H1_site_summary_macro_metrics.csv (00c's output), den_mat_fd.rds, bio_mat_fd.rds, the fuzzy and binary
# proxy trait files and trait_blocks.rds. Writes H1_site_summary_macro_metrics_twardochleb.csv and
# H1_trait_columns_provenance.csv beside the site table, which is never modified; 00b reads the first.
# Feeds every functional diversity metric screened by step 03 (Table 2, Table S4 part A) through 00b.
# Run:  Rscript R/00a_trait_columns.R   then 00b, 00d and run_all.R (./run.sh does all of this).
# =============================================================================
suppressPackageStartupMessages({ library(FD); library(dplyr); library(readr); library(tibble) })
args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
PIPE_DIR <- normalizePath(file.path(dirname(gsub("~\\+~", " ", f)), ".."))
source(file.path(PIPE_DIR, "R", "00_paths.R"))   # CH2_ROOT, DATA, RA, H1, REPO_LAYOUT (one definition for every script)
src <- file.path(H1, "H1_site_summary_macro_metrics.csv"); out <- file.path(H1, "H1_site_summary_macro_metrics_twardochleb.csv")
verify <- nzchar(Sys.getenv("CH2_00A_VERIFY"))
met <- read_csv(src, show_col_types = FALSE)
dm <- readRDS(file.path(DATA, "den_mat_fd.rds")); bm <- readRDS(file.path(DATA, "bio_mat_fd.rds"))
std <- function(x) { x <- as.data.frame(x); if (!"TAXANAME" %in% names(x)) x$TAXANAME <- x$Genus; x[!duplicated(x$TAXANAME), ] }
tf <- std(readRDS(file.path(DATA, "traits_fuzzy_proxy_main.rds"))); tb <- std(readRDS(file.path(DATA, "traits_binary_proxy_main.rds")))
tblk <- readRDS(file.path(DATA, "trait_blocks.rds"))
cl <- intersect(setdiff(names(tb), c("Genus", "TAXANAME")), setdiff(names(tf), c("Genus", "TAXANAME")))
cl <- cl[sapply(tb[cl], is.numeric) & sapply(tf[cl], is.numeric)]
tx <- intersect(union(colnames(dm), colnames(bm)), intersect(tb$TAXANAME, tf$TAXANAME))
build <- function(df) { m <- as.matrix(df[match(tx, df$TAXANAME), cl]); rownames(m) <- tx; storage.mode(m) <- "double"; m }  # NA kept
MF <- build(tf); MB <- build(tb)
keep <- (colSums(MB, na.rm = TRUE) > 0) | (colSums(MF, na.rm = TRUE) > 0)
MFk <- MF[, keep, drop = FALSE]; gw <- tblk$gowdis_w[colnames(MFk)]; stopifnot(all(is.finite(gw)))
sites <- met$SiteID
cat(sprintf("trait matrix: %d taxa x %d modalities (%d retained); NA cells %d on %d taxa\n",
            nrow(MF), ncol(MF), ncol(MFk), sum(is.na(MF)), sum(rowSums(is.na(MF)) > 0)))

# ---- functional diversity, as step 13 reproduces it, Cailliez throughout ---------------------------
run_fd <- function(a) {
  a <- a[, intersect(colnames(a), rownames(MFk)), drop = FALSE]
  suppressWarnings(FD::dbFD(x = MFk, a = a, w = gw, stand.x = TRUE, corr = "cailliez", w.abun = TRUE,
                            calc.FRic = TRUE, calc.FDiv = TRUE, stand.FRic = TRUE, calc.CWM = FALSE, m = 10, messages = FALSE))
}
fd <- list(den = run_fd(dm), bio = run_fd(bm))
pick <- function(o, nm) o[[nm]][match(sites, names(o[[nm]]))]
FDV <- list(RaoQ_den = pick(fd$den, "RaoQ"), FDis_den = pick(fd$den, "FDis"), FDiv_den = pick(fd$den, "FDiv"), FEve_den = pick(fd$den, "FEve"),
            RaoQ_bio = pick(fd$bio, "RaoQ"), FDis_bio = pick(fd$bio, "FDis"), FDiv_bio = pick(fd$bio, "FDiv"), FEve_bio = pick(fd$bio, "FEve"),
            FRic_den_fuzPool = pick(fd$den, "FRic"))
# binary geometry Q exactly as 00b computes it (unweighted Gower on the binary proxy matrix, all columns)
G_bin <- as.matrix(FD::gowdis(as.data.frame(MB)))
P <- function(M) { M <- M[sites, tx]; M / rowSums(M) }
q_of <- function(Pm, G) rowSums((Pm %*% G) * Pm)
FDV$RaoQ_gower_den <- unname(q_of(P(dm), G_bin)); FDV$RaoQ_gower_bio <- unname(q_of(P(bm), G_bin))

diffs <- sapply(names(FDV), function(k) max(abs(FDV[[k]] - met[[k]]), na.rm = TRUE))
tol <- ifelse(grepl("gower", names(FDV)), 1e-12, 1e-8)
prov <- tibble(column = names(FDV), max_abs_diff_from_upstream = signif(diffs, 3),
               spearman_with_upstream = sapply(names(FDV), function(k) round(cor(FDV[[k]], met[[k]], method = "spearman", use = "complete.obs"), 4)),
               n_nonfinite = sapply(FDV, function(v) sum(!is.finite(v))))
print(as.data.frame(prov), row.names = FALSE)
if (verify) {
  write_csv(prov %>% mutate(mode = "verify", tolerance = tol, reproduced = diffs < tol), file.path(H1, "H1_trait_columns_provenance.csv"))
  if (any(diffs >= tol)) stop("00a VERIFY FAILED: the recipe does not reproduce the upstream columns from these trait files; nothing written")
  cat("00a VERIFY: every functional diversity column reproduces the upstream table from these trait files\n"); quit(save = "no")
}
stopifnot(all(sapply(FDV, function(v) all(is.finite(v)))))

# ---- community-weighted means from the fuzzy affinities -------------------------------------------
cwm_of <- function(W) {
  W <- W[sites, tx, drop = FALSE]; W[is.na(W)] <- 0
  out <- sapply(cl, function(col) {
    a <- MF[, col]; rec <- !is.na(a); a[!rec] <- 0
    den <- as.vector(W %*% as.numeric(rec)); num <- as.vector(W %*% a)
    ifelse(den > 0, num / den, NA_real_)
  })
  rownames(out) <- sites; out
}
C_den <- cwm_of(dm); C_bio <- cwm_of(bm)
colnames(C_den) <- paste0(colnames(C_den), "_CWM_den"); colnames(C_bio) <- paste0(colnames(C_bio), "_CWM_bio")
grp <- sub("_.*", "", cl)
sums <- sapply(unique(grp), function(g) rowSums(C_den[, paste0(cl[grp == g], "_CWM_den"), drop = FALSE]))
stopifnot(all(abs(sums[is.finite(sums)] - 1) < 1e-9))                     # each group's shares sum to 1
sums_bio <- sapply(unique(grp), function(g) rowSums(C_bio[, paste0(cl[grp == g], "_CWM_bio"), drop = FALSE]))
stopifnot(all(abs(sums_bio[is.finite(sums_bio)] - 1) < 1e-9))             # and weighted by biomass (27 Sep 2026)

# ---- assemble --------------------------------------------------------------------------------------
drop_cols <- names(met)[grepl("_(binPool|binReg|fuzReg)$", names(met)) |
                        (grepl("_fuzPool$", names(met)) & names(met) != "FRic_den_fuzPool") |
                        grepl("_CWM_(den|bio)$", names(met)) |
                        grepl("^(CWM|REL)_(Response|Effect|Linking)", names(met))]
new <- met %>% select(-all_of(drop_cols))
for (k in names(FDV)) new[[k]] <- FDV[[k]]
new$FRic_den <- new$FRic_den_fuzPool
new <- bind_cols(new, as_tibble(C_den), as_tibble(C_bio))
write_csv(new, out)
prov <- bind_rows(prov %>% mutate(mode = "rebuild"),
                  tibble(column = c("dropped columns", "CWM columns written", "taxa with an NA trait group", "CWM group sums (all equal 1)"),
                         mode = "rebuild", note = c(paste(length(drop_cols), "stale or legacy-built"), as.character(ncol(C_den) + ncol(C_bio)),
                                                    as.character(sum(rowSums(is.na(MF)) > 0)), "checked")))
write_csv(prov, file.path(H1, "H1_trait_columns_provenance.csv"))
cat(sprintf("wrote %s: %d columns (%d dropped, %d CWM columns)\n", basename(out), ncol(new), length(drop_cols), ncol(C_den) + ncol(C_bio)))
