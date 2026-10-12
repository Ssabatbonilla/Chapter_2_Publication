# =============================================================================
# 00b_redundancy_geometry.R  (added 18 September 2026; run after 00a by ./run.sh, not part of run_all.R)
# =============================================================================
# Recomputes the redundancy family in the trait geometry Supplement S1.5 describes.
#
#   The input site table (00a's output; the binary geometry is computed in 00c and 00a) carries Rao's Q for the overlap redundancy index
#   1 - Q/D from FD::gowdis() on the BINARY trait matrix with no group weights. Every other functional
#   metric uses the fuzzy-coded matrix with the group weights in trait_blocks$gowdis_w (1 / group size).
#   This script computes Q on that fuzzy group-weighted Gower distance (uncorrected, so d is in [0, 1]
#   and Q <= D), rebuilds FR_RD_den, FR_RD_bio and FR_DB_RD, leaves every other column untouched, and
#   writes H1_site_summary_macro_metrics_fuzzyRD.csv beside the original, which step 01 reads under
#   CH2_RD_GEOMETRY=fuzzy (every redundancy row and FRic, FR_SF of Table 2 and Table S4). Since the later 18 September run it
#   also moves FRic_den to the fuzzy-pool value and recomputes FR_SF from it. The original file is not modified.
#
#   Provenance: the binary geometry reproduces the published RaoQ_gower to < 1e-12 and the fuzzy weighted
#   geometry reproduces the published dbFD RaoQ (d^2 / 2) to < 1e-12; both checks halt the script if they fail.
#   Output also includes H1_redundancy_geometry_provenance.csv with both checks and the site-level
#   agreement between the two versions of each index.
#
# Run:  Rscript R/00b_redundancy_geometry.R   then 00d, then   CH2_BIBI_MEDIAN=1 CH2_RD_GEOMETRY=fuzzy Rscript run_all.R
# =============================================================================
suppressPackageStartupMessages({ library(FD); library(dplyr); library(readr); library(tibble) })
args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
PIPE_DIR <- normalizePath(file.path(dirname(gsub("~\\+~", " ", f)), ".."))
source(file.path(PIPE_DIR, "R", "00_paths.R"))   # CH2_ROOT, DATA, RA, H1, REPO_LAYOUT (one definition for every script)
# Twardochleb fuzzy traits, 22 Sep 2026: read the table 00a rebuilt from the Twardochleb fuzzy file; stop if it is missing.
src <- file.path(H1, "H1_site_summary_macro_metrics_twardochleb.csv")
if (!file.exists(src)) stop("[00b] ", basename(src), " is missing: run R/00a_trait_columns.R first")
out <- file.path(H1, "H1_site_summary_macro_metrics_fuzzyRD.csv")
met <- read_csv(src, show_col_types = FALSE)
dm <- readRDS(file.path(DATA, "den_mat_fd.rds")); bm <- readRDS(file.path(DATA, "bio_mat_fd.rds"))
std <- function(x) { x <- as.data.frame(x); if (!"TAXANAME" %in% names(x)) x$TAXANAME <- x$Genus; x[!duplicated(x$TAXANAME), ] }
tb <- std(readRDS(file.path(DATA, "traits_binary_proxy_main.rds"))); tf <- std(readRDS(file.path(DATA, "traits_fuzzy_proxy_main.rds")))
tblk <- readRDS(file.path(DATA, "trait_blocks.rds"))
cl <- intersect(setdiff(names(tb), c("Genus", "TAXANAME")), setdiff(names(tf), c("Genus", "TAXANAME")))
cl <- cl[sapply(tb[cl], is.numeric) & sapply(tf[cl], is.numeric)]
tx <- intersect(intersect(tb$TAXANAME, tf$TAXANAME), union(colnames(dm), colnames(bm)))
build <- function(df) { m <- as.matrix(df[match(tx, df$TAXANAME), cl]); rownames(m) <- tx; storage.mode(m) <- "double"; m }   # NA kept: an empty trait group is missing, not zero
MB <- build(tb); MF <- build(tf); gw <- tblk$gowdis_w[colnames(MF)]; stopifnot(all(is.finite(gw)))
G_bin <- as.matrix(FD::gowdis(as.data.frame(MB)))
G_fuz <- as.matrix(suppressWarnings(FD::gowdis(as.data.frame(MF), w = gw)))
P <- function(M) { M <- M[met$SiteID, tx]; M / rowSums(M) }
P_den <- P(dm); P_bio <- P(bm); q_of <- function(P, G) rowSums((P %*% G) * P)
chk_bin <- max(abs(q_of(P_den, G_bin) - met$RaoQ_gower_den), abs(q_of(P_bio, G_bin) - met$RaoQ_gower_bio))
chk_fuz <- max(abs(q_of(P_den, G_fuz^2) / 2 - met$RaoQ_den), abs(q_of(P_bio, G_fuz^2) / 2 - met$RaoQ_bio))
cat(sprintf("provenance: binary geometry vs published RaoQ_gower max |diff| %.1e; fuzzy geometry vs published dbFD RaoQ max |diff| %.1e\n", chk_bin, chk_fuz))
if (chk_bin > 1e-12 || chk_fuz > 1e-12) stop("geometry provenance failed; nothing written")
qd <- q_of(P_den, G_fuz); qb <- q_of(P_bio, G_fuz)
stopifnot(all(qd <= met$Simpson_den + 1e-9), all(qb <= met$Simpson_bio + 1e-9))

# --- FRic and FR_SF into the same geometry (added later on 18 September 2026) ----------------------------------
# Since 00c (24 September 2026) the site table already carries FRic_den = FRic_den_fuzPool, the fuzzy-pool dbFD
# run (Cailliez, group weights, no collapsing); before 00c the upstream FRic_den was a binary-pool run. Provenance
# first: the fuzzy-pool value must reproduce from dbFD on the fuzzy matrix. Then FRic_den := FRic_den_fuzPool and
# FR_SF = Richness_FD / FRic_den, log10(x + 1), as upstream defined them.
a_fd <- dm[met$SiteID, tx]
fd_fuz <- suppressWarnings(FD::dbFD(x = as.data.frame(MF), a = a_fd, w = gw, stand.x = TRUE, corr = "cailliez",
                                    w.abun = TRUE, calc.FRic = TRUE, calc.FDiv = FALSE, calc.CWM = FALSE,
                                    stand.FRic = TRUE, m = min(nrow(MF) - 1, 10L), messages = FALSE))
fric_fuz <- fd_fuz$FRic[match(met$SiteID, names(fd_fuz$FRic))]
chk_fric <- max(abs(fric_fuz - met$FRic_den_fuzPool))
rich <- rowSums(a_fd > 0); chk_rich <- max(abs(rich - met$Richness_FD))
cat(sprintf("provenance: fuzzy-pool dbFD FRic vs published FRic_den_fuzPool max |diff| %.1e; richness max |diff| %.1e\n", chk_fric, chk_rich))
if (chk_fric > 1e-8 || chk_rich > 0) stop("FRic provenance failed; nothing written")
new <- met %>% mutate(RaoQ_gower_den = qd, RaoQ_gower_bio = qb,
                      FR_RD_den = 1 - qd / Simpson_den, FR_RD_bio = 1 - qb / Simpson_bio, FR_DB_RD = FR_RD_bio - FR_RD_den,
                      FRic_den = FRic_den_fuzPool, FR_SF_den = Richness_FD / FRic_den, FR_SF_den_log10 = log10(FR_SF_den + 1))
write_csv(new, out)
prov <- tibble(check = c("binary Gower Q reproduces published RaoQ_gower (max abs diff)", "fuzzy weighted Gower d^2/2 reproduces published dbFD RaoQ (max abs diff)",
                         "fuzzy-pool dbFD FRic reproduces published FRic_den_fuzPool (max abs diff)",
                         "Spearman, FR_RD_den binary vs fuzzy", "Spearman, FR_RD_bio binary vs fuzzy", "Spearman, FR_DB_RD binary vs fuzzy",
                         "taxa in trait matrix", "trait states", "group weights (unique)"),
               value = c(format(chk_bin, digits = 3), format(chk_fuz, digits = 3), format(chk_fric, digits = 3),
                         round(cor(met$FR_RD_den, new$FR_RD_den, method = "spearman"), 3), round(cor(met$FR_RD_bio, new$FR_RD_bio, method = "spearman"), 3),
                         round(cor(met$FR_DB_RD, new$FR_DB_RD, method = "spearman"), 3),
                         length(tx), length(cl), paste(round(sort(unique(gw)), 3), collapse = " ")))
write_csv(prov, file.path(H1, "H1_redundancy_geometry_provenance.csv"))
print(as.data.frame(prov), row.names = FALSE); cat("wrote", out, "\n")
