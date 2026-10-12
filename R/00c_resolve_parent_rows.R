# =============================================================================
# 00c_resolve_parent_rows.R -- parent rows resolved before any metric (24 September 2026)
# =============================================================================
# WHY THIS STEP EXISTS
#   Taxa enter the community table at the finest rank reached, so a family row (Baetidae) can sit beside genus rows
#   of the same family (Baetis, Acentrella) at the same site, and a chironomid family or subfamily row beside its
#   tribes. Kept both, the parent counts as another taxon and inflates richness, the chapter's structural control
#   (a median of 3 taxa per site before this step), and its melded trait profile enters trait space beside its
#   children's. Cuffney, Bilger & Haigler (2007, J. N. Am. Benthol. Soc. 26:286-307) show that none of the standard
#   resolutions keeps both rows, that for density- and biomass-weighted analyses "distribute parent among children"
#   (DPAC) and "remove parent or merge child" rank highest, and that one method must serve every metric.
#   Rule adopted 24 September 2026: resolve by DPAC at the sample level, before any metric.
#
# WHAT IT DOES
#   1. Reads the UNRESOLVED community table, matrices and site table: the *_unresolved_<CH2_00C_STAMP> record copies
#      (default 2026-10-06), never the files this step writes, so it can run again without resolving twice.
#   2. At each site, from the highest rank down (order or class, superfamily, family, subfamily, tribe), a parent
#      row whose children are present at that site gives its abundance, density and biomass to those children in
#      proportion to the children's median density there, and its row leaves that site. A parent with no child at
#      the site stays as its own taxon (DPAC-S). Children are found by lineage: the community table's Family column
#      for family rows; ITIS subfamily and tribe (Trait_Data/closest_kin_itis.csv) and the chironomid ranks for
#      subfamily and tribe rows; the superfamily's families for Tipuloidea; and a short lineage table below for the
#      rows above family (Gastropoda, Acari, Platyhelminthes, Megaloptera, Lepidoptera).
#   3. Rebuilds every input the pipeline reads from the resolved table, by the same definitions as upstream:
#      SiteTaxa_Density_Biomass_Transformed.csv; den_mat_fd.rds and bio_mat_fd.rds (median density and biomass of
#      the taxa with a trait record and a nonzero total); and the site table H1_site_summary_macro_metrics.csv:
#      Total_Density and Total_Biomass (sums over every taxon), Richness_alpha (taxa with abundance), Shannon (all
#      taxa, median density), Simpson_den and Simpson_bio and Richness_FD (the trait set), the tolerance metrics
#      (family tolerance values, bibi_tv_family_lookup.csv), the binary-geometry RaoQ_gower and FR_RD columns that
#      00b then recomputes in the fuzzy geometry, and the dbFD family (FRic, FEve, FDiv, FDis, RaoQ, both weightings;
#      Cailliez correction, group weights, ten axes, as step 13 verifies), with FRic_den_fuzPool, FRic_den and FR_SF
#      as 00b defines them. The community-weighted-mean columns are left for 00a, which rebuilds them.
#   4. Writes H1_parent_rows_resolution.csv (one row per parent resolved at a site: what moved to whom) and
#      H1_parent_rows_provenance.csv (the gate below and the before-and-after summary).
#
# PROVENANCE GATE (halts the script on failure; CH2_00C_VERIFY=1 runs the gate and writes nothing else)
#   Before anything is resolved, every column above is recomputed from the UNRESOLVED table and compared with the
#   unresolved site table (counts, sums, Shannon, Simpson, HBI_rel_abund and its coverage 1e-9; HBI_den, HBI_bio and
#   their coverage 1e-6), with 00a's output on the unresolved data (RaoQ_gower 1e-12, the dbFD family 1e-8) and with the
#   unresolved matrices (1e-9). So the definitions used on the resolved data are the ones that reproduce the record.
#
# FEEDS every site-level metric the pipeline screens (Table 2, Table S4 part A); H1_parent_rows_provenance.csv is
# also read by step 11 and by run_all.R's banner.
#
# Run:  Rscript R/00c_resolve_parent_rows.R   then 00a, 00b, 00d, then run_all.R (./run.sh; see README.md).
# =============================================================================
suppressPackageStartupMessages({ library(FD); library(dplyr); library(readr); library(tibble); library(tidyr) })
args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
PIPE_DIR <- normalizePath(file.path(dirname(gsub("~\\+~", " ", f)), ".."))
source(file.path(PIPE_DIR, "R", "00_paths.R"))   # CH2_ROOT, DATA, RA, H1, REPO_LAYOUT (one definition for every script)
VERIFY_ONLY <- nzchar(Sys.getenv("CH2_00C_VERIFY"))
SKIP_GATE <- nzchar(Sys.getenv("CH2_00C_SKIP_GATE"))   # development only: the gate takes minutes; skip it only right after a verify-only run passed on the same inputs
# 6 Oct 2026: the default is the unresolved record of the published run (the corrected 6 October density
# method). Earlier records (2026-09-24, 2026-10-03) stay reproducible with CH2_00C_STAMP set to them, where shipped.
# Until this date the default was 2026-09-24 and the record run set the stamp by hand, so a plain rerun would have
# resolved from the stale record under the same filenames.
STAMP <- Sys.getenv("CH2_00C_STAMP", "2026-10-06")
U <- list(long = file.path(DATA, sprintf("SiteTaxa_Density_Biomass_Transformed_unresolved_%s.csv", STAMP)),
          den  = file.path(DATA, sprintf("den_mat_fd_unresolved_%s.rds", STAMP)),
          bio  = file.path(DATA, sprintf("bio_mat_fd_unresolved_%s.rds", STAMP)),
          site = file.path(H1, sprintf("H1_site_summary_macro_metrics_unresolved_%s.csv", STAMP)),
          twa  = file.path(H1, sprintf("H1_site_summary_macro_metrics_twardochleb_unresolved_%s.csv", STAMP)))   # 00a's output on the unresolved data
W <- list(long = file.path(DATA, "SiteTaxa_Density_Biomass_Transformed.csv"),
          den  = file.path(DATA, "den_mat_fd.rds"), bio = file.path(DATA, "bio_mat_fd.rds"),
          site = file.path(H1, "H1_site_summary_macro_metrics.csv"),
          log  = file.path(H1, "H1_parent_rows_resolution.csv"), prov = file.path(H1, "H1_parent_rows_provenance.csv"))
# The unresolved copies are the record this step resolves from. They were made from the files as they stood before
# this step ran (first on 24 September 2026; the default record on 6 October), and they ship with the data. They are
# never recreated here: copying the resolved files into their place would pass the gate and resolve nothing.
miss <- unlist(U)[!file.exists(unlist(U))]
if (length(miss)) stop("[00c] unresolved record copies missing: ", paste(basename(miss), collapse = ", "))
std_site <- function(x) toupper(trimws(as.character(x)))

# ---- 1. Inputs ---------------------------------------------------------------------------------------------------
raw0 <- read_csv(U$long, show_col_types = FALSE) %>% mutate(SiteID = std_site(SiteID))
site0 <- read_csv(U$site, show_col_types = FALSE) %>% mutate(SiteID = std_site(SiteID))
twa0  <- read_csv(U$twa, show_col_types = FALSE) %>% mutate(SiteID = std_site(SiteID))
den0 <- readRDS(U$den); bio0 <- readRDS(U$bio)
tx  <- read.csv(file.path(DATA, "CONUS_Ancillary_Taxonomy.csv"), check.names = FALSE, stringsAsFactors = FALSE) %>% distinct(Submitted_name, .keep_all = TRUE)
kin <- read.csv(file.path(DATA, "Trait_Data", "closest_kin_itis.csv"), stringsAsFactors = FALSE)
tvl <- read_csv(file.path(DATA, "bibi_tv_family_lookup.csv"), show_col_types = FALSE) %>%
  transmute(Family = toupper(trimws(family)), TV = as.numeric(bibi_tv)) %>% distinct(Family, .keep_all = TRUE)
std_tr <- function(x) { x <- as.data.frame(x); if (!"TAXANAME" %in% names(x)) x$TAXANAME <- x$Genus; x[!duplicated(x$TAXANAME), ] }
tf <- std_tr(readRDS(file.path(DATA, "traits_fuzzy_proxy_main.rds"))); tb <- std_tr(readRDS(file.path(DATA, "traits_binary_proxy_main.rds")))
tblk <- readRDS(file.path(DATA, "trait_blocks.rds"))
ABUND_COLS <- c("Mean_Corrected_Abundance", "Median_Corrected_Abundance", "Total_Corrected_Abundance",
                "Mean_Density_per_m2", "Median_Density_per_m2", "Total_Density_per_m2",
                "Mean_Biomass_per_m2", "Median_Biomass_per_m2", "Total_Biomass_per_m2")
stopifnot(all(ABUND_COLS %in% names(raw0)))

# ---- 2. Lineage: the rank of every row and the children of every parent ------------------------------------------
HIGHER <- c("Acari", "Hydracarina", "Oligochaeta", "Gastropoda", "Turbellaria", "Platyhelminthes", "Collembola",
            "Lepidoptera", "Megaloptera", "Thysanoptera", "Chilopoda")
GASTROPOD_FAMS <- c("ANCYLIDAE", "HYDROBIIDAE", "LYMNAEIDAE", "PHYSIDAE", "PLANORBIDAE", "PLEUROCERIDAE", "VIVIPARIDAE", "BITHYNIIDAE")
TIPULOIDEA_FAMS <- c("TIPULIDAE", "LIMONIIDAE", "PEDICIIDAE", "CYLINDROTOMIDAE")
CHIRONOMINAE_TRIBES <- c("Chironomini", "Tanytarsini")
taxa <- raw0 %>% distinct(TAXANAME, Family) %>% mutate(Family = toupper(Family)) %>% distinct(TAXANAME, .keep_all = TRUE) %>%
  left_join(tx %>% transmute(TAXANAME = Submitted_name, Order, Genus_conus = Genus), by = "TAXANAME") %>%
  left_join(kin %>% transmute(TAXANAME = Genus, Subfamily, Tribe), by = "TAXANAME") %>%
  mutate(rank = case_when(
    !is.na(Genus_conus) & Genus_conus != "" ~ "genus",
    TAXANAME %in% HIGHER ~ "order_or_higher",
    grepl("oidea$|oidia$", TAXANAME) ~ "superfamily",
    grepl("idae$", TAXANAME) ~ "family",
    grepl("inae$", TAXANAME) ~ "subfamily",
    grepl("ini$", TAXANAME) ~ "tribe",
    TRUE ~ "genus"))                                   # a genus the CONUS file lacks (Asellus, Corbicula, Lirceus ...)
subfam_of <- function(t) { s <- taxa$Subfamily[match(t, taxa$TAXANAME)]; ifelse(is.na(s), "", s) }
tribe_of  <- function(t) { s <- taxa$Tribe[match(t, taxa$TAXANAME)];     ifelse(is.na(s), "", s) }
children_of <- function(parent) {
  r <- taxa$rank[taxa$TAXANAME == parent]; fam <- taxa$Family[taxa$TAXANAME == parent]; other <- taxa$TAXANAME != parent
  k <- switch(r,
    "family"    = taxa$TAXANAME[other & taxa$Family == fam & taxa$rank %in% c("genus", "subfamily", "tribe")],
    "subfamily" = if (parent == "Chironominae") taxa$TAXANAME[other & taxa$Family == "CHIRONOMIDAE" & (taxa$TAXANAME %in% CHIRONOMINAE_TRIBES | subfam_of(taxa$TAXANAME) == "Chironominae")]
                  else taxa$TAXANAME[other & taxa$Family == fam & tolower(subfam_of(taxa$TAXANAME)) == tolower(sub("Ceraptogoninae", "Ceratopogoninae", parent))],
    "tribe"     = taxa$TAXANAME[other & taxa$Family == fam & tribe_of(taxa$TAXANAME) == parent],
    "superfamily" = taxa$TAXANAME[other & taxa$Family %in% TIPULOIDEA_FAMS],
    "order_or_higher" = switch(parent,
        "Gastropoda" = taxa$TAXANAME[other & taxa$Family %in% GASTROPOD_FAMS],
        "Acari" = intersect("Hydracarina", taxa$TAXANAME),
        "Platyhelminthes" = intersect("Turbellaria", taxa$TAXANAME),
        "Megaloptera" = taxa$TAXANAME[other & (taxa$Order %in% "Megaloptera" | taxa$Family %in% c("SIALIDAE", "CORYDALIDAE"))],
        "Lepidoptera" = taxa$TAXANAME[other & (taxa$Order %in% "Lepidoptera" | taxa$Family %in% "CRAMBIDAE")],
        character(0)),
    character(0))
  k[!is.na(k)]
}
RANK_ORDER <- c("order_or_higher", "superfamily", "family", "subfamily", "tribe")
parents <- taxa %>% filter(rank %in% RANK_ORDER) %>% mutate(children = lapply(TAXANAME, children_of), n_children = lengths(children)) %>%
  arrange(match(rank, RANK_ORDER), TAXANAME)
cat(sprintf("[00c] %d rows above genus rank; %d have at least one child in the community table\n", nrow(parents), sum(parents$n_children > 0)))

# ---- 3. DPAC-S: resolve each parent at each site, from the top rank down ---------------------------------------------
resolve <- function(raw) {
  raw <- raw %>% mutate(.row = row_number()); logs <- list()
  for (i in seq_len(nrow(parents))) {
    p <- parents$TAXANAME[i]; kids <- parents$children[[i]]; if (!length(kids)) next
    psites <- raw %>% filter(TAXANAME == p, Median_Density_per_m2 > 0) %>% pull(SiteID)
    for (s in psites) {
      ch <- raw %>% filter(SiteID == s, TAXANAME %in% kids, Median_Density_per_m2 > 0)
      if (!nrow(ch)) next
      pr <- raw %>% filter(SiteID == s, TAXANAME == p); stopifnot(nrow(pr) == 1)
      share <- ch$Median_Density_per_m2 / sum(ch$Median_Density_per_m2)
      pos <- match(ch$.row, raw$.row)                                         # positions, not row numbers: rows leave the table as parents resolve
      for (cc in ABUND_COLS) raw[[cc]][pos] <- raw[[cc]][pos] + share * pr[[cc]]
      logs[[length(logs) + 1]] <- tibble(SiteID = s, Typology = pr$Typology, parent = p, rank = parents$rank[i], n_children = nrow(ch),
                                         children = paste(ch$TAXANAME, collapse = "; "),
                                         abundance_moved = pr$Total_Corrected_Abundance, density_moved = pr$Median_Density_per_m2, biomass_moved = pr$Median_Biomass_per_m2)
      raw <- raw %>% filter(.row != pr$.row)                                # the parent leaves that site
    }
  }
  # transformed columns follow the resolved values (log1p and sqrt, as upstream wrote them)
  for (b in c("Mean_Corrected_Abundance", "Median_Corrected_Abundance", "Mean_Density", "Mean_Biomass", "Median_Density", "Median_Biomass")) {
    src <- if (grepl("Abundance", b)) b else paste0(b, "_per_m2")
    raw[[paste0(b, "_log")]] <- log1p(raw[[src]]); raw[[paste0(b, "_sqrt")]] <- sqrt(raw[[src]])
  }
  list(raw = raw %>% select(-.row), log = bind_rows(logs))
}

# ---- 4. Every metric the pipeline reads, recomputed from a community table --------------------------------------------
metrics_from <- function(raw, site_ref) {
  sites <- site_ref$SiteID
  wide <- function(col) { m <- raw %>% select(SiteID, TAXANAME, v = all_of(col)) %>% pivot_wider(names_from = TAXANAME, values_from = v, values_fill = 0)
                          M <- as.matrix(m[, -1]); rownames(M) <- m$SiteID; M[sites, , drop = FALSE] }
  Mden <- wide("Median_Density_per_m2"); Mbio <- wide("Median_Biomass_per_m2")
  present <- colnames(Mden)[colSums(Mden) > 0]
  trait_set <- sort(intersect(intersect(tf$TAXANAME, tb$TAXANAME), present))
  den <- Mden[, trait_set, drop = FALSE]; bio <- Mbio[, trait_set, drop = FALSE]
  # structure and totals, every taxon
  st <- raw %>% group_by(SiteID) %>% summarise(Total_Density = sum(Median_Density_per_m2, na.rm = TRUE), Total_Biomass = sum(Median_Biomass_per_m2, na.rm = TRUE),
                                               Richness_alpha = sum(Total_Corrected_Abundance > 0, na.rm = TRUE), .groups = "drop")
  st <- st[match(sites, st$SiteID), ]
  st$Shannon <- as.numeric(vegan::diversity(Mden, index = "shannon"))
  st$Simpson_den <- as.numeric(vegan::diversity(den, index = "simpson")); st$Simpson_bio <- as.numeric(vegan::diversity(bio, index = "simpson"))
  st$Richness_FD <- as.numeric(rowSums(den > 0))
  # tolerance, family values by the community table's Family column
  wt <- function(w, TV) { ok <- is.finite(w) & is.finite(TV) & w > 0; if (!any(ok)) return(NA_real_); sum(w[ok] * TV[ok]) / sum(w[ok]) }
  h <- raw %>% mutate(Family = toupper(Family)) %>% left_join(tvl, by = "Family") %>% group_by(SiteID) %>%
    summarise(HBI_rel_abund = wt(Total_Corrected_Abundance, TV), HBI_den = wt(.data[[HBI_DEN_W]], TV), HBI_bio = wt(.data[[HBI_BIO_W]], TV),
              HBI_coverage_pct = 100 * sum(Total_Corrected_Abundance[is.finite(TV)], na.rm = TRUE) / sum(Total_Corrected_Abundance, na.rm = TRUE),
              HBI_density_coverage_pct = 100 * sum(.data[[HBI_DEN_W]][is.finite(TV)], na.rm = TRUE) / sum(.data[[HBI_DEN_W]], na.rm = TRUE),
              HBI_biomass_coverage_pct = 100 * sum(.data[[HBI_BIO_W]][is.finite(TV)], na.rm = TRUE) / sum(.data[[HBI_BIO_W]], na.rm = TRUE),
              HBI_taxa_with_TV = sum(is.finite(TV) & Total_Corrected_Abundance > 0), HBI_taxa_total = sum(Total_Corrected_Abundance > 0), .groups = "drop")
  st <- st %>% left_join(h, by = "SiteID")
  # trait geometry: the fuzzy file with group weights (dbFD family), and the binary file without weights (RaoQ_gower, FR_RD)
  cl <- intersect(setdiff(names(tb), c("Genus", "TAXANAME")), setdiff(names(tf), c("Genus", "TAXANAME")))
  cl <- cl[sapply(tb[cl], is.numeric) & sapply(tf[cl], is.numeric)]
  build <- function(df) { m <- as.matrix(df[match(trait_set, df$TAXANAME), cl]); rownames(m) <- trait_set; storage.mode(m) <- "double"; m }
  MB <- build(tb); MF <- build(tf)
  keep <- (colSums(MB, na.rm = TRUE) > 0) | (colSums(MF, na.rm = TRUE) > 0)                    # as 00a: a state no taxon has is not a column
  MFk <- MF[, keep, drop = FALSE]; gw <- tblk$gowdis_w[colnames(MFk)]; stopifnot(all(is.finite(gw)))
  G_bin <- as.matrix(FD::gowdis(as.data.frame(MB)))                                               # as 00a and 00b: unweighted, every column
  P <- function(M) M / rowSums(M); q_of <- function(P, G) rowSums((P %*% G) * P)
  st$RaoQ_gower_den <- as.numeric(q_of(P(den), G_bin)); st$RaoQ_gower_bio <- as.numeric(q_of(P(bio), G_bin))
  st$FR_RD_den <- 1 - st$RaoQ_gower_den / st$Simpson_den; st$FR_RD_bio <- 1 - st$RaoQ_gower_bio / st$Simpson_bio; st$FR_DB_RD <- st$FR_RD_bio - st$FR_RD_den
  run_fd <- function(a) suppressWarnings(FD::dbFD(x = MFk, a = a, w = gw, stand.x = TRUE, corr = "cailliez", w.abun = TRUE,
                                                  calc.FRic = TRUE, calc.FDiv = TRUE, stand.FRic = TRUE, calc.CWM = FALSE, m = 10, messages = FALSE))
  fd_den <- run_fd(den); fd_bio <- run_fd(bio)
  pick <- function(o, nm) as.numeric(o[[nm]][match(sites, names(o[[nm]]))])
  st$FRic_den_fuzPool <- pick(fd_den, "FRic"); st$FRic_den <- st$FRic_den_fuzPool
  st$FEve_den <- pick(fd_den, "FEve"); st$FDiv_den <- pick(fd_den, "FDiv"); st$FDis_den <- pick(fd_den, "FDis"); st$RaoQ_den <- pick(fd_den, "RaoQ")
  st$FEve_bio <- pick(fd_bio, "FEve"); st$FDiv_bio <- pick(fd_bio, "FDiv"); st$FDis_bio <- pick(fd_bio, "FDis"); st$RaoQ_bio <- pick(fd_bio, "RaoQ")
  st$FR_SF_den <- st$Richness_FD / st$FRic_den; st$FR_SF_den_log10 <- log10(st$FR_SF_den + 1)
  list(site = st, den = den, bio = bio, n_taxa = length(present), n_trait = length(trait_set))
}
# the weights that reproduce the upstream tolerance metrics exactly (checked 24 Sep 2026): total density and total biomass per m2
HBI_DEN_W <- "Total_Density_per_m2"; HBI_BIO_W <- "Total_Biomass_per_m2"

# ---- 5. Provenance gate on the unresolved table ------------------------------------------------------------------------
if (SKIP_GATE && !VERIFY_ONLY) {
  prev <- if (file.exists(W$prov)) read_csv(W$prov, show_col_types = FALSE) else NULL
  stopifnot(!is.null(prev), any(grepl("^gate", prev$mode)), all(prev$reproduced[grepl("^gate", prev$mode)] %in% TRUE))   # a passed gate must be on record
  cat("[00c] gate skipped on this run (CH2_00C_SKIP_GATE); the passed gate on record is kept in the provenance table\n")
  gate <- prev %>% filter(grepl("^gate", mode)) %>% mutate(mode = NULL, stamp = NULL)
  m0 <- list(n_taxa = length(unique(raw0$TAXANAME[raw0$Median_Density_per_m2 > 0])), n_trait = ncol(den0))
} else {
m0 <- metrics_from(raw0, site0)
EXACT <- c("Total_Density", "Total_Biomass", "Richness_alpha", "Shannon", "Simpson_den", "Simpson_bio", "Richness_FD", "HBI_rel_abund",
           "HBI_coverage_pct", "HBI_taxa_with_TV", "HBI_taxa_total")
GEOM  <- c("RaoQ_gower_den", "RaoQ_gower_bio")          # 00a's binary geometry, 1e-12
NEAR  <- c("HBI_den", "HBI_bio", "HBI_density_coverage_pct", "HBI_biomass_coverage_pct")
FDCOL <- c("FRic_den_fuzPool", "FEve_den", "FDiv_den", "FDis_den", "RaoQ_den", "FEve_bio", "FDiv_bio", "FDis_bio", "RaoQ_bio")
tol <- c(setNames(rep(1e-9, length(EXACT)), EXACT), setNames(rep(1e-6, length(NEAR)), NEAR), setNames(rep(1e-12, length(GEOM)), GEOM), setNames(rep(1e-8, length(FDCOL)), FDCOL))
REF <- c(setNames(rep("upstream", length(c(EXACT, NEAR))), c(EXACT, NEAR)), setNames(rep("00a", length(c(GEOM, FDCOL))), c(GEOM, FDCOL)))
gate <- bind_rows(lapply(names(tol), function(cc) {
  ref <- if (REF[[cc]] == "00a") twa0 else site0
  a <- m0$site[[cc]][match(ref$SiteID, m0$site$SiteID)]; b <- ref[[cc]]
  tibble(column = cc, reference = REF[[cc]], max_abs_diff = max(abs(a - b), na.rm = TRUE), n_na = sum(is.na(a) != is.na(b)), tolerance = tol[[cc]], reproduced = max(abs(a - b), na.rm = TRUE) <= tol[[cc]] && sum(is.na(a) != is.na(b)) == 0)
}))
gate <- bind_rows(gate, tibble(column = "den_mat_fd (median density, trait set)", reference = "upstream", max_abs_diff = max(abs(m0$den[rownames(den0), colnames(den0)] - den0)), n_na = 0L, tolerance = 1e-9, reproduced = identical(colnames(m0$den), colnames(den0)) && max(abs(m0$den[rownames(den0), colnames(den0)] - den0)) <= 1e-9),
                        tibble(column = "bio_mat_fd (median biomass, trait set)", reference = "upstream", max_abs_diff = max(abs(m0$bio[rownames(bio0), colnames(bio0)] - bio0)), n_na = 0L, tolerance = 1e-9, reproduced = identical(colnames(m0$bio), colnames(bio0)) && max(abs(m0$bio[rownames(bio0), colnames(bio0)] - bio0)) <= 1e-9))
cat("[00c] provenance gate on the unresolved table:\n"); print(as.data.frame(gate %>% mutate(max_abs_diff = signif(max_abs_diff, 3))), row.names = FALSE)
if (!all(gate$reproduced)) { write_csv(gate %>% mutate(mode = "gate", stamp = STAMP), W$prov); stop("[00c] provenance gate FAILED; nothing resolved, nothing written but the provenance table") }
if (VERIFY_ONLY) { write_csv(gate %>% mutate(mode = "gate (verify only)", stamp = STAMP), W$prov); cat("[00c] verify only: gate passed, nothing else written\n"); quit(status = 0) }
}

# ---- 6. Resolve, recompute, write --------------------------------------------------------------------------------------
res <- resolve(raw0); raw1 <- res$raw; lg <- res$log
m1 <- metrics_from(raw1, site0)
# the site table: every recomputed column replaces the unresolved value; everything else (BIBI, Typology, CWM columns for 00a) is carried
site1 <- site0
for (cc in intersect(names(m1$site), names(site1))) site1[[cc]] <- m1$site[[cc]][match(site1$SiteID, m1$site$SiteID)]
for (cc in setdiff(names(m1$site), names(site1))) site1[[cc]] <- m1$site[[cc]][match(site1$SiteID, m1$site$SiteID)]
write_csv(raw1, W$long); saveRDS(m1$den, W$den); saveRDS(m1$bio, W$bio); write_csv(site1, W$site); write_csv(lg, W$log)
rich0 <- site0$Richness_alpha; rich1 <- site1$Richness_alpha[match(site0$SiteID, site1$SiteID)]
summ <- tibble(column = c("parent-site cases resolved", "parents resolved at least once", "sites with a resolution", "taxa in the community table (before)", "taxa in the community table (after)",
                          "taxa in the trait set (before)", "taxa in the trait set (after)", "median taxa per site removed", "max taxa per site removed", "density moved (share of total)", "biomass moved (share of total)"),
               max_abs_diff = c(nrow(lg), n_distinct(lg$parent), n_distinct(lg$SiteID), m0$n_taxa, m1$n_taxa, m0$n_trait, m1$n_trait, median(rich0 - rich1), max(rich0 - rich1),
                                sum(lg$density_moved) / sum(raw0$Median_Density_per_m2, na.rm = TRUE), sum(lg$biomass_moved) / sum(raw0$Median_Biomass_per_m2, na.rm = TRUE)),
               reference = "resolution", n_na = NA_integer_, tolerance = NA_real_, reproduced = NA)
write_csv(bind_rows(gate %>% mutate(mode = "gate"), summ %>% mutate(mode = "resolution summary")) %>% mutate(stamp = STAMP), W$prov)
stopifnot(abs(sum(raw1$Median_Density_per_m2) - sum(raw0$Median_Density_per_m2)) < 1e-6, abs(sum(raw1$Median_Biomass_per_m2) - sum(raw0$Median_Biomass_per_m2)) < 1e-6)  # nothing lost, nothing made
cat(sprintf("[00c] resolved %d parent-site cases (%d parents, %d sites); taxa %d -> %d, trait set %d -> %d; richness per site fell by a median of %g (max %g)\n",
            nrow(lg), n_distinct(lg$parent), n_distinct(lg$SiteID), m0$n_taxa, m1$n_taxa, m0$n_trait, m1$n_trait, median(rich0 - rich1), max(rich0 - rich1)))
cat(sprintf("[00c] wrote %s, %s, %s, %s, %s, %s\n", basename(W$long), basename(W$den), basename(W$bio), basename(W$site), basename(W$log), basename(W$prov)))
