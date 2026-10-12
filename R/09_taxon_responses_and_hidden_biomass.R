# =============================================================================
# 09_taxon_responses_and_hidden_biomass.R
# =============================================================================
# THE TWO QUESTIONS THIS STEP ANSWERS
#
# (1) WHICH TAXA respond? Every taxon occurring at MIN_OCC or more sites in a
#     region is screened against the continuous B-IBI gradient, in both
#     weightings, within region (Table S8 part A; legacy Table S23), and the responders are named and
#     ranked. (Until 27 September 2026 this step also projected the responders
#     into a trait-space ordination read from a pre-pipeline file computed on
#     the retired binary traits before parent-row resolution; that section and
#     its figure were removed because nothing in the paper used them and the
#     coordinates no longer described the taxon set. The trait profile of the
#     responders, section D, is computed here from the current trait matrix.)
#
# (2) What is hiding inside a non-responding total biomass?
#     A flat total is not evidence of a static assemblage. Three different
#     things produce rho ~ 0 for total biomass:
#         (a) nothing changes;
#         (b) taxa are replaced but the total is conserved -- compensatory
#             turnover, which is an ECOLOGICAL RESULT, not a null;
#         (c) increases and decreases happen to cancel in the sum.
#     Where step 04 shows biomass-weighted structure tracking condition in a
#     region while total biomass does not, only (b) or (c) can hold. This step
#     decides which, by counting
#     the taxa that respond in biomass, measuring how much biomass they carry,
#     and testing whether the gainers offset the losers.
#
# OUTPUT
#   output/tables/09_taxon_response_screen.csv      supplement Table S8 part A (legacy Table S23); step 09b reads it
#   output/tables/09_hidden_biomass_summary.csv
#   output/tables/09_offset_test.csv
#   output/tables/09_responder_trait_profile.csv
#   output/tables/09_responding_taxa.csv
#   (the last four are diagnostic records; no later step or manuscript table reads them)
# =============================================================================

say("STEP 09 -- taxon-level responses and the biomass null")
set.seed(SEED)

MIN_OCC <- 8   # a taxon must occur at >= 8 of 30 sites in a region to be tested

# ---- A. Taxon-level response screen -----------------------------------------
screen_taxa_one <- function(rg, mat, currency) {
  m  <- dat %>% filter(Region == rg) %>% select(SiteID, BIBI)
  X  <- mat[m$SiteID, , drop = FALSE]
  occ <- colSums(X > 0)
  keep <- names(occ)[occ >= MIN_OCC]
  if (!length(keep)) return(NULL)
  map_dfr(keep, function(tx) {
    y <- X[, tx]
    ct <- suppressWarnings(cor.test(m$BIBI, y, method = "spearman", exact = FALSE))
    tibble(Region = rg, Currency = currency, TAXANAME = tx,
           n_occ = occ[[tx]],
           mean_share = mean(y / rowSums(X)),
           rho = unname(ct$estimate), p_raw = ct$p.value)
  })
}

taxa_screen <- bind_rows(
  map_dfr(REGIONS, screen_taxa_one, mat = den_mat, currency = "Density"),
  map_dfr(REGIONS, screen_taxa_one, mat = bio_mat, currency = "Biomass")
) %>%
  group_by(Region, Currency) %>%
  mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(Response = case_when(
    p_FDR <= 0.10 & rho > 0 ~ "increaser",
    p_FDR <= 0.10 & rho < 0 ~ "decreaser",
    TRUE                    ~ "no response"))

write_csv(taxa_screen, file.path(OUT_TAB, "09_taxon_response_screen.csv"))

sm <- taxa_screen %>% filter(Response != "no response") %>% count(Region, Currency, Response)
say("  responding taxa (FDR <= 0.10), of %d taxon x region x currency tests:", nrow(taxa_screen))
for (i in seq_len(nrow(sm)))
  say("    %-15s %-8s %-10s %d", sm$Region[i], sm$Currency[i], sm$Response[i], sm$n[i])

# ---- C. What is hiding inside the biomass null -------------------------------
hidden <- map_dfr(REGIONS, function(rg) {
  m <- dat %>% filter(Region == rg) %>% select(SiteID, BIBI)
  Xb <- bio_mat[m$SiteID, , drop = FALSE]
  s  <- taxa_screen %>% filter(Region == rg, Currency == "Biomass")
  inc <- s$TAXANAME[s$Response == "increaser"]
  dec <- s$TAXANAME[s$Response == "decreaser"]
  tot <- rowSums(Xb)
  b_inc <- if (length(inc)) rowSums(Xb[, inc, drop = FALSE]) else rep(0, nrow(Xb))
  b_dec <- if (length(dec)) rowSums(Xb[, dec, drop = FALSE]) else rep(0, nrow(Xb))
  sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
  tibble(
    Region = rg,
    total_biomass_rho     = sp(m$BIBI, tot),
    n_taxa_tested         = nrow(s),
    n_increasers          = length(inc),
    n_decreasers          = length(dec),
    pct_biomass_in_responders = 100 * mean((b_inc + b_dec) / tot),
    increaser_biomass_rho = if (length(inc)) sp(m$BIBI, b_inc) else NA_real_,
    decreaser_biomass_rho = if (length(dec)) sp(m$BIBI, b_dec) else NA_real_,
    # If gainers and losers offset, their site-level biomass should be negatively
    # related to each other (a plain Spearman correlation; the gradient is not partialled out).
    offset_rho            = if (length(inc) && length(dec)) sp(b_inc, b_dec) else NA_real_
  )
})
write_csv(hidden, file.path(OUT_TAB, "09_hidden_biomass_summary.csv"))

say("  the biomass null, decomposed:")
for (i in seq_len(nrow(hidden))) {
  h <- hidden[i, ]
  say("    %-15s total rho=%+.3f | %d increasers, %d decreasers of %d taxa | responders carry %.0f%% of biomass",
      h$Region, h$total_biomass_rho, h$n_increasers, h$n_decreasers, h$n_taxa_tested,
      h$pct_biomass_in_responders)
  say("                     increaser pool rho=%+.3f  decreaser pool rho=%+.3f",
      h$increaser_biomass_rho, h$decreaser_biomass_rho)
}

# ---- C3. The offset test: do gainers and losers cancel? ----------------------
# If a flat total conceals compensatory turnover, then the taxa going up and the
# taxa going down should each track the gradient while their SUM does not.
# FDR is too strict to define the pools where the signal is diffuse (a region
# may have no taxon surviving correction), so pools are defined at raw p <= 0.05.
# A region and weighting with no increaser or no decreaser at that level is skipped.
offset <- map_dfr(REGIONS, function(rg) map_dfr(c("Density","Biomass"), function(cu) {
  m <- dat %>% filter(Region == rg) %>% select(SiteID, BIBI)
  X <- (if (cu == "Density") den_mat else bio_mat)[m$SiteID, , drop = FALSE]
  s <- taxa_screen %>% filter(Region == rg, Currency == cu)
  inc <- s$TAXANAME[s$p_raw <= 0.05 & s$rho > 0]
  dec <- s$TAXANAME[s$p_raw <= 0.05 & s$rho < 0]
  if (!length(inc) || !length(dec)) return(NULL)
  b_inc <- rowSums(X[, inc, drop = FALSE]); b_dec <- rowSums(X[, dec, drop = FALSE])
  sp <- function(y) suppressWarnings(cor(m$BIBI, y, method = "spearman"))
  tibble(Region = rg, Currency = cu, n_inc = length(inc), n_dec = length(dec),
         rho_increaser_pool = sp(b_inc), rho_decreaser_pool = sp(b_dec),
         rho_pooled_sum = sp(b_inc + b_dec), rho_whole_assemblage = sp(rowSums(X)),
         offset_ratio = abs(sp(b_inc + b_dec)) /
                        max(abs(sp(b_inc)), abs(sp(b_dec))))
}))
write_csv(offset, file.path(OUT_TAB, "09_offset_test.csv"))

say("  offset test -- do the gainer and loser pools cancel in the sum?")
for (i in seq_len(nrow(offset))) {
  o <- offset[i, ]
  say("    %-15s %-8s increasers(%d) rho=%+.3f | decreasers(%d) rho=%+.3f | their SUM rho=%+.3f | whole assemblage rho=%+.3f",
      o$Region, o$Currency, o$n_inc, o$rho_increaser_pool, o$n_dec, o$rho_decreaser_pool,
      o$rho_pooled_sum, o$rho_whole_assemblage)
  say("                     offset ratio %.2f  (near 0 = the pools cancel; near 1 = no cancellation)",
      o$offset_ratio)
}

# ---- D. Trait profile: responders vs the non-responding background ----------
# The response turns out to be one-sided -- almost every taxon that survives FDR
# INCREASES with condition -- so a winner/loser contrast is not available. The
# informative contrast is responders against the background of taxa that were
# tested in the same region and did not respond. That is what turns "these taxa
# responded" into "these traits, therefore these processes".
traits <- readRDS(IN$traits) %>% filter(TAXANAME %in% colnames(den_mat))
mods   <- setdiff(names(traits), c("TAXANAME", "Genus"))   # Twardochleb fuzzy affinities since 22 Sep 2026

profile <- map_dfr(REGIONS, function(rg) map_dfr(c("Density","Biomass"), function(cu) {
  s <- taxa_screen %>% filter(Region == rg, Currency == cu) %>%
    inner_join(traits, by = "TAXANAME")
  resp <- s %>% filter(Response != "no response")
  back <- s %>% filter(Response == "no response")
  if (nrow(resp) < 3 || nrow(back) < 3) return(NULL)
  map_dfr(mods, function(md) {
    a <- resp[[md]]; b <- back[[md]]
    tibble(Region = rg, Currency = cu, Modality = md,
           responder_mean = mean(a, na.rm = TRUE), background_mean = mean(b, na.rm = TRUE),   # NA: no record in the family
           difference = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
           p_wilcox = suppressWarnings(wilcox.test(a, b)$p.value),
           n_resp = nrow(resp), n_back = nrow(back))
  }) %>% mutate(p_FDR = p.adjust(p_wilcox, method = "BH"))
}))
if (nrow(profile)) profile <- profile %>% arrange(Region, Currency, desc(abs(difference)))
write_csv(profile, file.path(OUT_TAB, "09_responder_trait_profile.csv"))

if (nrow(profile)) {
  say("  traits enriched in responders vs the tested background (FDR <= 0.10):")
  sig <- profile %>% filter(p_FDR <= 0.10)
  if (nrow(sig)) {
    for (i in seq_len(nrow(sig)))
      say("    %-15s %-8s %-42s responders %.2f vs %.2f  (diff %+.2f, FDR %.3f)",
          sig$Region[i], sig$Currency[i], sig$Modality[i], sig$responder_mean[i],
          sig$background_mean[i], sig$difference[i], sig$p_FDR[i])
  } else say("    none reached FDR <= 0.10")
}

# Name the responding taxa (FDR <= 0.10), ranked by |rho| within region and weighting.
responders <- taxa_screen %>% filter(Response != "no response") %>%
  arrange(Region, Currency, desc(abs(rho))) %>%
  select(Region, Currency, TAXANAME, n_occ, mean_share, rho, p_FDR, Response)
write_csv(responders, file.path(OUT_TAB, "09_responding_taxa.csv"))
say("  responding taxa written to 09_responding_taxa.csv (%d rows)", nrow(responders))

say("STEP 09 complete")
