# =============================================================================
# 03b_leverage_audit.R -- leave-few-out leverage audit of every responsive screen association (24 Sep 2026)
# =============================================================================
# WHY THIS STEP EXISTS
#   Methods 2.5 states that a leave-few-out leverage audit checked every responsive association. Until 24 September
#   2026 that audit lived only in an earlier script run on the observed index and the binary traits. This step is the
#   audit, on the published basis.
#
# WHAT IT DOES
#   For every screen association at the strong or exploratory tier (03_screen_all_metrics.csv), drop each site of
#   the region in turn (one run per site) and each pair of sites (one run per pair), recompute Spearman rho and its
#   raw p, and record the range of |rho| over the drops and the drop that weakens rho most. An association is flagged
#   where any single-site drop, or any pair drop, changes the sign of rho or takes it below RHO_MIN (the coefficient
#   half of the screening criterion; the FDR half is a family property and is not recomputed on n - 1 sites).
#   A flagged result is named with its flag where it is reported.
#
# OUTPUT
#   output/tables/03b_leverage_audit.csv   one row per responsive association; Table S3 part B (legacy Table S18)
# =============================================================================
say("STEP 03b -- leave-few-out leverage audit of the responsive screen associations")
screen_resp <- screen %>% filter(Tier %in% c("strong", "exploratory"))
dat_lev <- dat %>% left_join(conc, by = "SiteID")
drop_one <- function(x, y, i) suppressWarnings(cor.test(x[-i], y[-i], method = "spearman", exact = FALSE))
audit <- purrr::pmap_dfr(screen_resp %>% select(Family, Metric, Region, rho, Tier), function(Family, Metric, Region, rho, Tier) {
  s <- dat_lev %>% filter(Region == !!Region); x <- s$BIBI; y <- s[[Metric]]; n <- length(x); ids <- s$SiteID
  one <- purrr::map_dfr(seq_len(n), function(i) { ct <- drop_one(x, y, i); tibble(dropped = ids[i], rho_drop = unname(ct$estimate), p_drop = ct$p.value) })
  pairs <- utils::combn(n, 2)
  two <- purrr::map_dfr(seq_len(ncol(pairs)), function(k) { i <- pairs[, k]; ct <- drop_one(x, y, i)
                                                             tibble(dropped = paste(ids[i], collapse = "+"), rho_drop = unname(ct$estimate), p_drop = ct$p.value) })
  worst <- function(d) d[which.min(abs(d$rho_drop) * sign(d$rho_drop) * sign(rho)), ]     # the drop that weakens rho most (or flips it)
  w1 <- worst(one); w2 <- worst(two)
  tibble(Family, Metric, Region, Tier, rho_full = rho, n_sites = n,
         rho_min_abs_one = min(abs(one$rho_drop)), rho_max_abs_one = max(abs(one$rho_drop)), worst_site = w1$dropped, rho_worst_site = w1$rho_drop, p_worst_site = w1$p_drop,
         rho_min_abs_two = min(abs(two$rho_drop)), worst_pair = w2$dropped, rho_worst_pair = w2$rho_drop, p_worst_pair = w2$p_drop,
         sign_change_one = any(sign(one$rho_drop) != sign(rho)), below_030_one = any(abs(one$rho_drop) < RHO_MIN),
         sign_change_two = any(sign(two$rho_drop) != sign(rho)), below_030_two = any(abs(two$rho_drop) < RHO_MIN)) %>%
    mutate(flag = case_when(sign_change_one ~ "sign changes with one site dropped",
                            below_030_one ~ sprintf("falls below |rho| %.2f with one site dropped", RHO_MIN),
                            sign_change_two ~ "sign changes with two sites dropped",
                            below_030_two ~ sprintf("falls below |rho| %.2f with two sites dropped", RHO_MIN),
                            TRUE ~ "stable"))
})
write_csv(audit, file.path(OUT_TAB, "03b_leverage_audit.csv"))
say("  %d responsive associations audited: %d stable, %d flagged with one site dropped, %d flagged only with two",
    nrow(audit), sum(audit$flag == "stable"), sum(grepl("one site", audit$flag)), sum(grepl("two sites", audit$flag)))
for (i in which(audit$flag != "stable")) say("    %-18s %-15s %-12s rho %+.3f -> %+.3f (%s): %s", audit$Metric[i], audit$Region[i], audit$Tier[i], audit$rho_full[i],
                                              ifelse(grepl("one", audit$flag[i]), audit$rho_worst_site[i], audit$rho_worst_pair[i]),
                                              ifelse(grepl("one", audit$flag[i]), audit$worst_site[i], audit$worst_pair[i]), audit$flag[i])
say("STEP 03b complete")
