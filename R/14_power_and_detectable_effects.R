# =============================================================================
# 14_power_and_detectable_effects.R
# =============================================================================
# WHY THIS STEP EXISTS
#
#   Several of the paper's claims are absences: a region or a weighting in
#   which a metric does not track condition. A reader's first question about
#   any of them is whether 30 streams per region could have detected anything.
#
#   This step bounds them. It reports the smallest correlation the design could
#   reliably detect, and the achieved power behind every screened coefficient.
#   That converts "we found nothing" into "we can exclude effects larger than
#   X", which is a result rather than an absence.
#
# THE BOUND IS NOT THE WHOLE ARGUMENT
#
#   n = 30 cannot exclude moderate effects near rho = 0.3, and this step says
#   so. The comparison across regions carries the rest: the same n, design and
#   metrics, with the largest |rho| and the number of strong tests per region
#   logged at the end of this step.
#
# METHOD
#   Fisher z approximation for Spearman's rho. For n = 30 the z-transform is
#   accurate enough for a design bound; this is not a precise power analysis
#   and is not presented as one.
#
# OUTPUT
#   output/tables/14_detectable_effect_bounds.csv    log record; no later step reads it
#   output/tables/14_achieved_power_by_metric.csv    supplement Table S4 part B (legacy Table S14)
# =============================================================================

say("STEP 14 -- minimum detectable effect and achieved power")

# ---- 1. Design bound --------------------------------------------------------
mde <- function(n, alpha = 0.05, power = 0.80)
  tanh((qnorm(1 - alpha/2) + qnorm(power)) / sqrt(n - 3))

power_at <- function(rho, n, alpha = 0.05)
  pnorm(atanh(abs(rho)) * sqrt(n - 3) - qnorm(1 - alpha/2))

n_reg  <- unique(screen$n); stopifnot(length(n_reg) == 1)          # sites per region, from the screen, not typed
bounds <- tidyr::crossing(n = n_reg, alpha = c(0.05, 0.10), power = c(0.80, 0.90)) %>%
  mutate(min_detectable_rho = mde(n, alpha, power))

write_csv(bounds, file.path(OUT_TAB, "14_detectable_effect_bounds.csv"))

say("  minimum detectable |rho| at n = %d:", n_reg)
for (i in seq_len(nrow(bounds)))
  say("    alpha %.2f, power %.2f : |rho| >= %.3f",
      bounds$alpha[i], bounds$power[i], bounds$min_detectable_rho[i])

MDE_80 <- mde(n_reg, 0.05, 0.80)

# ---- 2. Achieved power behind every screened coefficient --------------------
pw <- screen %>%
  select(Family, Metric, Weighting, Region, n, rho, p_raw, p_FDR, Tier) %>%
  mutate(
    achieved_power = power_at(rho, n),
    # A null is INFORMATIVE only if the design could have caught a large effect.
    null_status = case_when(
      Tier != "not responsive"      ~ "responsive",
      abs(rho) >= MDE_80            ~ "null despite adequate power",
      TRUE                          ~ sprintf("null, bounded: excludes |rho| >= %.2f", MDE_80)
    ))

write_csv(pw, file.path(OUT_TAB, "14_achieved_power_by_metric.csv"))

# ---- 3. Two absences bounded: Valley & Ridge, and biomass-weighted geometry ----
vr <- pw %>% filter(Region == "Valley & Ridge")
say("")
say("  Valley & Ridge (%d metric x weighting tests):", nrow(vr))
say("    largest |rho| observed        : %.3f (%s)", max(abs(vr$rho), na.rm = TRUE),
    vr$Metric[which.max(abs(vr$rho))])
say("    achieved power at that value  : %.2f", max(vr$achieved_power, na.rm = TRUE))
say("    tests reaching |rho| >= %.2f   : %d", MDE_80, sum(abs(vr$rho) >= MDE_80, na.rm = TRUE))
say("    => the design excludes |rho| >= %.2f in Valley & Ridge with 80%% power.", MDE_80)
say("       It cannot exclude moderate effects near 0.3.")

geo_bio <- pw %>%
  filter(Weighting == "Biomass",
         Metric %in% c("FDiv_bio", "FDis_bio", "RaoQ_bio", "FEve_bio"))
say("")
say("  biomass-weighted FDiv, FDis, RaoQ and FEve (%d tests):", nrow(geo_bio))
say("    largest |rho| observed : %.3f", max(abs(geo_bio$rho), na.rm = TRUE))
say("    all below the %.2f detection bound : %s", MDE_80,
    all(abs(geo_bio$rho) < MDE_80, na.rm = TRUE))

# The comparison that carries the argument.
say("")
say("  THE INTERNAL COMPARISON across regions (same n, design and metrics):")
for (rg in REGIONS) {
  r <- pw %>% filter(Region == rg)
  say("    %-15s largest |rho| = %.3f | tests at 'strong' = %d of %d",
      rg, max(abs(r$rho), na.rm = TRUE), sum(r$Tier == "strong"), nrow(r))
}
say("STEP 14 complete")
