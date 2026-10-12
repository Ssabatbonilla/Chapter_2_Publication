#!/usr/bin/env bash
# Reproduce the paper's analysis from the shipped data (see README.md).
#   ./run.sh              resolve parent rows, rebuild the trait columns and the redundancy geometry, put the
#                         replicate-median site totals into the site tables, then run every analysis step
#   ./run.sh analysis     run the analysis steps only, on the derived data files as shipped
#   ./run.sh analysis 03  run one step (and the steps it depends on)
# The two environment flags select the published analysis; run_all.R refuses to run in this layout without them.
set -euo pipefail
cd "$(dirname "$0")"
if [ "${1:-}" != "analysis" ]; then
  Rscript R/00c_resolve_parent_rows.R
  Rscript R/00a_trait_columns.R
  Rscript R/00b_redundancy_geometry.R
  Rscript R/00d_site_replicate_medians.R
else
  shift
fi
CH2_BIBI_MEDIAN=1 CH2_RD_GEOMETRY=fuzzy Rscript run_all.R "$@"
