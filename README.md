# Stream macroinvertebrate structure and traits along a biological condition gradient, weighted by density and by biomass

Code and data for the paper by Sergio A. Sabat-Bonilla and co-authors (Department of Entomology, Virginia
Polytechnic Institute and State University), submitted to Freshwater Biology. This repository holds the R scripts,
the site-by-taxon density and biomass matrices, the trait matrix, the biological index values and a data dictionary,
and it reproduces every analysis, table and figure the paper reports. The macroinvertebrate collections themselves
are published, with full metadata, as U.S. Geological Survey data releases, which the paper cites; this repository
holds the processed tables derived from them.

## What the paper tests

Ninety wadeable streams were sampled in three physiographic regions of the Chesapeake Bay watershed (Valley and
Ridge, Piedmont, Coastal Plain; 30 streams each). Each stream's macroinvertebrate community was described by its
structure (richness, diversity, tolerance) and by the functional traits of its taxa (functional diversity,
redundancy and inequality), each trait metric weighted by areal density and by biomass, and by the difference
between the two weightings. Within each region, every metric was screened against a gradient of biological
condition, the Chessie Benthic Index of Biotic Integrity (B-IBI) of Smith et al. (2017), by Spearman rank correlation
with a Benjamini-Hochberg correction within each metric family. The supporting analyses are a leave-few-out leverage
audit, multivariate tests of community structure, statistical power, a trait-shuffling null model, a taxon-level
screen with tolerance values, and a sensitivity analysis that removes the four chironomid taxa (the tribes
Chironomini and Tanytarsini and the subfamilies Orthocladiinae and Tanypodinae) that make up a third to two thirds of
the trait-set density.

**Analyses the pipeline runs that this paper does not report.** Step 15 (RLQ and fourth-corner analysis of single
trait states) and step 15b (the density-against-biomass difference of each trait state) belong to a separate
manuscript on trait states. They still run here because step 20 checks its fourth-corner rerun against step 15's
tables and steps 06 and 18 build that manuscript's figures and tables from them; their manuscript assets are written to
`output/manuscript/trait_state_manuscript/` and are not part of this paper. Step 16 (beta diversity partitioned into
turnover and nestedness) is likewise not reported; it runs because step 18 builds the repository table
`output/manuscript/repository_tables/regional_pools_shared_taxa.csv` from one of its tables.

## Reproducing the paper

Requirements: R version 4.5.1 (2025-06-13) or later and the packages in the table below. `cairosvg` (a command-line tool,
`pip install cairosvg`) is optional: with it the PDF figures are rebuilt from the SVG with true minus signs and
en dashes; without it the PDF falls back to R's base device, which substitutes hyphens, and the run says so. The
PNG and SVG figures are correct either way. Times New Roman is used when the system has it, otherwise the
default serif face.

```bash
git clone https://github.com/Ssabatbonilla/Chapter_2_Publication.git
cd Chapter_2_Publication
./run.sh              # everything: about 37 minutes on an Apple silicon Mac (macOS, R 4.5.1), most of it R/00c's provenance gate
./run.sh analysis     # the analysis steps only, on the derived data files as shipped: about 9 minutes
```

`run.sh` runs, in order:

1. `R/00c_resolve_parent_rows.R`: resolves taxa recorded above genus beside their children at the same site
   (distribute parent among children, Cuffney et al. 2007), rebuilding the community table, the two site-by-taxon
   matrices and the site table from the `*_unresolved_2026-10-06` copies in `data/`. A provenance gate first
   recomputes every metric from the unresolved data and stops unless each reproduces the shipped site table.
2. `R/00a_trait_columns.R`: rebuilds the functional diversity and community-weighted mean columns of the site
   table from the trait matrix.
3. `R/00b_redundancy_geometry.R`: recomputes the redundancy family in the fuzzy group-weighted Gower geometry.
4. `R/00d_site_replicate_medians.R`: puts the site totals (`Total_Density`, `Total_Biomass`) and the site Hilsenhoff
   Biotic Index into the two site tables as the median over the processed replicates of each replicate's value, read
   from `data/Step2_SiteReplicateMedians.csv`, and keeps the mean versions beside them in the columns ending `_mean`.
   It first checks that the mean totals equal the site tables' totals, and a second run leaves the tables as the
   first left them.
5. `CH2_BIBI_MEDIAN=1 CH2_RD_GEOMETRY=fuzzy Rscript run_all.R`: the analysis steps listed below, in dependency order.

Steps 1 to 4 rewrite the derived files in `data/` (listed as "derived" in `data/DATA_DICTIONARY.md`); they come back
byte-identical to the shipped copies, so `./run.sh analysis` on the shipped files gives the same results.

**The two environment flags are required.** They select the published analysis: the rarefied-median B-IBI as
the condition gradient, and the fuzzy group-weighted geometry for the redundancy metrics. Without them the
runner would select two earlier bases under the same file names, so in this repository `run_all.R` refuses to start
unless both are set, and prints the full command. The run log (`output/logs/run_<date>_<time>.log`) opens with these
banner lines, which must read:

```
=== B-IBI basis: RAREFIED MEDIAN  (matches the manuscript) ===
=== Redundancy geometry: FUZZY group-weighted Gower  (matches the manuscript since 18 September 2026) ===
=== Step 15 weightings: density AND biomass  (both written) ===
=== Parent rows: resolved by DPAC before any metric (R/00c; 301 parent-site cases)  (matches the manuscript since 24 September 2026) ===
```

A single step and the steps it depends on: `./run.sh analysis 03`. The fourth-corner reruns in steps 15 and 20 run
their five seeds in parallel on Unix (`parallel::mclapply`); each seed sets its own random seed first, so the results
do not depend on the number of cores. Do not run other R jobs beside them. Step 20 calls `mclapply` with more than one core,
which R does not support on Windows; run the pipeline on macOS or Linux (on Windows, through the Windows Subsystem for
Linux).

`output/` as shipped is what `./run.sh` wrote in this repository (log `output/logs/run_20261012_180101.log`, with the session
record beside it). It reproduces the author's own run of the same scripts on the project data: every table, every
manuscript table and every PNG figure identical; the PDF figures carry a creation date and are not compared byte for
byte.

### How density and biomass are computed

The analysis follows the laboratory's density method (James et al. 2022). Each replicate sample's count, after the
subsample expansion, is divided by the sampled area of 0.25 square meter, and biomass is the sum of the length-mass
estimates of the counted macroinvertebrates over the same area. A site-by-taxon density or biomass is the mean over the
site's processed replicates, zeros included (the matrices, Bray-Curtis, the functional diversity metrics and the
inequality metrics all use these). The site totals and the site Hilsenhoff Biotic Index are the median over the
processed replicates of each replicate's value (R/00d). Replicates that were lost are treated as missing, not as
zeros: the divisor is the number of replicates processed (five at 86 of the 90 sites, four at three sites, three at
one site). The data chain that applies the method (the raw count and length sheets, the dated correction files that
record every correction to them, the length-mass estimates and the replicate roster) is upstream of this repository
and is not shipped; the shipped `*_unresolved_2026-10-06` records carry its result, and `data/DATA_DICTIONARY.md`
describes each file.

### Packages

| Package | Version used |
|---|---|
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| readr | 2.1.5 |
| tibble | 3.2.1 |
| purrr | 1.2.2 |
| stringr | 1.5.1 |
| vegan | 2.7.2 |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.2 |
| ragg | 1.4.0 |
| svglite | 2.1.3 |
| systemfonts | 1.2.3 |
| scales | 1.4.0 |
| FD | 1.0.12.3 |
| ade4 | 1.7.24 |
| betapart | 1.6.1 |
| mgcv | 1.9.3 |
| rlang | 1.3.0 |
| sf | 1.0.20 |
| maps | 3.4.2.1 |
| ggspatial | 1.1.9 |
| cowplot | 1.1.3 |
| parallel | 4.5.1 |

`parallel` is part of R. `FD` pulls in `geometry` and `ape`; `sf` needs GDAL, GEOS and PROJ on the system
(the CRAN binary carries them on macOS and Windows). Step 21 writes the versions of the run, and a citation for R and
for every package, to `output/tables/R_packages_and_versions.csv` and `output/manuscript/software/`.

## The steps

Each script opens with a header saying what the step does, what it reads and writes, and which table or figure of the
paper it feeds. `R/00_paths.R` decides where the data live (`data/` beside `R/` here); `R/00_config.R` holds every
threshold, method choice and label (the metric dictionary that every figure and table reads); `run_all.R` runs the
steps in dependency order and prints the banners above. Paper numbers are given first; the pipeline's working number,
which step 18 writes and step 18b maps, is in parentheses.

| Step | Script | What it does | In the paper |
|---|---|---|---|
| 01 | `01_load_and_verify.R` | Loads the inputs; stops unless every site's B-IBI equals its ICPRB release value and the structural guards hold | (gate) |
| 02 | `02_concentration_metrics.R` | Gini coefficient and N50 (taxa making up half the total), both weightings, and their divergence | Table 2, Figure 4 |
| 03 | `03_screen_metrics.R` | Spearman screen of every metric against the B-IBI within region; Benjamini-Hochberg within each metric family and region; the tier rule | Table 2, Table S4 part A, Figures 3 to 5 |
| 03b | `03b_leverage_audit.R` | Leave-few-out audit of every responsive association | Table S3 part B (S18) |
| 04 | `04_multivariate_within_region.R` | PERMANOVA, ANOSIM and dispersion tests on Bray-Curtis dissimilarity, within and across regions | Table S5 (S9) |
| 07 | `07_tables.R` | Screening and multivariate tables that step 18 numbers | Tables 2, S4 part A, S5 |
| 09 | `09_taxon_responses_and_hidden_biomass.R` | Taxon-level screen (taxa at eight or more sites per region), both weightings | Table S8 (S23) |
| 09b | `09b_taxon_tolerance.R` | Tolerance values of the responding taxa | Table S8 (S24) |
| 10 | `10_taxon_set_accounting.R` | Trait coverage by site; the excluded taxa; sensitivity to the 80 percent coverage criterion; family shares of density and biomass | Table S3 part A (S5) |
| 11 | `11_taxon_convention_and_traitspace_quality.R` | Community set against trait set; trait-space dimensionality (Maire et al. 2015); functional richness against richness | Table S4 part C (S25) |
| 12 | `12_evidence_sufficiency.R` | Every responsive association partialled on trait coverage | Table S3 part A (S5) |
| 13 | `13_fd_provenance.R` | Recomputes every functional diversity and redundancy metric from the matrices and trait files and stops unless each reproduces | (gate) |
| 14 | `14_power_and_detectable_effects.R` | Minimum detectable effect and achieved power per test | Table S4 part B (S14) |
| 15 | `15_rlq_fourthcorner.R` | RLQ and fourth-corner analysis per region-by-trait-group block, both weightings, five seeds | not in this paper (see above) |
| 15b | `15b_weighting_difference.R` | Density-against-biomass difference of each trait state's fourth-corner statistic | not in this paper |
| 06 | `06_figures.R` | Figures 3, 4 and 5 (and the trait-state figures of the separate manuscript) | Figures 3 to 5 |
| 16 | `16_beta_turnover_nestedness.R` | Beta diversity partition and cross-region pool comparison | not in this paper |
| 17 | `17_null_model_ses.R` | Null model: trait labels shuffled, standardized effect sizes of the functional metrics | Table S6 (S15) |
| 19 | `19_study_area_map.R` | The study area map | Figure 2 |
| 20 | `20_chironomid_sensitivity.R` | The four chironomid taxa: their share, the functional metrics without them against controls, the fourth-corner states without them, and a genus-level split | Table S7 parts A to C (S26, S27, S29) |
| 20b | `20b_figS1_from_record.R` | Draws Figure S1 from step 20's tables and the site matrices; recomputes nothing | Figure S1 |
| 18 | `18_manuscript_assets.R` | Numbers the tables and figures under the pipeline's working numbers, writes the asset manifest and generated captions | all assets (working numbers) |
| 18b | `18b_supplement_renumber.R` | Applies `R/supplement_asset_map.csv`: the paper's supplement numbers, and Table 2 (Tables 3a and 3b merged) | Table 2, Tables S1 to S8, Figure S1 |
| 21 | `21_session_record.R` | R version, platform, flags, package versions and citations of this run | software citations |

## Where each table and figure of the paper lands

| In the paper | File | Written by |
|---|---|---|
| Figure 1 (conceptual figure) | `output/manuscript/figures/Figure_1.*` | drawn by hand; step 18 copies `data/hand_drawn_figures/Figure_1.*` |
| Figure 2 (study area) | `output/manuscript/figures/Figure_2.*` | step 19 |
| Figures 3, 4, 5 | `output/manuscript/figures/Figure_3.*`, `Figure_4.*`, `Figure_5.*` | step 06 |
| Table 1 (predictions and outcomes) | authored in the paper; its outcomes are the screen results of Table 2 and Table S4 part A | |
| Table 2 (metrics meeting the screening criterion) | `output/manuscript/tables_final/Table_2.csv` | step 18b |
| Tables S1 to S8 (pipeline parts) | `output/manuscript/supplement_final/Table_S<number><part>.csv` | step 18b |
| Figure S1 (the four chironomid taxa) | `output/manuscript/supplement_final/Figure_S1.*` | step 20b, copied by step 18b |
| Which pipeline file each supplement table comes from, and the parts authored in the supplement | `output/manuscript/supplement_final/SUPPLEMENT_MANIFEST.csv` | step 18b |
| Graphical abstract | `output/manuscript/figures/Graphical_abstract.*` | see below |

Tables S1 part B and S2 part A have no pipeline source; they are authored in the supplement. Table S8 lays out, in the
supplement, the rows of `output/tables/09b_taxon_tolerance_screen.csv` and `output/tables/09b_tolerance_summary.csv`;
`supplement_final/Table_S8A.csv` and `Table_S8B.csv` hold the full taxon screen and the tolerance table behind it.

Other outputs:

| What | Where |
|---|---|
| Every number behind every table and figure | `output/tables/<step>_*.csv` |
| Assets under the pipeline's working numbers, the asset manifest, and captions generated with every coefficient read from the tables | `output/manuscript/tables/`, `output/manuscript/figures/`, `output/manuscript/ASSET_MANIFEST.csv`, `output/manuscript/CAPTIONS.md` |
| Supporting tables the paper does not number (trait coverage by site, trait-space dimensionality, metric provenance, regional pools, metric distributions, parent rows resolved, family shares) | `output/manuscript/repository_tables/` (see its `INDEX.csv`) |
| Assets of the separate trait-state manuscript | `output/manuscript/trait_state_manuscript/` |
| Working figures, including ones not in the paper | `output/figures/` |
| Software versions and citations | `output/tables/R_packages_and_versions.csv`, `output/manuscript/software/`, `output/logs/session_record_*.txt` |
| Run log | `output/logs/run_<date>_<time>.log` |

Figure 1 is drawn by hand and has no pipeline source. The graphical abstract is drawn by a separate script, outside
this repository, from this pipeline's screen table; step 18 skips it in this repository and says so in the log, and
the shipped copy is the file that script wrote. The captions in `CAPTIONS.md` are generated with every coefficient read
from the tables, as a check on the authored captions in the paper; they use the working numbers.

## Data

`data/` holds every file the pipeline reads. It mirrors the layout of the project's own data folder, so the
scripts run unchanged in both places. `data/DATA_DICTIONARY.md` documents every file and every column: what it
is, where its values come from, which script reads or writes it, and its units. In brief:

- the community table (site by taxon: abundance, density and biomass as the mean over the processed replicate
  samples) and the two site-by-taxon matrices (90 sites by 221 taxa before parent-row resolution, 215 after), before
  and after R/00c;
- the trait affinity matrix (576 taxa by 41 states in 11 trait groups, from the CONUS trait database of
  Twardochleb et al. 2021, with each proxy taxon's source genera and tier documented under `data/Trait_Data/`);
- the biological condition gradient: each site's Chessie B-IBI (median of 100 rarefied subsamples) and its
  rating, as released by the Interstate Commission on the Potomac River Basin, with the class cut-points;
- the site totals and site tolerance index as medians of the replicate values (`Step2_SiteReplicateMedians.csv`);
- the site table of screened metrics, before and after resolution, with the provenance tables the data steps write;
- the tolerance values behind the Hilsenhoff Biotic Index, the taxonomy used to place taxa, the site coordinates for
  the map, and the exported Figure 1 (`data/hand_drawn_figures/`).

## Author and citation

Sergio A. Sabat-Bonilla, Department of Entomology, Virginia Polytechnic Institute and State University,
Blacksburg, Virginia, USA. Cite the paper and this repository (see `CITATION.cff`). Code is released under the
MIT license and data under CC BY 4.0 (see `LICENSE`); values that originate with a third party keep their
source's terms.

## References

- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B (Methodological), 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Cuffney, T. F., Bilger, M. D., & Haigler, A. M. (2007). Ambiguous taxa: Effects on the characterization and interpretation of invertebrate assemblages. Journal of the North American Benthological Society, 26(2), 286–307. https://doi.org/10.1899/0887-3593(2007)26[286:ATEOTC]2.0.CO;2
- James, A., Pence, R., Pond, G., Schoenholtz, S., Timpano, A., Zipper, C., & Entrekin, S. (2022). Taxon and trait-based sampling curves can be used as a tool for assessing impairment in salinized headwater streams. Ecological Indicators, 139, 108942. https://doi.org/10.1016/j.ecolind.2022.108942
- Maire, E., Grenouillet, G., Brosse, S., & Villéger, S. (2015). How many dimensions are needed to accurately assess functional diversity? A pragmatic approach for assessing the quality of functional spaces. Global Ecology and Biogeography, 24(6), 728–740. https://doi.org/10.1111/geb.12299
- Smith, Z. M., Buchanan, C., & Nagel, A. (2017). Refinement of the basin-wide index of biotic integrity for non-tidal streams and wadeable rivers in the Chesapeake Bay watershed. Interstate Commission on the Potomac River Basin.
- Twardochleb, L., Hiltner, E., Pyne, M., & Zarnetske, P. (2021). Freshwater insects CONUS: A database of freshwater insect occurrences and traits for the contiguous United States. Global Ecology and Biogeography, 30(4), 826–841. https://doi.org/10.1111/geb.13257
