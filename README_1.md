# criminal-history-expiration-date

This repository provides replication materials for the manuscript:

**Criminal history's expiration date: Evidence of prediction decay**
Hamilton, Z., Kigerl, A., Tostlebe, J. J., & Ursino, J. — *Criminology*, the flagship, peer-reviewed journal of the American Society of Criminology.

---

## Overview

This project examines how the predictive value of prior criminal charges changes over time. Using large-scale administrative criminal history data, the analysis evaluates three core questions:

- Does recidivism risk decline as time since a prior charge increases? (H1)
- Is the rate of decline consistent across jurisdictions, demographic groups, offense type, incarceration exposure, and historical period? (H2)
- At what point does predicted recidivism risk converge toward the general population base rate? (H3)

The scripts in this repository reproduce all descriptive analyses, statistical models, tables, figures, and appendices reported in the manuscript, including the historical-period comparison and the one-event-per-person-year robustness check added during the revise-and-resubmit (R&R) process.

---

## Repository Structure

```
.
├── README.md
├── requirements.R
├── run_all.R
├── 01_register_subsamples.R
├── 02_descriptives.R
├── 03_H1_decay.R
├── 04_H2_invariance.R
├── 05_H3_convergence.R
├── 06_H2_historical_period.R
├── 07_register_sensitivity_subsamples.R
├── 08_H1_sensitivity_one_event_per_year.R
├── data/
│   ├── subsamples/
│   │   ├── subsample_00_raw.csv
│   │   ├── ...
│   │   └── subsample_19_raw.csv
│   ├── wa/
│   │   └── WA_decay_dataset.csv
│   ├── base_rates/
│   │   └── Baserate_for_modeling.csv
│   └── sensitivity/
│       ├── cohort_after2000.csv.gz
│       └── merged_state_events.csv.gz
└── output/
```

---

## Data

Replication datasets are available at:

https://uofnebraska-my.sharepoint.com/:f:/r/personal/49639061_nebraska_edu/Documents/Washington%20DOC/Decay/Manuscripts/Decay%20LifeCourse%20Manuscript/Documents/Data%20and%20scripts%20for%20Git%20Hub/data?csf=1&web=1&e=J8cr1a

After downloading, place the files in the following folders:

```
data/subsamples/
data/wa/
data/base_rates/
data/sensitivity/
```

File names must match those listed above for the scripts to run correctly.

### Core files (required for 01-05)

| File | Required columns |
|---|---|
| `subsamples/subsample_00_raw.csv` ... `subsample_19_raw.csv` | `ResearchID`, `State`, `LookBack`, `Recid` |
| `wa/WA_decay_dataset.csv` | `ResearchID`, `Recid`, `LookBack`, `Male`, `RaceEthnicity2`, `AgeCurveGrouped`, `CharType`, `PriorIncarYearsCount` |
| `base_rates/Baserate_for_modeling.csv` | `State` + either `base_rate`, or `TotalOffenders`/`TotalOffenders (NumeratorData)`, or `Denominator_Annual_State_Population`/`Numerator_Annual_Recidivism_Events` |

### Post-R&R sensitivity files (required for 06-08)

> **Note to maintainers:** these two files did not exist in the original data drop. They're exported from the SPSS extracts used during the R&R analysis via `export_sensitivity_data_from_spss.R` (a one-time local utility, not part of the numbered pipeline -- see that script's header) and added to the SharePoint folder linked above. Until they're added, `06_H2_historical_period.R` and `07_register_sensitivity_subsamples.R` will fail their file-existence check with a clear error message pointing at the missing path.

| File | Required columns | Used by |
|---|---|---|
| `sensitivity/cohort_after2000.csv.gz` | `ResearchID`, `State`, `LookBack`, `Recid`, `After2000` (0 = prior event before 2001, 1 = prior event in 2001+) | `06_H2_historical_period.R` |
| `sensitivity/merged_state_events.csv.gz` | `ResearchID`, `State`, `OffDate`, `LookBack`, `Recid` | `07_register_sensitivity_subsamples.R` |

`merged_state_events.csv.gz` spans the same five states as the primary analysis (AK, AZ, FL, TX, WA). `07_register_sensitivity_subsamples.R` validates this on read and will stop with a clear error if an unexpected state value appears.

**Why `.gz`:** a real run of the export script showed `cohort_after2000.csv` alone is on the order of hundreds of millions of rows -- comparable in scale to the manuscript's full person-year observation count -- and an uncompressed export ran out of local disk space partway through writing. Both files are written gzip-compressed instead. `06_H2_historical_period.R`, `07_register_sensitivity_subsamples.R`, and `check_data_schema.R` all read `.csv.gz` automatically (falling back to a plain `.csv` of the same name if that's what's present instead) via `data.table::fread()`, which requires the `R.utils` package (included in `requirements.R`).

One caveat worth knowing before running any of those three scripts against the real files: `fread()` has no partial/streaming read for `.gz` -- every read, even a small preview, fully decompresses the file to a scratch temp file first, which briefly needs about as much free space as the *uncompressed* file would, wherever R's `tempdir()` points. Each of those three scripts exposes a `FREAD_TMPDIR` variable near the top that defaults to `tempdir()` and can be pointed at a drive with more headroom if needed.

---

## Getting Started

1. Install required packages:

```r
source("requirements.R")
```

2. Run the full analysis:

```r
source("run_all.R")
```

Successful execution will produce output files in the `/output/` directory.

---

## Run Order (Manual)

Primary analysis (original submission):

```r
source("01_register_subsamples.R")
source("02_descriptives.R")
source("03_H1_decay.R")
source("04_H2_invariance.R")
source("05_H3_convergence.R")
```

Post-R&R additions (historical period + one-event-per-year robustness check):

```r
source("06_H2_historical_period.R")
source("07_register_sensitivity_subsamples.R")
source("08_H1_sensitivity_one_event_per_year.R")
```

`06` depends only on `data/sensitivity/cohort_after2000.csv`. `08` depends on `07`'s output and must run after it.

---

## Output Mapping

| Script | Main outputs | Manuscript components |
|------|-------------|----------------------|
| 02_descriptives.R | Descriptive tables | Table 1, Appendix A |
| 03_H1_decay.R | Decay models and figures | Table 2, Figures 3-4, Appendix C |
| 04_H2_invariance.R | State and subgroup models | Tables 3-5, Figures 5-7, Appendix D/E |
| 05_H3_convergence.R | Crossing-year summaries | Table 7, Figure 9 |
| 06_H2_historical_period.R | Historical-period GLMM estimates | Table 6, Figure 8 |
| 07_register_sensitivity_subsamples.R | One-event-per-person-year subsamples | Data prep for sensitivity check (Appendix, R&R response) |
| 08_H1_sensitivity_one_event_per_year.R | Pooled decay estimate, one event/person-year (sharded/parallel-fit — see Computational Notes) | Sensitivity check vs. Table 2 (R&R response) |

**Table/figure numbering note:** During the R&R, the editor requested that the historical-period comparison be elevated from a supplementary sensitivity analysis into a primary H2 subgroup analysis. This moved it into the main text as Table 6 / Figure 8, which pushed the convergence analysis (originally Table 6 / Figure 8 in the first-submission numbering this README used to follow) to its current position as **Table 7 / Figure 9**. The mapping above reflects the accepted-manuscript numbering; double-check against your final proofs before publishing this repo.

---

## Computational Notes

The distributed subsample files are the starting point for the replication workflow. Models are estimated within subsamples/partitions and pooled via random-effects meta-analysis (`metafor::rma.uni`) where appropriate. External data hosting is used due to file size constraints.

`06_H2_historical_period.R` follows the same 20-partition, meta-analytic-pooling design as the primary H1/H2 models, for the same computational-tractability reasons, and caches partition-level model results to disk so an interrupted run can resume without refitting completed partitions.

`07_register_sensitivity_subsamples.R` also fixes a cross-state `ResearchID` collision present in the merged five-state extract used for the one-event-per-year sensitivity check (IDs are prefixed with `State` before any person-level grouping). This collision is isolated to that sensitivity dataset — the distributed `subsample_00-19_raw.csv` files used by the primary pipeline (01-05) were already state-safe and are unaffected.

**`OffDate` format:** confirmed against a real extract that date-like columns in this data family (`OffDate`, `BackBegin`, `TrackBegin`, `TrackEnd`, `BackEnd`) are not consistently full parseable dates — some are stored as a bare 4-digit year (e.g. `2015`). `as.Date("2015")` throws a hard parse error rather than silently returning `NA`, so `07_register_sensitivity_subsamples.R` detects and handles the bare-year case explicitly before falling back to `as.Date()` parsing. If `merged_state_events.csv`'s `OffDate` turns out to already be a full date string, both forms are handled by the same code path — no reconfiguration needed either way.

**`08_H1_sensitivity_one_event_per_year.R` — sharding, parallel fitting, and checkpointing:** each of the 20 subsamples this script processes is itself on the order of ~597M rows / ~16M individuals — too large for a single `glmer()` call — so 08 shards each subsample into `N_SHARDS_PER_SUBSAMPLE` (default 20) smaller pieces, fits each shard's model independently, and pools shard-level estimates via `metafor::rma.uni` into a subsample-level estimate before that feeds into the same subsample-level pooling step used elsewhere in this repo. Shard fits run in parallel across up to `N_PARALLEL_WORKERS` R worker processes (default `min(2, cores - 1)`, base-R `parallel`/PSOCK clusters, Windows-safe), with each worker capped to a single BLAS thread so the workers don't oversubscribe the machine's cores. Checkpointing happens at both the shard level (`output/08_H1_sensitivity_one_event_per_year/shard_checkpoints/`) and the subsample level (`.../subsample_pooled_checkpoints/`), so an interrupted run — including a hard crash — resumes without refitting already-finished work; a truncated/corrupt checkpoint file is detected and deleted automatically rather than crashing the run. On a real full-data run this completed in a total of several days of wall-clock time across all 20 subsamples (individual shard fits ran roughly 60–110 minutes each). Full diagnostic detail is written per shard (`sensitivity_shard_coefficients.csv` / `_fit_stats.csv` / `_random_effects.csv`), per subsample (`sensitivity_subsample_coefficients.csv`), and as the final pooled result (`sensitivity_pooled_coefficients.csv`, `sensitivity_vs_primary_comparison.csv`). Only raise `N_PARALLEL_WORKERS` above the default after confirming a run completes without memory pressure at the default, and drop `N_SHARDS_PER_SUBSAMPLE`/`N_PARALLEL_WORKERS` toward 1 to fall back to fully serial, more conservative fitting if instability persists.

**Real-data file size:** a single `subsample_XX_raw.csv` partition was measured at ~23 GB in the actual data (vs. the illustrative small files used to smoke-test this pipeline's logic). `cohort_after2000.csv` was confirmed via a real export run to be on the order of hundreds of millions of rows uncompressed — comparable in scale to the manuscript's full person-year observation count, not a small summary file. Do not attempt to run this pipeline against real data anywhere with limited disk/memory — this is exactly why the 20-partition, meta-analytic-pooling design exists in the first place, and why the sensitivity exports are gzip-compressed (see "Post-R&R sensitivity files" above, including the `fread()`/scratch-space caveat for reading `.gz` files back).

**`CharType` blanks in `WA_decay_dataset.csv`:** a real extract of this file confirmed a substantial share of rows have a blank (not missing/`NA`) `CharType` value. `04_H2_invariance.R` now normalizes `""` to `NA` and codes it as **Non-Violent** in `CrimeType2`, per author confirmation, rather than silently dropping those rows from the offense-type H2 model as the original filter (`is.na(CharType) | CharType %in% valid_types`, combined with a strict `CharType == "Violent"` check) would have done — blank was neither `NA` nor one of the five named categories, so it fell through the filter and those individuals were dropped from Table 5 / Figure 7 without any error or warning. Confirmed via a synthetic regression test with ~38% blank `CharType` rows: all are now retained and classified Non-Violent.

---

## Contact

Zachary Hamilton
zhamilton@unomaha.edu
