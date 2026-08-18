# criminal-history-expiration-date

This repository provides replication materials for the manuscript:

**Criminal history’s expiration date: Universal evidence of decay**

---

## Overview

This project examines how the predictive value of prior criminal charges changes over time. Using large-scale administrative criminal history data, the analysis evaluates three core questions:

- Does recidivism risk decline as time since a prior charge increases?
- Is the rate of decline consistent across jurisdictions and subgroups?
- At what point does predicted recidivism risk converge toward the general population base rate?

The scripts in this repository reproduce all descriptive analyses, statistical models, tables, figures, and appendices reported in the manuscript.

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
├── data/
│   ├── subsamples/
│   │   ├── subsample_00_raw.csv
│   │   ├── ...
│   │   └── subsample_19_raw.csv
│   ├── wa/
│   │   └── WA_decay_dataset.csv
│   └── base_rates/
│       └── Baserate_for_modeling.csv
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
```

File names must match those listed above for the scripts to run correctly.

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

```r
source("01_register_subsamples.R")
source("02_descriptives.R")
source("03_H1_decay.R")
source("04_H2_invariance.R")
source("05_H3_convergence.R")
```

---

## Output Mapping

| Script | Main outputs | Manuscript components |
|------|-------------|----------------------|
| 02_descriptives.R | Descriptive tables | Table 1, Appendix A |
| 03_H1_decay.R | Decay models and figures | Table 2, Figures 3–4, Appendix C |
| 04_H2_invariance.R | State and subgroup models | Tables 3–5, Figures 5–7, Appendix D/E |
| 05_H3_convergence.R | Crossing-year summaries | Table 6, Figure 8 |

---

## Computational Notes

The distributed subsample files are the starting point for the replication workflow. Models are estimated within subsamples and pooled where appropriate. External data hosting is used due to file size constraints.

---

## Contact

Zachary Hamilton  
zhamilton@unomaha.edu
