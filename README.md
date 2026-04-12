[README.md](https://github.com/user-attachments/files/26649491/README.md)
# criminal-history-expiration-date

This repository provides complete replication materials for the manuscript, including data structure, analytic code, and scripts to reproduce all tables and figures:

**Decay in Recidivism Risk**

# criminal-history-expiration-date

This repository provides complete replication materials for the manuscript:

**Decay in Recidivism Risk**

---

## Overview

This project examines how the predictive value of prior criminal charges changes over time. Using large-scale administrative criminal history data, the analysis evaluates three core questions:

- Does recidivism risk decline as time since a prior charge increases?
- Is the rate of decline consistent across jurisdictions and subgroups?
- At what point does predicted recidivism risk converge toward the general population base rate?

These scripts reproduce all descriptive analyses, statistical models, tables, figures, and appendices reported in the manuscript.

---

## Repository Structure


```text
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
│   │   └── 2026.03.30 WA decay dataset.csv
│   └── base_rates/
│       └── state_base_rates.csv
└── output/
```


---

## Data Availability

Due to the size and administrative nature of the original dataset, the full raw data are not publicly distributed. This repository includes processed subsamples and supporting files necessary to reproduce all analyses reported in the manuscript.

Researchers seeking access to the original data should contact the corresponding author.

---

## Inputs

### Distributed subsamples

Place the 20 distributed subsample files in:

`data/subsamples/`

Expected file names:

- `subsample_00_raw.csv`
- `subsample_01_raw.csv`
- ...
- `subsample_19_raw.csv`

---

### Washington dataset

Place the Washington-specific dataset in:

`data/wa/`


Expected file name:

- `2026.03.30 WA decay dataset.csv`

---

### Base-rate file

Place the base-rate file in:

`data/base_rates/`


Expected file name:

- `state_base_rates.csv`

Required columns:

- `State`
- `base_rate`

---

## Getting Started

1. Install required packages:

```r
source("requirements.R")

- `state_base_rates.csv`

Expected columns:

- `State`
- `base_rate`

##Run the full analysis:

source("run_all.R")
```

##Run Order (Manual)

source("01_register_subsamples.R")
source("02_descriptives.R")
source("03_H1_decay.R")
source("04_H2_invariance.R")
source("05_H3_convergence.R")

## Output Mapping

| Script              | Main outputs              | Manuscript components                 |
| ------------------- | ------------------------- | ------------------------------------- |
| 02_descriptives.R   | Descriptive tables        | Table 1, Appendix A                   |
| 03_H1_decay.R       | Decay models and figures  | Table 2, Figures 3–4, Appendix C      |
| 04_H2_invariance.R  | State and subgroup models | Tables 3–5, Figures 5–7, Appendix D/E |
| 05_H3_convergence.R | Crossing-year summaries   | Table 6, Figure 8                     |


## Computational Notes

The distributed subsample files are the starting point for the public replication package. Models are estimated within subsamples and pooled where appropriate. The code is structured to prioritize clarity and reproducibility for reviewers and readers.

## Contact

Zachary Hamilton  
zhamilton@unomaha.edu


---


