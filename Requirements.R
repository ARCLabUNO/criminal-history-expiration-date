
---
  
  # requirements.R
  
  ```r
packages <- c(
  "data.table",
  "lme4",
  "broom.mixed",
  "metafor",
  "ggplot2",
  "scales",
  "patchwork",
  "openxlsx"
)

to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0L) {
  install.packages(to_install)
}