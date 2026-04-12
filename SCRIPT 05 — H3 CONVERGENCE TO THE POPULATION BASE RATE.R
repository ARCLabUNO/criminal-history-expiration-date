############################################################
# SCRIPT 05 — H3 CONVERGENCE TO THE POPULATION BASE RATE
#
# PURPOSE
#   Identify the lookback year at which predicted recidivism
#   probability converges toward the state-specific population
#   base rate.
#
# INPUT
#   /output/04_H2_invariance/H2_state_interaction_curves.csv
#   /data/base_rates/Baserate_for_modeling.csv
#
# OUTPUT
#   /output/05_H3_convergence/
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(openxlsx)
})

REPO_DIR <- "."
DATA_DIR <- file.path(REPO_DIR, "data")
H2_DIR <- file.path(REPO_DIR, "output", "04_H2_invariance")
BASE_RATE_DIR <- file.path(DATA_DIR, "base_rates")
OUT_DIR <- file.path(REPO_DIR, "output", "05_H3_convergence")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CURVES_FILE <- file.path(H2_DIR, "H2_state_interaction_curves.csv")
BASE_RATE_FILE <- file.path(BASE_RATE_DIR, "Baserate_for_modeling.csv")

REQUIRED_CURVE_COLUMNS <- c("subset", "LookBack", "State", "p_hat", "RR", "model")
REQUIRED_BASE_RATE_COLUMNS <- c("State", "base_rate")

EXPECTED_STATES <- c("AK", "AZ", "FL", "TX", "WA")

LB_MIN <- 0
LB_MAX <- 40
DELTAS <- c(1.0, 2.5)
BASE_FAMILY <- "Times New Roman"

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush.console()
}

validate_columns <- function(dt, required_cols, file_path) {
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in file:\n", file_path, "\nMissing columns: ", paste(missing_cols, collapse = ", "))
  }
}

find_crossing_true <- function(LB, p, base_rate, lb_min = 0, lb_max = Inf) {
  base_rate <- base_rate[1]
  o <- order(LB)
  LB <- LB[o]
  p <- p[o]
  
  keep <- LB >= lb_min & LB <= lb_max
  LB <- LB[keep]
  p <- p[keep]
  
  if (length(LB) < 2L) return(NA_real_)
  if (!is.finite(base_rate)) return(NA_real_)
  if (!is.finite(p[1])) return(NA_real_)
  if (p[1] <= base_rate) return(NA_real_)
  
  j <- which(p <= base_rate)[1]
  if (is.na(j) || j == 1L) return(NA_real_)
  
  x0 <- LB[j - 1L]
  y0 <- p[j - 1L]
  x1 <- LB[j]
  y1 <- p[j]
  
  if (!is.finite(y0) || !is.finite(y1)) return(NA_real_)
  if (isTRUE(all.equal(y0, y1))) return(as.numeric(x1))
  
  t <- (base_rate - y0) / (y1 - y0)
  x <- x0 + (x1 - x0) * t
  if (!is.finite(x)) return(NA_real_)
  as.numeric(x)
}

fmt_ci <- function(mid, lo, hi, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), mid, lo, hi)
}

if (!file.exists(CURVES_FILE)) stop("Missing H2 state interaction curves file:\n", CURVES_FILE)
if (!file.exists(BASE_RATE_FILE)) stop("Missing state base-rate file:\n", BASE_RATE_FILE)

curves <- fread(CURVES_FILE, showProgress = FALSE)
validate_columns(curves, REQUIRED_CURVE_COLUMNS, CURVES_FILE)

base_rates <- fread(BASE_RATE_FILE, showProgress = FALSE)
validate_columns(base_rates, REQUIRED_BASE_RATE_COLUMNS, BASE_RATE_FILE)

curves[, `:=`(
  subset = as.character(subset),
  State = trimws(as.character(State)),
  LookBack = as.numeric(LookBack),
  p_hat = as.numeric(p_hat),
  RR = as.numeric(RR),
  model = as.character(model)
)]

base_rates[, `:=`(
  State = trimws(as.character(State)),
  base_rate = as.numeric(base_rate)
)]

curves <- curves[State %in% EXPECTED_STATES & LookBack >= LB_MIN & LookBack <= LB_MAX]
base_rates <- base_rates[State %in% EXPECTED_STATES]

curves <- merge(curves, base_rates, by = "State", all.x = TRUE)
if (curves[is.na(base_rate), .N] > 0L) stop("One or more states in the curves file do not have a matching base rate.")

crossings_dt <- curves[
  ,
  .(
    crossing_year = find_crossing_true(LookBack, p_hat, base_rate[1], LB_MIN, LB_MAX),
    base_rate = base_rate[1]
  ),
  by = .(subset, State)
][order(subset, State)]

crossing_state_summary_dt <- crossings_dt[
  ,
  .(
    n_subsets_with_crossing = sum(!is.na(crossing_year)),
    crossing_median = if (all(is.na(crossing_year))) NA_real_ else median(crossing_year, na.rm = TRUE),
    crossing_q025 = if (all(is.na(crossing_year))) NA_real_ else quantile(crossing_year, 0.025, na.rm = TRUE),
    crossing_q975 = if (all(is.na(crossing_year))) NA_real_ else quantile(crossing_year, 0.975, na.rm = TRUE),
    crossing_mean = if (all(is.na(crossing_year))) NA_real_ else mean(crossing_year, na.rm = TRUE),
    crossing_sd = if (all(is.na(crossing_year))) NA_real_ else sd(crossing_year, na.rm = TRUE),
    base_rate = unique(base_rate)[1]
  ),
  by = State
][order(State)]

within_mean_metrics_dt <- rbindlist(
  lapply(DELTAS, function(delta) {
    crossings_dt[
      ,
      {
        x <- crossing_year[!is.na(crossing_year)]
        if (length(x) < 2L) {
          .(mean_crossing = NA_real_, prop_within = NA_real_, max_deviation = NA_real_, n_states = length(x))
        } else {
          mu <- mean(x)
          .(
            mean_crossing = mu,
            prop_within = mean(abs(x - mu) <= delta),
            max_deviation = max(abs(x - mu))
          )
        }
      },
      by = subset
    ][, delta_years := delta]
  }),
  use.names = TRUE, fill = TRUE
)

within_mean_summary_dt <- within_mean_metrics_dt[
  ,
  .(
    mean_crossing_median = median(mean_crossing, na.rm = TRUE),
    mean_crossing_q025 = quantile(mean_crossing, 0.025, na.rm = TRUE),
    mean_crossing_q975 = quantile(mean_crossing, 0.975, na.rm = TRUE),
    prop_within_median = median(prop_within, na.rm = TRUE),
    prop_within_q025 = quantile(prop_within, 0.025, na.rm = TRUE),
    prop_within_q975 = quantile(prop_within, 0.975, na.rm = TRUE),
    max_deviation_median = median(max_deviation, na.rm = TRUE),
    max_deviation_q025 = quantile(max_deviation, 0.025, na.rm = TRUE),
    max_deviation_q975 = quantile(max_deviation, 0.975, na.rm = TRUE)
  ),
  by = delta_years
][order(delta_years)]

table6_dt <- crossing_state_summary_dt[
  ,
  .(
    State,
    BaseRate = sprintf("%.3f", base_rate),
    CrossingYear = fmt_ci(crossing_median, crossing_q025, crossing_q975, 2),
    MeanSD = sprintf("%.2f (%.2f)", crossing_mean, crossing_sd),
    N_Subsets = n_subsets_with_crossing
  )
]

overall_crossing_summary_dt <- within_mean_summary_dt[
  ,
  .(
    delta_years,
    mean_crossing_median,
    mean_crossing_q025,
    mean_crossing_q975,
    prop_within_median,
    prop_within_q025,
    prop_within_q975,
    max_deviation_median,
    max_deviation_q025,
    max_deviation_q975
  )
]

curve_mean_dt <- curves[
  ,
  .(
    p_mean = mean(p_hat, na.rm = TRUE),
    base_rate = unique(base_rate)[1]
  ),
  by = .(State, LookBack)
]

curve_mean_dt <- merge(curve_mean_dt, crossing_state_summary_dt[, .(State, crossing_median)], by = "State", all.x = TRUE)

fig8 <- ggplot(curve_mean_dt, aes(x = LookBack, y = p_mean)) +
  geom_line(linewidth = 0.9) +
  geom_hline(aes(yintercept = base_rate), linetype = "dashed") +
  geom_vline(aes(xintercept = crossing_median), linetype = "dotted", na.rm = TRUE) +
  facet_wrap(~ State, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "State-specific convergence of recidivism probability to the base rate",
    x = "Lookback year",
    y = "Predicted probability"
  ) +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

ggsave(file.path(OUT_DIR, "H3_figure8_state_crossing.png"), fig8, width = 10, height = 7, dpi = 300)

fwrite(curves, file.path(OUT_DIR, "H3_curves_with_base_rates.csv"))
fwrite(crossings_dt, file.path(OUT_DIR, "H3_crossing_years_by_subset_state.csv"))
fwrite(crossing_state_summary_dt, file.path(OUT_DIR, "H3_crossing_years_state_summary.csv"))
fwrite(within_mean_metrics_dt, file.path(OUT_DIR, "H3_within_mean_metrics_by_subset.csv"))
fwrite(within_mean_summary_dt, file.path(OUT_DIR, "H3_within_mean_equivalence_summary.csv"))
fwrite(table6_dt, file.path(OUT_DIR, "H3_table6_crossing_summary.csv"))
fwrite(overall_crossing_summary_dt, file.path(OUT_DIR, "H3_overall_convergence_summary.csv"))

wb <- createWorkbook()
addWorksheet(wb, "CrossingsBySubsetState"); writeDataTable(wb, "CrossingsBySubsetState", crossings_dt)
addWorksheet(wb, "StateSummary"); writeDataTable(wb, "StateSummary", crossing_state_summary_dt)
addWorksheet(wb, "WithinMeanBySubset"); writeDataTable(wb, "WithinMeanBySubset", within_mean_metrics_dt)
addWorksheet(wb, "WithinMeanSummary"); writeDataTable(wb, "WithinMeanSummary", within_mean_summary_dt)
addWorksheet(wb, "Table6"); writeDataTable(wb, "Table6", table6_dt)
saveWorkbook(wb, file.path(OUT_DIR, "H3_results_workbook.xlsx"), overwrite = TRUE)

log_line("H3 analyses complete.")
print(crossing_state_summary_dt)
print(within_mean_summary_dt)
