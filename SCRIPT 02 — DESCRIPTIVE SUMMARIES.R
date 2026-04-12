############################################################
# SCRIPT 02 — DESCRIPTIVE SUMMARIES
#
# PURPOSE
#   Produce descriptive summaries from the distributed
#   subsample files.
#
# INPUT
#   /data/subsamples/subsample_00_raw.csv ... subsample_19_raw.csv
#
# REQUIRED COLUMNS
#   ResearchID, State, LookBack, Recid
#
# OPTIONAL COLUMNS
#   BackBegin, OffDate
#
# OUTPUT
#   /output/02_descriptives/
#     - table1_id_level_summary.csv
#     - appendixA_person_year_example.csv
#     - appendixA_population_style_descriptives.csv
#     - sequences_per_person_distribution_raw.csv
#     - sequences_per_person_grouped.csv
#     - recid_events_per_person_grouped.csv
#     - total_sequences_by_state.csv
#     - person_year_recap_by_state.csv
#     - rate_group_summary.csv
#     - diagnostics_descriptives.csv
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
})

REPO_DIR <- "."
SUBSAMPLE_DIR <- file.path(REPO_DIR, "data", "subsamples")
OUT_DIR <- file.path(REPO_DIR, "output", "02_descriptives")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SUBSAMPLE_IDS <- sprintf("%02d", 0:19)
INPUT_FILES <- file.path(SUBSAMPLE_DIR, sprintf("subsample_%s_raw.csv", SUBSAMPLE_IDS))

REQUIRED_COLUMNS <- c("ResearchID", "State", "LookBack", "Recid")
EXPECTED_STATES <- c("AK", "AZ", "FL", "TX", "WA")

LOW_RATE_MAX <- 0.125
MODERATE_RATE_MAX <- 0.35

data.table::setDTthreads(max(1L, parallel::detectCores(logical = FALSE) - 1L))

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush.console()
}

stop_if_missing_files <- function(paths) {
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files) > 0L) {
    stop("The following required files are missing:\n", paste(missing_files, collapse = "\n"))
  }
}

validate_columns <- function(dt, file_path) {
  missing_cols <- setdiff(REQUIRED_COLUMNS, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in file:\n", file_path, "\nMissing columns: ", paste(missing_cols, collapse = ", "))
  }
}

group_0_1to4_5to14_15plus <- function(x) {
  fcase(
    is.na(x), NA_character_,
    x == 0, "0",
    x >= 1 & x <= 4, "1-4",
    x >= 5 & x <= 14, "5-14",
    x >= 15, "15+",
    default = NA_character_
  )
}

summarize_id_block <- function(D) {
  data.table(
    N_IDs = as.integer(uniqueN(D$ResearchID)),
    events_mean = as.numeric(mean(D$recid_events, na.rm = TRUE)),
    events_sd = as.numeric(sd(D$recid_events, na.rm = TRUE)),
    events_median = as.numeric(median(D$recid_events, na.rm = TRUE)),
    lookback_mean = as.numeric(mean(D$lookback_median, na.rm = TRUE)),
    lookback_sd = as.numeric(sd(D$lookback_median, na.rm = TRUE)),
    lookback_median = as.numeric(median(D$lookback_median, na.rm = TRUE))
  )
}

log_line("Checking descriptive input files.")
stop_if_missing_files(INPUT_FILES)

id_parts <- vector("list", length(INPUT_FILES))
seq_parts <- vector("list", length(INPUT_FILES))
diag_parts <- vector("list", length(INPUT_FILES))
py_parts <- vector("list", length(INPUT_FILES))
example_parts <- vector("list", length(INPUT_FILES))

for (i in seq_along(INPUT_FILES)) {
  f <- INPUT_FILES[i]
  sid <- SUBSAMPLE_IDS[i]
  
  log_line("Reading subsample %s: %s", sid, basename(f))
  
  dt <- fread(f, showProgress = FALSE)
  validate_columns(dt, f)
  
  has_backbegin <- "BackBegin" %in% names(dt)
  has_offdate <- "OffDate" %in% names(dt)
  event_date_var <- if (has_backbegin) "BackBegin" else if (has_offdate) "OffDate" else NA_character_
  
  dt[, ResearchID := as.character(ResearchID)]
  dt[, State := trimws(as.character(State))]
  dt[, LookBack := suppressWarnings(as.integer(LookBack))]
  dt[, Recid := suppressWarnings(as.integer(Recid))]
  
  if (!is.na(event_date_var)) {
    dt[, (event_date_var) := as.character(get(event_date_var))]
  }
  
  dt <- dt[
    !is.na(ResearchID) &
      !is.na(State) &
      !is.na(LookBack) &
      !is.na(Recid)
  ]
  
  dt <- dt[State %chin% EXPECTED_STATES]
  
  diag_parts[[i]] <- data.table(
    subsample = sid,
    rows_after_filter = nrow(dt),
    n_unique_ids = uniqueN(dt$ResearchID),
    lookback_min = suppressWarnings(min(dt$LookBack, na.rm = TRUE)),
    lookback_max = suppressWarnings(max(dt$LookBack, na.rm = TRUE)),
    event_date_field = ifelse(is.na(event_date_var), "none", event_date_var)
  )
  
  py_parts[[i]] <- dt[
    ,
    .(
      person_year_obs = .N,
      total_recid_events = sum(Recid, na.rm = TRUE),
      py_with_recid = sum(Recid > 0L, na.rm = TRUE)
    ),
    by = .(State)
  ]
  
  if (!is.na(event_date_var)) {
    seq_parts[[i]] <- dt[
      LookBack == 1L,
      .(n_sequences = uniqueN(get(event_date_var))),
      by = .(State, ResearchID)
    ]
  } else {
    seq_parts[[i]] <- dt[
      LookBack == 1L,
      .(n_sequences = .N),
      by = .(State, ResearchID)
    ]
  }
  
  id_parts[[i]] <- dt[
    ,
    .(
      recid_events = sum(Recid, na.rm = TRUE),
      lookback_median = as.numeric(median(LookBack, na.rm = TRUE)),
      exposure_years = suppressWarnings(max(LookBack, na.rm = TRUE))
    ),
    by = .(State, ResearchID)
  ]
  
  id_parts[[i]][, exposure_years := fifelse(is.finite(exposure_years) & exposure_years > 0, exposure_years, NA_real_)]
  id_parts[[i]][, recid_rate := recid_events / exposure_years]
  
  example_ids <- head(unique(dt$ResearchID), 2L)
  
  example_dt <- dt[
    ResearchID %in% example_ids,
    .(ResearchID, State, LookBack, Recid)
  ][order(ResearchID, LookBack)]
  
  if (!is.na(event_date_var)) {
    example_dt <- dt[
      ResearchID %in% example_ids,
      .(ResearchID, State, LookBack, Recid, EventDate = get(event_date_var))
    ][order(ResearchID, EventDate, LookBack)]
    
    example_dt[, Sequence := rleidv(.SD, cols = c("ResearchID", "EventDate")), by = ResearchID]
    example_dt[, EventDate := NULL]
    setcolorder(example_dt, c("ResearchID", "State", "Sequence", "LookBack", "Recid"))
  } else {
    example_dt[, Sequence := 1L]
    setcolorder(example_dt, c("ResearchID", "State", "Sequence", "LookBack", "Recid"))
  }
  
  example_parts[[i]] <- example_dt
  
  rm(dt, example_dt)
  gc()
}

id_all <- rbindlist(id_parts, use.names = TRUE, fill = TRUE)
seq_all <- rbindlist(seq_parts, use.names = TRUE, fill = TRUE)
diag_all <- rbindlist(diag_parts, use.names = TRUE, fill = TRUE)
py_all <- rbindlist(py_parts, use.names = TRUE, fill = TRUE)
example_all <- rbindlist(example_parts, use.names = TRUE, fill = TRUE)

setkey(id_all, State, ResearchID)
setkey(seq_all, State, ResearchID)
id_all <- seq_all[id_all]
id_all[is.na(n_sequences), n_sequences := 0L]

table1_overall <- cbind(State = "ALL", summarize_id_block(id_all))
table1_state <- id_all[, summarize_id_block(.SD), by = State]
table1_dt <- rbindlist(list(table1_overall, table1_state), use.names = TRUE, fill = TRUE)
table1_dt[, events_mean_sd := sprintf("%.3f (%.3f)", events_mean, events_sd)]
table1_dt[, lookback_mean_sd := sprintf("%.3f (%.3f)", lookback_mean, lookback_sd)]

setcolorder(
  table1_dt,
  c("State", "N_IDs", "events_mean_sd", "events_median", "lookback_mean_sd", "lookback_median",
    "events_mean", "events_sd", "lookback_mean", "lookback_sd")
)

seq_raw_dt <- id_all[, .N, by = .(State, n_sequences)][order(State, n_sequences)]
seq_raw_overall <- id_all[, .N, by = .(n_sequences)][order(n_sequences)]
seq_raw_overall[, State := "ALL"]
seq_raw_dt <- rbindlist(list(seq_raw_overall[, .(State, n_sequences, N)], seq_raw_dt), use.names = TRUE, fill = TRUE)

seq_grouped_dt <- copy(id_all)
seq_grouped_dt[, seq_group := group_0_1to4_5to14_15plus(n_sequences)]
seq_grouped_summary_dt <- rbindlist(
  list(
    seq_grouped_dt[, .(n_persons = .N), by = .(State, seq_group)],
    seq_grouped_dt[, .(n_persons = .N), by = .(seq_group)][, State := "ALL"]
  ),
  use.names = TRUE,
  fill = TRUE
)
seq_grouped_summary_dt[, pct_persons := round(100 * n_persons / sum(n_persons), 3), by = State]
setorderv(seq_grouped_summary_dt, c("State", "seq_group"))

total_sequences_dt <- id_all[, .(total_sequences = sum(n_sequences, na.rm = TRUE)), by = State][order(State)]
total_sequences_dt <- rbindlist(
  list(
    data.table(State = "ALL", total_sequences = sum(total_sequences_dt$total_sequences, na.rm = TRUE)),
    total_sequences_dt
  ),
  use.names = TRUE,
  fill = TRUE
)

recid_grouped_dt <- copy(id_all)
recid_grouped_dt[, recid_group := group_0_1to4_5to14_15plus(recid_events)]
recid_grouped_summary_dt <- rbindlist(
  list(
    recid_grouped_dt[, .(n_persons = .N), by = .(State, recid_group)],
    recid_grouped_dt[, .(n_persons = .N), by = .(recid_group)][, State := "ALL"]
  ),
  use.names = TRUE,
  fill = TRUE
)
recid_grouped_summary_dt[, pct_persons := round(100 * n_persons / sum(n_persons), 3), by = State]
setorderv(recid_grouped_summary_dt, c("State", "recid_group"))

person_year_state_dt <- py_all[
  ,
  .(
    person_year_obs = sum(person_year_obs, na.rm = TRUE),
    total_recid_events = sum(total_recid_events, na.rm = TRUE),
    py_with_recid = sum(py_with_recid, na.rm = TRUE)
  ),
  by = State
][order(State)]
person_year_state_dt[, pct_py_with_recid := round(100 * py_with_recid / person_year_obs, 3)]

person_year_overall_dt <- person_year_state_dt[
  ,
  .(
    State = "ALL",
    person_year_obs = sum(person_year_obs, na.rm = TRUE),
    total_recid_events = sum(total_recid_events, na.rm = TRUE),
    py_with_recid = sum(py_with_recid, na.rm = TRUE)
  )
]
person_year_overall_dt[, pct_py_with_recid := round(100 * py_with_recid / person_year_obs, 3)]

person_year_recap_dt <- rbindlist(list(person_year_overall_dt, person_year_state_dt), use.names = TRUE, fill = TRUE)
person_year_recap_dt[, pct_py_with_recid_str := sprintf("%.3f%%", pct_py_with_recid)]

rate_group_dt <- copy(id_all)
rate_group_dt[
  ,
  rate_group := fcase(
    recid_events == 0L, "0 = Never",
    is.finite(recid_rate) & recid_rate > 0 & recid_rate <= LOW_RATE_MAX, "Rate: Low",
    is.finite(recid_rate) & recid_rate > LOW_RATE_MAX & recid_rate <= MODERATE_RATE_MAX, "Rate: Moderate",
    is.finite(recid_rate) & recid_rate > MODERATE_RATE_MAX, "Rate: High",
    default = NA_character_
  )
]

rate_group_summary_dt <- rbindlist(
  list(
    rate_group_dt[, .(n = .N), by = .(State, rate_group)],
    rate_group_dt[, .(n = .N), by = .(rate_group)][, State := "ALL"]
  ),
  use.names = TRUE,
  fill = TRUE
)
rate_group_summary_dt[, pct := round(100 * n / sum(n), 2), by = State]
setorderv(rate_group_summary_dt, c("State", "rate_group"))

appendixA_example_dt <- copy(example_all)
appendixA_example_dt <- appendixA_example_dt[order(ResearchID, Sequence, LookBack)]
first_ids <- head(unique(appendixA_example_dt$ResearchID), 2L)
appendixA_example_dt <- appendixA_example_dt[ResearchID %in% first_ids]

fwrite(table1_dt, file.path(OUT_DIR, "table1_id_level_summary.csv"))
fwrite(appendixA_example_dt, file.path(OUT_DIR, "appendixA_person_year_example.csv"))
fwrite(person_year_recap_dt, file.path(OUT_DIR, "appendixA_population_style_descriptives.csv"))
fwrite(seq_raw_dt, file.path(OUT_DIR, "sequences_per_person_distribution_raw.csv"))
fwrite(seq_grouped_summary_dt, file.path(OUT_DIR, "sequences_per_person_grouped.csv"))
fwrite(recid_grouped_summary_dt, file.path(OUT_DIR, "recid_events_per_person_grouped.csv"))
fwrite(total_sequences_dt, file.path(OUT_DIR, "total_sequences_by_state.csv"))
fwrite(person_year_recap_dt, file.path(OUT_DIR, "person_year_recap_by_state.csv"))
fwrite(rate_group_summary_dt, file.path(OUT_DIR, "rate_group_summary.csv"))
fwrite(diag_all, file.path(OUT_DIR, "diagnostics_descriptives.csv"))

log_line("Descriptive summaries complete.")
print(table1_dt)
print(person_year_recap_dt)