############################################################
# SCRIPT 01 — REGISTER DISTRIBUTED SUBSAMPLES
#
# PURPOSE
#   Verify the distributed subsample files, summarize their
#   contents, and write a manifest for downstream scripts.
#
# INPUT
#   /data/subsamples/subsample_00_raw.csv ... subsample_19_raw.csv
#
# REQUIRED COLUMNS
#   ResearchID, State, LookBack, Recid
#
# OUTPUT
#   /output/01_register_subsamples/
#     - subsample_manifest.csv
#     - subsample_validation_summary.csv
#     - state_counts_by_subsample.csv
#     - state_counts_overall.csv
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
})

REPO_DIR <- "."
DATA_DIR <- file.path(REPO_DIR, "data")
SUBSAMPLE_DIR <- file.path(DATA_DIR, "subsamples")
OUT_DIR <- file.path(REPO_DIR, "output", "01_register_subsamples")

SUBSAMPLE_IDS <- sprintf("%02d", 0:19)
INPUT_FILES <- file.path(SUBSAMPLE_DIR, sprintf("subsample_%s_raw.csv", SUBSAMPLE_IDS))

REQUIRED_COLUMNS <- c("ResearchID", "State", "LookBack", "Recid")
EXPECTED_STATES <- c("AK", "AZ", "FL", "TX", "WA")

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

log_line("Checking distributed subsample files.")
stop_if_missing_files(INPUT_FILES)

manifest_list <- vector("list", length(INPUT_FILES))
state_count_list <- vector("list", length(INPUT_FILES))

for (i in seq_along(INPUT_FILES)) {
  f <- INPUT_FILES[i]
  sid <- SUBSAMPLE_IDS[i]
  
  log_line("Reading subsample %s: %s", sid, basename(f))
  
  dt <- fread(f, showProgress = FALSE)
  validate_columns(dt, f)
  
  dt[, ResearchID := as.character(ResearchID)]
  dt[, State := trimws(as.character(State))]
  dt[, LookBack := suppressWarnings(as.integer(LookBack))]
  dt[, Recid := suppressWarnings(as.integer(Recid))]
  
  dt_valid <- dt[
    !is.na(ResearchID) &
      !is.na(State) &
      !is.na(LookBack) &
      !is.na(Recid)
  ]
  
  observed_states <- sort(unique(dt_valid$State))
  unexpected_states <- setdiff(observed_states, EXPECTED_STATES)
  
  manifest_list[[i]] <- data.table(
    subsample = sid,
    file_name = basename(f),
    file_path = normalizePath(f, mustWork = TRUE),
    n_rows = nrow(dt),
    n_valid_rows = nrow(dt_valid),
    n_unique_ids = uniqueN(dt_valid$ResearchID),
    lookback_min = suppressWarnings(min(dt_valid$LookBack, na.rm = TRUE)),
    lookback_max = suppressWarnings(max(dt_valid$LookBack, na.rm = TRUE)),
    recid_values = paste(sort(unique(dt_valid$Recid)), collapse = ","),
    unexpected_states = paste(unexpected_states, collapse = ", ")
  )
  
  state_count_list[[i]] <- dt_valid[
    ,
    .(
      n_rows = .N,
      n_unique_ids = uniqueN(ResearchID)
    ),
    by = .(subsample = sid, State)
  ][order(State)]
  
  rm(dt, dt_valid)
  gc()
}

manifest_dt <- rbindlist(manifest_list, use.names = TRUE, fill = TRUE)
state_counts_dt <- rbindlist(state_count_list, use.names = TRUE, fill = TRUE)

if (any(nzchar(manifest_dt$unexpected_states))) {
  stop("Unexpected state values were detected in one or more subsamples.")
}

state_summary_dt <- state_counts_dt[
  ,
  .(
    total_rows = sum(n_rows, na.rm = TRUE),
    total_unique_id_entries = sum(n_unique_ids, na.rm = TRUE),
    n_subsamples_present = .N
  ),
  by = State
][order(State)]

validation_summary_dt <- data.table(
  metric = c(
    "Number of expected subsamples",
    "Number of detected subsamples",
    "All files present",
    "Minimum lookback across subsamples",
    "Maximum lookback across subsamples",
    "Total valid rows across subsamples",
    "Total unique ID entries across subsamples"
  ),
  value = c(
    length(SUBSAMPLE_IDS),
    nrow(manifest_dt),
    "Yes",
    min(manifest_dt$lookback_min, na.rm = TRUE),
    max(manifest_dt$lookback_max, na.rm = TRUE),
    sum(manifest_dt$n_valid_rows, na.rm = TRUE),
    sum(manifest_dt$n_unique_ids, na.rm = TRUE)
  )
)

fwrite(manifest_dt, file.path(OUT_DIR, "subsample_manifest.csv"))
fwrite(validation_summary_dt, file.path(OUT_DIR, "subsample_validation_summary.csv"))
fwrite(state_counts_dt, file.path(OUT_DIR, "state_counts_by_subsample.csv"))
fwrite(state_summary_dt, file.path(OUT_DIR, "state_counts_overall.csv"))

log_line("Subsample registration complete.")
print(validation_summary_dt)
print(state_summary_dt)
