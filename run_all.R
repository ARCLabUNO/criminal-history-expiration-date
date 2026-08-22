############################################################
# run_all.R
#
# Reconstructed for this repository (not provided in the original
# upload) — sources every script in the pipeline in order. Verify
# this against your working copy before publishing.
#
# Each script begins with `rm(list = ls())`, which would wipe this
# script's own variables (e.g. the loop below) if sourced into the
# global environment directly. Each call below uses
# `source(..., local = TRUE)` so every script runs in its own throw-
# away environment and only `rm(list = ls())`s itself.
#
# 06-08 are the post-R&R additions (historical period + one-event-
# per-person-year sensitivity check). They depend on
# data/sensitivity/*.csv, which is not part of the original data
# drop (see README's "Post-R&R sensitivity files" section) and may
# not be present in every checkout. If those inputs are missing,
# this script logs a warning and skips 06/07 (08 is skipped
# automatically since it depends on 07's output) rather than
# aborting the whole run.
############################################################

# Auto-detect the repo root (the folder containing 01_*.R ... 08_*.R and a
# "data" folder) instead of assuming the R working directory is already
# pointed at it. Tries, in order: (1) how Rscript was invoked from a
# command line, (2) RStudio's "Source" button / Ctrl+Shift+Enter, (3)
# walking up from the current working directory looking for this repo's
# own data/ layout.
get_repo_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  hit <- grep(file_flag, cmd_args)
  if (length(hit) > 0L) {
    return(dirname(normalizePath(sub(file_flag, "", cmd_args[hit[1]]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path)))
    }
  }
  cur <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in 0:6) {
    if (dir.exists(file.path(cur, "data")) &&
        (dir.exists(file.path(cur, "data", "subsamples")) ||
         dir.exists(file.path(cur, "data", "sensitivity")) ||
         dir.exists(file.path(cur, "data", "wa")))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop(
    "Could not automatically determine the repo root (the folder containing\n",
    "01_register_subsamples.R ... 08_*.R and a 'data' folder).\n",
    "Fix: open this script directly in RStudio and click 'Source' (Session >\n",
    "Set Working Directory > To Source File Location also works), or run\n",
    "setwd(\"<path to that folder>\") yourself before sourcing this script."
  )
}

REPO_DIR <- get_repo_root()

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush.console()
}

run_script <- function(script_name) {
  # `local = TRUE` would source into THIS function's own call frame -- and
  # every script starts with `rm(list = ls())`, which would then wipe out
  # this function's own `script_name` argument before the line below runs.
  # `local = new.env(...)` gives each script a genuinely fresh, disposable
  # environment instead, so it can only ever clear its own variables.
  path <- file.path(REPO_DIR, script_name)
  log_line("Running %s", script_name)
  source(path, local = new.env(parent = globalenv()))
  log_line("Finished %s", script_name)
}

## --- Primary pipeline (original submission) ---------------------------

run_script("01_register_subsamples.R")
run_script("02_descriptives.R")
run_script("03_H1_decay.R")
run_script("04_H2_invariance.R")
run_script("05_H3_convergence.R")

## --- Post-R&R additions ------------------------------------------------

# Sensitivity inputs are distributed as .csv.gz (see README) -- they're
# large enough in the real data that an uncompressed export can run out of
# local disk space -- but either extension is accepted here.
sensitivity_dir <- file.path(REPO_DIR, "data", "sensitivity")
cohort_input <- Filter(file.exists, file.path(sensitivity_dir, c("cohort_after2000.csv.gz", "cohort_after2000.csv")))[1]
merged_states_input <- Filter(file.exists, file.path(sensitivity_dir, c("merged_state_events.csv.gz", "merged_state_events.csv")))[1]

if (!is.na(cohort_input)) {
  run_script("06_H2_historical_period.R")
} else {
  log_line("Skipping 06_H2_historical_period.R -- missing input: %s", file.path(sensitivity_dir, "cohort_after2000.csv[.gz]"))
}

if (!is.na(merged_states_input)) {
  run_script("07_register_sensitivity_subsamples.R")
  run_script("08_H1_sensitivity_one_event_per_year.R")
} else {
  log_line("Skipping 07/08 (one-event-per-year sensitivity) -- missing input: %s", file.path(sensitivity_dir, "merged_state_events.csv[.gz]"))
}

log_line("run_all.R complete. See /output/ for results.")
