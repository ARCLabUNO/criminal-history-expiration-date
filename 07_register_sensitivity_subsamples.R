############################################################
# SCRIPT 07 — REGISTER "ONE EVENT PER PERSON-YEAR" SENSITIVITY
#             SUBSAMPLES
#
# STATUS
#   Added post-R&R, as data preparation for the one-charge-per-
#   year robustness check requested during review. Plays the same
#   role for the sensitivity pipeline (07 -> 08) that
#   01_register_subsamples.R plays for the primary pipeline
#   (01 -> 03/04/05): it builds and persists the subsample files
#   that the modeling script consumes.
#
# PURPOSE
#   The primary analysis (01-05) can include multiple charge
#   events for the same person within the same calendar year.
#   This script tests whether that clustering drives the decay
#   pattern by collapsing to a single (first) event per
#   ResearchID-year, then draws 20 independent 80% subsamples of
#   individuals for downstream modeling in
#   08_H1_sensitivity_one_event_per_year.R.
#
#   It also fixes a cross-state ResearchID collision: the raw
#   merged-state extract does not guarantee ResearchID is unique
#   ACROSS states (only within-state), so IDs are prefixed with
#   State before any person-level grouping. This collision was
#   confirmed to be isolated to this sensitivity dataset; the
#   distributed subsample_00-19_raw.csv files used by the primary
#   pipeline (01-05) were already state-safe and are unaffected.
#
# INPUT
#   /data/sensitivity/merged_state_events.csv
#
# REQUIRED COLUMNS
#   ResearchID, State, OffDate, LookBack, Recid
#
# OUTPUT
#   /output/07_register_sensitivity_subsamples/
#     - sensitivity_subsample_manifest.csv
#     - collapsed_events.csv.gz  (the full one-event-per-person-year panel,
#       all individuals, written ONCE)
#     - subsample_ids/sensitivity_subsample_00_ids.csv.gz ... _19_ids.csv.gz
#       (just the sampled ResearchIDs for each of the 20 subsamples --
#       08_H1_sensitivity_one_event_per_year.R reads collapsed_events.csv.gz
#       once and filters it per subsample using these ID lists, rather than
#       each subsample being a fully separate materialized file -- see the
#       "Disk-space redesign" note below.)
#
# DISK-SPACE REDESIGN (2026-08-13)
#   An earlier version of this script wrote each of the 20 subsamples as
#   its own fully materialized CSV (sensitivity_subsample_00.csv ... _19.csv,
#   one row per (chosen event, LookBack horizon) per sampled individual).
#   That failed on a real run with "No space left on device" partway
#   through subsample 02: since the collapse (see below) barely reduces
#   row count -- collapsing removes only individuals' EXTRA same-year
#   charge events, and most person-years have just one event to begin
#   with -- each 80% subsample materializes to roughly 80% of the FULL
#   panel's row count (confirmed against the real run's own log: subsample
#   02's row count was exactly 80.0% of the full input file's row count).
#   Writing 20 such subsamples separately means writing ~16x
#   (20 subsamples x 80%) the size of the underlying panel to disk, which
#   is already comparable in scale to the entire primary pipeline's data.
#   Fixed by writing the full collapsed panel ONCE (compressed), plus one
#   small per-subsample file listing just that subsample's sampled
#   ResearchIDs -- 08 does the actual per-subsample filtering in memory
#   from the one shared file instead of needing 20 separate huge ones.
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
})

# Auto-detect the repo root (the folder containing 01_*.R ... 08_*.R and a
# "data" folder) instead of assuming the R working directory is already
# pointed at it. Tries, in order: (1) how Rscript was invoked from a
# command line, (2) RStudio's "Source" button / Ctrl+Shift+Enter, (3)
# walking up from the current working directory looking for this repo's
# own data/ layout (covers running via source() from a console whose
# working directory is already inside, or at, the repo root).
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
DATA_DIR <- file.path(REPO_DIR, "data")
SENSITIVITY_DIR <- file.path(DATA_DIR, "sensitivity")
OUT_DIR <- file.path(REPO_DIR, "output", "07_register_sensitivity_subsamples")
SUBSAMPLE_ID_OUT_DIR <- file.path(OUT_DIR, "subsample_ids")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUBSAMPLE_ID_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# merged_state_events.csv is exported as .csv.gz for the same disk-space
# reason as cohort_after2000.csv.gz (see export_sensitivity_data_from_spss.R
# and 06_H2_historical_period.R) -- fread() reads .csv.gz natively (requires
# the R.utils package). Falls back to a plain .csv for anyone who
# exported/hosted it uncompressed instead.
INPUT_FILE_GZ <- file.path(SENSITIVITY_DIR, "merged_state_events.csv.gz")
INPUT_FILE_PLAIN <- file.path(SENSITIVITY_DIR, "merged_state_events.csv")
INPUT_FILE <- if (file.exists(INPUT_FILE_GZ)) INPUT_FILE_GZ else INPUT_FILE_PLAIN

# IMPORTANT if INPUT_FILE is .gz: fread() always fully decompresses a .gz
# file to a scratch temp file before parsing it (there is no partial/
# streaming read), so reading this file needs roughly as much free space
# in FREAD_TMPDIR as the UNCOMPRESSED file would take, however briefly.
# Defaults to R's normal tempdir(); override with a path on a drive that
# has more room if needed (e.g. "D:/Rtemp").
FREAD_TMPDIR <- tempdir()

REQUIRED_COLUMNS <- c("ResearchID", "State", "OffDate", "LookBack", "Recid")

N_SUBSAMPLES <- 20L
SUBSAMPLE_FRACTION <- 0.80
RNG_SEED <- 12345

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

validate_columns <- function(dt, required_cols, file_path) {
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in file:\n", file_path, "\nMissing columns: ", paste(missing_cols, collapse = ", "))
  }
}

log_line("Checking sensitivity input file.")
stop_if_missing_files(INPUT_FILE)

log_line("Reading merged-state event file: %s", basename(INPUT_FILE))
dat <- fread(INPUT_FILE, showProgress = FALSE, tmpdir = FREAD_TMPDIR)
validate_columns(dat, REQUIRED_COLUMNS, INPUT_FILE)

log_line("Total raw rows loaded: %s", format(nrow(dat), big.mark = ","))

# --- Fix cross-state ID collisions ---
log_line("Enforcing global uniqueness by prefixing ResearchID with State.")
dat[, State := trimws(as.character(State))]
dat[, ResearchID := paste0(State, "_", as.character(ResearchID))]
log_line("Unique cross-state prefixed individuals: %s", format(uniqueN(dat$ResearchID), big.mark = ","))

EXPECTED_STATES <- c("AK", "AZ", "FL", "TX", "WA")
observed_states <- sort(unique(dat$State))
unexpected_states <- setdiff(observed_states, EXPECTED_STATES)
if (length(unexpected_states) > 0L) {
  stop(
    "Unexpected state value(s) in input file: ", paste(unexpected_states, collapse = ", "),
    ". Expected only: ", paste(EXPECTED_STATES, collapse = ", ")
  )
}
log_line("State coverage confirmed: %s", paste(observed_states, collapse = ", "))

# --- Collapse to one event per person-year ---
# OffDate in the underlying extracts is inconsistent: some pipeline files
# store it as a full parseable date, but sibling date-like columns in the
# same extract family (BackBegin, TrackBegin, TrackEnd, BackEnd) have been
# observed to hold a bare 4-digit year (e.g. 2015) instead. as.Date("2015")
# throws a hard error ("character string is not in a standard unambiguous
# format"), so a bare year would crash this script rather than silently
# misparse. Handle both forms explicitly instead of assuming one.
derive_year <- function(x) {
  x_chr <- trimws(as.character(x))
  is_bare_year <- grepl("^[12][0-9]{3}$", x_chr)

  year_out <- rep(NA_integer_, length(x_chr))
  year_out[is_bare_year] <- as.integer(x_chr[is_bare_year])

  needs_date_parse <- !is_bare_year & !is.na(x_chr) & nzchar(x_chr)
  if (any(needs_date_parse)) {
    parsed <- suppressWarnings(as.Date(x_chr[needs_date_parse]))
    year_out[needs_date_parse] <- as.integer(format(parsed, "%Y"))
  }
  year_out
}

log_line("Collapsing to 1 event per person per calendar year.")
dat[, Year := derive_year(OffDate)]
if (dat[, mean(is.na(Year))] > 0) {
  log_line("NOTE: %s of %s OffDate values could not be parsed into a year and were dropped.",
           format(sum(is.na(dat$Year)), big.mark = ","), format(nrow(dat), big.mark = ","))
}
dat <- dat[!is.na(ResearchID) & !is.na(Year) & !is.na(LookBack) & !is.na(Recid)]
dat[, LookBack := suppressWarnings(as.integer(LookBack))]
dat[, Recid := suppressWarnings(as.integer(Recid))]

# IMPORTANT: each row here is one (charge event, LookBack horizon)
# observation -- a single charge repeats across many rows (LookBack = 1 up
# to 40), all sharing that one charge's own OffDate/Year. A naive
# `dat[, .SD[1], by = .(ResearchID, Year)]` collapse (an earlier version of
# this script) groups on (ResearchID, Year) alone, which does NOT separate
# "two distinct charges in the same year" (what we actually want to
# collapse) from "the same one charge's 40 different LookBack rows" (which
# all share that year and must NOT be collapsed away). That bug was caught
# via a real-data run: every one of the 20 subsamples' lookback_min/
# lookback_max in the manifest came back as exactly 1 -- i.e. every
# LookBack-horizon row past LookBack=1 had silently been discarded for
# every person-year, which would have made 08's decay model unfittable.
#
# Fix: collapse at the EVENT (OffDate) level, not the row level. For each
# person-year, choose ONE representative charge event (the earliest
# OffDate that year), then keep EVERY LookBack row belonging to that
# chosen event.
setorder(dat, ResearchID, Year, OffDate, LookBack)
chosen_events <- unique(dat[, .(ResearchID, Year, OffDate)])[, .SD[1], by = .(ResearchID, Year)]
collapsed_dt <- dat[chosen_events, on = .(ResearchID, Year, OffDate)]
log_line("Collapsed row count (all LookBack rows of 1 chosen event per person-year): %s", format(nrow(collapsed_dt), big.mark = ","))
log_line("Sanity check -- LookBack range after collapse (should span more than just 1 if the fix is working): %s to %s",
         suppressWarnings(min(collapsed_dt$LookBack, na.rm = TRUE)), suppressWarnings(max(collapsed_dt$LookBack, na.rm = TRUE)))
rm(dat)
gc()

# Write the full collapsed panel ONCE, compressed -- this is the one large
# write this script performs; the 20 subsamples below are persisted as
# small ID lists instead of 20 separate copies of this same panel (see the
# "Disk-space redesign" note in the header).
COLLAPSED_OUT_PATH <- file.path(OUT_DIR, "collapsed_events.csv.gz")
log_line("Writing full collapsed panel (once) to: %s", COLLAPSED_OUT_PATH)
fwrite(collapsed_dt, COLLAPSED_OUT_PATH, compress = "gzip")

# --- Generate deterministic 80% subsamples ---
log_line("Generating %d independent %.0f%% sensitivity subsamples (seed = %d).", N_SUBSAMPLES, SUBSAMPLE_FRACTION * 100, RNG_SEED)
set.seed(RNG_SEED)
unique_ids <- unique(collapsed_dt$ResearchID)

manifest_list <- vector("list", N_SUBSAMPLES)

for (i in 0:(N_SUBSAMPLES - 1)) {
  sid <- sprintf("%02d", i)
  log_line("Generating subsample %s of %d.", sid, N_SUBSAMPLES)

  sampled_ids <- sample(unique_ids, size = floor(length(unique_ids) * SUBSAMPLE_FRACTION), replace = FALSE)

  # Persist just the sampled ID list (small) instead of a full materialized
  # copy of this subsample's rows (which would be ~80% of the entire
  # collapsed panel -- see "Disk-space redesign" in the header). The
  # manifest stats below are still computed from the real filtered rows,
  # just without writing those rows back out to disk.
  ids_file_path <- file.path(SUBSAMPLE_ID_OUT_DIR, sprintf("sensitivity_subsample_%s_ids.csv.gz", sid))
  fwrite(data.table(ResearchID = sampled_ids), ids_file_path, compress = "gzip")

  sub_dt <- collapsed_dt[ResearchID %in% sampled_ids]

  manifest_list[[i + 1L]] <- data.table(
    subsample = sid,
    ids_file_name = basename(ids_file_path),
    n_rows = nrow(sub_dt),
    n_unique_ids = uniqueN(sub_dt$ResearchID),
    lookback_min = suppressWarnings(min(sub_dt$LookBack, na.rm = TRUE)),
    lookback_max = suppressWarnings(max(sub_dt$LookBack, na.rm = TRUE)),
    recid_events = sum(sub_dt$Recid, na.rm = TRUE)
  )
  rm(sub_dt)
}

manifest_dt <- rbindlist(manifest_list, use.names = TRUE, fill = TRUE)
fwrite(manifest_dt, file.path(OUT_DIR, "sensitivity_subsample_manifest.csv"))

log_line("Sensitivity subsample registration complete.")
log_line("Collapsed panel written to: %s", COLLAPSED_OUT_PATH)
log_line("Subsample ID lists written to: %s", SUBSAMPLE_ID_OUT_DIR)
print(manifest_dt)
