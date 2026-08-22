############################################################
# check_data_schema.R  (diagnostic only — not part of the pipeline)
#
# PURPOSE
#   Inspect column names, types, and value encodings in the real
#   data files against what 01-08 expect, WITHOUT printing any
#   individual-level rows or ResearchID values. Safe to paste the
#   console output back into chat.
#
# USAGE ("repo root" = the folder you got by unzipping
#   decay-github-repo.zip -- the one that directly contains
#   01_register_subsamples.R, 02_descriptives.R, ... 08_*.R,
#   README.md, and a "data" folder.)
#
#   1. Save this file into that same folder, alongside 01_*.R ... 08_*.R.
#   2. Open it in RStudio and click the "Source" button at the top of
#      the editor pane (or press Ctrl+Shift+Enter / Cmd+Shift+Enter on
#      Mac). This script figures out its own folder automatically --
#      you do NOT need to set a working directory by hand.
#   3. Copy everything printed in the Console pane and paste it back.
#
#   (If you instead run it via source("check_data_schema.R") from a
#   console whose working directory is already the repo root, that
#   also works -- the auto-detect below falls back to the working
#   directory in that case.)
############################################################

suppressPackageStartupMessages(library(data.table))

# Figure out which folder this script itself is sitting in, so the user
# never has to think about "working directory" or "repo root" by hand.
# Tries, in order: (1) how Rscript was invoked from a command line,
# (2) RStudio's "Source" button / Ctrl+Shift+Enter, (3) falls back to
# the current working directory (correct if the user already cd'd /
# setwd()'d into the repo root before sourcing this file).
get_script_dir <- function() {
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
  getwd()
}

REPO_DIR <- get_script_dir()
DATA_DIR <- file.path(REPO_DIR, "data")

# fread() always fully decompresses a .gz file to a scratch temp file
# before parsing it -- there is no partial/streaming read, even for a
# small `nrows` preview like the ones below -- so previewing a .gz
# sensitivity file needs roughly as much free space in FREAD_TMPDIR as
# the UNCOMPRESSED file would take, however briefly. The real
# cohort_after2000.csv.gz / merged_state_events.csv.gz have been observed
# to be on the order of hundreds of millions of rows uncompressed, so this
# is not a small ask. Defaults to R's normal tempdir(); override with a
# path on a drive with more room if needed (e.g. "D:/Rtemp").
FREAD_TMPDIR <- tempdir()

hr <- function() cat(strrep("-", 70), "\n")

summarize_column <- function(x, name) {
  cls <- paste(class(x), collapse = "/")
  n_na <- sum(is.na(x))
  n_unique <- uniqueN(x)

  if (n_unique <= 12) {
    vals <- sort(unique(x))
    val_str <- paste(vals, collapse = ", ")
    cat(sprintf("  %-24s %-12s NA=%-6d unique=%-4d values=[%s]\n", name, cls, n_na, n_unique, val_str))
  } else if (is.numeric(x)) {
    cat(sprintf("  %-24s %-12s NA=%-6d unique=%-4d range=[%s, %s]\n", name, cls, n_na, n_unique,
                format(suppressWarnings(min(x, na.rm = TRUE))), format(suppressWarnings(max(x, na.rm = TRUE)))))
  } else {
    ex <- head(unique(x), 3)
    cat(sprintf("  %-24s %-12s NA=%-6d unique=%-4d example_values=[%s ...]\n", name, cls, n_na, n_unique,
                paste(ex, collapse = ", ")))
  }
}

check_file <- function(path, required_cols, label, n_preview = 5000) {
  hr()
  cat(label, "\n")
  cat("  path:", path, "\n")
  if (!file.exists(path)) {
    cat("  STATUS: FILE NOT FOUND -- skipping.\n")
    return(invisible(NULL))
  }

  if (grepl("\\.gz$", path)) {
    cat("  NOTE: this is a .gz file -- fread() will fully decompress it to a scratch\n")
    cat("        temp file first (in", FREAD_TMPDIR, "), which briefly needs about as much\n")
    cat("        free space as the file would take uncompressed, even for this preview.\n")
  }

  dt <- tryCatch(fread(path, nrows = n_preview, showProgress = FALSE, tmpdir = FREAD_TMPDIR), error = function(e) {
    cat("  STATUS: FAILED TO READ --", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(dt)) return(invisible(NULL))

  cat(sprintf("  STATUS: read OK (first %d rows previewed)\n", nrow(dt)))
  cat("  all column names found:", paste(names(dt), collapse = ", "), "\n")

  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    cat("  *** MISSING REQUIRED COLUMN(S):", paste(missing_cols, collapse = ", "), "***\n")
  } else {
    cat("  all required columns present.\n")
  }

  cat("  column summaries (no individual values shown for ResearchID):\n")
  for (col in intersect(c(required_cols, names(dt)), names(dt))) {
    if (col == "ResearchID") {
      cat(sprintf("  %-24s %-12s NA=%-6d unique=%d (values withheld)\n",
                   col, paste(class(dt[[col]]), collapse = "/"), sum(is.na(dt[[col]])), uniqueN(dt[[col]])))
    } else {
      summarize_column(dt[[col]], col)
    }
  }
  invisible(NULL)
}

cat("======================================================================\n")
cat("SCHEMA AUDIT --", format(Sys.time()), "\n")
cat("Detected repo root:", normalizePath(REPO_DIR, mustWork = FALSE), "\n")
cat("Looking for data in:", normalizePath(DATA_DIR, mustWork = FALSE), "\n")
if (!dir.exists(DATA_DIR)) {
  cat("\n")
  cat("*** No 'data' folder found at that location. ***\n")
  cat("This usually means one of two things:\n")
  cat("  1. This script isn't saved in the repo root -- move it into the same\n")
  cat("     folder as 01_register_subsamples.R, 02_descriptives.R, etc.\n")
  cat("     (that folder should also contain a 'data' subfolder), then re-run.\n")
  cat("  2. You're running this a way where auto-detection couldn't find its\n")
  cat("     own location -- run setwd(\"<path to that folder>\") first (with your\n")
  cat("     actual folder path in place of <path to that folder>), then\n")
  cat("     source(\"check_data_schema.R\") again.\n")
  cat("Continuing anyway -- every check below will report FILE NOT FOUND.\n")
}
cat("======================================================================\n")

## --- Core files (01-05) ---
check_file(
  file.path(DATA_DIR, "subsamples", "subsample_00_raw.csv"),
  c("ResearchID", "State", "LookBack", "Recid"),
  "subsamples/subsample_00_raw.csv"
)

check_file(
  file.path(DATA_DIR, "wa", "WA_decay_dataset.csv"),
  c("ResearchID", "Recid", "LookBack", "Male", "RaceEthnicity2", "AgeCurveGrouped", "CharType", "PriorIncarYearsCount"),
  "wa/WA_decay_dataset.csv"
)

check_file(
  file.path(DATA_DIR, "base_rates", "Baserate_for_modeling.csv"),
  c("State"),
  "base_rates/Baserate_for_modeling.csv (needs State + one of: base_rate | TotalOffenders+'TotalOffenders (NumeratorData)' | Denominator_Annual_State_Population+Numerator_Annual_Recidivism_Events)",
  n_preview = 1000
)

## --- Post-R&R sensitivity files (06-08), if they exist yet ---
## Distributed as .csv.gz (real data is large enough that an uncompressed
## export can run out of local disk space) -- fread() reads .csv.gz
## natively (requires the R.utils package) -- but either extension works.
pick_existing <- function(paths) {
  found <- paths[file.exists(paths)]
  if (length(found) > 0L) found[1] else paths[1]
}

check_file(
  pick_existing(file.path(DATA_DIR, "sensitivity", c("cohort_after2000.csv.gz", "cohort_after2000.csv"))),
  c("ResearchID", "State", "LookBack", "Recid", "After2000"),
  "sensitivity/cohort_after2000.csv[.gz]",
  n_preview = 5000
)

check_file(
  pick_existing(file.path(DATA_DIR, "sensitivity", c("merged_state_events.csv.gz", "merged_state_events.csv"))),
  c("ResearchID", "State", "OffDate", "LookBack", "Recid"),
  "sensitivity/merged_state_events.csv[.gz]",
  n_preview = 5000
)

hr()
cat("Also useful: how many subsample_XX_raw.csv files exist, and their sizes.\n")
subsample_files <- list.files(file.path(DATA_DIR, "subsamples"), pattern = "subsample_.*\\.csv$", full.names = TRUE)
if (length(subsample_files) == 0L) {
  cat("  none found in", file.path(DATA_DIR, "subsamples"), "\n")
} else {
  sizes_mb <- round(file.info(subsample_files)$size / 1e6, 1)
  cat(sprintf("  %d files found, size range: %.1f MB - %.1f MB, total: %.1f MB\n",
              length(subsample_files), min(sizes_mb), max(sizes_mb), sum(sizes_mb)))
}

wa_file <- file.path(DATA_DIR, "wa", "WA_decay_dataset.csv")
if (file.exists(wa_file)) cat(sprintf("  WA_decay_dataset.csv size: %.1f MB\n", file.info(wa_file)$size / 1e6))

hr()
cat("DONE. Copy everything above this line back into chat.\n")
