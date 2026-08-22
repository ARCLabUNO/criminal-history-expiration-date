############################################################
# SCRIPT 08 — H1 SENSITIVITY: ONE EVENT PER PERSON-YEAR
#
# STATUS
#   Added post-R&R, as a robustness check requested during review.
#   Fits the primary decay model (H1, cf. 03_H1_decay.R) on the
#   one-event-per-person-year subsamples registered by
#   07_register_sensitivity_subsamples.R, to test whether the
#   decay pattern reported in Table 2 is an artifact of people
#   contributing multiple charge events within the same year.
#
# PURPOSE
#   For each of the 20 subsamples produced by Script 07, fit:
#     Recid ~ LookBack_z + (1 | ResearchID)
#   then pool the 20 subsample-level estimates via random-effects
#   meta-analysis (metafor::rma.uni), the same pooling approach
#   used throughout the primary pipeline (03/04/06). This produces
#   a single comparable pooled OR for the Lookback slope that can
#   be checked against the primary model's Table 2 estimate
#   (OR = 0.39 GLMM / 0.44 GLM in the primary analysis).
#
# INPUT
#   /output/07_register_sensitivity_subsamples/
#     - collapsed_events.csv.gz  (the full one-event-per-person-year panel,
#       all individuals -- streamed and filtered per subsample, never
#       fully loaded into memory at once -- see "Memory redesign" below)
#     - subsample_ids/sensitivity_subsample_00_ids.csv.gz ... _19_ids.csv.gz
#       (per-subsample ID lists)
#
# REQUIRED COLUMNS
#   ResearchID, LookBack, Recid
#     (LookBack_z is computed here if not already present. Only these 3
#     columns are ever materialized in R -- see "Memory redesign" below.)
#
# OUTPUT
#   /output/08_H1_sensitivity_one_event_per_year/
#     - subsample_shards/subsample_XX_shard_YY.csv.gz  (scratch -- each
#       subsample's rows partitioned into N_SHARDS_PER_SUBSAMPLE smaller
#       files; deleted after that subsample's shards are all fit unless
#       KEEP_SHARD_FILES is set)
#     - shard_checkpoints/subsample_XX_shard_YY.rds  (per-SHARD fitted
#       results -- a re-run skips any shard that already has a checkpoint)
#     - subsample_pooled_checkpoints/subsample_XX.rds  (per-SUBSAMPLE
#       already-shard-pooled results -- a re-run skips a subsample entirely,
#       without even re-reading its shard checkpoints, once this exists)
#     - sensitivity_shard_coefficients.csv / _fit_stats.csv / _random_effects.csv
#       (one row per shard -- the fine-grained diagnostic detail)
#     - sensitivity_subsample_coefficients.csv  (one row per subsample --
#       each subsample's shards pooled via meta-analysis first)
#     - sensitivity_pooled_coefficients.csv  (the 20 subsample estimates
#       pooled again, into the final overall result)
#     - sensitivity_vs_primary_comparison.csv
#
# MEMORY REDESIGN (2026-08-13, two rounds)
#   Round 1: a real run crashed from an out-of-memory error partway through
#   subsample 2 of 20. Fixed two inefficiencies: the script had been
#   holding the ENTIRE collapsed panel (~746M rows) in memory for the whole
#   run in addition to each subsample's ~80%-sized filtered copy, and that
#   copy carried 3 unused columns (State/OffDate/Year) alongside the 3 the
#   model actually needs. Fixed by streaming collapsed_events.csv.gz
#   directly per subsample (never holding the full panel) and dropping
#   unused columns immediately -- see the streaming code below.
#
#   Round 2: even with that fix, a single subsample is still ~597M rows /
#   ~16M individuals -- one `glmer()` call on that much data was still
#   stalling out real machines. This mirrors exactly why the PRIMARY
#   pipeline (01/03/04/06) never fits one monster model either: it splits
#   into 20 partitions and pools via meta-analysis. The same idea is
#   applied one level deeper here: each of the 20 subsamples is further
#   split into N_SHARDS_PER_SUBSAMPLE smaller shards (default 8, so 160
#   shards total -- in the same ballpark as the ~156-shard scale used for
#   the primary analyses), `glmer` fits each shard individually (far
#   smaller and far less memory-hungry), and shard-level estimates are
#   pooled via `metafor::rma.uni` into a subsample-level estimate -- which
#   then feeds into the EXISTING subsample-level pooling step unchanged, so
#   the final output is produced the same way as before, just built from
#   smaller pieces. Checkpointing now happens at the shard level, so a
#   stall loses at most one shard's progress, not a whole subsample's.
#
#   Note on "data compression": 03_H1_decay.R and 06_H2_historical_period.R
#   use a binomial success/failure COUNT compression (collapsing rows that
#   share the same ResearchID + LookBack into one row with `y`/`n` counts)
#   before fitting -- that's a real, different technique from sharding, and
#   it's NOT ported here on purpose: it only helps when multiple rows
#   legitimately share the same (ResearchID, LookBack) pair, which is
#   exactly the clustering 07's one-event-per-person-year collapse already
#   removed. Applying it to 08's data would be a no-op (every group already
#   has exactly 1 row) -- sharding is the applicable lever here, not count
#   compression. Also ported over from 03/06: `nAGQ = 0` and
#   `calc.derivs = FALSE` in glmerControl, a real, direct memory/speed win
#   (cheaper Laplace approximation, skips a post-fit derivative check) that
#   was already established as safe in this repo but had been missed here.
#
#   Round 3 (2026-08-14): a real run showed sharding alone wasn't enough --
#   bucketing subsample 00 took ~54 minutes, and each of its shards then
#   took 60-107 minutes to FIT (one real shard took nearly 2 hours). At
#   that rate, 160 shards fit one at a time is on the order of a week of
#   wall-clock time. Bucketing (one streaming read+write pass) is a small
#   fraction of that total -- the `glmer()` FIT of each shard is the
#   dominant cost by roughly an order of magnitude. Since shard fits are
#   fully independent of each other, they now run in parallel across up to
#   N_PARALLEL_WORKERS R worker processes (base-R `parallel` package, PSOCK
#   cluster, no new dependencies, Windows-safe) instead of one at a time --
#   see `fit_shards_parallel()` below. This is a pure engineering change
#   (same model, same data, same settings per shard); it does not touch
#   bucketing, sharding, or pooling logic, and every existing checkpoint
#   file (shard-level or subsample-level) from a prior serial run is still
#   recognized and skipped exactly as before, so a resume picks up right
#   where it left off, and any shards already fit under the old serial
#   behavior are not refit.
#
#   Round 4 (2026-08-15): the first parallel run still looked slow -- a real
#   log showed a 3-shard batch dispatched with no further output for a long
#   stretch, which turned out to be TWO separate problems, not one. (1) The
#   original design dispatched shards in FIXED batches of N_PARALLEL_WORKERS
#   via parLapply(), which only reports progress once the ENTIRE batch
#   returns -- with real per-shard fit times varying 60-107 minutes, that
#   meant long stretches with no log output at all even though work was
#   progressing normally, which reads exactly like a hang. Fixed by
#   switching to a dynamic queue (fit_shards_parallel() now uses
#   parallel:::sendCall()/recvOneResult(), the same low-level mechanism
#   parLapply() itself is built on): the moment ANY worker finishes, that
#   shard is logged immediately and the freed-up worker is handed the next
#   pending shard right away, instead of waiting for its batch-mates. (2) A
#   real, likely-more-serious problem: R's BLAS library (OpenBLAS/MKL, used
#   internally by glmer()'s matrix algebra) multi-threads by default on many
#   systems -- with N_PARALLEL_WORKERS separate R processes EACH ALSO trying
#   to spawn several BLAS threads, a machine can end up with, e.g., 3
#   workers x 4 BLAS threads = 12 threads competing for 4 physical cores.
#   That oversubscription can make "parallel" fitting SLOWER than plain
#   serial fitting, not faster -- a very plausible explanation for a
#   parallelized run still feeling slow. Fixed by setting
#   OMP_NUM_THREADS/OPENBLAS_NUM_THREADS/MKL_NUM_THREADS=1 inside each
#   worker before any BLAS-touching library loads, so each of the
#   N_PARALLEL_WORKERS worker processes uses exactly one thread and the
#   worker count lines up with physical cores as intended. Neither change
#   touches the model, the data, or per-shard settings -- both are pure
#   scheduling/threading fixes, and checkpoint/resume behavior (including
#   checkpoints from earlier serial OR fixed-batch-parallel runs) is
#   unaffected.
#
#   Unresolved: the exact "156 shards across twenty sample data sets"
#   figure mentioned for the primary analyses doesn't appear anywhere in
#   this repo's R code (checked 01/03/04/05/06) -- it likely comes from an
#   upstream step (e.g. how subsample_00-19_raw.csv were originally
#   partitioned before this repo ever sees them). N_SHARDS_PER_SUBSAMPLE=8
#   is a reasonable default landing close to that same total shard count,
#   not a reproduction of a specific documented method -- raise/lower it
#   based on what your machine actually handles.
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(broom.mixed)
  library(metafor)
  library(parallel)  # base R -- no install needed, works on Windows via PSOCK clusters
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
SENSITIVITY_07_DIR <- file.path(REPO_DIR, "output", "07_register_sensitivity_subsamples")
COLLAPSED_FILE <- file.path(SENSITIVITY_07_DIR, "collapsed_events.csv.gz")
SUBSAMPLE_ID_DIR <- file.path(SENSITIVITY_07_DIR, "subsample_ids")
OUT_DIR <- file.path(REPO_DIR, "output", "08_H1_sensitivity_one_event_per_year")
SHARD_DIR <- file.path(OUT_DIR, "subsample_shards")
SHARD_CHECKPOINT_DIR <- file.path(OUT_DIR, "shard_checkpoints")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SHARD_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SHARD_CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

# Rows are streamed out of COLLAPSED_FILE this many lines at a time while
# bucketing into shards (see bucket_subsample_into_shards() below) -- lower
# this if a single chunk's in-memory data.table is still too large for your
# machine; raise it for fewer, larger chunks (a bit faster, more transient
# memory per chunk).
CHUNK_LINES <- 5e6

# Each of the 20 subsamples (~597M rows / ~16M individuals) is further
# split into this many shards before any model is fit -- see "MEMORY
# REDESIGN, round 2" above for why. Raised from an earlier default of 8 to
# 20 (400 total shards) after a real run caused a full system crash (blank
# screen, forced restart) with multiple ~75M-row shards fitting at once --
# smaller shards mean each PARALLEL WORKER's memory footprint is smaller
# too, which matters more than raw shard count once you're running several
# shards concurrently (see N_PARALLEL_WORKERS below). This does NOT reduce
# total computation -- glmer's cost scales roughly linearly with rows, so
# more/smaller shards mostly trade "fewer, slower fits" for "more, faster
# fits" at roughly the same total CPU time -- its real value here is
# capping peak memory per shard and giving much finer checkpointing (a
# crash now loses a few minutes of one shard's work, not up to two hours).
# Raise further (e.g. towards the ~156/subsample scale used upstream in
# the primary analyses, per the user) if memory is still an issue at 20;
# lower it only once you've confirmed the machine is stable and want less
# per-shard/pooling overhead.
N_SHARDS_PER_SUBSAMPLE <- 20L

# Shard scratch files are deleted after that subsample's shards are all
# successfully fit and pooled (they're regenerable by re-running the
# bucketing step). Set TRUE to keep them around for inspection/debugging.
KEEP_SHARD_FILES <- FALSE

# Number of R worker processes used to FIT shards in parallel -- this is
# the dominant real-world cost in this script. Shards are fully
# independent model fits, so they parallelize cleanly, BUT each worker
# needs roughly the same RAM as one serial shard fit, so more workers
# means more simultaneous RAM use, not just more speed -- if that exceeds
# available RAM, the OS starts swapping/thrashing, which can make things
# far slower than serial and, on a real run, caused a full system crash
# (blank screen, forced restart) with several ~75M-row shards fitting at
# once. Dropped from an earlier default of min(4, cores-1) to a much more
# conservative min(2, cores-1) for that reason -- running only 2 shards at
# once (now much smaller too, see N_SHARDS_PER_SUBSAMPLE above) leaves far
# more memory headroom. Set to 1 to fall back to the original fully-serial
# behavior if instability persists even at 2 -- slow and correct beats
# fast and crashing. Only raise this once a run has completed several
# subsamples without incident, and even then raise it one step at a time
# rather than jumping back up to 3-4. Checkpointing is unaffected either
# way -- each worker saves its own shard's .rds the moment that shard
# finishes, so an interrupted run still loses at most one shard's worth of
# progress (now a few minutes at N_SHARDS_PER_SUBSAMPLE=20, not up to two
# hours), not a whole batch's or subsample's.
N_PARALLEL_WORKERS <- max(1L, min(2L, parallel::detectCores() - 1L))

N_SUBSAMPLES <- 20L
SUBSAMPLE_IDS <- sprintf("%02d", 0:(N_SUBSAMPLES - 1L))
ID_LIST_FILES <- file.path(SUBSAMPLE_ID_DIR, sprintf("sensitivity_subsample_%s_ids.csv.gz", SUBSAMPLE_IDS))

REQUIRED_COLUMNS <- c("ResearchID", "LookBack", "Recid")
RESID_VAR <- (pi^2) / 3

# Primary-model estimates from Table 2 (03_H1_decay.R), for side-by-side comparison.
PRIMARY_LOOKBACK_OR_GLMM <- 0.39
PRIMARY_LOOKBACK_OR_GLM <- 0.44

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush.console()
}

# A checkpoint file that EXISTS is normally trusted as "already done" (see
# every `file.exists(checkpoint_path)` check throughout this script) -- but
# saveRDS() does not write atomically, so a hard kill mid-write (a real
# system crash, as opposed to a clean R-level interrupt) can leave a
# truncated, unreadable .rds file on disk that still passes file.exists().
# Reading it back with a plain readRDS() would then throw partway through
# final aggregation, potentially crashing the whole run after it had
# otherwise finished. This wrapper treats a read failure as "not actually
# done": it logs a clear warning, deletes the corrupt file so future runs
# don't keep tripping over it, and returns NULL so the caller's existing
# "missing checkpoint -- will retry on next run" logic picks it back up
# normally, rather than the script dying on a bad file.
safe_read_rds <- function(path) {
  result <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(result, "error")) {
    log_line("  WARNING: checkpoint file appears corrupt/truncated (likely from a hard crash) -- deleting and will retry: %s (%s)",
              path, conditionMessage(result))
    unlink(path)
    return(NULL)
  }
  result
}

stop_if_missing_files <- function(paths) {
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files) > 0L) {
    stop(
      "The following required files are missing:\n", paste(missing_files, collapse = "\n"),
      "\n\nRun 07_register_sensitivity_subsamples.R first."
    )
  }
}

extract_random_intercept_variance <- function(model, id_var = "ResearchID") {
  vc <- tryCatch(as.data.frame(VarCorr(model)), error = function(e) NULL)
  if (is.null(vc)) return(NA_real_)
  vv <- vc$vcov[vc$grp == id_var][1]
  if (length(vv) == 0L || is.na(vv)) return(NA_real_)
  as.numeric(vv)
}

calc_icc <- function(var_u0) {
  if (!is.finite(var_u0) || is.na(var_u0)) return(NA_real_)
  var_u0 / (var_u0 + RESID_VAR)
}

meta_pool_coefficients <- function(coef_dt) {
  coef_dt <- copy(coef_dt)
  coef_dt[, vi := std.error^2]
  combos <- unique(coef_dt[, .(term)])

  pooled_list <- lapply(seq_len(nrow(combos)), function(i) {
    t <- combos$term[i]
    sub <- coef_dt[term == t]

    fit <- tryCatch(
      metafor::rma.uni(yi = sub$estimate, vi = sub$vi, method = "REML"),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(data.table(
        term = t, pooled_est = NA_real_, pooled_se = NA_real_,
        z = NA_real_, p = NA_real_, ci_low = NA_real_, ci_high = NA_real_
      ))
    }

    data.table(
      term = t,
      pooled_est = as.numeric(fit$b), pooled_se = as.numeric(fit$se),
      z = as.numeric(fit$zval), p = as.numeric(fit$pval),
      ci_low = as.numeric(fit$ci.lb), ci_high = as.numeric(fit$ci.ub)
    )
  })

  pooled_dt <- rbindlist(pooled_list, use.names = TRUE, fill = TRUE)
  pooled_dt[, `:=`(OR = exp(pooled_est), OR_low = exp(ci_low), OR_high = exp(ci_high))]
  pooled_dt[]
}

shard_file_path <- function(sid, shard_no) file.path(SHARD_DIR, sprintf("subsample_%s_shard_%02d.csv.gz", sid, shard_no))
shard_checkpoint_path <- function(sid, shard_no) file.path(SHARD_CHECKPOINT_DIR, sprintf("subsample_%s_shard_%02d.rds", sid, shard_no))

# Streams collapsed_gz_path ONCE via a gzfile() connection (pure base R +
# data.table -- no external CLI tools, so this stays portable to Windows),
# and fans this subsample's rows out into N_SHARDS_PER_SUBSAMPLE separate
# small files (round-robin assignment over sampled_ids, which are already
# in random order -- see 07's `sample()` call -- so this is an unbiased
# split, not a systematic one). At most one chunk (CHUNK_LINES rows) is
# held in memory at a time; each chunk's shard pieces are written and
# discarded immediately, so peak memory here is roughly one chunk, not one
# subsample.
bucket_subsample_into_shards <- function(collapsed_gz_path, sampled_ids, cols_keep,
                                          sid, n_shards, chunk_lines = CHUNK_LINES) {
  shard_paths <- vapply(seq_len(n_shards), function(s) shard_file_path(sid, s), character(1))
  unlink(shard_paths[file.exists(shard_paths)])  # clear any stale partial files before starting

  id_shard_dt <- data.table(ResearchID = sampled_ids,
                             shard = (seq_along(sampled_ids) - 1L) %% n_shards + 1L)
  setkey(id_shard_dt, ResearchID)

  con <- gzfile(collapsed_gz_path, open = "rt")
  on.exit(close(con), add = TRUE)

  header_line <- readLines(con, n = 1L)
  header_cols <- strsplit(header_line, ",", fixed = TRUE)[[1]]
  if (!all(cols_keep %in% header_cols)) {
    stop("Missing required column(s) in ", collapsed_gz_path, ": ",
         paste(setdiff(cols_keep, header_cols), collapse = ", "))
  }

  wrote_any <- rep(FALSE, n_shards)
  repeat {
    lines <- readLines(con, n = chunk_lines)
    if (length(lines) == 0L) break
    chunk_dt <- fread(text = lines, header = FALSE, col.names = header_cols, showProgress = FALSE)
    chunk_dt <- chunk_dt[id_shard_dt, on = "ResearchID", nomatch = NULL]  # inner join: filters + tags shard in one step
    if (nrow(chunk_dt) == 0L) next
    for (s in seq_len(n_shards)) {
      piece <- chunk_dt[shard == s, ..cols_keep]
      if (nrow(piece) == 0L) next
      fwrite(piece, shard_paths[s], compress = "gzip", append = wrote_any[s])
      wrote_any[s] <- TRUE
    }
  }
  shard_paths
}

fit_one_shard <- function(shard_gz_path, sid, shard_no) {
  label <- sprintf("%s/shard %02d", sid, shard_no)
  log_line("  Fitting %s: %s", label, basename(shard_gz_path))

  dt <- fread(shard_gz_path, showProgress = FALSE)
  dt[, ResearchID := as.character(ResearchID)]
  dt[, LookBack := suppressWarnings(as.integer(LookBack))]
  dt[, Recid := suppressWarnings(as.integer(Recid))]
  dt <- dt[!is.na(ResearchID) & !is.na(LookBack) & !is.na(Recid)]
  dt <- dt[Recid %in% c(0L, 1L)]

  dt[, ResearchID := droplevels(factor(ResearchID))]
  dt[, LookBack_z := as.numeric(scale(LookBack))]

  # nAGQ=0 (cheaper Laplace approximation) and calc.derivs=FALSE (skips a
  # post-fit derivative check) match 03_H1_decay.R / 06_H2_historical_period.R
  # exactly -- both a real speed/memory win, and consistency with the rest
  # of this repo's established, already-real-data-tested glmer settings.
  fit <- tryCatch(
    glmer(
      Recid ~ LookBack_z + (1 | ResearchID),
      data = dt, family = binomial(link = "logit"), nAGQ = 0L,
      control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE, optCtrl = list(maxfun = 2e5))
    ),
    error = function(e) {
      log_line("  Model failed to converge for %s: %s", label, conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) return(NULL)

  coef_dt <- as.data.table(broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE))
  coef_dt[, `:=`(subset = sid, shard = shard_no)]

  ll <- as.numeric(logLik(fit))
  p_hat <- fitted(fit)

  fit_dt <- data.table(
    subset = sid, shard = shard_no,
    minus2LL = -2 * ll,
    AIC = AIC(fit),
    BIC = BIC(fit),
    RMSE = sqrt(mean((dt$Recid - p_hat)^2)),
    Brier = mean((dt$Recid - p_hat)^2),
    n_obs = nobs(fit),
    n_ids = uniqueN(dt$ResearchID)
  )

  var_u0 <- extract_random_intercept_variance(fit)
  re_dt <- data.table(
    subset = sid, shard = shard_no,
    var_u0 = var_u0,
    sd_u0 = sqrt(var_u0),
    ICC = calc_icc(var_u0)
  )

  list(coef = coef_dt, fit = fit_dt, re = re_dt)
}

# Fits every PENDING shard in `jobs` (a data.table with columns sid,
# shard_no, shard_path, checkpoint_path -- already filtered to exclude
# anything with an existing checkpoint) using up to `n_workers` R worker
# processes running concurrently instead of one at a time. This is the
# single biggest speed lever in this script -- a real run showed one
# shard's glmer() fit alone takes 60-107 minutes, and shards are otherwise
# fit strictly one after another. Each worker fits its shard completely
# independently and saves ITS OWN checkpoint file directly to disk the
# moment it finishes (no shared state between workers, no write
# conflicts), so this is exactly as resumable as the serial version: if
# the whole run is killed mid-batch, every shard that had already finished
# on its worker still has a valid checkpoint on disk, and a re-run skips
# it. Falls back to a plain serial loop when n_workers <= 1 (e.g. on a
# single-core machine, or if you set N_PARALLEL_WORKERS <- 1L deliberately
# to reproduce the original one-at-a-time behavior).
fit_shards_parallel <- function(jobs, n_workers) {
  if (nrow(jobs) == 0L) return(invisible(NULL))
  n_workers <- max(1L, min(n_workers, nrow(jobs)))

  if (n_workers <= 1L) {
    for (i in seq_len(nrow(jobs))) {
      result <- fit_one_shard(jobs$shard_path[i], jobs$sid[i], jobs$shard_no[i])
      if (!is.null(result)) {
        saveRDS(result, jobs$checkpoint_path[i])
        log_line("  %s/shard %02d checkpointed to: %s", jobs$sid[i], jobs$shard_no[i], jobs$checkpoint_path[i])
      } else {
        log_line("  %s/shard %02d did not converge -- not checkpointed, will retry on next run.", jobs$sid[i], jobs$shard_no[i])
      }
    }
    return(invisible(NULL))
  }

  log_line("Fitting %d pending shard(s) using %d parallel worker(s) (dynamic queue -- shards checkpoint and log as each one finishes, not in dispatch order)...",
            nrow(jobs), n_workers)
  cl <- makePSOCKcluster(n_workers)
  on.exit(stopCluster(cl), add = TRUE)

  # Cap BLAS-level threading INSIDE each worker to 1 thread, before any
  # library that touches linear algebra loads. Without this, each of the
  # N_PARALLEL_WORKERS separate R processes may ALSO try to multi-thread
  # its own matrix operations (OpenBLAS/MKL do this by default on many
  # systems) -- N workers x several BLAS threads each can vastly
  # OVERSUBSCRIBE the machine's actual core count (e.g. 3 workers x 4
  # threads each = 12 threads competing for 4 physical cores), which can
  # make "parallel" fitting SLOWER than plain serial fitting, not faster.
  # This keeps each worker to exactly one thread, so N_PARALLEL_WORKERS
  # lines up with N_PARALLEL_WORKERS physical cores as intended.
  clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
               NUMEXPR_NUM_THREADS = "1")
    suppressPackageStartupMessages({
      library(data.table); library(lme4); library(broom.mixed)
    })
  })
  clusterExport(cl, c("fit_one_shard", "extract_random_intercept_variance", "calc_icc", "RESID_VAR", "log_line"))

  # Each worker fits its shard, saves ITS OWN checkpoint directly to disk
  # (no shared state, no write conflicts), and returns only a small status
  # list -- the (already tiny, tidy-coefficient-sized) result never has to
  # travel back over the cluster socket. Jobs are dispatched as a DYNAMIC
  # queue (not fixed batches): the moment any worker finishes, it's handed
  # the next pending shard immediately, and that shard's outcome is logged
  # right away -- so a slow shard (real per-shard times varied from ~60 to
  # ~107 minutes in testing) no longer blocks the whole batch's log output,
  # and no worker sits idle waiting for batch-mates to finish.
  worker_fit <- function(shard_path, sid, shard_no, checkpoint_path) {
    tryCatch({
      result <- fit_one_shard(shard_path, sid, shard_no)
      if (!is.null(result)) saveRDS(result, checkpoint_path)
      list(sid = sid, shard_no = shard_no, ok = !is.null(result))
    }, error = function(e) {
      list(sid = sid, shard_no = shard_no, ok = FALSE, error = conditionMessage(e))
    })
  }

  n_jobs <- nrow(jobs)
  n_start <- min(n_workers, n_jobs)
  for (w in seq_len(n_start)) {
    parallel:::sendCall(cl[[w]], worker_fit,
                         list(jobs$shard_path[w], jobs$sid[w], jobs$shard_no[w], jobs$checkpoint_path[w]))
  }
  next_job <- n_start + 1L
  completed <- 0L
  while (completed < n_jobs) {
    res <- parallel:::recvOneResult(cl)
    completed <- completed + 1L
    r <- res$value
    if (isTRUE(r$ok)) {
      log_line("  %s/shard %02d checkpointed. (%d/%d shards done)", r$sid, r$shard_no, completed, n_jobs)
    } else {
      log_line("  %s/shard %02d did not converge%s -- not checkpointed, will retry on next run. (%d/%d shards done)",
                r$sid, r$shard_no, if (!is.null(r$error)) paste0(" (", r$error, ")") else "", completed, n_jobs)
    }
    if (next_job <= n_jobs) {
      parallel:::sendCall(cl[[res$node]], worker_fit,
                           list(jobs$shard_path[next_job], jobs$sid[next_job], jobs$shard_no[next_job], jobs$checkpoint_path[next_job]))
      next_job <- next_job + 1L
    }
  }
  invisible(NULL)
}

# Orchestrates one subsample end-to-end: bucket into shards (skipped if all
# shard files already exist), fit (or skip, if already checkpointed) each
# shard one at a time, then pool that subsample's shard-level estimates
# into a single subsample-level estimate via the same meta-analysis used
# for the final pooling step -- so from here on, this subsample behaves
# exactly like the pre-sharding version of this script did.
process_subsample <- function(sid, id_list_path) {
  shard_paths <- vapply(seq_len(N_SHARDS_PER_SUBSAMPLE), function(s) shard_file_path(sid, s), character(1))
  checkpoint_paths <- vapply(seq_len(N_SHARDS_PER_SUBSAMPLE), function(s) shard_checkpoint_path(sid, s), character(1))

  if (all(file.exists(checkpoint_paths))) {
    # Every shard already fit on a previous run (shard scratch files were
    # likely already cleaned up too) -- no need to touch the source panel
    # at all, just load the checkpoints below.
    log_line("Subsample %s: all %d shards already checkpointed -- skipping bucketing entirely.", sid, N_SHARDS_PER_SUBSAMPLE)
  } else if (all(file.exists(shard_paths))) {
    log_line("Subsample %s: shard files already present -- skipping bucketing, fitting remaining shards.", sid)
  } else {
    ids_dt <- fread(id_list_path, showProgress = FALSE)
    if (!"ResearchID" %in% names(ids_dt)) {
      stop("Missing ResearchID column in ID list file:\n", id_list_path)
    }
    sampled_ids <- as.character(ids_dt$ResearchID)
    rm(ids_dt)
    log_line("Bucketing subsample %s into %d shards from %s (reads the full panel once; ~%d individuals expected)",
              sid, N_SHARDS_PER_SUBSAMPLE, basename(COLLAPSED_FILE), length(sampled_ids))
    bucket_subsample_into_shards(COLLAPSED_FILE, sampled_ids, REQUIRED_COLUMNS, sid, N_SHARDS_PER_SUBSAMPLE)
    log_line("Subsample %s: bucketing complete.", sid)
  }

  # Skip anything already checkpointed (from this run's own bucketing above,
  # or from an earlier -- possibly serial, pre-parallel -- run); fit the
  # rest as a batch, in parallel where possible (see fit_shards_parallel()).
  already_done <- file.exists(checkpoint_paths)
  for (s in which(already_done)) {
    log_line("  %s/shard %02d already checkpointed -- skipping.", sid, s)
  }
  pending_idx <- which(!already_done)
  # NOTE: `sid` must be built via rep(sid, length(pending_idx)), NOT passed
  # as the bare scalar `sid`. When pending_idx is empty (every shard for
  # this subsample already checkpointed -- a real case, e.g. resuming a run
  # that was interrupted after its last shard finished but before this
  # subsample's own pooled checkpoint was written), data.table() silently
  # recycles/pads a length-1 scalar column against length-0 vector columns
  # into a bogus 1-row table full of NAs instead of a proper 0-row table --
  # caught during testing, where it produced a "Fitting shard NA: NA" crash.
  pending_jobs <- data.table(
    sid = rep(sid, length(pending_idx)),
    shard_no = pending_idx,
    shard_path = shard_paths[pending_idx],
    checkpoint_path = checkpoint_paths[pending_idx]
  )
  fit_shards_parallel(pending_jobs, N_PARALLEL_WORKERS)
  gc()

  shard_results <- lapply(seq_len(N_SHARDS_PER_SUBSAMPLE), function(s) {
    if (file.exists(checkpoint_paths[s])) safe_read_rds(checkpoint_paths[s]) else NULL
  })
  shard_results <- Filter(Negate(is.null), shard_results)

  if (length(shard_results) == 0L) {
    log_line("Subsample %s: no shards converged -- cannot pool this subsample.", sid)
    return(NULL)
  }
  is_complete <- length(shard_results) == N_SHARDS_PER_SUBSAMPLE
  if (!is_complete) {
    log_line("Subsample %s: pooling %d of %d shards -- the rest have not converged/checkpointed yet (this subsample will NOT be marked done -- it will be retried on the next run until all shards succeed).",
              sid, length(shard_results), N_SHARDS_PER_SUBSAMPLE)
  }

  shard_coef_all <- rbindlist(lapply(shard_results, `[[`, "coef"), use.names = TRUE, fill = TRUE)
  shard_fit_all <- rbindlist(lapply(shard_results, `[[`, "fit"), use.names = TRUE, fill = TRUE)
  shard_re_all <- rbindlist(lapply(shard_results, `[[`, "re"), use.names = TRUE, fill = TRUE)

  # Inner pooling: this subsample's N shard-level estimates -> one
  # subsample-level estimate (same rma.uni machinery as the outer pool).
  subsample_pooled <- meta_pool_coefficients(shard_coef_all)
  subsample_coef_dt <- data.table(
    term = subsample_pooled$term,
    estimate = subsample_pooled$pooled_est,
    std.error = subsample_pooled$pooled_se,
    subset = sid
  )

  # Only clean up shard scratch files once EVERY shard has genuinely
  # succeeded -- if we're incomplete (a shard hasn't converged yet, or its
  # checkpoint was just found corrupt and deleted by safe_read_rds()),
  # deleting its still-good bucket file too would force a wasteful full
  # re-bucket of the whole subsample on the next run for no reason.
  if (is_complete && !KEEP_SHARD_FILES) {
    unlink(shard_paths[file.exists(shard_paths)])
  }

  list(coef = subsample_coef_dt, shard_coef = shard_coef_all, shard_fit = shard_fit_all, shard_re = shard_re_all,
       complete = is_complete)
}

log_line("Checking sensitivity input files.")
stop_if_missing_files(c(COLLAPSED_FILE, ID_LIST_FILES))

subsample_checkpoint_path <- function(subset_id) file.path(OUT_DIR, "subsample_pooled_checkpoints", sprintf("subsample_%s.rds", subset_id))
dir.create(dirname(subsample_checkpoint_path("00")), recursive = TRUE, showWarnings = FALSE)

# Process (or skip, if already checkpointed) one subsample at a time --
# each subsample internally shards, fits, and inner-pools (see
# process_subsample() above) before its own result is checkpointed here.
for (i in seq_along(ID_LIST_FILES)) {
  sid <- SUBSAMPLE_IDS[i]
  cp_path <- subsample_checkpoint_path(sid)

  if (file.exists(cp_path)) {
    log_line("Subsample %s already fully pooled -- skipping (delete %s to force a redo).", sid, basename(cp_path))
    next
  }

  result <- process_subsample(sid, ID_LIST_FILES[i])
  if (!is.null(result) && isTRUE(result$complete)) {
    # Only persist the subsample-level "done" checkpoint once every shard
    # genuinely succeeded -- a partial result (some shards still pending,
    # or a checkpoint dropped for being corrupt) is deliberately NOT saved
    # here, so this subsample is picked back up and finished on the next
    # run instead of being silently, permanently marked complete while
    # quietly missing a shard's contribution.
    saveRDS(result, cp_path)
    log_line("Subsample %s pooled and checkpointed to: %s", sid, cp_path)
  } else if (!is.null(result)) {
    log_line("Subsample %s: partially pooled but NOT checkpointed (some shards still pending) -- will retry remaining shards on next run.", sid)
  } else {
    log_line("Subsample %s: no shards converged -- not checkpointed, will retry on next run.", sid)
  }
  rm(result)
  gc()
}

# Load every available subsample-level checkpoint for the final pooling.
results_list <- lapply(SUBSAMPLE_IDS, function(sid) {
  cp_path <- subsample_checkpoint_path(sid)
  if (file.exists(cp_path)) safe_read_rds(cp_path) else NULL
})
results_list <- Filter(Negate(is.null), results_list)

if (length(results_list) == 0L) {
  stop("No subsample models converged; nothing to pool.")
}
if (length(results_list) < N_SUBSAMPLES) {
  log_line("NOTE: pooling %d of %d subsamples -- %d have not converged/checkpointed yet.",
            length(results_list), N_SUBSAMPLES, N_SUBSAMPLES - length(results_list))
}

coef_all <- rbindlist(lapply(results_list, `[[`, "coef"), use.names = TRUE, fill = TRUE)
shard_coef_all <- rbindlist(lapply(results_list, `[[`, "shard_coef"), use.names = TRUE, fill = TRUE)
shard_fit_all <- rbindlist(lapply(results_list, `[[`, "shard_fit"), use.names = TRUE, fill = TRUE)
shard_re_all <- rbindlist(lapply(results_list, `[[`, "shard_re"), use.names = TRUE, fill = TRUE)

# Shard-level detail (one row per shard -- ~160 rows by default) for full
# transparency/diagnostics, plus the subsample-level coefficients (one row
# per subsample per term, each already pooled across that subsample's shards).
fwrite(shard_coef_all, file.path(OUT_DIR, "sensitivity_shard_coefficients.csv"))
fwrite(shard_fit_all, file.path(OUT_DIR, "sensitivity_shard_fit_stats.csv"))
fwrite(shard_re_all, file.path(OUT_DIR, "sensitivity_shard_random_effects.csv"))
fwrite(coef_all, file.path(OUT_DIR, "sensitivity_subsample_coefficients.csv"))

log_line("Pooling %d subsample-level estimates (each already pooled across its own shards) via random-effects meta-analysis.", length(results_list))
pooled_dt <- meta_pool_coefficients(coef_all)
fwrite(pooled_dt, file.path(OUT_DIR, "sensitivity_pooled_coefficients.csv"))

lookback_row <- pooled_dt[term == "LookBack_z"]
comparison_dt <- data.table(
  model = c("Primary (GLMM, multiple events/person-year allowed)", "Primary (GLM)", "Sensitivity (GLMM, one event/person-year)"),
  lookback_OR = c(PRIMARY_LOOKBACK_OR_GLMM, PRIMARY_LOOKBACK_OR_GLM, lookback_row$OR),
  lookback_OR_low = c(NA_real_, NA_real_, lookback_row$OR_low),
  lookback_OR_high = c(NA_real_, NA_real_, lookback_row$OR_high)
)
fwrite(comparison_dt, file.path(OUT_DIR, "sensitivity_vs_primary_comparison.csv"))

log_line("H1 one-event-per-year sensitivity analysis complete.")
print(pooled_dt)
print(comparison_dt)
