############################################################
# SCRIPT 06 — H2 HISTORICAL PERIOD (PRE-2001 vs. 2001+)
#
# STATUS
#   Added post-R&R. Originally developed as a reviewer-requested
#   sensitivity check on historical period, this comparison was
#   elevated by the editor to a PRIMARY H2 subgroup analysis
#   (see manuscript response letter). It now reproduces the
#   manuscript's Table 6 and Figure 8, alongside the state,
#   demographic, and offense/incarceration comparisons in
#   04_H2_invariance.R (Tables 3-5, Figures 5-7).
#
# PURPOSE
#   Test whether the rate of recidivism prediction decay differs
#   for prior charge events occurring before 2001 versus in 2001
#   or later, using the same additive-vs-interaction GLMM logic
#   as 04_H2_invariance.R.
#
# MODEL LOGIC
#   Additive:
#     cbind(y, fail) ~ After2000 + LookBack_z + (1 | ResearchID)
#
#   Interaction:
#     cbind(y, fail) ~ After2000 * LookBack_z + (1 | ResearchID)
#
# INTERPRETATION
#   After2000 main effect:
#     Difference in baseline recidivism risk between period groups.
#
#   After2000 x LookBack_z:
#     Difference in the decay rate between period groups.
#     If OR is close to 1.00, the decay slope is substantively similar.
#
# INPUT
#   /data/sensitivity/cohort_after2000.csv
#
# REQUIRED COLUMNS
#   ResearchID, State, LookBack, Recid, After2000
#     After2000: 0 = prior event before 2001, 1 = prior event in 2001+
#
# COMPUTATIONAL NOTE
#   As with the primary H1/H2 models, individuals are split into
#   20 deterministic partitions and pooled via random-effects
#   meta-analysis (metafor::rma.uni) for computational tractability
#   on the full population. Partition results are cached to disk
#   so an interrupted run can resume without refitting completed
#   partitions.
#
# OUTPUT
#   /output/06_H2_historical_period/
#     - after2000_partition_manifest.csv
#     - after2000_fixed_effects_by_partition.csv
#     - after2000_pooled_fixed_effects.csv
#     - after2000_table6_publication_table.csv
#     - after2000_fit_stats_by_partition.csv
#     - after2000_fit_summary.csv
#     - after2000_random_effects_by_partition.csv
#     - after2000_random_effects_summary.csv
#     - after2000_lrt_by_partition.csv
#     - after2000_curve_predictions.csv
#     - H2_figure8_historical_period_probability.png
#     - H2_figure8_historical_period_riskratio.png
#     - H2_historical_period_results_workbook.xlsx
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)
options(scipen = 999)
options(contrasts = c("contr.treatment", "contr.poly"))

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(broom.mixed)
  library(metafor)
  library(ggplot2)
  library(scales)
  library(openxlsx)
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
OUT_DIR <- file.path(REPO_DIR, "output", "06_H2_historical_period")
PARTITION_DIR <- file.path(OUT_DIR, "partitions")
RESULT_DIR <- file.path(OUT_DIR, "partition_results")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PARTITION_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# cohort_after2000.csv turned out to be ~757 million rows in the real data
# (see export_sensitivity_data_from_spss.R) -- too large to write
# uncompressed without running out of local disk space, so it's exported
# as .csv.gz. fread() reads .csv.gz natively (requires the R.utils
# package). Falls back to a plain .csv for anyone who exported/hosted it
# uncompressed instead.
INPUT_FILE_GZ <- file.path(SENSITIVITY_DIR, "cohort_after2000.csv.gz")
INPUT_FILE_PLAIN <- file.path(SENSITIVITY_DIR, "cohort_after2000.csv")
INPUT_FILE <- if (file.exists(INPUT_FILE_GZ)) INPUT_FILE_GZ else INPUT_FILE_PLAIN

# IMPORTANT if INPUT_FILE is .gz: fread() always fully decompresses a .gz
# file to a scratch temp file before parsing it (this is how data.table's
# R.utils fallback works -- there is no partial/streaming read, even when
# only previewing a few rows), so reading this file needs roughly as much
# free space in FREAD_TMPDIR as the UNCOMPRESSED file would take, however
# briefly, on top of whatever's already used by the compressed copy in
# data/sensitivity/. Defaults to R's normal tempdir(); override with a
# path on a drive that has more room if needed (e.g. "D:/Rtemp").
FREAD_TMPDIR <- tempdir()

REQUIRED_COLUMNS <- c("ResearchID", "State", "LookBack", "Recid", "After2000")

N_PARTITIONS <- 20L
WRITE_PARTITION_FILES <- TRUE

# Common observable lookback window for both historical periods, so pre-2001
# and post-2000 events are compared over the same follow-up horizon.
MAX_LOOKBACK <- 20L

# Reference group for historical-period comparison.
PERIOD_REFERENCE <- "Pre-2001"

BASE_FAMILY <- "Times New Roman"
PLOT_WIDTH <- 8.5
PLOT_HEIGHT <- 5.5
PLOT_DPI <- 300

USE_COMPRESSION <- TRUE
RESID_VAR <- (pi^2) / 3

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

validate_columns <- function(dt, required_cols, file_path = "data") {
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }
}

format_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "< .001", sprintf("%.3f", p)))
}

format_or_ci <- function(or, low, high) {
  ifelse(is.na(or), "", sprintf("%.2f (%.2f, %.2f)", or, low, high))
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

compress_binomial_counts <- function(dt, by_vars) {
  comp <- dt[, .(y = sum(Recid, na.rm = TRUE), n = .N), by = by_vars]
  comp[, fail := n - y]
  comp[]
}

extract_fixed_effects <- function(model, model_name) {
  tt <- as.data.table(broom.mixed::tidy(model, effects = "fixed"))
  tt[, `:=`(
    model = model_name,
    OR = exp(estimate),
    OR_low = exp(estimate - 1.96 * std.error),
    OR_high = exp(estimate + 1.96 * std.error)
  )]
  tt[]
}

compute_fit_stats <- function(model, model_name, dt_uncompressed) {
  p_hat <- predict(model, newdata = dt_uncompressed, type = "response", re.form = NA)
  y_obs <- dt_uncompressed$Recid
  brier <- mean((y_obs - p_hat)^2, na.rm = TRUE)
  rmse <- sqrt(brier)
  ll <- as.numeric(logLik(model))
  data.table(
    model = model_name,
    minus2LL = -2 * ll,
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = ll,
    RMSE = rmse,
    Brier = brier,
    n_obs = nrow(dt_uncompressed),
    n_ids = uniqueN(dt_uncompressed$ResearchID)
  )
}

fit_glmm_pair <- function(dt_comp) {
  f_add <- cbind(y, fail) ~ After2000 + LookBack_z + (1 | ResearchID)
  f_int <- cbind(y, fail) ~ After2000 * LookBack_z + (1 | ResearchID)

  m_add <- glmer(
    f_add, data = dt_comp, family = binomial(link = "logit"), nAGQ = 0L,
    control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE, optCtrl = list(maxfun = 2e5))
  )
  m_int <- glmer(
    f_int, data = dt_comp, family = binomial(link = "logit"), nAGQ = 0L,
    control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE, optCtrl = list(maxfun = 2e5))
  )

  list(additive = m_add, interaction = m_int)
}

meta_pool_coefficients <- function(coef_dt) {
  coef_dt <- copy(coef_dt)
  coef_dt[, vi := std.error^2]
  combos <- unique(coef_dt[, .(model, term)])

  pooled_list <- lapply(seq_len(nrow(combos)), function(i) {
    m <- combos$model[i]
    t <- combos$term[i]
    sub <- coef_dt[model == m & term == t]

    fit <- tryCatch(
      metafor::rma.uni(yi = sub$estimate, vi = sub$vi, method = "REML"),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(data.table(
        model = m, term = t, pooled_est = NA_real_, pooled_se = NA_real_,
        z = NA_real_, p = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        tau2 = NA_real_, I2 = NA_real_, Q = NA_real_, Q_p = NA_real_
      ))
    }

    data.table(
      model = m, term = t,
      pooled_est = as.numeric(fit$b), pooled_se = as.numeric(fit$se),
      z = as.numeric(fit$zval), p = as.numeric(fit$pval),
      ci_low = as.numeric(fit$ci.lb), ci_high = as.numeric(fit$ci.ub),
      tau2 = as.numeric(fit$tau2), I2 = as.numeric(fit$I2),
      Q = as.numeric(fit$QE), Q_p = as.numeric(fit$QEp)
    )
  })

  pooled_dt <- rbindlist(pooled_list, use.names = TRUE, fill = TRUE)
  pooled_dt[, `:=`(OR = exp(pooled_est), OR_low = exp(ci_low), OR_high = exp(ci_high))]
  pooled_dt[]
}

make_term_labels <- function(x) {
  x <- gsub("\\(Intercept\\)", "Intercept", x)
  x <- gsub("^LookBack_z$", "Lookback", x)
  x <- gsub("LookBack_z", "Lookback", x, fixed = TRUE)
  # The fitted factor level is literally "Post-2000 event", so glmer's
  # coefficient name is "After2000Post-2000 event" -- replace the WHOLE
  # thing before the generic fallback below, or the fallback doubles it
  # into "Post-2000 eventPost-2000 event".
  x <- gsub("After2000Post-2000 event", "Post-2000 event", x, fixed = TRUE)
  x <- gsub("After2000Pre-2001 event", "Pre-2001 event", x, fixed = TRUE)
  x <- gsub("After2000", "Post-2000 event", x, fixed = TRUE)
  x <- gsub(":", " x ", x, fixed = TRUE)
  x
}

save_plot <- function(p, filename) {
  ggsave(filename = file.path(OUT_DIR, filename), plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
}

############################################################
# 1) READ INPUT AND CREATE PARTITIONS
############################################################

PARTITION_IDS <- sprintf("%02d", 0:(N_PARTITIONS - 1L))
PARTITION_FILES <- file.path(PARTITION_DIR, sprintf("after2000_partition_%s.csv", PARTITION_IDS))
RESULT_FILES <- file.path(RESULT_DIR, sprintf("after2000_result_partition_%s.rds", PARTITION_IDS))

REUSE_EXISTING_PARTITIONS <- all(file.exists(PARTITION_FILES))
if (REUSE_EXISTING_PARTITIONS) {
  log_line("Found all %d partition files. Reusing existing partitions.", N_PARTITIONS)
} else {
  log_line("Reading historical-period input file: %s", INPUT_FILE)
  stop_if_missing_files(INPUT_FILE)

  dt <- fread(INPUT_FILE, showProgress = FALSE, tmpdir = FREAD_TMPDIR)
  validate_columns(dt, REQUIRED_COLUMNS, INPUT_FILE)

  dt[, ResearchID := as.character(ResearchID)]
  dt[, State := trimws(as.character(State))]
  dt[, LookBack := suppressWarnings(as.integer(LookBack))]
  dt[, Recid := suppressWarnings(as.integer(Recid))]
  dt[, After2000 := suppressWarnings(as.integer(After2000))]

  dt <- dt[!is.na(ResearchID) & !is.na(LookBack) & !is.na(Recid) & !is.na(After2000)]
  dt <- dt[Recid %in% c(0L, 1L)]
  dt <- dt[After2000 %in% c(0L, 1L)]
  if (!is.na(MAX_LOOKBACK)) dt <- dt[LookBack <= MAX_LOOKBACK]

  dt[, After2000 := factor(
    fifelse(After2000 == 1L, "Post-2000 event", "Pre-2001 event"),
    levels = c("Pre-2001 event", "Post-2000 event")
  )]
  if (PERIOD_REFERENCE == "Post-2000") {
    dt[, After2000 := relevel(After2000, ref = "Post-2000 event")]
  } else {
    dt[, After2000 := relevel(After2000, ref = "Pre-2001 event")]
  }

  log_line("Assigning deterministic partitions by ResearchID.")
  id_dt <- unique(dt[, .(ResearchID)])
  setorder(id_dt, ResearchID)
  id_dt[, partition_num := ((seq_len(.N) - 1L) %% N_PARTITIONS)]
  id_dt[, partition := sprintf("%02d", partition_num)]
  id_dt[, partition_num := NULL]

  dt <- merge(dt, id_dt, by = "ResearchID", all.x = TRUE, sort = FALSE)

  manifest <- dt[, .(
    n_rows = .N,
    n_unique_ids = uniqueN(ResearchID),
    lookback_min = min(LookBack, na.rm = TRUE),
    lookback_max = max(LookBack, na.rm = TRUE),
    recid_events = sum(Recid, na.rm = TRUE),
    pre2001_rows = sum(After2000 == "Pre-2001 event", na.rm = TRUE),
    post2000_rows = sum(After2000 == "Post-2000 event", na.rm = TRUE)
  ), by = partition][order(partition)]

  fwrite(manifest, file.path(OUT_DIR, "after2000_partition_manifest.csv"))

  if (WRITE_PARTITION_FILES) {
    log_line("Writing %d partition files.", N_PARTITIONS)
    for (sid in PARTITION_IDS) {
      f <- file.path(PARTITION_DIR, sprintf("after2000_partition_%s.csv", sid))
      fwrite(dt[partition == sid, .(ResearchID, State, LookBack, Recid, After2000)], f)
    }
  }
  rm(dt, id_dt)
  gc()
}

stop_if_missing_files(PARTITION_FILES)

############################################################
# 2) FIT ONE PARTITION
############################################################

fit_one_partition <- function(file_path, partition_id) {
  log_line("Fitting partition %s", partition_id)

  dt <- fread(file_path, showProgress = FALSE)

  dt[, ResearchID := as.character(ResearchID)]
  dt[, LookBack := suppressWarnings(as.integer(LookBack))]
  dt[, Recid := suppressWarnings(as.integer(Recid))]

  if (is.numeric(dt$After2000) || is.integer(dt$After2000)) {
    dt[, After2000 := factor(
      fifelse(as.integer(After2000) == 1L, "Post-2000 event", "Pre-2001 event"),
      levels = c("Pre-2001 event", "Post-2000 event")
    )]
  } else {
    dt[, After2000 := trimws(as.character(After2000))]
    dt[After2000 %in% c("0", "Pre2001", "Pre-2001", "Pre-2001 event"), After2000 := "Pre-2001 event"]
    dt[After2000 %in% c("1", "After2000", "Post2000", "Post-2000", "Post-2000 event"), After2000 := "Post-2000 event"]
    dt[, After2000 := factor(After2000, levels = c("Pre-2001 event", "Post-2000 event"))]
  }
  if (PERIOD_REFERENCE == "Post-2000") {
    dt[, After2000 := relevel(After2000, ref = "Post-2000 event")]
  } else {
    dt[, After2000 := relevel(After2000, ref = "Pre-2001 event")]
  }

  dt <- dt[!is.na(ResearchID) & !is.na(LookBack) & !is.na(Recid) & !is.na(After2000)]
  dt <- dt[Recid %in% c(0L, 1L)]
  if (!is.na(MAX_LOOKBACK)) dt <- dt[LookBack <= MAX_LOOKBACK]

  mu_lb <- mean(dt$LookBack, na.rm = TRUE)
  sd_lb <- sd(dt$LookBack, na.rm = TRUE)
  if (!is.finite(sd_lb) || sd_lb == 0) stop("LookBack SD is invalid in partition ", partition_id)
  dt[, LookBack_z := (LookBack - mu_lb) / sd_lb]

  dt_comp <- if (USE_COMPRESSION) {
    compress_binomial_counts(dt, c("ResearchID", "After2000", "LookBack", "LookBack_z"))
  } else {
    tmp <- copy(dt)
    tmp[, y := Recid]
    tmp[, n := 1L]
    tmp[, fail := 1L - Recid]
    tmp[, .(ResearchID, After2000, LookBack, LookBack_z, y, n, fail)]
  }
  dt_comp[, ResearchID := factor(ResearchID)]
  dt_comp[, After2000 := factor(After2000, levels = levels(dt$After2000))]

  mods <- fit_glmm_pair(dt_comp)

  coef_dt <- rbindlist(
    list(extract_fixed_effects(mods$additive, "additive"), extract_fixed_effects(mods$interaction, "interaction")),
    use.names = TRUE, fill = TRUE
  )
  coef_dt[, partition := partition_id]

  fit_dt <- rbindlist(
    list(compute_fit_stats(mods$additive, "additive", dt), compute_fit_stats(mods$interaction, "interaction", dt)),
    use.names = TRUE, fill = TRUE
  )
  fit_dt[, partition := partition_id]

  var_add <- extract_random_intercept_variance(mods$additive)
  var_int <- extract_random_intercept_variance(mods$interaction)
  re_dt <- rbindlist(
    list(
      data.table(partition = partition_id, model = "additive", var_u0 = var_add, sd_u0 = sqrt(var_add), ICC = calc_icc(var_add)),
      data.table(partition = partition_id, model = "interaction", var_u0 = var_int, sd_u0 = sqrt(var_int), ICC = calc_icc(var_int))
    ),
    use.names = TRUE, fill = TRUE
  )

  lrt_dt <- as.data.table(anova(mods$additive, mods$interaction, test = "Chisq"))
  lrt_dt[, partition := partition_id]

  curve_grid <- unique(dt[, .(After2000, LookBack)])
  setorder(curve_grid, After2000, LookBack)
  curve_grid[, LookBack_z := (LookBack - mu_lb) / sd_lb]
  curve_grid[, ResearchID := dt$ResearchID[1]]
  curve_grid[, p_additive := as.numeric(predict(mods$additive, newdata = curve_grid, type = "response", re.form = NA))]
  curve_grid[, p_interaction := as.numeric(predict(mods$interaction, newdata = curve_grid, type = "response", re.form = NA))]
  curve_grid[, rr_additive := p_additive / p_additive[LookBack == min(LookBack, na.rm = TRUE)][1], by = After2000]
  curve_grid[, rr_interaction := p_interaction / p_interaction[LookBack == min(LookBack, na.rm = TRUE)][1], by = After2000]
  curve_grid[, partition := partition_id]

  list(coef = coef_dt, fit = fit_dt, re = re_dt, lrt = lrt_dt, curves = curve_grid)
}

############################################################
# 3) RUN MODELS ACROSS PARTITIONS WITH RESTART PROTECTION
############################################################

log_line("Existing completed result files: %d of %d.", sum(file.exists(RESULT_FILES)), length(RESULT_FILES))

for (i in seq_along(PARTITION_FILES)) {
  partition_id <- PARTITION_IDS[i]
  result_file <- RESULT_FILES[i]

  if (file.exists(result_file)) {
    log_line("Skipping partition %s; saved result already exists.", partition_id)
    next
  }

  res_i <- fit_one_partition(PARTITION_FILES[i], partition_id)
  saveRDS(res_i, result_file)
  rm(res_i)
  gc()
}

missing_results <- RESULT_FILES[!file.exists(RESULT_FILES)]
if (length(missing_results) > 0L) {
  stop("The following partition result files are still missing:\n", paste(missing_results, collapse = "\n"))
}

log_line("Loading all saved partition results.")
results <- lapply(RESULT_FILES, readRDS)

coef_all <- rbindlist(lapply(results, `[[`, "coef"), use.names = TRUE, fill = TRUE)
fit_all <- rbindlist(lapply(results, `[[`, "fit"), use.names = TRUE, fill = TRUE)
re_all <- rbindlist(lapply(results, `[[`, "re"), use.names = TRUE, fill = TRUE)
lrt_all <- rbindlist(lapply(results, `[[`, "lrt"), use.names = TRUE, fill = TRUE)
curves_all <- rbindlist(lapply(results, `[[`, "curves"), use.names = TRUE, fill = TRUE)

fwrite(coef_all, file.path(OUT_DIR, "after2000_fixed_effects_by_partition.csv"))
fwrite(fit_all, file.path(OUT_DIR, "after2000_fit_stats_by_partition.csv"))
fwrite(re_all, file.path(OUT_DIR, "after2000_random_effects_by_partition.csv"))
fwrite(lrt_all, file.path(OUT_DIR, "after2000_lrt_by_partition.csv"))

############################################################
# 4) POOL AND SUMMARIZE
############################################################

pooled_coef <- meta_pool_coefficients(coef_all)
pooled_coef[, term_label := make_term_labels(term)]
fwrite(pooled_coef, file.path(OUT_DIR, "after2000_pooled_fixed_effects.csv"))

fit_summary <- fit_all[, .(
  mean_minus2LL = mean(minus2LL, na.rm = TRUE),
  mean_AIC = mean(AIC, na.rm = TRUE),
  mean_BIC = mean(BIC, na.rm = TRUE),
  mean_RMSE = mean(RMSE, na.rm = TRUE),
  mean_Brier = mean(Brier, na.rm = TRUE),
  total_n_obs = sum(n_obs, na.rm = TRUE),
  total_n_id_entries = sum(n_ids, na.rm = TRUE)
), by = model][order(model)]
fwrite(fit_summary, file.path(OUT_DIR, "after2000_fit_summary.csv"))

re_summary <- re_all[, .(
  mean_var_u0 = mean(var_u0, na.rm = TRUE),
  mean_sd_u0 = mean(sd_u0, na.rm = TRUE),
  mean_ICC = mean(ICC, na.rm = TRUE)
), by = model][order(model)]
fwrite(re_summary, file.path(OUT_DIR, "after2000_random_effects_summary.csv"))

# Table 6 in the published manuscript (historical-period GLMM estimates).
table6_dt <- pooled_coef[, .(
  Model = fifelse(model == "additive", "Additive", "Interaction"),
  Predictor = term_label,
  `Logit (SE)` = sprintf("%.2f (%.2f)", pooled_est, pooled_se),
  `OR (95% CI)` = format_or_ci(OR, OR_low, OR_high),
  `p-value` = format_p(p)
)]
predictor_order <- c("Intercept", "Lookback", "Post-2000 event", "Post-2000 event x Lookback")
table6_dt[, predictor_sort := match(Predictor, predictor_order)]
setorder(table6_dt, Model, predictor_sort)
table6_dt[, predictor_sort := NULL]
fwrite(table6_dt, file.path(OUT_DIR, "after2000_table6_publication_table.csv"))

partition_weights <- fit_all[, .(n_obs = unique(n_obs)[1]), by = partition]
curves_long <- rbindlist(
  list(
    curves_all[, .(partition, LookBack, After2000, model = "additive", p_hat = p_additive, RR = rr_additive)],
    curves_all[, .(partition, LookBack, After2000, model = "interaction", p_hat = p_interaction, RR = rr_interaction)]
  ),
  use.names = TRUE, fill = TRUE
)
curves_long <- merge(curves_long, partition_weights, by = "partition", all.x = TRUE)
curve_summary <- curves_long[, .(
  p_hat = weighted.mean(p_hat, w = n_obs, na.rm = TRUE),
  RR = weighted.mean(RR, w = n_obs, na.rm = TRUE)
), by = .(LookBack, After2000, model)][order(model, After2000, LookBack)]
fwrite(curve_summary, file.path(OUT_DIR, "after2000_curve_predictions.csv"))

############################################################
# 5) FIGURE 8 (PUBLISHED MANUSCRIPT FIGURE)
############################################################

plot_dt <- curve_summary[model == "interaction"]
plot_dt[, After2000 := factor(as.character(After2000), levels = c("Pre-2001 event", "Post-2000 event"))]

p_prob <- ggplot(plot_dt, aes(x = LookBack, y = p_hat, linetype = After2000)) +
  geom_line(linewidth = 1.0) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  labs(
    x = "Lookback year", y = "Predicted recidivism probability", linetype = NULL,
    title = "Figure 8A. Predicted recidivism probability decay by historical period"
  ) +
  theme_classic(base_family = BASE_FAMILY) +
  theme(plot.title = element_text(face = "bold", size = 13), legend.position = "bottom")

p_rr <- ggplot(plot_dt, aes(x = LookBack, y = RR, linetype = After2000)) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_line(linewidth = 1.0) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  labs(
    x = "Lookback year", y = "Risk ratio relative to first observed lookback year", linetype = NULL,
    title = "Figure 8B. Risk-ratio decay by historical period"
  ) +
  theme_classic(base_family = BASE_FAMILY) +
  theme(plot.title = element_text(face = "bold", size = 13), legend.position = "bottom")

save_plot(p_prob, "H2_figure8_historical_period_probability.png")
save_plot(p_rr, "H2_figure8_historical_period_riskratio.png")

############################################################
# 6) EXCEL WORKBOOK
############################################################

wb <- createWorkbook()
addWorksheet(wb, "Table6_HistoricalPeriod"); writeDataTable(wb, "Table6_HistoricalPeriod", table6_dt)
addWorksheet(wb, "PooledFixedEffects"); writeDataTable(wb, "PooledFixedEffects", pooled_coef)
addWorksheet(wb, "FitSummary"); writeDataTable(wb, "FitSummary", fit_summary)
addWorksheet(wb, "RandomEffects"); writeDataTable(wb, "RandomEffects", re_summary)
addWorksheet(wb, "LRTByPartition"); writeDataTable(wb, "LRTByPartition", lrt_all)
addWorksheet(wb, "CurvePredictions"); writeDataTable(wb, "CurvePredictions", curve_summary)

manifest_file <- file.path(OUT_DIR, "after2000_partition_manifest.csv")
if (file.exists(manifest_file)) {
  addWorksheet(wb, "PartitionManifest")
  writeDataTable(wb, "PartitionManifest", fread(manifest_file))
}

saveWorkbook(wb, file.path(OUT_DIR, "H2_historical_period_results_workbook.xlsx"), overwrite = TRUE)

log_line("H2 historical-period analysis (Table 6 / Figure 8) complete.")
print(table6_dt)
print(fit_summary)
