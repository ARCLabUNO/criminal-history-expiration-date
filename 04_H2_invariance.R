############################################################
# SCRIPT 04 — H2 INVARIANCE OF DECAY
#
# PURPOSE
#   Evaluate whether the rate of recidivism decay is consistent
#   across:
#     1. states
#     2. sex
#     3. race/ethnicity
#     4. age group
#     5. offense type
#     6. incarceration exposure
#
# INPUT
#   /data/subsamples/subsample_00_raw.csv ... subsample_19_raw.csv
#   /data/wa/WA_decay_dataset.csv
#
# OUTPUT
#   /output/04_H2_invariance/
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
  library(patchwork)
  library(openxlsx)
})

############################################################
# PATH SETUP — AUTO-DETECT REPO ROOT
############################################################

# Tries, in order: (1) how Rscript was invoked from a command line, (2)
# RStudio's "Source" button / Ctrl+Shift+Enter (both of these are exact --
# they don't depend on the working directory at all), (3) the original
# walk-up-from-cwd strategy as a fallback for other invocation styles.
# NOTE: strategy (3) alone (the original implementation) required an
# "output" folder to already exist, which fails on a first-ever run of a
# freshly cloned repo before anything has written output yet -- relaxed
# below to only require "data" + one of its known subfolders.
find_repo_root <- function(start_dir = getwd(), max_up = 6L) {
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

  cur <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

  for (i in 0:max_up) {
    has_data <- dir.exists(file.path(cur, "data"))
    has_subsamples <- dir.exists(file.path(cur, "data", "subsamples"))
    has_wa <- dir.exists(file.path(cur, "data", "wa"))
    has_sensitivity <- dir.exists(file.path(cur, "data", "sensitivity"))

    if (has_data && (has_subsamples || has_wa || has_sensitivity)) {
      return(cur)
    }

    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }

  stop(
    "Could not locate project root.\n",
    "Fix: open this script directly in RStudio and click 'Source' (Session >\n",
    "Set Working Directory > To Source File Location also works), or run\n",
    "setwd(\"<path to that folder>\") yourself before sourcing this script.\n",
    "Current directory:\n",
    normalizePath(start_dir)
  )
}

REPO_DIR <- find_repo_root()

DATA_DIR <- file.path(REPO_DIR, "data")
SUBSAMPLE_DIR <- file.path(DATA_DIR, "subsamples")
WA_DIR <- file.path(DATA_DIR, "wa")

OUT_DIR <- file.path(REPO_DIR, "output", "04_H2_invariance")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SUBSAMPLE_IDS <- sprintf("%02d", 0:19)
SUBSAMPLE_FILES <- file.path(SUBSAMPLE_DIR, sprintf("subsample_%s_raw.csv", SUBSAMPLE_IDS))
WA_FILE <- file.path(WA_DIR, "WA_decay_dataset.csv")

REQUIRED_STATE_COLUMNS <- c("ResearchID", "State", "LookBack", "Recid")
REQUIRED_WA_COLUMNS <- c("ResearchID", "Recid", "LookBack", "Male", "RaceEthnicity2", "AgeCurveGrouped", "CharType", "PriorIncarYearsCount")

EXPECTED_STATES <- c("AK", "AZ", "FL", "TX", "WA")
STATE_REFERENCE <- "FL"

ID_VAR <- "ResearchID"
Y_VAR <- "Recid"
LB_VAR <- "LookBack"

BASE_FAMILY <- "Times New Roman"
RESID_VAR <- (pi^2) / 3
USE_COMPRESSION <- TRUE

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

format_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "< .001", sprintf("%.3f", p)))
}

extract_random_intercept_variance <- function(model, id_var = ID_VAR) {
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
  comp <- dt[, .(y = sum(get(Y_VAR), na.rm = TRUE), n = .N), by = by_vars]
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
  y_obs <- dt_uncompressed[[Y_VAR]]
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
    n_ids = uniqueN(dt_uncompressed[[ID_VAR]])
  )
}

fit_glmm_pair <- function(dt_comp, formula_add, formula_int) {
  m_add <- glmer(
    formula_add,
    data = dt_comp,
    family = binomial(link = "logit"),
    nAGQ = 0L,
    control = glmerControl(
      optimizer = "bobyqa",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 2e5)
    )
  )
  
  m_int <- glmer(
    formula_int,
    data = dt_comp,
    family = binomial(link = "logit"),
    nAGQ = 0L,
    control = glmerControl(
      optimizer = "bobyqa",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 2e5)
    )
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
        model = m,
        term = t,
        pooled_est = NA_real_,
        pooled_se = NA_real_,
        z = NA_real_,
        p = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      ))
    }
    
    data.table(
      model = m,
      term = t,
      pooled_est = as.numeric(fit$b),
      pooled_se = as.numeric(fit$se),
      z = as.numeric(fit$zval),
      p = as.numeric(fit$pval),
      ci_low = as.numeric(fit$ci.lb),
      ci_high = as.numeric(fit$ci.ub)
    )
  })
  
  pooled_dt <- rbindlist(pooled_list, use.names = TRUE, fill = TRUE)
  pooled_dt[, `:=`(
    OR = exp(pooled_est),
    OR_low = exp(ci_low),
    OR_high = exp(ci_high)
  )]
  pooled_dt[]
}

make_term_labels <- function(x) {
  x <- gsub("^State", "", x)
  x <- gsub("^LookBack_z$", "LookBack", x)
  x <- gsub("LookBack_z", "LookBack", x, fixed = TRUE)
  x <- gsub(":", " × ", x, fixed = TRUE)
  x <- gsub("\\(Intercept\\)", "Intercept", x)
  x
}

stop_if_missing_files(SUBSAMPLE_FILES)
if (!file.exists(WA_FILE)) stop("Washington file not found:\n", WA_FILE)

fit_one_state_subsample <- function(file_path, subset_id) {
  dt <- fread(file_path, showProgress = FALSE)
  validate_columns(dt, REQUIRED_STATE_COLUMNS, file_path)
  
  dt[, ResearchID := as.character(ResearchID)]
  dt[, State := trimws(as.character(State))]
  dt[, LookBack := suppressWarnings(as.integer(LookBack))]
  dt[, Recid := suppressWarnings(as.integer(Recid))]
  
  dt <- dt[!is.na(ResearchID) & !is.na(State) & !is.na(LookBack) & !is.na(Recid)]
  dt <- dt[State %chin% EXPECTED_STATES]
  dt <- dt[Recid %in% c(0L, 1L)]
  
  dt[, State := factor(State, levels = EXPECTED_STATES)]
  dt[, State := relevel(State, ref = STATE_REFERENCE)]
  
  mu_lb <- mean(dt$LookBack, na.rm = TRUE)
  sd_lb <- sd(dt$LookBack, na.rm = TRUE)
  if (!is.finite(sd_lb) || sd_lb == 0) {
    stop("LookBack SD is invalid in state subsample ", subset_id, ".")
  }
  
  dt[, LookBack_z := (LookBack - mu_lb) / sd_lb]
  
  dt_comp <- if (USE_COMPRESSION) {
    compress_binomial_counts(dt, c("ResearchID", "State", "LookBack", "LookBack_z"))
  } else {
    tmp <- copy(dt)
    tmp[, y := Recid]
    tmp[, n := 1L]
    tmp[, fail := 1L - Recid]
    tmp[, .(ResearchID, State, LookBack, LookBack_z, y, n, fail)]
  }
  
  dt_comp[, ResearchID := factor(ResearchID)]
  dt_comp[, State := factor(State, levels = EXPECTED_STATES)]
  dt_comp[, State := relevel(State, ref = STATE_REFERENCE)]
  
  f_add <- cbind(y, fail) ~ State + LookBack_z + (1 | ResearchID)
  f_int <- cbind(y, fail) ~ State * LookBack_z + (1 | ResearchID)
  
  mods <- fit_glmm_pair(dt_comp, f_add, f_int)
  
  coef_dt <- rbindlist(
    list(
      extract_fixed_effects(mods$additive, "additive"),
      extract_fixed_effects(mods$interaction, "interaction")
    ),
    use.names = TRUE,
    fill = TRUE
  )
  coef_dt[, subset := subset_id]
  
  fit_dt <- rbindlist(
    list(
      compute_fit_stats(mods$additive, "additive", dt),
      compute_fit_stats(mods$interaction, "interaction", dt)
    ),
    use.names = TRUE,
    fill = TRUE
  )
  fit_dt[, subset := subset_id]
  
  var_add <- extract_random_intercept_variance(mods$additive)
  var_int <- extract_random_intercept_variance(mods$interaction)
  
  re_dt <- rbindlist(
    list(
      data.table(
        subset = subset_id,
        model = "additive",
        var_u0 = var_add,
        sd_u0 = sqrt(var_add),
        ICC = calc_icc(var_add)
      ),
      data.table(
        subset = subset_id,
        model = "interaction",
        var_u0 = var_int,
        sd_u0 = sqrt(var_int),
        ICC = calc_icc(var_int)
      )
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
  lb_seq <- sort(unique(dt$LookBack))
  curve_grid <- CJ(LookBack = lb_seq, State = factor(EXPECTED_STATES, levels = EXPECTED_STATES))
  curve_grid[, State := relevel(State, ref = STATE_REFERENCE)]
  curve_grid[, LookBack_z := (LookBack - mu_lb) / sd_lb]
  curve_grid[, p_additive := as.numeric(predict(mods$additive, newdata = curve_grid, type = "response", re.form = NA))]
  curve_grid[, p_interaction := as.numeric(predict(mods$interaction, newdata = curve_grid, type = "response", re.form = NA))]
  curve_grid[, rr_additive := p_additive / p_additive[LookBack == 1L][1], by = State]
  curve_grid[, rr_interaction := p_interaction / p_interaction[LookBack == 1L][1], by = State]
  curve_grid[, subset := subset_id]
  
  list(
    coef = coef_dt,
    fit = fit_dt,
    re = re_dt,
    curves = curve_grid
  )
}

state_results <- lapply(
  seq_along(SUBSAMPLE_FILES),
  function(i) fit_one_state_subsample(SUBSAMPLE_FILES[i], SUBSAMPLE_IDS[i])
)

state_coef_all <- rbindlist(lapply(state_results, `[[`, "coef"), use.names = TRUE, fill = TRUE)
state_fit_all <- rbindlist(lapply(state_results, `[[`, "fit"), use.names = TRUE, fill = TRUE)
state_re_all <- rbindlist(lapply(state_results, `[[`, "re"), use.names = TRUE, fill = TRUE)
state_curves_all <- rbindlist(lapply(state_results, `[[`, "curves"), use.names = TRUE, fill = TRUE)

state_meta_coef <- meta_pool_coefficients(state_coef_all)

state_fit_summary <- state_fit_all[
  ,
  .(
    mean_minus2LL = mean(minus2LL, na.rm = TRUE),
    mean_AIC = mean(AIC, na.rm = TRUE),
    mean_BIC = mean(BIC, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    mean_Brier = mean(Brier, na.rm = TRUE)
  ),
  by = model
][order(model)]

state_re_summary <- state_re_all[
  ,
  .(
    mean_var_u0 = mean(var_u0, na.rm = TRUE),
    mean_sd_u0 = mean(sd_u0, na.rm = TRUE),
    mean_ICC = mean(ICC, na.rm = TRUE)
  ),
  by = model
][order(model)]

state_subset_weights <- state_fit_all[, .(n_obs = unique(n_obs)[1]), by = subset]

state_curves_long <- rbindlist(
  list(
    state_curves_all[, .(subset, LookBack, State = as.character(State), model = "additive", p_hat = p_additive, RR = rr_additive)],
    state_curves_all[, .(subset, LookBack, State = as.character(State), model = "interaction", p_hat = p_interaction, RR = rr_interaction)]
  ),
  use.names = TRUE,
  fill = TRUE
)

state_curves_long <- merge(
  state_curves_long,
  state_subset_weights,
  by = "subset",
  all.x = TRUE
)

state_curve_summary <- state_curves_long[
  ,
  .(
    p_hat = weighted.mean(p_hat, w = n_obs, na.rm = TRUE),
    RR = weighted.mean(RR, w = n_obs, na.rm = TRUE)
  ),
  by = .(LookBack, State, model)
][order(model, State, LookBack)]

state_gap_dt <- state_curve_summary[
  model == "interaction",
  .(
    p_max = max(p_hat, na.rm = TRUE),
    p_min = min(p_hat, na.rm = TRUE),
    gap = max(p_hat, na.rm = TRUE) - min(p_hat, na.rm = TRUE),
    state_max = State[which.max(p_hat)],
    state_min = State[which.min(p_hat)]
  ),
  by = LookBack
][order(LookBack)]

wa_dt <- fread(WA_FILE, showProgress = FALSE)
validate_columns(wa_dt, REQUIRED_WA_COLUMNS, WA_FILE)

wa_dt[, ResearchID := as.character(ResearchID)]
wa_dt[, Recid := suppressWarnings(as.integer(Recid))]
wa_dt[, LookBack := suppressWarnings(as.integer(LookBack))]
wa_dt[, PriorIncarYearsCount := suppressWarnings(as.numeric(PriorIncarYearsCount))]
wa_dt[, CharType := trimws(as.character(CharType))]
# Blank CharType is common in the real extract (confirmed via schema audit)
# and is coded as Non-Violent below, not dropped -- normalize "" to NA here
# so it's treated identically to a true missing value.
wa_dt[CharType == "", CharType := NA_character_]
wa_dt <- wa_dt[!is.na(ResearchID) & !is.na(Recid) & !is.na(LookBack)]
wa_dt <- wa_dt[Recid %in% c(0L, 1L)]
if ("State" %in% names(wa_dt)) wa_dt <- wa_dt[State == "WA"]

if (is.numeric(wa_dt$Male) || is.integer(wa_dt$Male)) {
  wa_dt[, Male := factor(Male, levels = c(0, 1), labels = c("Female", "Male"))]
} else {
  wa_dt[, Male := trimws(as.character(Male))]
  wa_dt[Male %in% c("0", "Female", "female"), Male := "Female"]
  wa_dt[Male %in% c("1", "Male", "male"), Male := "Male"]
  wa_dt[, Male := factor(Male, levels = c("Female", "Male"))]
}
wa_dt[, Male := relevel(Male, ref = "Female")]

if (is.numeric(wa_dt$RaceEthnicity2) || is.integer(wa_dt$RaceEthnicity2)) {
  wa_dt[, RaceEthnicity2 := factor(RaceEthnicity2, levels = c(0, 1, 2, 3), labels = c("White", "Black", "Hispanic", "Other"))]
} else {
  wa_dt[, RaceEthnicity2 := trimws(as.character(RaceEthnicity2))]
  wa_dt[RaceEthnicity2 %in% c("White", "white", "1"), RaceEthnicity2 := "White"]
  wa_dt[RaceEthnicity2 %in% c("Black", "black", "2"), RaceEthnicity2 := "Black"]
  wa_dt[RaceEthnicity2 %in% c("Hispanic", "hispanic", "3"), RaceEthnicity2 := "Hispanic"]
  wa_dt[RaceEthnicity2 %in% c("Other", "other", "4"), RaceEthnicity2 := "Other"]
  wa_dt[, RaceEthnicity2 := factor(RaceEthnicity2, levels = c("White", "Black", "Hispanic", "Other"))]
}
wa_dt[, RaceEthnicity2 := relevel(RaceEthnicity2, ref = "White")]

if (is.numeric(wa_dt$AgeCurveGrouped) || is.integer(wa_dt$AgeCurveGrouped)) {
  wa_dt[, AgeCurveGrouped := factor(AgeCurveGrouped, levels = c(1, 2, 3), labels = c("18-24", "25-34", "35+"))]
} else {
  wa_dt[, AgeCurveGrouped := trimws(as.character(AgeCurveGrouped))]
  wa_dt[AgeCurveGrouped %in% c("18-24", "18–24"), AgeCurveGrouped := "18-24"]
  wa_dt[AgeCurveGrouped %in% c("25-34", "25–34"), AgeCurveGrouped := "25-34"]
  wa_dt[AgeCurveGrouped %in% c("35+"), AgeCurveGrouped := "35+"]
  wa_dt[, AgeCurveGrouped := factor(AgeCurveGrouped, levels = c("18-24", "25-34", "35+"))]
}
wa_dt[, AgeCurveGrouped := relevel(AgeCurveGrouped, ref = "18-24")]

valid_types <- c("Drug", "Other", "Property", "Sex", "Violent")
wa_dt <- wa_dt[is.na(CharType) | CharType %in% valid_types]
# Blank/missing CharType is coded as Non-Violent (not dropped): NA no longer
# excludes a row from CrimeType2, it just falls into the Non-Violent bucket
# along with every other non-"Violent" value.
wa_dt[, CrimeType2 := fifelse(!is.na(CharType) & CharType == "Violent", "Violent", "Non-Violent")]
wa_dt[, CrimeType2 := factor(CrimeType2, levels = c("Non-Violent", "Violent"))]
wa_dt[is.na(PriorIncarYearsCount), PriorIncarYearsCount := 0]
wa_dt <- wa_dt[!is.na(Male) & !is.na(RaceEthnicity2) & !is.na(AgeCurveGrouped) & !is.na(CrimeType2)]

wa_mu_lb <- mean(wa_dt$LookBack, na.rm = TRUE)
wa_sd_lb <- sd(wa_dt$LookBack, na.rm = TRUE)
if (!is.finite(wa_sd_lb) || wa_sd_lb == 0) stop("LookBack SD is invalid in WA data.")
wa_dt[, LookBack_z := (LookBack - wa_mu_lb) / wa_sd_lb]

fit_wa_group_models <- function(dt, group_var, model_label) {
  dt_local <- copy(dt)
  by_vars <- c("ResearchID", "LookBack", "LookBack_z", group_var)
  
  dt_comp <- if (USE_COMPRESSION) {
    compress_binomial_counts(dt_local, by_vars = by_vars)
  } else {
    tmp <- copy(dt_local)
    tmp[, y := Recid]
    tmp[, n := 1L]
    tmp[, fail := 1L - Recid]
    tmp[, c("ResearchID", "LookBack", "LookBack_z", group_var, "y", "n", "fail"), with = FALSE]
  }
  
  if (is.factor(dt_local[[group_var]])) {
    dt_comp[, (group_var) := factor(get(group_var), levels = levels(dt_local[[group_var]]))]
  }
  dt_comp[, ResearchID := factor(ResearchID)]
  
  f_add <- as.formula(paste0("cbind(y, fail) ~ ", group_var, " + LookBack_z + (1 | ResearchID)"))
  f_int <- as.formula(paste0("cbind(y, fail) ~ ", group_var, " * LookBack_z + (1 | ResearchID)"))
  
  mods <- fit_glmm_pair(dt_comp, f_add, f_int)
  
  coef_dt <- rbindlist(
    list(
      extract_fixed_effects(mods$additive, "additive"),
      extract_fixed_effects(mods$interaction, "interaction")
    ),
    use.names = TRUE,
    fill = TRUE
  )
  coef_dt[, domain := model_label]
  
  fit_dt <- rbindlist(
    list(
      compute_fit_stats(mods$additive, "additive", dt_local),
      compute_fit_stats(mods$interaction, "interaction", dt_local)
    ),
    use.names = TRUE,
    fill = TRUE
  )
  fit_dt[, domain := model_label]
  
  var_add <- extract_random_intercept_variance(mods$additive)
  var_int <- extract_random_intercept_variance(mods$interaction)
  
  re_dt <- rbindlist(
    list(
      data.table(domain = model_label, model = "additive", var_u0 = var_add, sd_u0 = sqrt(var_add), ICC = calc_icc(var_add)),
      data.table(domain = model_label, model = "interaction", var_u0 = var_int, sd_u0 = sqrt(var_int), ICC = calc_icc(var_int))
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
  lrt_dt <- as.data.table(anova(mods$additive, mods$interaction, test = "Chisq"))
  lrt_dt[, domain := model_label]
  
  lb_seq <- sort(unique(dt_local$LookBack))
  group_levels <- if (is.factor(dt_local[[group_var]])) levels(dt_local[[group_var]]) else sort(unique(dt_local[[group_var]]))
  
  curve_grid <- CJ(LookBack = lb_seq, GroupLevel = group_levels)
  setnames(curve_grid, "GroupLevel", group_var)
  if (is.factor(dt_local[[group_var]])) {
    curve_grid[, (group_var) := factor(get(group_var), levels = group_levels)]
  }
  curve_grid[, LookBack_z := (LookBack - wa_mu_lb) / wa_sd_lb]
  curve_grid[, p_hat := as.numeric(predict(mods$interaction, newdata = curve_grid, type = "response", re.form = NA))]
  curve_grid[, RR := p_hat / p_hat[LookBack == 1L][1], by = group_var]
  curve_grid[, domain := model_label]
  
  list(
    coefficients = coef_dt,
    fit = fit_dt,
    re = re_dt,
    lrt = lrt_dt,
    curves = curve_grid
  )
}

fit_wa_incarceration_models <- function(dt) {
  dt_local <- copy(dt)
  by_vars <- c("ResearchID", "LookBack", "LookBack_z", "PriorIncarYearsCount")
  
  dt_comp <- if (USE_COMPRESSION) {
    compress_binomial_counts(dt_local, by_vars = by_vars)
  } else {
    tmp <- copy(dt_local)
    tmp[, y := Recid]
    tmp[, n := 1L]
    tmp[, fail := 1L - Recid]
    tmp[, .(ResearchID, LookBack, LookBack_z, PriorIncarYearsCount, y, n, fail)]
  }
  
  dt_comp[, ResearchID := factor(ResearchID)]
  
  f_add <- cbind(y, fail) ~ LookBack_z + PriorIncarYearsCount + (1 | ResearchID)
  f_int <- cbind(y, fail) ~ LookBack_z * PriorIncarYearsCount + (1 | ResearchID)
  
  mods <- fit_glmm_pair(dt_comp, f_add, f_int)
  
  coef_dt <- rbindlist(
    list(
      extract_fixed_effects(mods$additive, "additive"),
      extract_fixed_effects(mods$interaction, "interaction")
    ),
    use.names = TRUE,
    fill = TRUE
  )
  coef_dt[, domain := "Incarceration"]
  
  fit_dt <- rbindlist(
    list(
      compute_fit_stats(mods$additive, "additive", dt_local),
      compute_fit_stats(mods$interaction, "interaction", dt_local)
    ),
    use.names = TRUE,
    fill = TRUE
  )
  fit_dt[, domain := "Incarceration"]
  
  var_add <- extract_random_intercept_variance(mods$additive)
  var_int <- extract_random_intercept_variance(mods$interaction)
  
  re_dt <- rbindlist(
    list(
      data.table(domain = "Incarceration", model = "additive", var_u0 = var_add, sd_u0 = sqrt(var_add), ICC = calc_icc(var_add)),
      data.table(domain = "Incarceration", model = "interaction", var_u0 = var_int, sd_u0 = sqrt(var_int), ICC = calc_icc(var_int))
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
  lrt_dt <- as.data.table(anova(mods$additive, mods$interaction, test = "Chisq"))
  lrt_dt[, domain := "Incarceration"]
  
  lb_seq <- sort(unique(dt_local$LookBack))
  incar_grid <- data.table(LookBack = rep(lb_seq, 2L))
  incar_grid[, PriorIncarYearsCount := rep(c(0, 1), each = length(lb_seq))]
  incar_grid[, LookBack_z := (LookBack - wa_mu_lb) / wa_sd_lb]
  incar_grid[, GroupLabel := factor(ifelse(PriorIncarYearsCount == 0, "0 years", "1 year"), levels = c("0 years", "1 year"))]
  incar_grid[, p_hat := as.numeric(predict(mods$interaction, newdata = incar_grid, type = "response", re.form = NA))]
  incar_grid[, RR := p_hat / p_hat[LookBack == 1L][1], by = GroupLabel]
  incar_grid[, domain := "Incarceration"]
  
  list(
    coefficients = coef_dt,
    fit = fit_dt,
    re = re_dt,
    lrt = lrt_dt,
    curves = incar_grid
  )
}

sex_results <- fit_wa_group_models(wa_dt, "Male", "Sex")
race_results <- fit_wa_group_models(wa_dt, "RaceEthnicity2", "RaceEthnicity")
age_results <- fit_wa_group_models(wa_dt, "AgeCurveGrouped", "Age")
crime_results <- fit_wa_group_models(wa_dt, "CrimeType2", "CrimeType")
incar_results <- fit_wa_incarceration_models(wa_dt)

wa_coef_all <- rbindlist(
  list(
    sex_results$coefficients,
    race_results$coefficients,
    age_results$coefficients,
    crime_results$coefficients,
    incar_results$coefficients
  ),
  use.names = TRUE,
  fill = TRUE
)

wa_fit_all <- rbindlist(
  list(
    sex_results$fit,
    race_results$fit,
    age_results$fit,
    crime_results$fit,
    incar_results$fit
  ),
  use.names = TRUE,
  fill = TRUE
)

wa_re_all <- rbindlist(
  list(
    sex_results$re,
    race_results$re,
    age_results$re,
    crime_results$re,
    incar_results$re
  ),
  use.names = TRUE,
  fill = TRUE
)

wa_lrt_all <- rbindlist(
  list(
    sex_results$lrt,
    race_results$lrt,
    age_results$lrt,
    crime_results$lrt,
    incar_results$lrt
  ),
  use.names = TRUE,
  fill = TRUE
)

table3_dt <- dcast(
  state_meta_coef[, .(
    model,
    term,
    logit_se = sprintf("%.2f (%.2f)", pooled_est, pooled_se),
    or_ci = sprintf("%.2f (%.2f, %.2f)", OR, OR_low, OR_high),
    p_fmt = format_p(p)
  )],
  term ~ model,
  value.var = c("logit_se", "or_ci", "p_fmt")
)
table3_dt[, Predictor := make_term_labels(term)]
table3_dt[, term := NULL]

table4_dt <- wa_coef_all[
  domain %in% c("Sex", "RaceEthnicity", "Age"),
  .(
    domain,
    model,
    Predictor = make_term_labels(term),
    logit_se = sprintf("%.2f (%.2f)", estimate, std.error),
    or_ci = sprintf("%.2f (%.2f, %.2f)", OR, OR_low, OR_high),
    p_fmt = format_p(p.value)
  )
]

table5_dt <- wa_coef_all[
  domain %in% c("CrimeType", "Incarceration"),
  .(
    domain,
    model,
    Predictor = make_term_labels(term),
    logit_se = sprintf("%.2f (%.2f)", estimate, std.error),
    or_ci = sprintf("%.2f (%.2f, %.2f)", OR, OR_low, OR_high),
    p_fmt = format_p(p.value)
  )
]

appendixE_dt <- copy(wa_coef_all)[
  ,
  .(
    domain,
    model,
    Predictor = make_term_labels(term),
    estimate,
    std.error,
    OR,
    OR_low,
    OR_high,
    p_value = p.value
  )
]

state_fig_dt <- state_curve_summary[model == "interaction"]

fig5_prob <- ggplot(state_fig_dt, aes(x = LookBack, y = p_hat, linetype = State, color = State)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "A. State-specific predicted recidivism probability",
    x = "Lookback year",
    y = "Predicted probability"
  ) +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig5_rr <- ggplot(state_fig_dt, aes(x = LookBack, y = RR, linetype = State, color = State)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  labs(
    title = "B. State-specific risk ratio decay",
    x = "Lookback year",
    y = "Risk ratio"
  ) +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig5 <- fig5_prob + fig5_rr + plot_layout(ncol = 2)
ggsave(file.path(OUT_DIR, "H2_figure5_state_decay.png"), fig5, width = 12, height = 5, dpi = 300)

sex_curve_dt <- copy(sex_results$curves)
race_curve_dt <- copy(race_results$curves)
age_curve_dt <- copy(age_results$curves)

fig6a <- ggplot(sex_curve_dt, aes(x = LookBack, y = p_hat, linetype = Male, color = Male)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "A. Sex", x = "Lookback year", y = "Predicted probability") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig6b <- ggplot(race_curve_dt, aes(x = LookBack, y = p_hat, linetype = RaceEthnicity2, color = RaceEthnicity2)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "B. Race/ethnicity", x = "Lookback year", y = "Predicted probability") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig6c <- ggplot(age_curve_dt, aes(x = LookBack, y = p_hat, linetype = AgeCurveGrouped, color = AgeCurveGrouped)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "C. Age", x = "Lookback year", y = "Predicted probability") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig6 <- fig6a / fig6b / fig6c
ggsave(file.path(OUT_DIR, "H2_figure6_demographic_decay.png"), fig6, width = 8, height = 12, dpi = 300)

crime_curve_dt <- copy(crime_results$curves)
incar_curve_dt <- copy(incar_results$curves)

fig7a <- ggplot(crime_curve_dt, aes(x = LookBack, y = p_hat, linetype = CrimeType2, color = CrimeType2)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "A. Offense type", x = "Lookback year", y = "Predicted probability") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig7b <- ggplot(incar_curve_dt, aes(x = LookBack, y = p_hat, linetype = GroupLabel, color = GroupLabel)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "B. Incarceration exposure", x = "Lookback year", y = "Predicted probability") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig7 <- fig7a + fig7b + plot_layout(ncol = 2)
ggsave(file.path(OUT_DIR, "H2_figure7_context_decay.png"), fig7, width = 12, height = 5, dpi = 300)

############################################################
# SAVE FILES
############################################################

# Pooled H2 state interaction curves (used for H2 figures and summaries)
state_interaction_curves_pooled <- state_curve_summary[model == "interaction"][order(State, LookBack)]

# Subset-level H2 state interaction curves (used for H3)
state_interaction_curves_by_subset <- state_curves_long[model == "interaction"][
  ,
  .(subset, LookBack, State, p_hat, RR, model)
][order(subset, State, LookBack)]

fwrite(state_meta_coef, file.path(OUT_DIR, "H2_state_meta_coefficients.csv"))
fwrite(state_fit_summary, file.path(OUT_DIR, "H2_state_fit_summary.csv"))
fwrite(table3_dt, file.path(OUT_DIR, "H2_table3_state.csv"))
fwrite(table4_dt, file.path(OUT_DIR, "H2_table4_demographic_models.csv"))
fwrite(table5_dt, file.path(OUT_DIR, "H2_table5_context_models.csv"))
fwrite(state_gap_dt, file.path(OUT_DIR, "H2_appendixD_state_gap_by_lookback.csv"))
fwrite(appendixE_dt, file.path(OUT_DIR, "H2_appendixE_alternate_models.csv"))
fwrite(wa_lrt_all, file.path(OUT_DIR, "H2_wa_model_comparisons.csv"))
fwrite(wa_fit_all, file.path(OUT_DIR, "H2_wa_fit_stats.csv"))
fwrite(wa_re_all, file.path(OUT_DIR, "H2_wa_random_effects.csv"))
fwrite(state_interaction_curves_pooled, file.path(OUT_DIR, "H2_state_interaction_curves.csv"))
fwrite(state_interaction_curves_by_subset, file.path(OUT_DIR, "H2_state_interaction_curves_by_subset.csv"))

wb <- createWorkbook()
addWorksheet(wb, "Table3_State"); writeDataTable(wb, "Table3_State", table3_dt)
addWorksheet(wb, "Table4_Demographics"); writeDataTable(wb, "Table4_Demographics", table4_dt)
addWorksheet(wb, "Table5_Context"); writeDataTable(wb, "Table5_Context", table5_dt)
addWorksheet(wb, "AppendixD_Gap"); writeDataTable(wb, "AppendixD_Gap", state_gap_dt)
addWorksheet(wb, "AppendixE_Alternate"); writeDataTable(wb, "AppendixE_Alternate", appendixE_dt)
addWorksheet(wb, "StateMetaCoefficients"); writeDataTable(wb, "StateMetaCoefficients", state_meta_coef)
addWorksheet(wb, "WAFit"); writeDataTable(wb, "WAFit", wa_fit_all)
saveWorkbook(wb, file.path(OUT_DIR, "H2_results_workbook.xlsx"), overwrite = TRUE)

log_line("H2 analyses complete.")
print(state_fit_summary)
print(wa_fit_all[, .(domain, model, RMSE, Brier)])