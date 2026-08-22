############################################################
# SCRIPT 03 — H1 DECAY: DESCRIPTIVE PATTERN AND MODEL FORM
#
# PURPOSE
#   Estimate the functional form of recidivism decay using the
#   distributed subsample files.
#
# INPUT
#   /data/subsamples/subsample_00_raw.csv ... subsample_19_raw.csv
#
# REQUIRED COLUMNS
#   ResearchID, State, LookBack, Recid
#
# OUTPUT
#   /output/03_H1_decay/
#     - H1_yearly_probabilities_and_rr.csv
#     - H1_subset_level_coefficients.csv
#     - H1_subset_level_fit_stats.csv
#     - H1_subset_level_random_effects.csv
#     - H1_subset_level_model_comparisons.csv
#     - H1_meta_pooled_coefficients.csv
#     - H1_meta_fit_summary.csv
#     - H1_stouffer_model_comparisons.csv
#     - H1_table2_glm_vs_glmm.csv
#     - H1_appendixC_quadratic_table.csv
#     - H1_figure3_descriptive_decay.png
#     - H1_figure4_glm_vs_glmm_linear.png
#     - H1_appendixC_linear_vs_quadratic.png
#     - H1_results_workbook.xlsx
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
SUBSAMPLE_DIR <- file.path(DATA_DIR, "subsamples")
OUT_DIR <- file.path(REPO_DIR, "output", "03_H1_decay")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SUBSAMPLE_IDS <- sprintf("%02d", 0:19)
INPUT_FILES <- file.path(SUBSAMPLE_DIR, sprintf("subsample_%s_raw.csv", SUBSAMPLE_IDS))

REQUIRED_COLUMNS <- c("ResearchID", "State", "LookBack", "Recid")
EXPECTED_STATES <- c("AK", "AZ", "FL", "TX", "WA")

ID_VAR <- "ResearchID"
Y_VAR <- "Recid"
LB_VAR <- "LookBack"

BASE_FAMILY <- "Times New Roman"
RESID_VAR <- (pi^2) / 3

USE_COMPRESSION <- TRUE
MAX_IDS_PER_SUBSAMPLE <- Inf

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

compress_binomial_counts <- function(dt) {
  comp <- dt[
    ,
    .(
      y = sum(get(Y_VAR), na.rm = TRUE),
      n = .N
    ),
    by = .(ResearchID, LookBack, LookBack_z, LookBack_z2)
  ]
  comp[, fail := n - y]
  comp[]
}

compute_descriptive_yearly_stats <- function(dt) {
  yearly <- dt[
    ,
    .(
      person_year = .N,
      recid_events = sum(get(Y_VAR) == 1L, na.rm = TRUE)
    ),
    by = .(LookBack = get(LB_VAR))
  ][order(LookBack)]

  yearly[, recid_probability := recid_events / person_year]

  p1 <- yearly[LookBack == 1L, recid_probability]
  if (length(p1) == 0L || is.na(p1) || p1 == 0) {
    stop("LookBack == 1 is missing or has zero probability.")
  }

  yearly[, risk_ratio_vs_year1 := recid_probability / p1]
  yearly[, proportionate_reduction_vs_year1 := 1 - risk_ratio_vs_year1]
  yearly[]
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
        model = m, term = t,
        pooled_est = NA_real_, pooled_se = NA_real_,
        z = NA_real_, p = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_
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

extract_lrt_row <- function(anova_obj, subset_id, comparison_name, simpler_model, complex_model, n_obs) {
  tab <- as.data.table(anova_obj)
  if (nrow(tab) < 2L) {
    return(data.table(
      subset = subset_id,
      comparison = comparison_name,
      simpler_model = simpler_model,
      complex_model = complex_model,
      llr = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      direction = 0L,
      preferred_model = NA_character_,
      n_obs = n_obs
    ))
  }

  chisq_col <- if ("Chisq" %in% names(tab)) "Chisq" else if ("LRT" %in% names(tab)) "LRT" else NA_character_
  p_col <- grep("^Pr\\(", names(tab), value = TRUE)[1]

  llr_val <- if (!is.na(chisq_col)) suppressWarnings(as.numeric(tab[[chisq_col]][2])) else NA_real_
  p_val <- if (!is.na(p_col)) suppressWarnings(as.numeric(tab[[p_col]][2])) else NA_real_

  if ("Df" %in% names(tab)) {
    df_vals <- suppressWarnings(as.numeric(tab[["Df"]]))
    df_diff <- if (length(df_vals) >= 2L) df_vals[2] - df_vals[1] else NA_real_
  } else {
    df_diff <- NA_real_
  }

  if (!is.finite(df_diff) || is.na(df_diff)) df_diff <- 1

  direction <- if (is.finite(llr_val) && !is.na(llr_val) && llr_val > 0) 1L else 0L
  preferred_model <- if (direction == 1L) complex_model else simpler_model

  data.table(
    subset = subset_id,
    comparison = comparison_name,
    simpler_model = simpler_model,
    complex_model = complex_model,
    llr = llr_val,
    df_diff = df_diff,
    p_value = p_val,
    direction = direction,
    preferred_model = preferred_model,
    n_obs = n_obs
  )
}

stouffer_combine <- function(comp_dt) {
  if (nrow(comp_dt) == 0L) return(data.table())

  comp_dt <- copy(comp_dt)
  comp_dt <- comp_dt[!is.na(p_value) & !is.na(n_obs)]
  if (nrow(comp_dt) == 0L) return(data.table())

  comp_dt[, p_clipped := pmin(pmax(p_value, 1e-300), 1 - 1e-16)]
  comp_dt[, z_i := qnorm(1 - (p_clipped / 2)) * fifelse(direction > 0, 1, -1)]
  comp_dt[, weight := sqrt(n_obs)]

  comp_dt[
    ,
    .(
      k_subsamples = .N,
      stouffer_z = sum(weight * z_i, na.rm = TRUE) / sqrt(sum(weight^2, na.rm = TRUE)),
      stouffer_p = 2 * pnorm(-abs(sum(weight * z_i, na.rm = TRUE) / sqrt(sum(weight^2, na.rm = TRUE)))),
      n_favor_complex = sum(direction > 0, na.rm = TRUE),
      n_favor_simple = sum(direction <= 0, na.rm = TRUE),
      favored_model = fifelse(
        (sum(weight * z_i, na.rm = TRUE) / sqrt(sum(weight^2, na.rm = TRUE))) > 0,
        unique(complex_model)[1],
        unique(simpler_model)[1]
      )
    ),
    by = .(comparison, simpler_model, complex_model)
  ][order(comparison)]
}

fit_models_one_subsample <- function(file_path, subset_id) {
  log_line("Reading subsample %s: %s", subset_id, basename(file_path))

  dt <- fread(file_path, showProgress = FALSE)
  validate_columns(dt, file_path)

  dt[, ResearchID := as.character(ResearchID)]
  dt[, State := trimws(as.character(State))]
  dt[, LookBack := suppressWarnings(as.integer(LookBack))]
  dt[, Recid := suppressWarnings(as.integer(Recid))]

  dt <- dt[
    !is.na(ResearchID) &
      !is.na(State) &
      !is.na(LookBack) &
      !is.na(Recid)
  ]
  dt <- dt[State %chin% EXPECTED_STATES]
  dt <- dt[Recid %in% c(0L, 1L)]

  if (is.finite(MAX_IDS_PER_SUBSAMPLE)) {
    keep_ids <- head(unique(dt$ResearchID), MAX_IDS_PER_SUBSAMPLE)
    dt <- dt[ResearchID %in% keep_ids]
  }

  mu_lb <- mean(dt$LookBack, na.rm = TRUE)
  sd_lb <- sd(dt$LookBack, na.rm = TRUE)
  if (!is.finite(sd_lb) || sd_lb == 0) stop("LookBack SD is invalid in subsample ", subset_id, ".")

  dt[, LookBack_z := (LookBack - mu_lb) / sd_lb]
  dt[, LookBack_z2 := LookBack_z^2]

  yearly <- compute_descriptive_yearly_stats(dt)
  yearly[, subset := subset_id]

  dt_comp <- if (USE_COMPRESSION) compress_binomial_counts(dt) else {
    tmp <- copy(dt)
    tmp[, y := Recid]
    tmp[, fail := 1L - Recid]
    tmp[, .(ResearchID, LookBack, LookBack_z, LookBack_z2, y, fail)]
  }

  dt_comp[, ResearchID := factor(ResearchID)]

  f_glm_lin <- cbind(y, fail) ~ LookBack_z
  f_glm_quad <- cbind(y, fail) ~ LookBack_z + LookBack_z2
  f_glmm_lin <- cbind(y, fail) ~ LookBack_z + (1 | ResearchID)
  f_glmm_quad <- cbind(y, fail) ~ LookBack_z + LookBack_z2 + (1 | ResearchID)

  m_glm_lin <- glm(f_glm_lin, data = dt_comp, family = binomial(link = "logit"))
  m_glm_quad <- glm(f_glm_quad, data = dt_comp, family = binomial(link = "logit"))

  m_glmm_lin <- glmer(
    f_glmm_lin, data = dt_comp, family = binomial(link = "logit"), nAGQ = 0L,
    control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE, optCtrl = list(maxfun = 2e5))
  )

  m_glmm_quad <- glmer(
    f_glmm_quad, data = dt_comp, family = binomial(link = "logit"), nAGQ = 0L,
    control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE, optCtrl = list(maxfun = 2e5))
  )

  extract_fixed <- function(model, model_name) {
    tt <- as.data.table(broom.mixed::tidy(model, effects = "fixed"))
    tt[, `:=`(
      subset = subset_id,
      model = model_name,
      OR = exp(estimate),
      OR_low = exp(estimate - 1.96 * std.error),
      OR_high = exp(estimate + 1.96 * std.error)
    )]
    tt[]
  }

  fixed_dt <- rbindlist(
    list(
      extract_fixed(m_glm_lin, "glm_linear"),
      extract_fixed(m_glm_quad, "glm_quadratic"),
      extract_fixed(m_glmm_lin, "glmm_linear"),
      extract_fixed(m_glmm_quad, "glmm_quadratic")
    ),
    use.names = TRUE, fill = TRUE
  )

  compute_fit_stats <- function(model, model_name, dt_uncompressed, terms_n) {
    p_hat <- predict(model, newdata = dt_uncompressed, type = "response", re.form = NA)
    y_obs <- dt_uncompressed[[Y_VAR]]
    brier <- mean((y_obs - p_hat)^2, na.rm = TRUE)
    rmse <- sqrt(brier)
    ll <- as.numeric(logLik(model))
    data.table(
      subset = subset_id,
      model = model_name,
      minus2LL = -2 * ll,
      AIC = AIC(model),
      BIC = BIC(model),
      logLik = ll,
      RMSE = rmse,
      Brier = brier,
      n_terms = terms_n,
      n_obs = nrow(dt_uncompressed),
      n_ids = uniqueN(dt_uncompressed[[ID_VAR]])
    )
  }

  fit_dt <- rbindlist(
    list(
      compute_fit_stats(m_glm_lin, "glm_linear", dt, 2L),
      compute_fit_stats(m_glm_quad, "glm_quadratic", dt, 3L),
      compute_fit_stats(m_glmm_lin, "glmm_linear", dt, 2L),
      compute_fit_stats(m_glmm_quad, "glmm_quadratic", dt, 3L)
    ),
    use.names = TRUE, fill = TRUE
  )

  re_dt <- rbindlist(
    list(
      data.table(
        subset = subset_id, model = "glmm_linear",
        var_u0 = extract_random_intercept_variance(m_glmm_lin),
        sd_u0 = sqrt(extract_random_intercept_variance(m_glmm_lin)),
        ICC = calc_icc(extract_random_intercept_variance(m_glmm_lin))
      ),
      data.table(
        subset = subset_id, model = "glmm_quadratic",
        var_u0 = extract_random_intercept_variance(m_glmm_quad),
        sd_u0 = sqrt(extract_random_intercept_variance(m_glmm_quad)),
        ICC = calc_icc(extract_random_intercept_variance(m_glmm_quad))
      )
    ),
    use.names = TRUE, fill = TRUE
  )

  ll_glm_lin <- as.numeric(logLik(m_glm_lin))
  ll_glm_quad <- as.numeric(logLik(m_glm_quad))
  ll_glmm_lin <- as.numeric(logLik(m_glmm_lin))
  ll_glmm_quad <- as.numeric(logLik(m_glmm_quad))

  comp_dt <- rbindlist(
    list(
      data.table(
        subset = subset_id,
        comparison = "glm_linear_vs_glmm_linear",
        simpler_model = "glm_linear",
        complex_model = "glmm_linear",
        llr = 2 * (ll_glmm_lin - ll_glm_lin),
        df_diff = 1,
        p_value = ifelse(
          is.finite(2 * (ll_glmm_lin - ll_glm_lin)) && (2 * (ll_glmm_lin - ll_glm_lin)) > 0,
          pchisq(2 * (ll_glmm_lin - ll_glm_lin), df = 1, lower.tail = FALSE),
          1
        ),
        direction = ifelse(is.finite(ll_glmm_lin - ll_glm_lin) && (ll_glmm_lin - ll_glm_lin) > 0, 1L, 0L),
        preferred_model = ifelse(is.finite(ll_glmm_lin - ll_glm_lin) && (ll_glmm_lin - ll_glm_lin) > 0, "glmm_linear", "glm_linear"),
        n_obs = nrow(dt)
      ),
      data.table(
        subset = subset_id,
        comparison = "glm_quadratic_vs_glmm_quadratic",
        simpler_model = "glm_quadratic",
        complex_model = "glmm_quadratic",
        llr = 2 * (ll_glmm_quad - ll_glm_quad),
        df_diff = 1,
        p_value = ifelse(
          is.finite(2 * (ll_glmm_quad - ll_glm_quad)) && (2 * (ll_glmm_quad - ll_glm_quad)) > 0,
          pchisq(2 * (ll_glmm_quad - ll_glm_quad), df = 1, lower.tail = FALSE),
          1
        ),
        direction = ifelse(is.finite(ll_glmm_quad - ll_glm_quad) && (ll_glmm_quad - ll_glm_quad) > 0, 1L, 0L),
        preferred_model = ifelse(is.finite(ll_glmm_quad - ll_glm_quad) && (ll_glmm_quad - ll_glm_quad) > 0, "glmm_quadratic", "glm_quadratic"),
        n_obs = nrow(dt)
      ),
      extract_lrt_row(anova(m_glm_lin, m_glm_quad, test = "Chisq"), subset_id, "glm_linear_vs_glm_quadratic", "glm_linear", "glm_quadratic", nrow(dt)),
      extract_lrt_row(anova(m_glmm_lin, m_glmm_quad, test = "Chisq"), subset_id, "glmm_linear_vs_glmm_quadratic", "glmm_linear", "glmm_quadratic", nrow(dt))
    ),
    use.names = TRUE, fill = TRUE
  )

  lb_seq <- sort(unique(dt$LookBack))
  curve_grid <- data.table(LookBack = lb_seq)
  curve_grid[, LookBack_z := (LookBack - mu_lb) / sd_lb]
  curve_grid[, LookBack_z2 := LookBack_z^2]
  curve_grid[, p_glm_linear := as.numeric(predict(m_glm_lin, newdata = curve_grid, type = "response"))]
  curve_grid[, p_glm_quadratic := as.numeric(predict(m_glm_quad, newdata = curve_grid, type = "response"))]
  curve_grid[, p_glmm_linear := as.numeric(predict(m_glmm_lin, newdata = curve_grid, type = "response", re.form = NA))]
  curve_grid[, p_glmm_quadratic := as.numeric(predict(m_glmm_quad, newdata = curve_grid, type = "response", re.form = NA))]
  curve_grid[, rr_glm_linear := p_glm_linear / p_glm_linear[LookBack == 1L][1]]
  curve_grid[, rr_glm_quadratic := p_glm_quadratic / p_glm_quadratic[LookBack == 1L][1]]
  curve_grid[, rr_glmm_linear := p_glmm_linear / p_glmm_linear[LookBack == 1L][1]]
  curve_grid[, rr_glmm_quadratic := p_glmm_quadratic / p_glmm_quadratic[LookBack == 1L][1]]
  curve_grid[, subset := subset_id]

  list(yearly = yearly, fixed = fixed_dt, fit = fit_dt, re = re_dt, comp = comp_dt, curves = curve_grid)
}

log_line("Checking H1 input files.")
stop_if_missing_files(INPUT_FILES)

results_list <- lapply(seq_along(INPUT_FILES), function(i) {
  fit_models_one_subsample(INPUT_FILES[i], SUBSAMPLE_IDS[i])
})

yearly_all <- rbindlist(lapply(results_list, `[[`, "yearly"), use.names = TRUE, fill = TRUE)
fixed_all <- rbindlist(lapply(results_list, `[[`, "fixed"), use.names = TRUE, fill = TRUE)
fit_all <- rbindlist(lapply(results_list, `[[`, "fit"), use.names = TRUE, fill = TRUE)
re_all <- rbindlist(lapply(results_list, `[[`, "re"), use.names = TRUE, fill = TRUE)
comp_all <- rbindlist(lapply(results_list, `[[`, "comp"), use.names = TRUE, fill = TRUE)
curves_all <- rbindlist(lapply(results_list, `[[`, "curves"), use.names = TRUE, fill = TRUE)

yearly_summary <- yearly_all[
  ,
  .(
    person_year = sum(person_year, na.rm = TRUE),
    recid_events = sum(recid_events, na.rm = TRUE)
  ),
  by = LookBack
][order(LookBack)]
yearly_summary[, recid_probability := recid_events / person_year]
p1 <- yearly_summary[LookBack == 1L, recid_probability]
yearly_summary[, risk_ratio_vs_year1 := recid_probability / p1]
yearly_summary[, proportionate_reduction_vs_year1 := 1 - risk_ratio_vs_year1]

meta_coef <- meta_pool_coefficients(fixed_all)

fit_summary <- fit_all[
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

re_summary <- re_all[
  ,
  .(
    mean_var_u0 = mean(var_u0, na.rm = TRUE),
    mean_sd_u0 = mean(sd_u0, na.rm = TRUE),
    mean_ICC = mean(ICC, na.rm = TRUE)
  ),
  by = model
][order(model)]

stouffer_summary <- stouffer_combine(comp_all)
stouffer_summary[, stouffer_p_fmt := format_p(stouffer_p)]

subset_weights <- fit_all[, .(n_obs = unique(n_obs)[1]), by = subset]

curves_long_main <- rbindlist(
  list(
    curves_all[, .(subset, LookBack, model = "GLM linear", p_hat = p_glm_linear, RR = rr_glm_linear)],
    curves_all[, .(subset, LookBack, model = "GLMM linear", p_hat = p_glmm_linear, RR = rr_glmm_linear)]
  ),
  use.names = TRUE, fill = TRUE
)
curves_long_main <- merge(curves_long_main, subset_weights, by = "subset", all.x = TRUE)
curves_main_summary <- curves_long_main[
  ,
  .(
    p_hat = weighted.mean(p_hat, w = n_obs, na.rm = TRUE),
    RR = weighted.mean(RR, w = n_obs, na.rm = TRUE)
  ),
  by = .(LookBack, model)
][order(model, LookBack)]

curves_long_quad <- rbindlist(
  list(
    curves_all[, .(subset, LookBack, model = "Linear", p_hat = p_glmm_linear)],
    curves_all[, .(subset, LookBack, model = "Quadratic", p_hat = p_glmm_quadratic)]
  ),
  use.names = TRUE, fill = TRUE
)
curves_long_quad <- merge(curves_long_quad, subset_weights, by = "subset", all.x = TRUE)
curves_quad_summary <- curves_long_quad[
  ,
  .(p_hat = weighted.mean(p_hat, w = n_obs, na.rm = TRUE)),
  by = .(LookBack, model)
][order(model, LookBack)]

curves_quad_wide <- dcast(curves_quad_summary, LookBack ~ model, value.var = "p_hat")
curves_quad_wide[, prob_diff := Quadratic - Linear]
curves_quad_wide[, logodds_linear := log(Linear / (1 - Linear))]
curves_quad_wide[, logodds_quadratic := log(Quadratic / (1 - Quadratic))]

table2_dt <- merge(
  meta_coef[model == "glm_linear", .(
    term,
    glm_logit_se = sprintf("%.2f (%.2f)", pooled_est, pooled_se),
    glm_or_ci = sprintf("%.2f (%.2f, %.2f)", OR, OR_low, OR_high),
    glm_p = format_p(p)
  )],
  meta_coef[model == "glmm_linear", .(
    term,
    glmm_logit_se = sprintf("%.2f (%.2f)", pooled_est, pooled_se),
    glmm_or_ci = sprintf("%.2f (%.2f, %.2f)", OR, OR_low, OR_high),
    glmm_p = format_p(p)
  )],
  by = "term", all = TRUE
)
table2_dt[, Predictor := fifelse(term == "(Intercept)", "Intercept", fifelse(term == "LookBack_z", "LookBack (z)", term))]
table2_dt[, term := NULL]

appendixC_quad_dt <- merge(
  meta_coef[model == "glm_quadratic", .(
    term,
    glm_logit_se = sprintf("%.2f (%.2f)", pooled_est, pooled_se),
    glm_or_ci = sprintf("%.2f (%.2f, %.2f)", OR, OR_low, OR_high),
    glm_p = format_p(p)
  )],
  meta_coef[model == "glmm_quadratic", .(
    term,
    glmm_logit_se = sprintf("%.2f (%.2f)", pooled_est, pooled_se),
    glmm_or_ci = sprintf("%.2f (%.2f, %.2f)", OR, OR_low, OR_high),
    glmm_p = format_p(p)
  )],
  by = "term", all = TRUE
)
appendixC_quad_dt[, Predictor := fifelse(
  term == "(Intercept)", "Intercept",
  fifelse(term == "LookBack_z", "LookBack (z)", fifelse(term == "LookBack_z2", "LookBack (z)^2", term))
)]
appendixC_quad_dt[, term := NULL]

p_fig3a <- ggplot(yearly_summary, aes(x = LookBack, y = recid_probability)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "A. Observed recidivism probability by lookback year", x = "Lookback year", y = "Recidivism probability") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

p_fig3b <- ggplot(yearly_summary, aes(x = LookBack, y = risk_ratio_vs_year1)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  labs(title = "B. Risk ratio relative to lookback year 1", x = "Lookback year", y = "Risk ratio") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig3 <- p_fig3a + p_fig3b + plot_layout(ncol = 2)
ggsave(file.path(OUT_DIR, "H1_figure3_descriptive_decay.png"), fig3, width = 12, height = 5, dpi = 300)

p_fig4a <- ggplot(curves_main_summary, aes(x = LookBack, y = p_hat, linetype = model)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "A. Predicted recidivism probability", x = "Lookback year", y = "Predicted probability", linetype = "Model") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

p_fig4b <- ggplot(curves_main_summary, aes(x = LookBack, y = RR, linetype = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  labs(title = "B. Risk ratio decay", x = "Lookback year", y = "Risk ratio", linetype = "Model") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig4 <- p_fig4a + p_fig4b + plot_layout(ncol = 2)
ggsave(file.path(OUT_DIR, "H1_figure4_glm_vs_glmm_linear.png"), fig4, width = 12, height = 5, dpi = 300)

logodds_dt <- rbindlist(
  list(
    curves_quad_wide[, .(LookBack, model = "Linear", logodds = logodds_linear)],
    curves_quad_wide[, .(LookBack, model = "Quadratic", logodds = logodds_quadratic)]
  )
)
prob_dt <- rbindlist(
  list(
    curves_quad_wide[, .(LookBack, model = "Linear", p_hat = Linear)],
    curves_quad_wide[, .(LookBack, model = "Quadratic", p_hat = Quadratic)]
  )
)

p_appC_a <- ggplot(logodds_dt, aes(x = LookBack, y = logodds, linetype = model)) +
  geom_line(linewidth = 0.9) +
  labs(title = "A. Log-odds", x = "Lookback year", y = "Predicted log-odds", linetype = "Model") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

p_appC_b <- ggplot(prob_dt, aes(x = LookBack, y = p_hat, linetype = model)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "B. Probability", x = "Lookback year", y = "Predicted probability", linetype = "Model") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

p_appC_c <- ggplot(curves_quad_wide, aes(x = LookBack, y = prob_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.01)) +
  labs(title = "C. Quadratic minus linear probability", x = "Lookback year", y = "Difference") +
  theme_classic(base_family = BASE_FAMILY, base_size = 12)

fig_appC <- p_appC_a + p_appC_b + p_appC_c + plot_layout(ncol = 1)
ggsave(file.path(OUT_DIR, "H1_appendixC_linear_vs_quadratic.png"), fig_appC, width = 8, height = 12, dpi = 300)

fwrite(yearly_summary, file.path(OUT_DIR, "H1_yearly_probabilities_and_rr.csv"))
fwrite(fixed_all, file.path(OUT_DIR, "H1_subset_level_coefficients.csv"))
fwrite(fit_all, file.path(OUT_DIR, "H1_subset_level_fit_stats.csv"))
fwrite(re_all, file.path(OUT_DIR, "H1_subset_level_random_effects.csv"))
fwrite(comp_all, file.path(OUT_DIR, "H1_subset_level_model_comparisons.csv"))
fwrite(meta_coef, file.path(OUT_DIR, "H1_meta_pooled_coefficients.csv"))
fwrite(fit_summary, file.path(OUT_DIR, "H1_meta_fit_summary.csv"))
fwrite(stouffer_summary, file.path(OUT_DIR, "H1_stouffer_model_comparisons.csv"))
fwrite(table2_dt, file.path(OUT_DIR, "H1_table2_glm_vs_glmm.csv"))
fwrite(appendixC_quad_dt, file.path(OUT_DIR, "H1_appendixC_quadratic_table.csv"))

wb <- createWorkbook()
addWorksheet(wb, "YearlyDecay"); writeDataTable(wb, "YearlyDecay", yearly_summary)
addWorksheet(wb, "SubsetCoefficients"); writeDataTable(wb, "SubsetCoefficients", fixed_all)
addWorksheet(wb, "SubsetFit"); writeDataTable(wb, "SubsetFit", fit_all)
addWorksheet(wb, "SubsetRE"); writeDataTable(wb, "SubsetRE", re_all)
addWorksheet(wb, "SubsetModelComparisons"); writeDataTable(wb, "SubsetModelComparisons", comp_all)
addWorksheet(wb, "MetaCoefficients"); writeDataTable(wb, "MetaCoefficients", meta_coef)
addWorksheet(wb, "StoufferComparisons"); writeDataTable(wb, "StoufferComparisons", stouffer_summary)
addWorksheet(wb, "Table2"); writeDataTable(wb, "Table2", table2_dt)
addWorksheet(wb, "AppendixC_Quadratic"); writeDataTable(wb, "AppendixC_Quadratic", appendixC_quad_dt)
saveWorkbook(wb, file.path(OUT_DIR, "H1_results_workbook.xlsx"), overwrite = TRUE)

log_line("H1 analyses complete.")
print(yearly_summary)
print(fit_summary)
print(stouffer_summary)
