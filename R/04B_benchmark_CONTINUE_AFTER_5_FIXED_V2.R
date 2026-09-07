# ============================================================================
# 04B_benchmark_CONTINUE_AFTER_5_FIXED.R
# Continue benchmark from the CURRENT R session after the [5/10] error.
# IMPORTANT: DO NOT clear the workspace and DO NOT rerun 04_clinical_score...
#
# Fixes:
# 1) Robust flexible calibration for duplicated/discrete predictions by first
#    aggregating identical predicted probabilities and fitting weighted LOESS.
# 2) Recomputes ONLY score-model B=1000 calibration-error optimism/curves.
#    It does NOT rerun MICE, LASSO, CV, or the B=1000 scalar benchmark metrics.
# 3) Reuses the already-created B=5000 patient bootstrap indices and recomputes
#    ONLY ICI/E50/E90; AUC/Brier/intercept/slope are reused from memory.
# 4) Paired bootstrap distributions preserve the TRUE bootstrap replicate ID.
# 5) Hard-stop checks prevent silent all-NA output.
# ============================================================================

options(stringsAsFactors = FALSE)

.required <- c(
  "dat", "y", "pred_list", "score_specs", "MODEL_LABELS", "model_names",
  "metric_names", "err_names", "score_corrected", "score_boot_scalar",
  "bootstrap_indices", "B_OPTIMISM", "B_UNCERTAINTY", "uncertainty_store",
  "uncertainty_indices", "full_corrected", "cal_metrics_old", "OUT_PREFIX",
  "OUTCOME", "DCA_THRESHOLDS", "imp", "coef_out", "recon_check",
  "candidate_vars", "FINAL_VARS", "SCORE_VARS", "M", "MICE_MAXIT",
  "SEED_APPARENT", "SEED_BOOT_BASE", "SEED_UNCERTAINTY", "EXPECTED_N",
  "EXPECTED_EVENTS", "full_app_metrics"
)
.missing <- .required[!vapply(.required, exists, logical(1), inherits = TRUE)]
if (length(.missing) > 0L) {
  stop(
    "Required objects are missing from the current R session: ",
    paste(.missing, collapse = ", "),
    "\nDo NOT clear the workspace. This continuation script must be run in the session that reached the [5/10] error."
  )
}

if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
if (!requireNamespace("mice", quietly = TRUE)) stop("Package 'mice' is required.")

# ----------------------------- robust helpers --------------------------------
clip_prob2 <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

# Robust flexible calibration using a binomial GLM with a natural cubic spline
# of the logit of predicted risk. This avoids LOESS recursion/stack failures with
# duplicated bootstrap rows and discrete clinical scores.
flex_cal_predict <- function(y, p, eval_p) {
  y <- as.numeric(y)
  p <- clip_prob2(p)
  eval_p <- clip_prob2(eval_p)

  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]

  if (length(y) < 20L || length(unique(y)) < 2L) {
    return(rep(NA_real_, length(eval_p)))
  }

  x <- stats::qlogis(p)
  xnew <- stats::qlogis(eval_p)
  ux <- length(unique(x))

  # Use up to 3 df, but never more flexibility than the support allows.
  # For very discrete scores, fall back to ordinary logistic recalibration.
  df_use <- min(3L, max(1L, ux - 1L))

  if (df_use >= 2L) {
    fit <- tryCatch(
      suppressWarnings(stats::glm(
        y ~ splines::ns(x, df = df_use),
        family = stats::binomial()
      )),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      pr <- tryCatch(
        suppressWarnings(as.numeric(stats::predict(
          fit,
          newdata = data.frame(x = xnew),
          type = "response"
        ))),
        error = function(e) rep(NA_real_, length(xnew))
      )
      if (sum(is.finite(pr)) >= max(5L, floor(0.80 * length(xnew)))) {
        return(pmin(pmax(pr, 0), 1))
      }
    }
  }

  # Stable fallback: logistic calibration with intercept + slope.
  fit2 <- tryCatch(
    suppressWarnings(stats::glm(y ~ x, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit2)) return(rep(NA_real_, length(xnew)))

  pr2 <- tryCatch(
    suppressWarnings(as.numeric(stats::predict(
      fit2,
      newdata = data.frame(x = xnew),
      type = "response"
    ))),
    error = function(e) rep(NA_real_, length(xnew))
  )
  pmin(pmax(pr2, 0), 1)
}

get_calibration_errors_robust <- function(y, p) {
  p0 <- clip_prob2(p)
  obs_hat <- flex_cal_predict(y = y, p = p0, eval_p = p0)
  err <- abs(obs_hat - p0)
  err <- err[is.finite(err)]
  if (length(err) < max(20L, floor(0.80 * length(p0)))) {
    return(c(ICI = NA_real_, E50 = NA_real_, E90 = NA_real_))
  }
  c(
    ICI = mean(err),
    E50 = unname(stats::quantile(err, 0.50, na.rm = TRUE, names = FALSE)),
    E90 = unname(stats::quantile(err, 0.90, na.rm = TRUE, names = FALSE))
  )
}

safe_curve_robust <- function(y, p, grid) {
  flex_cal_predict(y = y, p = p, eval_p = grid)
}

safe_quantile2 <- function(x, prob, label) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) stop("No finite bootstrap values for: ", label)
  unname(stats::quantile(x, prob, na.rm = TRUE, names = FALSE))
}

# ---------------- REPAIR B=1000 SCORE CALIBRATION ERROR + CURVES -------------
message("[PATCH A] Recomputing ONLY score-model B=1000 ICI/E50/E90 optimism and flexible curves...")

score_boot_errors <- list()
score_boot_curves <- list()
score_cal_grid <- list()

for (nm in names(score_specs)) {
  score_boot_errors[[nm]] <- matrix(
    NA_real_, nrow = B_OPTIMISM, ncol = 2L * length(err_names),
    dimnames = list(NULL, c(paste0(err_names, "_boot"), paste0(err_names, "_orig")))
  )

  p_app <- pred_list[[nm]]
  lo <- as.numeric(stats::quantile(p_app, 0.025, na.rm = TRUE, names = FALSE))
  hi <- as.numeric(stats::quantile(p_app, 0.975, na.rm = TRUE, names = FALSE))
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    lo <- min(p_app, na.rm = TRUE)
    hi <- max(p_app, na.rm = TRUE)
  }
  score_cal_grid[[nm]] <- seq(lo, hi, length.out = 101L)
  score_boot_curves[[nm]] <- list(
    grid = score_cal_grid[[nm]],
    boot = matrix(NA_real_, nrow = B_OPTIMISM, ncol = 101L),
    orig = matrix(NA_real_, nrow = B_OPTIMISM, ncol = 101L)
  )
}

for (b in seq_len(B_OPTIMISM)) {
  idx <- bootstrap_indices[[b]]
  if (is.null(idx) || length(idx) != nrow(dat)) {
    stop("Invalid saved B=1000 bootstrap index at replicate ", b)
  }
  d_boot <- dat[idx, , drop = FALSE]

  for (nm in names(score_specs)) {
    v <- score_specs[[nm]]
    fit_b <- tryCatch(
      suppressWarnings(stats::glm(
        stats::reformulate(v, response = OUTCOME),
        data = d_boot,
        family = stats::binomial()
      )),
      error = function(e) NULL
    )
    if (is.null(fit_b)) next

    p_boot <- tryCatch(
      as.numeric(stats::predict(fit_b, newdata = d_boot, type = "response")),
      error = function(e) rep(NA_real_, nrow(d_boot))
    )
    p_orig <- tryCatch(
      as.numeric(stats::predict(fit_b, newdata = dat, type = "response")),
      error = function(e) rep(NA_real_, nrow(dat))
    )
    if (any(!is.finite(p_boot)) || any(!is.finite(p_orig))) next

    score_boot_errors[[nm]][b, paste0(err_names, "_boot")] <-
      get_calibration_errors_robust(d_boot[[OUTCOME]], p_boot)
    score_boot_errors[[nm]][b, paste0(err_names, "_orig")] <-
      get_calibration_errors_robust(y, p_orig)

    g <- score_cal_grid[[nm]]
    score_boot_curves[[nm]]$boot[b, ] <- safe_curve_robust(d_boot[[OUTCOME]], p_boot, g)
    score_boot_curves[[nm]]$orig[b, ] <- safe_curve_robust(y, p_orig, g)
  }

  if (b %% 100L == 0L) message("  calibration repair completed ", b, "/", B_OPTIMISM)
}

score_err_corrected <- list()
score_curve_corrected <- list()
score_optimism_rows <- list()

for (nm in names(score_specs)) {
  app <- get_metrics(y, pred_list[[nm]])
  app_err <- get_calibration_errors_robust(y, pred_list[[nm]])

  mat <- score_boot_scalar[[nm]]
  opt_scalar <- sapply(metric_names, function(mt) {
    mean(mat[, paste0(mt, "_boot")] - mat[, paste0(mt, "_orig")], na.rm = TRUE)
  })

  emat <- score_boot_errors[[nm]]
  eopt <- sapply(err_names, function(mt) {
    z <- emat[, paste0(mt, "_boot")] - emat[, paste0(mt, "_orig")]
    z <- z[is.finite(z)]
    if (length(z) < floor(0.80 * B_OPTIMISM)) {
      stop("Too few valid B=1000 calibration-error replicates for ", nm, " / ", mt,
           ": ", length(z), "/", B_OPTIMISM)
    }
    mean(z)
  })
  score_err_corrected[[nm]] <- app_err - eopt

  dif_curve <- score_boot_curves[[nm]]$boot - score_boot_curves[[nm]]$orig
  valid_per_grid <- colSums(is.finite(dif_curve))
  if (any(valid_per_grid < floor(0.70 * B_OPTIMISM))) {
    warning("Some curve grid points have <70% valid B=1000 replicates for ", nm,
            "; affected points will remain NA in the benchmark curve.")
  }
  copt <- apply(dif_curve, 2L, function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0L) return(NA_real_)
    mean(z)
  })
  app_curve <- safe_curve_robust(y, pred_list[[nm]], score_cal_grid[[nm]])
  score_curve_corrected[[nm]] <- pmin(pmax(app_curve - copt, 0), 1)

  score_optimism_rows[[nm]] <- data.frame(
    model = MODEL_LABELS[[nm]],
    metric = c(metric_names, err_names),
    apparent = c(app, app_err),
    mean_optimism_B1000 = c(opt_scalar, eopt),
    optimism_corrected = c(score_corrected[[nm]], score_err_corrected[[nm]]),
    stringsAsFactors = FALSE
  )
}

score_optimism_table <- do.call(rbind, score_optimism_rows)
write.csv(score_optimism_table,
          paste0(OUT_PREFIX, "05_score_model_B1000_optimism_summary.csv"),
          row.names = FALSE)

# ------------------------ FULL MODEL ICI/E50/E90 ------------------------------
# Use the already-finalized calibration post-processing CSV; do NOT recompute
# the full-model complete-pipeline optimism here.
required_cal_cols <- c(
  "Metric", "Apparent", "Mean_optimism_complete_pipeline_B1000", "Optimism_corrected"
)
if (!all(required_cal_cols %in% names(cal_metrics_old))) {
  stop("390CAL_POSTV3_02 file is missing required columns: ",
       paste(setdiff(required_cal_cols, names(cal_metrics_old)), collapse = ", "))
}

full_err_app <- setNames(rep(NA_real_, 3L), err_names)
full_err_opt <- setNames(rep(NA_real_, 3L), err_names)
full_err_corr <- setNames(rep(NA_real_, 3L), err_names)
for (mt in err_names) {
  rr <- cal_metrics_old[toupper(cal_metrics_old$Metric) == toupper(mt), , drop = FALSE]
  if (nrow(rr) != 1L) stop("Could not identify exactly one ", mt, " row in calibration CSV.")
  full_err_app[mt] <- as.numeric(rr$Apparent)
  full_err_opt[mt] <- as.numeric(rr$Mean_optimism_complete_pipeline_B1000)
  full_err_corr[mt] <- as.numeric(rr$Optimism_corrected)
}

full_corrected[err_names] <- full_err_corr
full_mean_opt <- c(
  full_app_metrics[metric_names] - full_corrected[metric_names],
  full_err_opt
)
full_mean_opt <- full_mean_opt[c(metric_names, err_names)]

# ---------------- REPAIR B=5000 ICI/E50/E90 USING SAVED INDICES --------------
message("[PATCH B] Recomputing ONLY B=5000 ICI/E50/E90 with robust spline calibration...")

if (!is.matrix(uncertainty_indices) || nrow(uncertainty_indices) != B_UNCERTAINTY) {
  stop("uncertainty_indices is not the expected B=5000 index matrix.")
}

for (b in seq_len(B_UNCERTAINTY)) {
  idx <- uncertainty_indices[b, ]
  yb <- y[idx]
  for (nm in model_names) {
    pb <- pred_list[[nm]][idx]
    uncertainty_store[[nm]][b, err_names] <- get_calibration_errors_robust(yb, pb)
  }
  if (b %% 500L == 0L) message("  uncertainty repair completed ", b, "/", B_UNCERTAINTY)
}

# Hard audit before producing any table.
finite_audit <- do.call(rbind, lapply(model_names, function(nm) {
  data.frame(
    model = MODEL_LABELS[[nm]],
    metric = colnames(uncertainty_store[[nm]]),
    valid_B = colSums(is.finite(uncertainty_store[[nm]])),
    stringsAsFactors = FALSE
  )
}))
write.csv(finite_audit, paste0(OUT_PREFIX, "06B_B5000_finite_metric_audit.csv"), row.names = FALSE)

bad_audit <- finite_audit[finite_audit$valid_B < floor(0.80 * B_UNCERTAINTY), , drop = FALSE]
if (nrow(bad_audit) > 0L) {
  print(bad_audit)
  stop("At least one B=5000 metric has <80% valid replicates. See 06B_B5000_finite_metric_audit.csv")
}

# ------------------ REBUILD CORRECTED METRICS + UNCERTAINTY ------------------
all_metric_names <- c(metric_names, err_names)
model_corrected <- list(Full = full_corrected[all_metric_names])
model_mean_opt <- list(Full = full_mean_opt[all_metric_names])

for (nm in names(score_specs)) {
  app_all <- c(
    get_metrics(y, pred_list[[nm]]),
    get_calibration_errors_robust(y, pred_list[[nm]])
  )
  corr_all <- c(score_corrected[[nm]], score_err_corrected[[nm]])
  model_corrected[[nm]] <- corr_all[all_metric_names]
  model_mean_opt[[nm]] <- app_all[all_metric_names] - corr_all[all_metric_names]
}

for (nm in model_names) {
  if (any(!is.finite(model_corrected[[nm]][all_metric_names]))) {
    stop("Non-finite corrected estimate remains for model: ", nm)
  }
  if (any(!is.finite(model_mean_opt[[nm]][all_metric_names]))) {
    stop("Non-finite mean optimism remains for model: ", nm)
  }
}

uncertainty_summary_rows <- list()
shifted_boot <- list()
for (nm in model_names) {
  shifted <- sweep(
    uncertainty_store[[nm]],
    2L,
    model_mean_opt[[nm]][all_metric_names],
    FUN = "-"
  )
  shifted_boot[[nm]] <- shifted

  for (mt in all_metric_names) {
    vals <- shifted[, mt]
    good <- is.finite(vals)
    if (sum(good) < floor(0.80 * B_UNCERTAINTY)) {
      stop("Too few finite shifted bootstrap values for ", nm, " / ", mt,
           ": ", sum(good), "/", B_UNCERTAINTY)
    }
    uncertainty_summary_rows[[length(uncertainty_summary_rows) + 1L]] <- data.frame(
      model = MODEL_LABELS[[nm]],
      metric = mt,
      optimism_corrected_estimate = as.numeric(model_corrected[[nm]][mt]),
      lower_95_uncertainty = safe_quantile2(vals, 0.025, paste(nm, mt)),
      upper_95_uncertainty = safe_quantile2(vals, 0.975, paste(nm, mt)),
      valid_B = sum(good),
      method = "Patient-level bootstrap of fixed apparent predictions shifted by B=1000 mean optimism; not a two-stage bootstrap CI",
      stringsAsFactors = FALSE
    )
  }
}
uncertainty_summary <- do.call(rbind, uncertainty_summary_rows)
write.csv(uncertainty_summary,
          paste0(OUT_PREFIX, "07_all_models_corrected_metrics_95_uncertainty.csv"),
          row.names = FALSE)

# --------------------- PAIRED DIFFERENCES WITH TRUE IDs -----------------------
paired_rows <- list()
paired_dist_rows <- list()
for (nm in names(score_specs)) {
  for (mt in all_metric_names) {
    d_raw <- shifted_boot$Full[, mt] - shifted_boot[[nm]][, mt]
    good_id <- which(is.finite(d_raw))
    d <- d_raw[good_id]

    if (length(d) < floor(0.80 * B_UNCERTAINTY)) {
      stop("Too few finite paired differences for ", nm, " / ", mt,
           ": ", length(d), "/", B_UNCERTAINTY)
    }

    est <- as.numeric(model_corrected$Full[mt] - model_corrected[[nm]][mt])
    paired_rows[[length(paired_rows) + 1L]] <- data.frame(
      comparison = paste0("Final model - ", MODEL_LABELS[[nm]]),
      metric = mt,
      corrected_difference = est,
      lower_95_uncertainty = safe_quantile2(d, 0.025, paste("paired", nm, mt)),
      upper_95_uncertainty = safe_quantile2(d, 0.975, paste("paired", nm, mt)),
      valid_B = length(d),
      stringsAsFactors = FALSE
    )

    paired_dist_rows[[length(paired_dist_rows) + 1L]] <- data.frame(
      bootstrap_id = good_id,
      comparison = paste0("Final model - ", MODEL_LABELS[[nm]]),
      metric = mt,
      difference = d,
      stringsAsFactors = FALSE
    )
  }
}
paired_summary <- do.call(rbind, paired_rows)
paired_distribution <- do.call(rbind, paired_dist_rows)
write.csv(paired_summary,
          paste0(OUT_PREFIX, "08_paired_corrected_differences_95_uncertainty.csv"),
          row.names = FALSE)

gz_path <- paste0(OUT_PREFIX, "09_paired_bootstrap_difference_distributions.csv.gz")
gz_con <- gzfile(gz_path, "wt")
on.exit(try(close(gz_con), silent = TRUE), add = TRUE)
write.csv(paired_distribution, gz_con, row.names = FALSE)
close(gz_con)

# -------------------------- DELONG AUC TESTS ----------------------------------
delong_rows <- list()
roc_full <- pROC::roc(y, pred_list$Full, quiet = TRUE, direction = "<")
for (nm in names(score_specs)) {
  roc_score <- pROC::roc(y, pred_list[[nm]], quiet = TRUE, direction = "<")
  tst <- pROC::roc.test(roc_full, roc_score, paired = TRUE, method = "delong")
  delong_rows[[nm]] <- data.frame(
    comparison = paste0("Final model vs ", MODEL_LABELS[[nm]]),
    apparent_delta_AUC = safe_auc(y, pred_list$Full) - safe_auc(y, pred_list[[nm]]),
    p_value_delong = as.numeric(tst$p.value),
    stringsAsFactors = FALSE
  )
}
delong_table <- do.call(rbind, delong_rows)
delong_table$p_value_holm <- stats::p.adjust(delong_table$p_value_delong, method = "holm")
write.csv(delong_table, paste0(OUT_PREFIX, "10_paired_DeLong_AUC_tests_Holm.csv"), row.names = FALSE)

# ------------------- FORMAL INCREMENTAL CONTRIBUTION: D1 ---------------------
message("[6/10] Testing formal incremental contribution beyond SAPS II (pooled D1)...")
fit_full_mira <- with(
  imp,
  stats::glm(death_28days ~ sapsii + rr + pf_ratio + immunosuppressant,
             family = stats::binomial())
)
fit_sapsii_mira <- with(
  imp,
  stats::glm(death_28days ~ sapsii, family = stats::binomial())
)

d1_obj <- tryCatch(mice::D1(fit_full_mira, fit_sapsii_mira), error = function(e) e)
if (inherits(d1_obj, "error")) {
  d1_table <- data.frame(
    comparison = "SAPS II vs SAPS II + RR + P/F ratio + immunosuppressant",
    status = "D1_ERROR",
    message = conditionMessage(d1_obj),
    stringsAsFactors = FALSE
  )
} else {
  d1_table <- as.data.frame(d1_obj)
  d1_table$comparison <- "SAPS II vs SAPS II + RR + P/F ratio + immunosuppressant"
  d1_table$status <- "OK"
}
write.csv(d1_table, paste0(OUT_PREFIX, "11_incremental_value_D1_nested_test.csv"), row.names = FALSE)
saveRDS(
  list(full = fit_full_mira, sapsii = fit_sapsii_mira, D1 = d1_obj),
  paste0(OUT_PREFIX, "11_incremental_value_D1_objects.rds")
)

# ----------------------- DCA + CALIBRATION CURVES -----------------------------
message("[7/10] Computing DCA and benchmark calibration curves...")

dca_rows <- list()
for (nm in model_names) {
  dd <- calc_dca(y, pred_list[[nm]], DCA_THRESHOLDS)
  dd$model <- MODEL_LABELS[[nm]]
  dca_rows[[nm]] <- dd
}
dca_table <- do.call(rbind, dca_rows)
write.csv(dca_table, paste0(OUT_PREFIX, "12_decision_curve_net_benefit.csv"), row.names = FALSE)

full_curve_for_benchmark <- data.frame(
  predicted_probability = cal_curve_old$Predicted_probability,
  observed_probability = cal_curve_old$Optimism_corrected_observed,
  model = "Final model",
  stringsAsFactors = FALSE
)
score_curve_rows <- lapply(names(score_specs), function(nm) {
  data.frame(
    predicted_probability = score_cal_grid[[nm]],
    observed_probability = score_curve_corrected[[nm]],
    model = MODEL_LABELS[[nm]],
    stringsAsFactors = FALSE
  )
})
names(score_curve_rows) <- names(score_specs)
benchmark_cal_curves <- rbind(full_curve_for_benchmark, do.call(rbind, score_curve_rows))
write.csv(benchmark_cal_curves,
          paste0(OUT_PREFIX, "13_optimism_corrected_calibration_curves.csv"),
          row.names = FALSE)

# Save complete curve-bootstrap objects now so no future refit is needed.
saveRDS(
  list(
    score_cal_grid = score_cal_grid,
    score_boot_curves = score_boot_curves,
    score_curve_corrected = score_curve_corrected,
    score_boot_errors = score_boot_errors,
    score_error_corrected = score_err_corrected
  ),
  paste0(OUT_PREFIX, "13B_score_calibration_bootstrap_objects.rds"),
  compress = "xz"
)

# ------------------------------- FIGURES --------------------------------------
message("[8/10] Saving publication-ready benchmark figures...")

roc_plot_df <- do.call(rbind, lapply(model_names, function(nm) {
  r <- pROC::roc(y, pred_list[[nm]], quiet = TRUE, direction = "<")
  data.frame(
    specificity = as.numeric(r$specificities),
    sensitivity = as.numeric(r$sensitivities),
    model = MODEL_LABELS[[nm]],
    stringsAsFactors = FALSE
  )
}))

p_roc <- ggplot2::ggplot(
  roc_plot_df,
  ggplot2::aes(x = 1 - specificity, y = sensitivity, linetype = model)
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.5) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(x = "1 - Specificity", y = "Sensitivity", linetype = NULL) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(paste0(OUT_PREFIX, "14_ROC_benchmark.png"), p_roc, width = 7, height = 7, dpi = 600)
ggplot2::ggsave(paste0(OUT_PREFIX, "14_ROC_benchmark.pdf"), p_roc, width = 7, height = 7)

p_cal <- ggplot2::ggplot(
  benchmark_cal_curves,
  ggplot2::aes(x = predicted_probability, y = observed_probability, linetype = model)
) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.5) +
  ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    x = "Predicted 28-day mortality risk",
    y = "Observed 28-day mortality probability",
    linetype = NULL
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(paste0(OUT_PREFIX, "15_calibration_benchmark.png"), p_cal, width = 7.5, height = 7, dpi = 600)
ggplot2::ggsave(paste0(OUT_PREFIX, "15_calibration_benchmark.pdf"), p_cal, width = 7.5, height = 7)

model_dca <- dca_table[, c("threshold", "net_benefit", "model")]
base_dca <- unique(dca_table[, c("threshold", "treat_all", "treat_none")])
base_long <- rbind(
  data.frame(threshold = base_dca$threshold, net_benefit = base_dca$treat_all, model = "Treat all"),
  data.frame(threshold = base_dca$threshold, net_benefit = base_dca$treat_none, model = "Treat none")
)
dca_plot_df <- rbind(model_dca, base_long)

p_dca <- ggplot2::ggplot(
  dca_plot_df,
  ggplot2::aes(x = threshold, y = net_benefit, linetype = model)
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::labs(
    x = "Threshold probability",
    y = "Net benefit",
    linetype = NULL,
    caption = paste0(
      "Thresholds represent risk levels at which intensified prognostic assessment/monitoring would be considered; ",
      "they do not represent decisions to initiate mechanical ventilation."
    )
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    legend.position = "right",
    plot.caption = ggplot2::element_text(size = 9, hjust = 0)
  )
ggplot2::ggsave(paste0(OUT_PREFIX, "16_decision_curve_analysis.png"), p_dca, width = 8.5, height = 6.5, dpi = 600)
ggplot2::ggsave(paste0(OUT_PREFIX, "16_decision_curve_analysis.pdf"), p_dca, width = 8.5, height = 6.5)

auc_df <- uncertainty_summary[uncertainty_summary$metric == "AUC", , drop = FALSE]
auc_df$model <- factor(auc_df$model, levels = rev(unname(MODEL_LABELS)))
p_auc <- ggplot2::ggplot(
  auc_df,
  ggplot2::aes(x = optimism_corrected_estimate, y = model)
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = lower_95_uncertainty,
      xend = upper_95_uncertainty,
      y = model,
      yend = model
    ),
    linewidth = 0.7
  ) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(x = "Optimism-corrected AUC (95% uncertainty interval)", y = NULL) +
  ggplot2::theme_classic(base_size = 13)
ggplot2::ggsave(paste0(OUT_PREFIX, "17_AUC_forest.png"), p_auc, width = 7.5, height = 4.8, dpi = 600)
ggplot2::ggsave(paste0(OUT_PREFIX, "17_AUC_forest.pdf"), p_auc, width = 7.5, height = 4.8)

# -------------------------- MANUSCRIPT-READY TABLES ---------------------------
message("[9/10] Creating manuscript-ready summary tables...")

main_metrics <- c("AUC", "Brier", "Calibration_intercept", "Calibration_slope")
main_table <- uncertainty_summary[
  uncertainty_summary$metric %in% main_metrics,
  c("model", "metric", "optimism_corrected_estimate", "lower_95_uncertainty", "upper_95_uncertainty", "valid_B")
]
write.csv(main_table, paste0(OUT_PREFIX, "18_main_benchmark_table.csv"), row.names = FALSE)

cal_error_table <- uncertainty_summary[
  uncertainty_summary$metric %in% err_names,
  c("model", "metric", "optimism_corrected_estimate", "lower_95_uncertainty", "upper_95_uncertainty", "valid_B")
]
write.csv(cal_error_table, paste0(OUT_PREFIX, "19_calibration_error_metrics_table.csv"), row.names = FALSE)

run_status <- data.frame(
  item = c(
    "N", "Events", "Event_rate", "Initial_candidate_predictors", "Crude_EPV_44",
    "Final_model_variables", "Score_benchmarks", "B_score_optimism", "B_uncertainty",
    "Same_V3_bootstrap_indices_used_for_scores", "LASSO_CV_rerun", "Bootstrap_MICE_rerun",
    "Apparent_MICE_runs", "Calibration_error_patch", "Output_prefix"
  ),
  value = c(
    nrow(dat), sum(y), mean(y), length(candidate_vars), sum(y) / length(candidate_vars),
    paste(FINAL_VARS, collapse = " + "), paste(SCORE_VARS, collapse = ", "),
    B_OPTIMISM, B_UNCERTAINTY, TRUE, FALSE, FALSE, 1L,
    "natural-cubic-spline logistic calibration; B1000 score calibration + B5000 error metrics recomputed",
    OUT_PREFIX
  ),
  stringsAsFactors = FALSE
)
write.csv(run_status, paste0(OUT_PREFIX, "20_run_status.csv"), row.names = FALSE)

# Save EVERYTHING needed for future post-processing.
benchmark_object <- list(
  settings = list(
    expected_n = EXPECTED_N,
    expected_events = EXPECTED_EVENTS,
    outcome = OUTCOME,
    final_vars = FINAL_VARS,
    score_vars = SCORE_VARS,
    candidate_vars = candidate_vars,
    M = M,
    mice_maxit = MICE_MAXIT,
    seed_apparent = SEED_APPARENT,
    seed_boot_base = SEED_BOOT_BASE,
    B_optimism = B_OPTIMISM,
    seed_uncertainty = SEED_UNCERTAINTY,
    B_uncertainty = B_UNCERTAINTY,
    dca_thresholds = DCA_THRESHOLDS,
    calibration_patch = "natural cubic spline logistic calibration on logit predicted risk"
  ),
  patient_predictions = patient_predictions,
  bootstrap_indices_B1000 = bootstrap_indices,
  uncertainty_indices_B5000 = uncertainty_indices,
  full_model_reconstruction_check = recon_check,
  full_corrected = full_corrected,
  full_mean_optimism = full_mean_opt,
  score_corrected = score_corrected,
  score_error_corrected = score_err_corrected,
  score_boot_scalar = score_boot_scalar,
  score_boot_errors = score_boot_errors,
  score_cal_grid = score_cal_grid,
  score_boot_curves = score_boot_curves,
  score_curve_corrected = score_curve_corrected,
  uncertainty_store_apparent = uncertainty_store,
  uncertainty_store_shifted = shifted_boot,
  uncertainty_summary = uncertainty_summary,
  paired_summary = paired_summary,
  paired_distribution = paired_distribution,
  finite_metric_audit = finite_audit,
  DeLong = delong_table,
  D1 = d1_obj,
  DCA = dca_table,
  calibration_curves = benchmark_cal_curves,
  coefficients = coef_out,
  run_status = run_status
)
saveRDS(benchmark_object,
        paste0(OUT_PREFIX, "21_COMPLETE_BENCHMARK_OBJECT.rds"),
        compress = "xz")

writeLines(capture.output(sessionInfo()), paste0(OUT_PREFIX, "22_sessionInfo.txt"))
manifest <- data.frame(
  file = sort(list.files(pattern = paste0("^", OUT_PREFIX))),
  stringsAsFactors = FALSE
)
write.csv(manifest, paste0(OUT_PREFIX, "23_output_manifest.csv"), row.names = FALSE)

message("[10/10] DONE — benchmark continuation completed successfully.")
message("Key outputs:")
message("  ", OUT_PREFIX, "06B_B5000_finite_metric_audit.csv")
message("  ", OUT_PREFIX, "07_all_models_corrected_metrics_95_uncertainty.csv")
message("  ", OUT_PREFIX, "08_paired_corrected_differences_95_uncertainty.csv")
message("  ", OUT_PREFIX, "10_paired_DeLong_AUC_tests_Holm.csv")
message("  ", OUT_PREFIX, "11_incremental_value_D1_nested_test.csv")
message("  ", OUT_PREFIX, "12_decision_curve_net_benefit.csv")
message("  ", OUT_PREFIX, "18_main_benchmark_table.csv")
message("  ", OUT_PREFIX, "19_calibration_error_metrics_table.csv")
message("  ", OUT_PREFIX, "21_COMPLETE_BENCHMARK_OBJECT.rds")
message("No MICE/LASSO/CV was rerun in this continuation script.")
