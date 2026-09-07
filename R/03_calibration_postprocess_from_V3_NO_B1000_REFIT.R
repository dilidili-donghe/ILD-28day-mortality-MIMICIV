# ==============================================================================
# 03_calibration_postprocess_from_V3_NO_B1000_REFIT.R
# FINAL calibration post-processing for the frozen 390-patient V3 analysis
#
# PURPOSE
#   Build the remaining calibration outputs WITHOUT rerunning the B=1000
#   model-development pipeline (no repeated LASSO/CV and no repeated logistic
#   refitting).
#
# KEY IDEA
#   The V3 complete object already stores, for every bootstrap replicate:
#     - the selected variables
#     - the fitted coefficients for each of the 5 imputations
#     - deterministic seed settings
#     - the stored scalar bootstrap/original performance
#   Although patient-level predictions were not saved, they can be reconstructed
#   deterministically by replaying ONLY the bootstrap sampling + MICE imputation
#   and then applying the SAVED coefficients. Selection and model fitting are NOT
#   repeated.
#
# OUTPUTS
#   1) optimism-corrected flexible calibration curve based on the ORIGINAL V3 B=1000
#   2) optimism-corrected ICI / E50 / E90 based on the ORIGINAL V3 B=1000
#   3) optimism-shifted bootstrap 95% intervals for calibration intercept/slope
#   4) 95% pointwise bootstrap uncertainty band for the corrected flexible curve
#   5) publication-ready PNG/PDF calibration figure
#
# IMPORTANT TERMINOLOGY
#   - The scalar 95% intervals and pointwise band are post-processing uncertainty
#     intervals/bands. They are NOT a two-stage/nested bootstrap confidence interval.
#   - The pointwise MEAN optimism for the flexible curve, and the mean optimism for
#     ICI/E50/E90, DO come from reconstructing the ORIGINAL V3 complete-pipeline
#     bootstrap predictions using the stored coefficients and original seeds.
#
# REQUIREMENTS IN work_dir
#   - original 390-patient .xlsx workbook used by V3
#   - 390DD_BOOT1000_V3_11_complete_bootstrap_object.rds
#
# This script intentionally does NOT call glmnet::cv.glmnet() and does NOT refit
# bootstrap logistic models.
# ==============================================================================

rm(list = ls())
gc()

# ----------------------------- 0. USER SETTINGS -------------------------------
work_dir <- "D:/R代码/大修1_1"
setwd(work_dir)

OBJECT_FILE <- file.path(work_dir, "390DD_BOOT1000_V3_11_complete_bootstrap_object.rds")
OUT_PREFIX <- "390CAL_POSTV3_"
B_UNCERTAINTY <- 5000L
SEED_UNCERTAINTY <- 202608271L
GRID_N <- 101L
# Flexible curve is reported over the central prediction range to avoid unstable
# loess extrapolation at the extreme tails.
GRID_Q_LO <- 0.025
GRID_Q_HI <- 0.975
# Reconstruction must reproduce stored scalar performance very closely.
VERIFY_TOL <- 5e-6

# ----------------------------- 1. PACKAGES ------------------------------------
required_pkgs <- c("readxl", "mice")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Please install required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readxl)
  library(mice)
})

# ----------------------------- 2. LOAD V3 OBJECT ------------------------------
if (!file.exists(OBJECT_FILE)) {
  stop(
    "Cannot find: ", basename(OBJECT_FILE), "\n",
    "This file was written by the completed V3 B=1000 script. ",
    "Do NOT rerun B=1000. Locate/copy this RDS into work_dir and rerun this post-processing script."
  )
}

obj <- readRDS(OBJECT_FILE)
required_obj <- c("settings", "apparent", "results", "performance", "summary", "status")
missing_obj <- setdiff(required_obj, names(obj))
if (length(missing_obj) > 0L) {
  stop("V3 complete object is missing: ", paste(missing_obj, collapse = ", "))
}

S <- obj$settings
EXPECTED_N <- as.integer(S$EXPECTED_N)
OUTCOME <- as.character(S$OUTCOME)
M <- as.integer(S$M)
MICE_MAXIT <- as.integer(S$MICE_MAXIT)
B <- as.integer(S$B)
candidate_vars_saved <- as.character(S$candidate_vars)
SEED_APPARENT <- as.integer(S$seed_apparent)
SEED_BOOT_BASE <- as.integer(S$seed_boot_base)
SEED_MICE_BASE <- as.integer(S$seed_mice_base)

if (EXPECTED_N != 390L) stop("The V3 object is not the frozen N=390 analysis.")
if (B != 1000L) stop("The V3 object does not contain B=1000 settings.")
if (length(obj$results) < B) stop("The V3 object contains fewer than 1000 result slots.")

success_idx <- which(vapply(obj$results, function(x) !is.null(x), logical(1)))
if (length(success_idx) != 1000L) {
  stop("Expected 1000 successful V3 replicates in the complete object, found ", length(success_idx), ".")
}

# ----------------------------- 3. HELPERS -------------------------------------
normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "null")] <- NA_character_
  x
}

is_yes <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z %in% c("是", "yes", "y", "1", "true", "include", "included", "纳入")
}

find_excel_file <- function(dir) {
  ff <- list.files(dir, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  ff <- ff[!grepl("^~\\$", basename(ff))]
  ff390 <- ff[grepl("^390", basename(ff), ignore.case = TRUE)]
  if (length(ff390) == 1L) return(ff390)
  if (length(ff390) > 1L) {
    cand <- ff390[!grepl("BOOT|output|result", basename(ff390), ignore.case = TRUE)]
    if (length(cand) == 1L) return(cand)
    stop("More than one 390*.xlsx workbook found. Keep only the analysis workbook or set excel_file manually.\n",
         paste(basename(ff390), collapse = "\n"))
  }
  stop("No 390*.xlsx analysis workbook found in: ", dir)
}

extract_candidate_pool <- function(sheet2, data_names) {
  nm <- names(sheet2)
  nm_low <- tolower(trimws(nm))
  include_idx <- which(grepl("是否纳入|纳入|include|included", nm_low, perl = TRUE))
  if (length(include_idx) == 0L) stop("Cannot identify inclusion column in Sheet2.")
  include_idx <- include_idx[1L]

  overlap_n <- vapply(seq_along(sheet2), function(j) {
    vals <- normalize_chr(sheet2[[j]])
    sum(vals %in% data_names, na.rm = TRUE)
  }, numeric(1))
  overlap_n[include_idx] <- -Inf
  var_idx <- which.max(overlap_n)
  if (!is.finite(overlap_n[var_idx]) || overlap_n[var_idx] <= 0) {
    hits <- which(grepl("变量|variable|predictor|字段|field", nm_low, perl = TRUE))
    hits <- setdiff(hits, include_idx)
    if (length(hits) == 0L) stop("Cannot identify variable-name column in Sheet2.")
    var_idx <- hits[1L]
  }

  vars <- normalize_chr(sheet2[[var_idx]])
  keep <- is_yes(sheet2[[include_idx]])
  vars <- unique(vars[keep & !is.na(vars)])
  vars <- vars[vars %in% data_names]
  setdiff(vars, OUTCOME)
}

coerce_numeric_predictors <- function(dat, vars) {
  missing_tokens <- c("", "na", "n/a", "n.a.", "null", "none", "missing", ".", "-")
  true_tokens <- c("1", "yes", "y", "true", "是")
  false_tokens <- c("0", "no", "n", "false", "否")

  for (v in vars) {
    x0 <- dat[[v]]
    if (is.logical(x0)) x0 <- as.integer(x0)
    if (is.factor(x0)) x0 <- as.character(x0)

    if (is.character(x0)) {
      raw_chr <- trimws(x0)
      z <- tolower(raw_chr)
      token_missing <- is.na(raw_chr) | z %in% missing_tokens
      raw_chr[token_missing] <- NA_character_
      z[token_missing] <- NA_character_
      vals <- unique(z[!is.na(z)])

      if (length(vals) > 0L && all(vals %in% c(true_tokens, false_tokens))) {
        xx <- rep(NA_real_, length(z))
        xx[z %in% true_tokens] <- 1
        xx[z %in% false_tokens] <- 0
        x <- xx
      } else {
        numeric_chr <- gsub(",", "", raw_chr, fixed = TRUE)
        suppressWarnings(xx <- as.numeric(numeric_chr))
        bad <- !is.na(raw_chr) & is.na(xx)
        if (any(bad)) {
          ex <- unique(raw_chr[bad])
          stop("Predictor '", v, "' contains invalid nonnumeric value(s): ",
               paste(head(ex, 5L), collapse = ", "))
        }
        x <- xx
      }
    } else {
      suppressWarnings(x <- as.numeric(x0))
    }

    x[!is.finite(x)] <- NA_real_
    dat[[v]] <- x
  }
  dat
}

build_mice_setup <- function(dat, candidate_vars) {
  meth <- rep("", ncol(dat)); names(meth) <- names(dat)
  for (v in candidate_vars) if (anyNA(dat[[v]])) meth[v] <- "pmm"
  meth[OUTCOME] <- ""

  pred <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat),
                 dimnames = list(names(dat), names(dat)))
  for (v in candidate_vars) if (meth[v] != "") pred[v, candidate_vars] <- 1L
  diag(pred) <- 0L
  pred[, OUTCOME] <- 0L
  pred[OUTCOME, ] <- 0L

  pt_names <- intersect(c("pt", "PT"), names(dat))
  inr_names <- intersect(c("inr", "INR"), names(dat))
  if (length(pt_names) > 0L && length(inr_names) > 0L) {
    for (ptv in pt_names) for (inrv in inr_names) {
      pred[ptv, inrv] <- 0L
      pred[inrv, ptv] <- 0L
    }
  }
  list(method = meth, predictorMatrix = pred)
}

run_mice <- function(dat, candidate_vars, seed, ignore = NULL) {
  setup <- build_mice_setup(dat, candidate_vars)
  suppressWarnings(
    mice(
      dat,
      m = M,
      maxit = MICE_MAXIT,
      method = setup$method,
      predictorMatrix = setup$predictorMatrix,
      seed = seed,
      printFlag = FALSE,
      ignore = ignore
    )
  )
}

clip_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

auc_rank <- function(y, p) {
  ok <- is.finite(p) & !is.na(y)
  y <- as.integer(y[ok]); p <- as.numeric(p[ok])
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

brier_score <- function(y, p) {
  mean((as.numeric(y) - as.numeric(p))^2)
}

calibration_metrics <- function(y, p) {
  ok <- !is.na(y) & is.finite(p)
  y <- as.integer(y[ok]); p <- clip_prob(p[ok])
  if (length(y) < 10L || length(unique(y)) < 2L) {
    return(c(cal_intercept = NA_real_, cal_slope = NA_real_))
  }
  lp <- qlogis(p)
  cint <- tryCatch({
    fit_i <- suppressWarnings(glm(y ~ 1, family = binomial(), offset = lp))
    as.numeric(coef(fit_i)[1L])
  }, error = function(e) NA_real_)
  cslope <- tryCatch({
    fit_s <- suppressWarnings(glm(y ~ lp, family = binomial()))
    as.numeric(coef(fit_s)["lp"])
  }, error = function(e) NA_real_)
  c(cal_intercept = cint, cal_slope = cslope)
}

scalar_performance <- function(y, p) {
  cm <- calibration_metrics(y, p)
  c(
    AUC = auc_rank(y, p),
    Brier = brier_score(y, p),
    cal_intercept = unname(cm["cal_intercept"]),
    cal_slope = unname(cm["cal_slope"])
  )
}

# Apply coefficients SAVED by the original V3 run. No model fitting occurs here.
predict_from_saved_coefficients <- function(df, coef_df, imputation_id) {
  zz <- coef_df[coef_df$imputation == imputation_id, , drop = FALSE]
  if (nrow(zz) == 0L) stop("No saved coefficients for imputation ", imputation_id)

  beta <- setNames(as.numeric(zz$estimate), as.character(zz$term))
  if (!"(Intercept)" %in% names(beta)) stop("Saved coefficients lack an intercept.")

  lp <- rep(unname(beta["(Intercept)"]), nrow(df))
  terms <- setdiff(names(beta), "(Intercept)")
  if (length(terms) > 0L) {
    missing_terms <- setdiff(terms, names(df))
    if (length(missing_terms) > 0L) {
      stop("Saved coefficient term(s) absent from completed data: ",
           paste(missing_terms, collapse = ", "))
    }
    for (v in terms) lp <- lp + beta[v] * as.numeric(df[[v]])
  }
  plogis(lp)
}

average_saved_predictions <- function(completed_list, coef_df) {
  pmat <- matrix(NA_real_, nrow = nrow(completed_list[[1L]]), ncol = M)
  for (i in seq_len(M)) {
    pmat[, i] <- predict_from_saved_coefficients(completed_list[[i]], coef_df, i)
  }
  rowMeans(pmat)
}

safe_loess_predict <- function(y, p, new_p) {
  ok <- is.finite(p) & !is.na(y)
  y <- as.numeric(y[ok]); p <- as.numeric(p[ok])
  if (length(y) < 20L || length(unique(y)) < 2L || length(unique(p)) < 8L) {
    return(rep(NA_real_, length(new_p)))
  }
  fit <- try(
    loess(y ~ p, degree = 1, span = 0.75,
          control = loess.control(surface = "direct")),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(rep(NA_real_, length(new_p)))
  z <- try(predict(fit, newdata = data.frame(p = new_p)), silent = TRUE)
  if (inherits(z, "try-error")) return(rep(NA_real_, length(new_p)))
  as.numeric(z)
}

smooth_error_metrics <- function(y, p) {
  # Estimate the flexible observed probability at each subject's own predicted risk,
  # then summarize |observed_smooth - predicted|.
  obs_hat <- safe_loess_predict(y, p, p)
  ae <- abs(obs_hat - p)
  ae <- ae[is.finite(ae)]
  if (length(ae) == 0L) return(c(ICI = NA_real_, E50 = NA_real_, E90 = NA_real_))
  c(
    ICI = mean(ae),
    E50 = unname(quantile(ae, 0.50, na.rm = TRUE, names = FALSE, type = 6)),
    E90 = unname(quantile(ae, 0.90, na.rm = TRUE, names = FALSE, type = 6))
  )
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 6))
}

get_summary_value <- function(metric_name, col_name) {
  sm <- obj$summary
  row <- which(sm$metric == metric_name)
  if (length(row) != 1L) stop("Cannot uniquely locate metric in V3 summary: ", metric_name)
  as.numeric(sm[row, col_name])
}

# ----------------------------- 4. RELOAD ORIGINAL DATA ------------------------
excel_file <- find_excel_file(work_dir)
message("Analysis workbook: ", basename(excel_file))
sheets <- excel_sheets(excel_file)
if (length(sheets) < 2L) stop("Workbook must contain at least two sheets.")

raw <- as.data.frame(read_excel(excel_file, sheet = 1), check.names = FALSE)
sheet2 <- as.data.frame(read_excel(excel_file, sheet = 2), check.names = FALSE)
if (nrow(raw) != EXPECTED_N) stop("Expected ", EXPECTED_N, " rows; found ", nrow(raw), ".")
if (!OUTCOME %in% names(raw)) stop("Outcome missing: ", OUTCOME)

candidate_vars_now <- extract_candidate_pool(sheet2, names(raw))
if (!identical(candidate_vars_now, candidate_vars_saved)) {
  stop(
    "Current Sheet2 candidate pool/order does not exactly match the frozen V3 object.\n",
    "Do not continue until the original V3 workbook is restored."
  )
}

dat <- raw[, c(OUTCOME, candidate_vars_saved), drop = FALSE]
dat <- coerce_numeric_predictors(dat, candidate_vars_saved)

y <- dat[[OUTCOME]]
if (is.factor(y)) y <- as.character(y)
if (is.logical(y)) y <- as.integer(y)
if (is.character(y)) {
  z <- trimws(tolower(y))
  if (!all(z %in% c("0", "1"))) stop("Outcome must be complete and coded 0/1.")
  y <- as.numeric(z)
} else {
  suppressWarnings(y <- as.numeric(y))
}
if (anyNA(y) || !all(y %in% c(0, 1))) stop("Outcome must be complete and coded 0/1.")
dat[[OUTCOME]] <- as.integer(y)

# ----------------------------- 5. RECONSTRUCT APPARENT PREDICTIONS ------------
message("\nReconstructing apparent predictions from saved V3 coefficients...")
imp_app <- run_mice(dat, candidate_vars_saved, seed = SEED_APPARENT)
app_list <- lapply(seq_len(M), function(i) complete(imp_app, i))
p_app <- average_saved_predictions(app_list, obj$apparent$coefficients)
y_app <- dat[[OUTCOME]]

app_reconstructed <- scalar_performance(y_app, p_app)
app_stored <- obj$apparent$performance[names(app_reconstructed)]
app_diff <- abs(as.numeric(app_reconstructed) - as.numeric(app_stored))

verification_app <- data.frame(
  metric = names(app_reconstructed),
  reconstructed = as.numeric(app_reconstructed),
  stored_V3 = as.numeric(app_stored),
  abs_difference = app_diff,
  stringsAsFactors = FALSE
)
write.csv(verification_app,
          paste0(OUT_PREFIX, "00_apparent_reconstruction_check.csv"),
          row.names = FALSE)

if (any(app_diff > VERIFY_TOL, na.rm = TRUE)) {
  stop(
    "Apparent prediction reconstruction does not reproduce the frozen V3 metrics within tolerance. ",
    "Do not use the calibration outputs. Check workbook/package versions."
  )
}
message("Apparent reconstruction check PASSED.")

# Grid based on central 95% of the apparent prediction distribution.
grid_lo <- unname(quantile(p_app, GRID_Q_LO, na.rm = TRUE, names = FALSE))
grid_hi <- unname(quantile(p_app, GRID_Q_HI, na.rm = TRUE, names = FALSE))
cal_grid <- seq(grid_lo, grid_hi, length.out = GRID_N)
app_curve <- safe_loess_predict(y_app, p_app, cal_grid)
app_smooth_metrics <- smooth_error_metrics(y_app, p_app)

# ----------------------------- 6. RECONSTRUCT ORIGINAL B=1000 PREDICTIONS -----
# IMPORTANT: This section replays only sampling + MICE and applies SAVED V3
# coefficients. No LASSO/CV and no logistic model refitting.
message("\nReconstructing patient-level predictions for the frozen V3 B=1000...")

curve_opt_mat <- matrix(NA_real_, nrow = B, ncol = GRID_N)
smooth_opt_mat <- matrix(NA_real_, nrow = B, ncol = 3L,
                         dimnames = list(NULL, c("ICI", "E50", "E90")))
verify_rows <- vector("list", B)

for (b in success_idx) {
  set.seed(SEED_BOOT_BASE + b)
  idx <- sample.int(EXPECTED_N, size = EXPECTED_N, replace = TRUE)
  boot_raw <- dat[idx, , drop = FALSE]

  combined <- rbind(boot_raw, dat)
  ignore <- c(rep(FALSE, EXPECTED_N), rep(TRUE, EXPECTED_N))
  imp_b <- run_mice(
    combined,
    candidate_vars_saved,
    seed = SEED_MICE_BASE + b,
    ignore = ignore
  )

  train_list <- vector("list", M)
  test_list <- vector("list", M)
  for (i in seq_len(M)) {
    comp <- complete(imp_b, i)
    train_list[[i]] <- comp[seq_len(EXPECTED_N), , drop = FALSE]
    test_list[[i]] <- comp[EXPECTED_N + seq_len(EXPECTED_N), , drop = FALSE]
  }

  coef_b <- obj$results[[b]]$coefficients
  p_boot <- average_saved_predictions(train_list, coef_b)
  p_orig <- average_saved_predictions(test_list, coef_b)

  # Exact reconstruction audit against the scalar metrics already stored by V3.
  perf_boot_r <- scalar_performance(boot_raw[[OUTCOME]], p_boot)
  perf_orig_r <- scalar_performance(dat[[OUTCOME]], p_orig)
  stored <- obj$results[[b]]$performance[1, , drop = FALSE]

  max_diff <- max(c(
    abs(perf_boot_r["AUC"] - stored$AUC_boot),
    abs(perf_orig_r["AUC"] - stored$AUC_original),
    abs(perf_boot_r["Brier"] - stored$Brier_boot),
    abs(perf_orig_r["Brier"] - stored$Brier_original),
    abs(perf_boot_r["cal_intercept"] - stored$cal_intercept_boot),
    abs(perf_orig_r["cal_intercept"] - stored$cal_intercept_original),
    abs(perf_boot_r["cal_slope"] - stored$cal_slope_boot),
    abs(perf_orig_r["cal_slope"] - stored$cal_slope_original)
  ), na.rm = TRUE)

  verify_rows[[b]] <- data.frame(bootstrap = b, max_abs_difference = max_diff)
  if (!is.finite(max_diff) || max_diff > VERIFY_TOL) {
    stop(
      "Bootstrap reconstruction mismatch at replicate ", b,
      " (max abs difference = ", signif(max_diff, 6), ").\n",
      "Stop: do not use post-processing results. Check that the same workbook and package versions are being used."
    )
  }

  curve_boot <- safe_loess_predict(boot_raw[[OUTCOME]], p_boot, cal_grid)
  curve_orig <- safe_loess_predict(dat[[OUTCOME]], p_orig, cal_grid)
  curve_opt_mat[b, ] <- curve_boot - curve_orig

  sm_boot <- smooth_error_metrics(boot_raw[[OUTCOME]], p_boot)
  sm_orig <- smooth_error_metrics(dat[[OUTCOME]], p_orig)
  smooth_opt_mat[b, ] <- sm_boot - sm_orig

  if (b %% 25L == 0L || b == B) {
    cat(sprintf("\rReconstruction: %d/%d", b, B))
    flush.console()
  }
}
cat("\n")

verification_boot <- do.call(rbind, verify_rows[success_idx])
write.csv(verification_boot,
          paste0(OUT_PREFIX, "01_B1000_reconstruction_check.csv"),
          row.names = FALSE)
message("All B=1000 reconstruction checks PASSED. No LASSO/logistic refitting was performed.")

# ----------------------------- 7. COMPLETE-PIPELINE SMOOTH CALIBRATION --------
pointwise_mean_optimism <- apply(curve_opt_mat, 2, mean, na.rm = TRUE)
pointwise_valid_B <- colSums(is.finite(curve_opt_mat))
corrected_curve <- app_curve - pointwise_mean_optimism
corrected_curve <- pmin(pmax(corrected_curve, 0), 1)

smooth_mean_optimism <- colMeans(smooth_opt_mat, na.rm = TRUE)
smooth_corrected <- app_smooth_metrics - smooth_mean_optimism

smooth_summary <- data.frame(
  Metric = c("ICI", "E50", "E90"),
  Apparent = as.numeric(app_smooth_metrics[c("ICI", "E50", "E90")]),
  Mean_optimism_complete_pipeline_B1000 = as.numeric(smooth_mean_optimism[c("ICI", "E50", "E90")]),
  Optimism_corrected = as.numeric(smooth_corrected[c("ICI", "E50", "E90")]),
  B_with_defined_optimism = colSums(is.finite(smooth_opt_mat)),
  stringsAsFactors = FALSE
)
write.csv(smooth_summary,
          paste0(OUT_PREFIX, "02_optimism_corrected_ICI_E50_E90.csv"),
          row.names = FALSE)

curve_core <- data.frame(
  Predicted_probability = cal_grid,
  Apparent_observed = app_curve,
  Mean_pointwise_optimism_complete_pipeline_B1000 = pointwise_mean_optimism,
  Optimism_corrected_observed = corrected_curve,
  B_with_defined_pointwise_optimism = pointwise_valid_B,
  stringsAsFactors = FALSE
)
write.csv(curve_core,
          paste0(OUT_PREFIX, "03_optimism_corrected_flexible_curve_core.csv"),
          row.names = FALSE)

# ----------------------------- 8. UNCERTAINTY POST-PROCESSING -----------------
# Patient-level bootstrap of the fixed apparent predictions. Then shift the
# bootstrap distribution by the FORMAL V3 B=1000 mean optimism (scalar or pointwise).
message("\nRunning lightweight patient-level bootstrap for calibration uncertainty (B=5000)...")
set.seed(SEED_UNCERTAINTY)

intercept_boot <- rep(NA_real_, B_UNCERTAINTY)
slope_boot <- rep(NA_real_, B_UNCERTAINTY)
curve_boot_fixed <- matrix(NA_real_, nrow = B_UNCERTAINTY, ncol = GRID_N)

for (r in seq_len(B_UNCERTAINTY)) {
  idx <- sample.int(EXPECTED_N, size = EXPECTED_N, replace = TRUE)
  y_r <- y_app[idx]
  p_r <- p_app[idx]
  if (length(unique(y_r)) < 2L) next

  cm <- calibration_metrics(y_r, p_r)
  intercept_boot[r] <- cm["cal_intercept"]
  slope_boot[r] <- cm["cal_slope"]
  curve_boot_fixed[r, ] <- safe_loess_predict(y_r, p_r, cal_grid)

  if (r %% 500L == 0L || r == B_UNCERTAINTY) {
    cat(sprintf("\rUncertainty bootstrap: %d/%d", r, B_UNCERTAINTY))
    flush.console()
  }
}
cat("\n")

mean_opt_intercept <- get_summary_value("cal_intercept", "mean_optimism")
mean_opt_slope <- get_summary_value("cal_slope", "mean_optimism")
corr_intercept <- get_summary_value("cal_intercept", "optimism_corrected")
corr_slope <- get_summary_value("cal_slope", "optimism_corrected")

valid_i <- intercept_boot[is.finite(intercept_boot)]
valid_s <- slope_boot[is.finite(slope_boot)]
if (length(valid_i) < 0.95 * B_UNCERTAINTY) warning(">5% invalid intercept uncertainty replicates.")
if (length(valid_s) < 0.95 * B_UNCERTAINTY) warning(">5% invalid slope uncertainty replicates.")

ci_i_app <- safe_quantile(valid_i, c(0.025, 0.975))
ci_s_app <- safe_quantile(valid_s, c(0.025, 0.975))
ci_i_corr <- ci_i_app - mean_opt_intercept
ci_s_corr <- ci_s_app - mean_opt_slope

scalar_uncertainty <- data.frame(
  Metric = c("Calibration intercept", "Calibration slope"),
  Optimism_corrected_estimate = c(corr_intercept, corr_slope),
  Lower_95 = c(ci_i_corr[1], ci_s_corr[1]),
  Upper_95 = c(ci_i_corr[2], ci_s_corr[2]),
  Mean_optimism_from_formal_V3_B1000 = c(mean_opt_intercept, mean_opt_slope),
  Valid_patient_level_bootstrap_replicates = c(length(valid_i), length(valid_s)),
  Method = rep(
    "Patient-level bootstrap of fixed apparent predictions, shifted by formal V3 B=1000 mean optimism; not a two-stage bootstrap CI",
    2L
  ),
  stringsAsFactors = FALSE
)
write.csv(scalar_uncertainty,
          paste0(OUT_PREFIX, "04_calibration_intercept_slope_95_uncertainty_interval.csv"),
          row.names = FALSE)

# Pointwise band: patient-level apparent curves shifted by the ACTUAL pointwise
# mean optimism reconstructed from the frozen V3 complete-pipeline B=1000 run.
curve_corrected_uncertainty <- sweep(
  curve_boot_fixed,
  MARGIN = 2,
  STATS = pointwise_mean_optimism,
  FUN = "-"
)

point_lo <- rep(NA_real_, GRID_N)
point_hi <- rep(NA_real_, GRID_N)
valid_curve_n <- rep(0L, GRID_N)
for (j in seq_len(GRID_N)) {
  z <- curve_corrected_uncertainty[, j]
  z <- z[is.finite(z)]
  valid_curve_n[j] <- length(z)
  if (length(z) >= 100L) {
    qq <- safe_quantile(z, c(0.025, 0.975))
    point_lo[j] <- qq[1]
    point_hi[j] <- qq[2]
  }
}
point_lo <- pmax(0, point_lo)
point_hi <- pmin(1, point_hi)

curve_band <- data.frame(
  Predicted_probability = cal_grid,
  Apparent_observed = app_curve,
  Mean_pointwise_optimism_complete_pipeline_B1000 = pointwise_mean_optimism,
  Optimism_corrected_observed = corrected_curve,
  Lower_95_pointwise_uncertainty = point_lo,
  Upper_95_pointwise_uncertainty = point_hi,
  V3_B1000_valid_pointwise_optimism = pointwise_valid_B,
  Patient_bootstrap_valid_curves = valid_curve_n,
  stringsAsFactors = FALSE
)
write.csv(curve_band,
          paste0(OUT_PREFIX, "05_flexible_calibration_curve_with_95_pointwise_band.csv"),
          row.names = FALSE)

# ----------------------------- 9. PUBLICATION-READY FIGURE --------------------
make_plot <- function(device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(paste0(OUT_PREFIX, "06_publication_calibration.png"),
        width = 1800, height = 1800, res = 300)
  } else {
    pdf(paste0(OUT_PREFIX, "06_publication_calibration.pdf"),
        width = 6.2, height = 6.2, useDingbats = FALSE)
  }

  op <- par(mar = c(4.6, 4.8, 1.1, 1.1), las = 1,
            cex.axis = 0.9, cex.lab = 1.0)
  on.exit({par(op); dev.off()}, add = TRUE)

  plot(
    cal_grid, corrected_curve,
    type = "n",
    xlim = c(0, 1), ylim = c(0, 1),
    xaxs = "i", yaxs = "i",
    xlab = "Predicted 28-day mortality risk",
    ylab = "Observed 28-day mortality probability",
    bty = "l"
  )
  abline(a = 0, b = 1, lty = 2, lwd = 1)

  good <- which(is.finite(point_lo) & is.finite(point_hi) & is.finite(corrected_curve))
  if (length(good) > 1L) {
    polygon(
      c(cal_grid[good], rev(cal_grid[good])),
      c(point_lo[good], rev(point_hi[good])),
      border = NA,
      density = 18,
      angle = 45
    )
  }
  lines(cal_grid, corrected_curve, lwd = 2)

  # Small rug of the apparent predicted-risk distribution.
  rug(p_app, side = 1, ticksize = 0.025)

  legend(
    "topleft",
    legend = c("Ideal", "Optimism-corrected flexible calibration", "95% pointwise uncertainty band"),
    lty = c(2, 1, NA),
    lwd = c(1, 2, NA),
    pch = c(NA, NA, 15),
    pt.cex = c(NA, NA, 1.2),
    bty = "n",
    cex = 0.82
  )
}

make_plot("png")
make_plot("pdf")

# ----------------------------- 10. SAVE AUDIT OBJECT --------------------------
saveRDS(
  list(
    source_V3_object = basename(OBJECT_FILE),
    settings = list(
      N = EXPECTED_N,
      B_formal = B,
      B_uncertainty = B_UNCERTAINTY,
      seed_uncertainty = SEED_UNCERTAINTY,
      grid_quantiles = c(GRID_Q_LO, GRID_Q_HI),
      grid_n = GRID_N,
      verify_tolerance = VERIFY_TOL,
      note = paste(
        "Patient-level predictions for formal V3 B=1000 were reconstructed by",
        "replaying only original bootstrap sampling and MICE with frozen seeds,",
        "then applying saved coefficients. No LASSO/CV or logistic refitting."
      )
    ),
    apparent_predictions = p_app,
    apparent_reconstruction_check = verification_app,
    B1000_reconstruction_check = verification_boot,
    smooth_summary = smooth_summary,
    curve_core = curve_core,
    scalar_uncertainty = scalar_uncertainty,
    curve_band = curve_band
  ),
  paste0(OUT_PREFIX, "07_complete_calibration_postprocessing.rds")
)

# ----------------------------- 11. CONSOLE SUMMARY ----------------------------
cat("\n============================================================\n")
cat("V3 CALIBRATION POST-PROCESSING COMPLETED\n")
cat("============================================================\n")
cat("No B=1000 LASSO/CV rerun. No bootstrap logistic refitting.\n")
cat("Formal V3 B=1000 predictions were reconstructed using saved coefficients.\n\n")
cat("Optimism-corrected ICI / E50 / E90:\n")
print(smooth_summary, row.names = FALSE)
cat("\nCalibration intercept/slope uncertainty intervals:\n")
print(scalar_uncertainty, row.names = FALSE)
cat("\nOutputs:\n")
cat(paste0(OUT_PREFIX, "00_apparent_reconstruction_check.csv\n"))
cat(paste0(OUT_PREFIX, "01_B1000_reconstruction_check.csv\n"))
cat(paste0(OUT_PREFIX, "02_optimism_corrected_ICI_E50_E90.csv\n"))
cat(paste0(OUT_PREFIX, "03_optimism_corrected_flexible_curve_core.csv\n"))
cat(paste0(OUT_PREFIX, "04_calibration_intercept_slope_95_uncertainty_interval.csv\n"))
cat(paste0(OUT_PREFIX, "05_flexible_calibration_curve_with_95_pointwise_band.csv\n"))
cat(paste0(OUT_PREFIX, "06_publication_calibration.png\n"))
cat(paste0(OUT_PREFIX, "06_publication_calibration.pdf\n"))
cat(paste0(OUT_PREFIX, "07_complete_calibration_postprocessing.rds\n"))
cat("============================================================\n")
