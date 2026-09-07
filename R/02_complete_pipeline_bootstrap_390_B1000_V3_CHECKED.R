# ==============================================================================
# 02_complete_pipeline_bootstrap_390_B1000_V3_CHECKED.R
# Complete-pipeline internal validation for the FINAL 390-patient data-driven model
#
# Design frozen for the revised manuscript:
#   - N = 390 (patients dying before the 24-h landmark already excluded)
#   - candidate predictors = Sheet2 variables marked for inclusion
#   - NO forced-in ventilation variable
#   - multiple imputation: m = 5, PMM, maxit = 20
#   - 10-fold CV LASSO, alpha = 1, lambda.1se
#   - final predictor must be selected in 5/5 imputed datasets
#   - complete-pipeline bootstrap: B = 1000
#
# V3 checked version:
#   1) mice::loggedEvents are AUDITED but no longer treated as automatic failures.
#   2) if the prespecified lambda.1se + 5/5 rule selects no predictor,
#      the replicate is evaluated as an intercept-only model (no rule relaxation).
#   3) calibration-in-the-large and calibration slope are calculated robustly from
#      the averaged MI predictions for both bootstrap and original assessment data.
#   4) deterministic replicate-specific seeds allow safe checkpoint/resume.
#   5) all outputs use a V3 prefix so previous runs are never overwritten.
#
# IMPORTANT:
#   - Do NOT change 5/5 to 4/5 after seeing results.
#   - Do NOT force vent_mode / IMV_24h into the model.
#   - Do NOT pre-impute the full data before bootstrap resampling.
# ============================================================================== 

rm(list = ls())
gc()

# ----------------------------- 0. USER SETTINGS -------------------------------
work_dir <- "D:/R代码/大修1_1"
setwd(work_dir)

EXPECTED_N <- 390L
OUTCOME <- "death_28days"
M <- 5L
MICE_MAXIT <- 20L
K_FOLDS <- 10L
B <- 1000L
CHECKPOINT_EVERY <- 10L

# Set TRUE only for a short technical test. Test outputs use a different prefix.
TEST_MODE <- FALSE
if (TEST_MODE) B <- 5L

PREFIX <- if (TEST_MODE) "390DD_BOOTTEST5_V3_" else "390DD_BOOT1000_V3_"
CHECKPOINT_FILE <- file.path(work_dir, paste0(PREFIX, "checkpoint.rds"))

# Deterministic seed bases (do not change midway through a run)
SEED_APPARENT <- 202608290L
SEED_BOOT_BASE <- 202609000L
SEED_MICE_BASE <- 202610000L
SEED_CV_BASE <- 202611000L

# ----------------------------- 1. PACKAGES ------------------------------------
required_pkgs <- c("readxl", "mice", "glmnet")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Please install required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readxl)
  library(mice)
  library(glmnet)
})

# ----------------------------- 2. HELPERS -------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

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
    # Prefer names that do not look like generated outputs
    cand <- ff390[!grepl("BOOT|output|result", basename(ff390), ignore.case = TRUE)]
    if (length(cand) == 1L) return(cand)
    stop("More than one Excel file beginning with '390' was found. Keep only the analysis workbook in work_dir or set the file manually.\n",
         paste(basename(ff390), collapse = "\n"))
  }
  stop("No .xlsx analysis workbook beginning with '390' was found in: ", dir)
}

extract_candidate_pool <- function(sheet2, data_names) {
  nm <- names(sheet2)
  nm_low <- tolower(trimws(nm))

  include_idx <- which(
    grepl("是否纳入|纳入|include|included", nm_low, perl = TRUE)
  )
  if (length(include_idx) == 0L) {
    stop("Could not identify the inclusion-indicator column in Sheet2.")
  }
  include_idx <- include_idx[1L]

  # Identify the variable-name column by overlap with actual Sheet1 column names.
  overlap_n <- vapply(seq_along(sheet2), function(j) {
    vals <- normalize_chr(sheet2[[j]])
    sum(vals %in% data_names, na.rm = TRUE)
  }, numeric(1))
  overlap_n[include_idx] <- -Inf
  var_idx <- which.max(overlap_n)
  if (!is.finite(overlap_n[var_idx]) || overlap_n[var_idx] <= 0) {
    # Fallback to a likely heading.
    hits <- which(grepl("变量|variable|predictor|字段|field", nm_low, perl = TRUE))
    hits <- setdiff(hits, include_idx)
    if (length(hits) == 0L) stop("Could not identify the variable-name column in Sheet2.")
    var_idx <- hits[1L]
  }

  vars <- normalize_chr(sheet2[[var_idx]])
  keep <- is_yes(sheet2[[include_idx]])
  vars <- unique(vars[keep & !is.na(vars)])
  vars <- vars[vars %in% data_names]
  vars <- setdiff(vars, OUTCOME)

  if (length(vars) == 0L) stop("Candidate pool extracted from Sheet2 is empty.")
  vars
}

coerce_numeric_predictors <- function(dat, vars) {
  # readxl may read an otherwise numeric column as character when Excel cells
  # contain text such as "NA". Treat only explicit missing-value tokens as NA;
  # do not silently coerce genuinely non-numeric labels.
  missing_tokens <- c("", "na", "n/a", "n.a.", "null", "none", "missing", ".", "-")
  true_tokens <- c("1", "yes", "y", "true", "是")
  false_tokens <- c("0", "no", "n", "false", "否")

  audit <- vector("list", length(vars))
  names(audit) <- vars

  for (v in vars) {
    x0 <- dat[[v]]
    original_class <- paste(class(x0), collapse = "/")

    if (is.logical(x0)) x0 <- as.integer(x0)
    if (is.factor(x0)) x0 <- as.character(x0)

    n_token_to_na <- 0L
    if (is.character(x0)) {
      raw_chr <- trimws(x0)
      z <- tolower(raw_chr)

      token_missing <- is.na(raw_chr) | z %in% missing_tokens
      n_token_to_na <- sum(!is.na(raw_chr) & z %in% missing_tokens)
      raw_chr[token_missing] <- NA_character_
      z[token_missing] <- NA_character_

      vals <- unique(z[!is.na(z)])

      # Binary text coding, if and only if every observed token is recognizable.
      if (length(vals) > 0L && all(vals %in% c(true_tokens, false_tokens))) {
        xx <- rep(NA_real_, length(z))
        xx[z %in% true_tokens] <- 1
        xx[z %in% false_tokens] <- 0
        x <- xx
      } else {
        # Remove commas used as thousands separators, then require all remaining
        # observed cells to be numeric. This safely converts values such as
        # "72.5" while still stopping on unexpected category labels.
        numeric_chr <- gsub(",", "", raw_chr, fixed = TRUE)
        suppressWarnings(xx <- as.numeric(numeric_chr))
        bad <- !is.na(raw_chr) & is.na(xx)
        if (any(bad)) {
          ex <- unique(raw_chr[bad])
          stop(
            "Candidate predictor '", v,
            "' contains non-numeric/non-binary values after explicit missing tokens were cleaned. Examples: ",
            paste(head(ex, 5L), collapse = ", "),
            ". Please correct those cells in Excel rather than allowing silent recoding."
          )
        }
        x <- xx
      }
    } else {
      x <- suppressWarnings(as.numeric(x0))
    }

    x[!is.finite(x)] <- NA_real_
    dat[[v]] <- x

    audit[[v]] <- data.frame(
      variable = v,
      original_class = original_class,
      n = length(x),
      n_missing_after_cleaning = sum(is.na(x)),
      missing_percent = 100 * mean(is.na(x)),
      n_text_missing_tokens_cleaned = n_token_to_na,
      n_unique_nonmissing = length(unique(x[!is.na(x)])),
      min = if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE),
      max = if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  attr(dat, "coercion_audit") <- do.call(rbind, audit)
  dat
}

make_stratified_foldid <- function(y, k = 10L, seed = 1L) {
  set.seed(seed)
  y <- as.integer(y)
  if (length(unique(y)) < 2L) stop("Outcome has fewer than two classes.")
  foldid <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
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
  ok <- is.finite(p) & !is.na(y)
  if (!any(ok)) return(NA_real_)
  mean((as.numeric(y[ok]) - as.numeric(p[ok]))^2)
}

calibration_metrics <- function(y, p) {
  ok <- !is.na(y) & is.finite(p)
  y <- as.integer(y[ok]); p <- as.numeric(p[ok])
  if (length(y) < 10L || length(unique(y)) < 2L) {
    return(c(cal_intercept = NA_real_, cal_slope = NA_real_))
  }
  eps <- 1e-6
  p <- pmin(pmax(p, eps), 1 - eps)
  lp <- qlogis(p)

  cint <- tryCatch({
    fit_i <- suppressWarnings(glm(y ~ 1, family = binomial(), offset = lp))
    as.numeric(coef(fit_i)[1L])
  }, error = function(e) NA_real_)

  cslope <- tryCatch({
    fit_s <- suppressWarnings(glm(y ~ lp, family = binomial()))
    as.numeric(coef(fit_s)["lp"])
  }, error = function(e) NA_real_)

  if (!is.finite(cint)) cint <- NA_real_
  if (!is.finite(cslope)) cslope <- NA_real_
  c(cal_intercept = cint, cal_slope = cslope)
}

performance_metrics <- function(y, p) {
  cm <- calibration_metrics(y, p)
  c(
    AUC = auc_rank(y, p),
    Brier = brier_score(y, p),
    cal_intercept = unname(cm["cal_intercept"]),
    cal_slope = unname(cm["cal_slope"])
  )
}

safe_logistic_predict <- function(train_df, new_df, outcome, vars) {
  y <- as.integer(train_df[[outcome]])
  if (length(vars) == 0L) {
    pr <- mean(y == 1L)
    return(list(
      fit = NULL,
      pred_train = rep(pr, nrow(train_df)),
      pred_new = rep(pr, nrow(new_df)),
      coefficients = c(`(Intercept)` = qlogis(pmin(pmax(pr, 1e-6), 1 - 1e-6)))
    ))
  }

  f <- reformulate(vars, response = outcome)
  fit <- suppressWarnings(glm(f, data = train_df, family = binomial()))
  if (!isTRUE(fit$converged)) warning("glm did not report convergence")
  p_train <- suppressWarnings(predict(fit, newdata = train_df, type = "response"))
  p_new <- suppressWarnings(predict(fit, newdata = new_df, type = "response"))
  if (any(!is.finite(p_train)) || any(!is.finite(p_new))) {
    stop("Non-finite logistic predictions encountered.")
  }
  list(fit = fit, pred_train = p_train, pred_new = p_new, coefficients = coef(fit))
}

build_mice_setup <- function(dat, candidate_vars) {
  # PMM for all incompletely observed numeric candidate predictors.
  meth <- rep("", ncol(dat)); names(meth) <- names(dat)
  for (v in candidate_vars) {
    if (anyNA(dat[[v]])) meth[v] <- "pmm"
  }
  meth[OUTCOME] <- ""

  pred <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat),
                 dimnames = list(names(dat), names(dat)))
  # Each imputed candidate can use the other candidate predictors.
  for (v in candidate_vars) {
    if (meth[v] != "") pred[v, candidate_vars] <- 1L
  }
  diag(pred) <- 0L

  # Outcome is never used to impute predictors and is never itself imputed.
  pred[, OUTCOME] <- 0L
  pred[OUTCOME, ] <- 0L

  # Avoid mutual use of PT and INR if both are present.
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
  imp <- suppressWarnings(
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
  imp
}

extract_mice_events <- function(imp, replicate_id, stage) {
  ev <- imp$loggedEvents
  if (is.null(ev) || nrow(ev) == 0L) return(NULL)
  ev <- as.data.frame(ev, stringsAsFactors = FALSE)
  ev$bootstrap <- replicate_id
  ev$stage <- stage
  # Put identifiers first.
  ev[, c("bootstrap", "stage", setdiff(names(ev), c("bootstrap", "stage"))), drop = FALSE]
}

lasso_select_5of5 <- function(completed_train_list, candidate_vars, seed_base) {
  selected_by_imp <- vector("list", M)
  lambda_1se <- rep(NA_real_, M)

  for (i in seq_len(M)) {
    d <- completed_train_list[[i]]
    y <- as.integer(d[[OUTCOME]])

    # glmnet requires finite numeric matrix. Variables with zero variance in this
    # imputation are not eligible for selection in that imputation.
    usable <- candidate_vars[vapply(candidate_vars, function(v) {
      x <- d[[v]]
      all(is.finite(x)) && stats::sd(x) > 0
    }, logical(1))]

    if (length(usable) == 0L) {
      selected_by_imp[[i]] <- character(0)
      next
    }

    x <- as.matrix(d[, usable, drop = FALSE])
    storage.mode(x) <- "double"
    foldid <- make_stratified_foldid(y, k = K_FOLDS, seed = seed_base + i)

    cvfit <- cv.glmnet(
      x = x,
      y = y,
      family = "binomial",
      alpha = 1,
      foldid = foldid,
      standardize = TRUE,
      type.measure = "deviance"
    )
    lambda_1se[i] <- cvfit$lambda.1se
    cf <- as.matrix(coef(cvfit, s = "lambda.1se"))
    sel <- rownames(cf)[cf[, 1] != 0]
    sel <- setdiff(sel, "(Intercept)")
    selected_by_imp[[i]] <- intersect(sel, candidate_vars)
  }

  counts <- setNames(integer(length(candidate_vars)), candidate_vars)
  for (s in selected_by_imp) counts[s] <- counts[s] + 1L
  final_vars <- names(counts)[counts == M]

  list(
    selected_by_imp = selected_by_imp,
    counts = counts,
    final_vars = final_vars,
    lambda_1se = lambda_1se
  )
}

fit_and_average_predictions <- function(train_list, test_list, final_vars) {
  p_train_mat <- matrix(NA_real_, nrow = nrow(train_list[[1L]]), ncol = M)
  p_test_mat  <- matrix(NA_real_, nrow = nrow(test_list[[1L]]),  ncol = M)
  coef_rows <- vector("list", M)

  for (i in seq_len(M)) {
    z <- safe_logistic_predict(train_list[[i]], test_list[[i]], OUTCOME, final_vars)
    p_train_mat[, i] <- z$pred_train
    p_test_mat[, i] <- z$pred_new
    coef_rows[[i]] <- data.frame(
      imputation = i,
      term = names(z$coefficients),
      estimate = as.numeric(z$coefficients),
      stringsAsFactors = FALSE
    )
  }

  list(
    pred_train = rowMeans(p_train_mat),
    pred_test = rowMeans(p_test_mat),
    coefficients = do.call(rbind, coef_rows)
  )
}

run_apparent_pipeline <- function(dat, candidate_vars) {
  imp <- run_mice(dat, candidate_vars, seed = SEED_APPARENT)
  train_list <- lapply(seq_len(M), function(i) complete(imp, i))
  sel <- lasso_select_5of5(train_list, candidate_vars, seed_base = SEED_CV_BASE)
  fit <- fit_and_average_predictions(train_list, train_list, sel$final_vars)
  perf <- performance_metrics(dat[[OUTCOME]], fit$pred_train)
  list(
    performance = perf,
    final_vars = sel$final_vars,
    selection_counts = sel$counts,
    lambda_1se = sel$lambda_1se,
    mice_events = extract_mice_events(imp, 0L, "apparent"),
    coefficients = fit$coefficients
  )
}

run_one_bootstrap <- function(b, dat, candidate_vars) {
  n <- nrow(dat)

  # 1) Bootstrap resampling of original rows.
  set.seed(SEED_BOOT_BASE + b)
  idx <- sample.int(n, size = n, replace = TRUE)
  boot_raw <- dat[idx, , drop = FALSE]

  if (length(unique(boot_raw[[OUTCOME]])) < 2L) {
    stop("Bootstrap sample contains only one outcome class.")
  }

  # 2) Combine bootstrap training rows + original assessment rows.
  #    Assessment rows are ignored when estimating the imputation models, but
  #    are still imputed by the bootstrap-trained MICE process.
  combined <- rbind(boot_raw, dat)
  ignore <- c(rep(FALSE, n), rep(TRUE, n))

  imp <- run_mice(combined, candidate_vars,
                  seed = SEED_MICE_BASE + b,
                  ignore = ignore)
  mice_events <- extract_mice_events(imp, b, "bootstrap_train_plus_original_assessment")

  train_list <- vector("list", M)
  test_list <- vector("list", M)
  for (i in seq_len(M)) {
    comp <- complete(imp, i)
    train_list[[i]] <- comp[seq_len(n), , drop = FALSE]
    test_list[[i]]  <- comp[n + seq_len(n), , drop = FALSE]
  }

  # 3) Entire selection procedure is repeated inside this bootstrap sample.
  sel <- lasso_select_5of5(train_list, candidate_vars,
                           seed_base = SEED_CV_BASE + b * 100L)

  # 4) If 5/5 selects zero predictors, DO NOT relax the rule. Use the
  #    intercept-only model as the literal output of the prespecified pipeline.
  final_vars <- sel$final_vars
  zero_predictor <- length(final_vars) == 0L

  # 5) Fit the selected logistic model separately in each imputation, average
  #    probabilities across imputations, then evaluate in bootstrap and original.
  fit <- fit_and_average_predictions(train_list, test_list, final_vars)
  perf_boot <- performance_metrics(boot_raw[[OUTCOME]], fit$pred_train)
  perf_orig <- performance_metrics(dat[[OUTCOME]], fit$pred_test)
  optimism <- perf_boot - perf_orig

  list(
    performance = data.frame(
      bootstrap = b,
      n_unique = length(unique(idx)),
      events_boot = sum(boot_raw[[OUTCOME]] == 1L),
      model_size = length(final_vars),
      zero_predictor = zero_predictor,
      AUC_boot = unname(perf_boot["AUC"]),
      AUC_original = unname(perf_orig["AUC"]),
      optimism_AUC = unname(optimism["AUC"]),
      Brier_boot = unname(perf_boot["Brier"]),
      Brier_original = unname(perf_orig["Brier"]),
      optimism_Brier = unname(optimism["Brier"]),
      cal_intercept_boot = unname(perf_boot["cal_intercept"]),
      cal_intercept_original = unname(perf_orig["cal_intercept"]),
      optimism_cal_intercept = unname(optimism["cal_intercept"]),
      cal_slope_boot = unname(perf_boot["cal_slope"]),
      cal_slope_original = unname(perf_orig["cal_slope"]),
      optimism_cal_slope = unname(optimism["cal_slope"]),
      stringsAsFactors = FALSE
    ),
    final_vars = final_vars,
    selection_counts = sel$counts,
    lambda_1se = sel$lambda_1se,
    mice_events = mice_events,
    coefficients = fit$coefficients
  )
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 6))
}

# ----------------------------- 3. LOAD DATA -----------------------------------
excel_file <- find_excel_file(work_dir)
message("Analysis workbook: ", basename(excel_file))

sheet_names <- excel_sheets(excel_file)
if (length(sheet_names) < 2L) stop("Workbook must contain at least Sheet1 and Sheet2.")

raw <- as.data.frame(read_excel(excel_file, sheet = 1), check.names = FALSE)
sheet2 <- as.data.frame(read_excel(excel_file, sheet = 2), check.names = FALSE)

if (nrow(raw) != EXPECTED_N) {
  stop("Expected exactly ", EXPECTED_N, " rows in Sheet1, but found ", nrow(raw), ".")
}
if (!OUTCOME %in% names(raw)) stop("Outcome column not found: ", OUTCOME)

candidate_vars <- extract_candidate_pool(sheet2, names(raw))
message("Candidate pool size from Sheet2: ", length(candidate_vars))
message("Candidates: ", paste(candidate_vars, collapse = ", "))

dat <- raw[, c(OUTCOME, candidate_vars), drop = FALSE]
dat <- coerce_numeric_predictors(dat, candidate_vars)
coercion_audit <- attr(dat, "coercion_audit")
if (!is.null(coercion_audit)) {
  write.csv(coercion_audit, paste0(PREFIX, "PRECHECK_predictor_coercion_audit.csv"), row.names = FALSE)
  cleaned <- coercion_audit[coercion_audit$n_text_missing_tokens_cleaned > 0L, , drop = FALSE]
  if (nrow(cleaned) > 0L) {
    message("Explicit text missing-value tokens converted to NA in: ",
            paste(paste0(cleaned$variable, " (", cleaned$n_text_missing_tokens_cleaned, ")"), collapse = ", "))
  }
}

y <- dat[[OUTCOME]]
if (is.factor(y)) y <- as.character(y)
if (is.logical(y)) y <- as.integer(y)
if (is.character(y)) {
  y_chr <- trimws(tolower(y))
  if (any(is.na(y_chr) | y_chr %in% c("", "na", "n/a", "null", ".", "-"))) {
    stop("Outcome contains missing values; outcome must be complete and coded 0/1.")
  }
  if (all(y_chr %in% c("0", "1"))) {
    y_num <- as.numeric(y_chr)
  } else {
    stop("Outcome contains non-0/1 values. Examples: ", paste(head(unique(y[!y_chr %in% c("0", "1")]), 5L), collapse = ", "))
  }
} else {
  suppressWarnings(y_num <- as.numeric(y))
}
if (anyNA(y_num)) stop("Outcome contains missing/non-numeric values.")
if (!all(y_num %in% c(0, 1))) stop("Outcome must be coded 0/1.")
dat[[OUTCOME]] <- as.integer(y_num)

if (anyNA(dat[[OUTCOME]])) stop("Outcome has missing values.")
if (length(unique(dat[[OUTCOME]])) != 2L) stop("Outcome does not contain both classes.")

# Candidate sanity checks before any expensive computation.
if (anyDuplicated(candidate_vars)) stop("Duplicate candidate variable names were extracted from Sheet2.")
if (anyDuplicated(names(dat))) stop("Duplicate column names exist in the analysis data.")

message("Outcome counts: death=", sum(dat[[OUTCOME]] == 1L),
        ", survival=", sum(dat[[OUTCOME]] == 0L))

# Candidate sanity check: no predictor may be completely missing.
all_missing <- candidate_vars[vapply(candidate_vars, function(v) all(is.na(dat[[v]])), logical(1))]
if (length(all_missing) > 0L) stop("Completely missing candidate(s): ", paste(all_missing, collapse = ", "))

zero_var_observed <- candidate_vars[vapply(candidate_vars, function(v) {
  z <- dat[[v]][!is.na(dat[[v]])]
  length(z) > 0L && length(unique(z)) < 2L
}, logical(1))]
if (length(zero_var_observed) > 0L) {
  stop("Candidate(s) have zero variance in the original 390-patient dataset: ",
       paste(zero_var_observed, collapse = ", "),
       ". Remove/correct these in Sheet2 before running the bootstrap.")
}

message("Preflight data checks passed. Starting apparent pipeline...")

# ----------------------------- 4. APPARENT PIPELINE ---------------------------
message("\nRunning apparent development pipeline on the original 390 patients...")
app <- run_apparent_pipeline(dat, candidate_vars)
message("Apparent final variables (5/5): ",
        if (length(app$final_vars)) paste(app$final_vars, collapse = ", ") else "<none>")
message(sprintf("Apparent AUC = %.4f; Brier = %.4f; calibration intercept = %.4f; slope = %.4f",
                app$performance["AUC"], app$performance["Brier"],
                app$performance["cal_intercept"], app$performance["cal_slope"]))

write.csv(data.frame(metric = names(app$performance), value = as.numeric(app$performance)),
          paste0(PREFIX, "00_apparent_performance.csv"), row.names = FALSE)
write.csv(data.frame(variable = candidate_vars,
                     selected_n_of_5 = as.integer(app$selection_counts[candidate_vars]),
                     selected_5of5 = app$selection_counts[candidate_vars] == M),
          paste0(PREFIX, "00_apparent_selection.csv"), row.names = FALSE)

# ----------------------------- 5. INIT / RESUME -------------------------------
if (file.exists(CHECKPOINT_FILE)) {
  ck <- readRDS(CHECKPOINT_FILE)
  if (!identical(ck$B, B) || !identical(ck$EXPECTED_N, EXPECTED_N) ||
      !identical(ck$candidate_vars, candidate_vars)) {
    stop("Checkpoint settings do not match this run. Rename/delete the checkpoint before restarting.")
  }
  results <- ck$results
  failures <- ck$failures
  start_b <- ck$last_completed + 1L
  message("Resuming from checkpoint after bootstrap ", ck$last_completed, ".")
} else {
  results <- vector("list", B)
  failures <- vector("list", B)
  start_b <- 1L
}

if (start_b <= B) {
  for (b in seq.int(start_b, B)) {
    t0 <- proc.time()[3L]
    ans <- tryCatch(
      run_one_bootstrap(b, dat, candidate_vars),
      error = function(e) e
    )

    if (inherits(ans, "error")) {
      failures[[b]] <- data.frame(
        bootstrap = b,
        error = conditionMessage(ans),
        stringsAsFactors = FALSE
      )
      message(sprintf("[%d/%d] ERROR: %s", b, B, conditionMessage(ans)))
    } else {
      results[[b]] <- ans
      elapsed <- proc.time()[3L] - t0
      message(sprintf("[%d/%d] OK | model size=%d | zero=%s | AUCboot=%.3f | AUCorig=%.3f | %.1fs",
                      b, B, ans$performance$model_size,
                      ans$performance$zero_predictor,
                      ans$performance$AUC_boot,
                      ans$performance$AUC_original,
                      elapsed))
    }

    if (b %% CHECKPOINT_EVERY == 0L || b == B) {
      saveRDS(list(
        version = "V3",
        B = B,
        EXPECTED_N = EXPECTED_N,
        candidate_vars = candidate_vars,
        last_completed = b,
        results = results,
        failures = failures
      ), CHECKPOINT_FILE)
    }
  }
}

# ----------------------------- 6. COLLATE -------------------------------------
success_idx <- which(vapply(results, function(x) !is.null(x), logical(1)))
failure_idx <- which(vapply(failures, function(x) !is.null(x), logical(1)))

message("\nSuccessful replicates: ", length(success_idx), "/", B)
message("Hard failures: ", length(failure_idx), "/", B)

if (length(success_idx) == 0L) stop("No successful bootstrap replicates.")

perf_df <- do.call(rbind, lapply(results[success_idx], `[[`, "performance"))
write.csv(perf_df, paste0(PREFIX, "01_replicate_performance.csv"), row.names = FALSE)

# ---- Optimism-corrected performance summary ----
metric_map <- list(
  AUC = c("optimism_AUC", "AUC"),
  Brier = c("optimism_Brier", "Brier"),
  cal_intercept = c("optimism_cal_intercept", "cal_intercept"),
  cal_slope = c("optimism_cal_slope", "cal_slope")
)

summary_rows <- lapply(names(metric_map), function(metric) {
  opt_col <- metric_map[[metric]][1L]
  app_name <- metric_map[[metric]][2L]
  optimism <- perf_df[[opt_col]]
  finite <- is.finite(optimism)
  optimism_f <- optimism[finite]
  apparent <- as.numeric(app$performance[app_name])
  mean_opt <- if (length(optimism_f)) mean(optimism_f) else NA_real_
  corrected <- apparent - mean_opt

  # Distribution of replicate-specific corrected values (descriptive bootstrap
  # distribution, not labelled as a formal confidence interval).
  corr_dist <- apparent - optimism_f
  q_opt <- safe_quantile(optimism_f, c(0.025, 0.975))
  q_corr <- safe_quantile(corr_dist, c(0.025, 0.975))

  data.frame(
    metric = metric,
    apparent = apparent,
    n_with_defined_optimism = length(optimism_f),
    mean_optimism = mean_opt,
    optimism_corrected = corrected,
    optimism_p2.5 = q_opt[1L],
    optimism_p97.5 = q_opt[2L],
    corrected_distribution_p2.5 = q_corr[1L],
    corrected_distribution_p97.5 = q_corr[2L],
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, paste0(PREFIX, "02_optimism_corrected_summary.csv"), row.names = FALSE)

# ---- Complete-pipeline selection frequency ----
sel_mat <- matrix(0L, nrow = length(success_idx), ncol = length(candidate_vars),
                  dimnames = list(NULL, candidate_vars))
for (j in seq_along(success_idx)) {
  vv <- results[[success_idx[j]]]$final_vars
  sel_mat[j, vv] <- 1L
}
sel_freq <- data.frame(
  variable = candidate_vars,
  selected_n = colSums(sel_mat),
  denominator_successful = length(success_idx),
  selection_percent = 100 * colMeans(sel_mat),
  stringsAsFactors = FALSE
)
sel_freq <- sel_freq[order(-sel_freq$selection_percent, sel_freq$variable), ]
write.csv(sel_freq, paste0(PREFIX, "03_selection_frequency.csv"), row.names = FALSE)

# ---- Model size and zero-predictor frequency ----
size_tab <- as.data.frame(table(perf_df$model_size), stringsAsFactors = FALSE)
names(size_tab) <- c("model_size", "n")
size_tab$percent <- 100 * size_tab$n / nrow(perf_df)
write.csv(size_tab, paste0(PREFIX, "04_model_size_distribution.csv"), row.names = FALSE)

# ---- Coefficients across all MI-specific fitted models ----
coef_long <- do.call(rbind, lapply(success_idx, function(b) {
  z <- results[[b]]$coefficients
  z$bootstrap <- b
  z[, c("bootstrap", "imputation", "term", "estimate")]
}))
write.csv(coef_long, paste0(PREFIX, "05_bootstrap_coefficients_long.csv"), row.names = FALSE)

coef_stability <- do.call(rbind, lapply(split(coef_long, coef_long$term), function(z) {
  est <- z$estimate[is.finite(z$estimate)]
  data.frame(
    term = z$term[1L],
    n_MI_fits = length(est),
    mean = if (length(est)) mean(est) else NA_real_,
    sd = if (length(est) > 1L) sd(est) else NA_real_,
    median = if (length(est)) median(est) else NA_real_,
    p2.5 = safe_quantile(est, 0.025),
    p97.5 = safe_quantile(est, 0.975),
    stringsAsFactors = FALSE
  )
}))
coef_stability <- coef_stability[order(coef_stability$term), ]
write.csv(coef_stability, paste0(PREFIX, "06_coefficient_stability.csv"), row.names = FALSE)

# ---- Hard failures ----
if (length(failure_idx)) {
  fail_df <- do.call(rbind, failures[failure_idx])
} else {
  fail_df <- data.frame(bootstrap = integer(0), error = character(0))
}
write.csv(fail_df, paste0(PREFIX, "07_failures.csv"), row.names = FALSE)

# ---- MICE loggedEvents audit: events are NOT failures ----
event_list <- list()
if (!is.null(app$mice_events)) event_list[[length(event_list) + 1L]] <- app$mice_events
for (b in success_idx) {
  ev <- results[[b]]$mice_events
  if (!is.null(ev) && nrow(ev) > 0L) event_list[[length(event_list) + 1L]] <- ev
}
if (length(event_list)) {
  mice_events_df <- do.call(rbind, event_list)
} else {
  mice_events_df <- data.frame(
    bootstrap = integer(0), stage = character(0),
    stringsAsFactors = FALSE
  )
}
write.csv(mice_events_df, paste0(PREFIX, "08_mice_loggedEvents_audit.csv"), row.names = FALSE)

# Compact per-replicate event counts
if (nrow(mice_events_df) > 0L && "bootstrap" %in% names(mice_events_df)) {
  tmp <- aggregate(rep(1L, nrow(mice_events_df)),
                   by = list(bootstrap = mice_events_df$bootstrap), FUN = sum)
  names(tmp)[2L] <- "n_loggedEvents"
  event_counts <- merge(data.frame(bootstrap = 0:B), tmp, by = "bootstrap", all.x = TRUE)
  event_counts$n_loggedEvents[is.na(event_counts$n_loggedEvents)] <- 0L
} else {
  event_counts <- data.frame(bootstrap = 0:B, n_loggedEvents = 0L)
}
write.csv(event_counts, paste0(PREFIX, "09_mice_loggedEvents_counts.csv"), row.names = FALSE)

# ---- Overall run status ----
status <- data.frame(
  item = c(
    "script_version", "B_requested", "B_successful", "B_hard_failures",
    "zero_predictor_replicates", "zero_predictor_percent",
    "replicates_with_defined_AUC_optimism",
    "replicates_with_defined_Brier_optimism",
    "replicates_with_defined_calibration_intercept_optimism",
    "replicates_with_defined_calibration_slope_optimism",
    "candidate_pool_size", "apparent_model_size"
  ),
  value = c(
    "V3", B, length(success_idx), length(failure_idx),
    sum(perf_df$zero_predictor), 100 * mean(perf_df$zero_predictor),
    sum(is.finite(perf_df$optimism_AUC)),
    sum(is.finite(perf_df$optimism_Brier)),
    sum(is.finite(perf_df$optimism_cal_intercept)),
    sum(is.finite(perf_df$optimism_cal_slope)),
    length(candidate_vars), length(app$final_vars)
  ),
  stringsAsFactors = FALSE
)
write.csv(status, paste0(PREFIX, "10_run_status.csv"), row.names = FALSE)

# ---- Save complete object and session information ----
saveRDS(list(
  settings = list(
    EXPECTED_N = EXPECTED_N, OUTCOME = OUTCOME, M = M, MICE_MAXIT = MICE_MAXIT,
    K_FOLDS = K_FOLDS, B = B, candidate_vars = candidate_vars,
    seed_apparent = SEED_APPARENT, seed_boot_base = SEED_BOOT_BASE,
    seed_mice_base = SEED_MICE_BASE, seed_cv_base = SEED_CV_BASE
  ),
  apparent = app,
  results = results,
  failures = failures,
  performance = perf_df,
  summary = summary_df,
  selection_frequency = sel_freq,
  mice_loggedEvents = mice_events_df,
  status = status
), paste0(PREFIX, "11_complete_bootstrap_object.rds"))

sink(paste0(PREFIX, "12_sessionInfo.txt"))
print(sessionInfo())
sink()

# ----------------------------- 7. CONSOLE SUMMARY -----------------------------
cat("\n============================================================\n")
cat("COMPLETE-PIPELINE BOOTSTRAP V3 FINISHED\n")
cat("============================================================\n")
cat("Successful replicates:", length(success_idx), "/", B, "\n")
cat("Hard failures:", length(failure_idx), "\n")
cat("Zero-predictor (intercept-only) replicates:", sum(perf_df$zero_predictor),
    sprintf("(%.2f%%)\n", 100 * mean(perf_df$zero_predictor)))
cat("\nApparent final variables:\n  ",
    if (length(app$final_vars)) paste(app$final_vars, collapse = ", ") else "<none>", "\n", sep = "")
cat("\nOptimism-corrected summary:\n")
print(summary_df, row.names = FALSE)
cat("\nTop selection frequencies:\n")
print(utils::head(sel_freq, 15L), row.names = FALSE)
cat("\nCalibration optimism defined in:\n")
cat("  intercept:", sum(is.finite(perf_df$optimism_cal_intercept)), "/", nrow(perf_df), "\n")
cat("  slope:    ", sum(is.finite(perf_df$optimism_cal_slope)), "/", nrow(perf_df), "\n")
cat("\nMICE loggedEvents are stored for audit and were NOT automatically treated as failures.\n")
cat("Outputs written with prefix:", PREFIX, "\n")
cat("============================================================\n")
