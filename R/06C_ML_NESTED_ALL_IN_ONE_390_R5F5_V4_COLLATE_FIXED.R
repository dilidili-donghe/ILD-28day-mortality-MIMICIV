# ======================================================================
# 06C_ML_NESTED_ALL_IN_ONE_390.R
#
# ONE-SCRIPT, RESUMABLE, LEAKAGE-CONTROLLED NESTED ML COMPARISON
#
# Study cohort:
#   Original, UNIMPUTED N=390 cohort.
# Candidate pool:
#   The same prespecified 44 candidate predictors for every formally compared
#   algorithm.
#
# Formal algorithms:
#   1) Primary_LASSO_5of5
#      - alpha = 1
#      - lambda chosen by inner CV using a 1-SE rule separately within each
#        of the 5 imputation streams
#      - outer-training LASSO selection within each imputation
#      - predictor retained only if selected in all 5 imputations
#      - unpenalized logistic refit on the retained predictors
#   2) ElasticNet
#   3) RandomForest
#   4) XGBoost
#   5) SVM_RBF
#
# Validation:
#   - Frozen outer CV: 5 repeats x 5 folds (25 outer test folds)
#   - Frozen inner CV: 5 folds within every outer training set
#   - All algorithms use the SAME frozen outer and inner splits
#   - Imputation/preprocessing is learned from training rows only
#   - Inner validation rows NEVER contribute to their imputation models
#   - Outer test rows NEVER contribute to their imputation models
#   - Outcome is NOT used as a predictor in MICE, preventing test-outcome
#     leakage into predictor imputation.
#
# Missing data:
#   - MICE m=5, PMM, maxit=20
#   - NEW fold-specific imputations are REQUIRED for nested CV.
#   - Previously generated full-cohort / bootstrap imputations CANNOT be
#     reused because they were estimated under different training samples.
#
# Resumability / saving:
#   - Every inner and outer imputation cache is saved immediately.
#   - Every completed outer model result is saved immediately.
#   - Existing valid cache/result files are automatically reused.
#   - Patient-level OOF predictions, split IDs, seeds, tuning tables,
#     selected predictors, coefficients/importances, timing, failures,
#     complete RDS object, figures, sessionInfo, and manifest are saved.
#
# IMPORTANT:
#   Do NOT delete the frozen split files and do NOT regenerate them.
# ======================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(warn = 1)

# ----------------------------------------------------------------------
# 0. SETTINGS
# ----------------------------------------------------------------------
WORK_DIR <- "D:/R代码/大修1_1"

DATA_FILE <- file.path(WORK_DIR, "390run.xlsx")
OUTER_FILE <- file.path(WORK_DIR, "390ML_04_FROZEN_OUTER_SPLITS_R10_F5.csv")
INNER_FILE <- file.path(WORK_DIR, "390ML_06B_06_FROZEN_INNER_SPLITS_ALL_OUTER_FOLDS.csv")

OUTCOME <- "death_28days"

EXPECTED_N <- 390L
EXPECTED_EVENTS <- 111L
EXPECTED_CANDIDATES <- 44L
N_REPEATS <- 5L
N_OUTER_FOLDS <- 5L
N_INNER_FOLDS <- 5L

M <- 5L
MICE_MAXIT <- 20L

# PRE-SPECIFIED COMPUTATIONAL DESIGN: outer 5-fold CV repeated 5 times; inner 5-fold CV.
# This is a repeated nested-CV design chosen before the formal ML run.
# Deterministic seed bases.
SEED_MICE_INNER_BASE <- 202613000L
SEED_MICE_OUTER_BASE <- 202614000L
SEED_MODEL_BASE <- 202615000L

# Conservative, prespecified tuning grids to control computation.
# glmnet expects lambda in decreasing order.
LAMBDA_GRID <- exp(seq(log(1), log(1e-4), length.out = 25L))
ENET_ALPHA_GRID <- c(0.10, 0.50, 0.90)

RF_GRID <- expand.grid(
  mtry = c(7L, 15L, 30L),
  min_node_size = c(5L, 15L),
  sample_fraction = c(0.80, 1.00),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
RF_TREES <- 500L

XGB_GRID <- expand.grid(
  max_depth = c(2L, 4L),
  eta = c(0.03, 0.10),
  min_child_weight = c(1, 5),
  subsample = 0.80,
  colsample_bytree = 0.80,
  nrounds = c(200L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

SVM_GRID <- expand.grid(
  C = c(0.5, 2, 8),
  sigma = c(0.005, 0.02, 0.08),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

ALGORITHMS <- c(
  "Primary_LASSO_5of5",
  "ElasticNet",
  "RandomForest",
  "XGBoost",
  "SVM_RBF"
)

# CPU threads used INSIDE ranger/xgboost fits.
AVAILABLE_CORES <- parallel::detectCores(logical = TRUE)
N_THREADS <- max(1L, min(4L, AVAILABLE_CORES - 1L))

# Output structure.
ROOT_DIR <- file.path(WORK_DIR, "390ML_06C_R5F5_ALL")
CACHE_INNER_DIR <- file.path(ROOT_DIR, "cache_inner")
CACHE_OUTER_DIR <- file.path(ROOT_DIR, "cache_outer")
RESULT_DIR <- file.path(ROOT_DIR, "outer_results")
TUNING_DIR <- file.path(ROOT_DIR, "tuning_tables")
FIG_DIR <- file.path(ROOT_DIR, "figures")
LOG_DIR <- file.path(ROOT_DIR, "logs")

for (d in c(ROOT_DIR, CACHE_INNER_DIR, CACHE_OUTER_DIR, RESULT_DIR,
            TUNING_DIR, FIG_DIR, LOG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

RUN_LOG <- file.path(LOG_DIR, "390ML_06C_console_log.txt")

# ----------------------------------------------------------------------
# 1. PACKAGES
# ----------------------------------------------------------------------
required_pkgs <- c(
  "readxl", "mice", "glmnet", "ranger", "xgboost", "kernlab", "pROC"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

# ----------------------------------------------------------------------
# 2. LOGGING / GENERAL HELPERS
# ----------------------------------------------------------------------
log_msg <- function(...) {
  txt <- paste0(..., collapse = "")
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", txt)
  cat(line, "\n")
  cat(line, "\n", file = RUN_LOG, append = TRUE)
}

fmt_time <- function(seconds) {
  seconds <- as.numeric(seconds)
  if (!is.finite(seconds)) return(NA_character_)
  h <- floor(seconds / 3600)
  m <- floor((seconds %% 3600) / 60)
  s <- round(seconds %% 60, 1)
  sprintf("%02d:%02d:%04.1f", h, m, s)
}

seed_from <- function(base, repeat_id, outer_fold, inner_fold = 0L,
                      imp = 0L, param = 0L) {
  val <- as.numeric(base) +
    as.numeric(repeat_id) * 100000 +
    as.numeric(outer_fold) * 10000 +
    as.numeric(inner_fold) * 1000 +
    as.numeric(imp) * 100 +
    as.numeric(param)
  as.integer(val %% 2147483000)
}

clip_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

logloss <- function(y, p) {
  p <- clip_prob(p)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

brier <- function(y, p) {
  mean((as.numeric(y) - as.numeric(p))^2)
}

auc_value <- function(y, p) {
  if (length(unique(y)) < 2L) return(NA_real_)
  as.numeric(
    pROC::auc(
      pROC::roc(
        response = y,
        predictor = p,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      )
    )
  )
}

calibration_intercept <- function(y, p) {
  lp <- qlogis(clip_prob(p))
  fit <- try(
    stats::glm(y ~ 1 + offset(lp), family = stats::binomial()),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NA_real_)
  cf <- stats::coef(fit)
  if (length(cf) < 1L || !is.finite(cf[1L])) return(NA_real_)
  unname(cf[1L])
}

calibration_slope <- function(y, p) {
  lp <- qlogis(clip_prob(p))
  fit <- try(
    stats::glm(y ~ lp, family = stats::binomial()),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NA_real_)
  cf <- stats::coef(fit)
  if (length(cf) < 2L || !is.finite(cf[2L])) return(NA_real_)
  unname(cf[2L])
}

metric_row <- function(y, p) {
  data.frame(
    AUC = auc_value(y, p),
    Brier = brier(y, p),
    Cal_Intercept = calibration_intercept(y, p),
    Cal_Slope = calibration_slope(y, p),
    stringsAsFactors = FALSE
  )
}

safe_trim <- function(x) trimws(as.character(x))

coerce_missing_tokens <- function(x) {
  if (is.character(x)) {
    z <- trimws(x)
    z[toupper(z) %in% c("", "NA", "N/A", "NAN", "NULL", ".")] <- NA_character_
    return(z)
  }
  x
}

safe_numeric <- function(x, varname) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  z <- suppressWarnings(as.numeric(as.character(x)))
  bad <- !is.na(x) & is.na(z)
  if (any(bad)) {
    stop(
      "Variable '", varname, "' has non-numeric value(s): ",
      paste(utils::head(unique(as.character(x[bad])), 10L), collapse = ", ")
    )
  }
  z
}

get_meta_columns <- function(meta, dat_names) {
  incl_candidates <- c("是否纳入", "纳入", "include", "included", "Include", "Included")
  incl_col <- incl_candidates[incl_candidates %in% names(meta)][1L]
  if (length(incl_col) == 0L || is.na(incl_col)) {
    stop("Could not identify inclusion column in Sheet2.")
  }

  name_candidates <- c(
    "表格中原始名称", "变量名", "变量", "variable",
    "Variable", "var", "name", "变量名称"
  )
  name_col <- name_candidates[name_candidates %in% names(meta)][1L]

  if (length(name_col) == 0L || is.na(name_col)) {
    overlaps <- vapply(
      meta,
      function(z) sum(safe_trim(z) %in% dat_names, na.rm = TRUE),
      numeric(1)
    )
    if (max(overlaps) == 0L) stop("Could not identify variable-name column.")
    name_col <- names(which.max(overlaps))
  }

  type_candidates <- c("变量类型", "类型", "type", "Type", "variable_type")
  type_col <- type_candidates[type_candidates %in% names(meta)][1L]
  if (length(type_col) == 0L || is.na(type_col)) type_col <- NA_character_

  list(incl_col = incl_col, name_col = name_col, type_col = type_col)
}

# ----------------------------------------------------------------------
# 3. READ / VERIFY RAW DATA AND FROZEN SPLITS
# ----------------------------------------------------------------------
OVERALL_START <- Sys.time()
if (file.exists(RUN_LOG)) file.remove(RUN_LOG)

log_msg("06C ALL-IN-ONE nested ML started.")
log_msg("Threads for ranger/xgboost: ", N_THREADS)

for (f in c(DATA_FILE, OUTER_FILE, INNER_FILE)) {
  if (!file.exists(f)) stop("Required file not found: ", f)
}

dat <- as.data.frame(
  readxl::read_excel(
    DATA_FILE,
    sheet = "Sheet1",
    na = c("", "NA", "N/A", "NaN", "NULL", ".")
  )
)
meta <- as.data.frame(
  readxl::read_excel(
    DATA_FILE,
    sheet = "Sheet2",
    na = c("", "NA", "N/A", "NaN", "NULL", ".")
  )
)
outer_splits <- utils::read.csv(
  OUTER_FILE, stringsAsFactors = FALSE, check.names = FALSE
)
inner_splits <- utils::read.csv(
  INNER_FILE, stringsAsFactors = FALSE, check.names = FALSE
)

for (j in seq_along(dat)) dat[[j]] <- coerce_missing_tokens(dat[[j]])

if (nrow(dat) != EXPECTED_N) stop("N mismatch.")
if (!OUTCOME %in% names(dat)) stop("Outcome missing.")

y <- suppressWarnings(as.integer(as.character(dat[[OUTCOME]])))
if (anyNA(y) || !all(y %in% c(0L, 1L))) stop("Outcome must be complete 0/1.")
if (sum(y == 1L) != EXPECTED_EVENTS) stop("Event count mismatch.")

mc <- get_meta_columns(meta, names(dat))
incl_raw <- toupper(safe_trim(meta[[mc$incl_col]]))
incl_yes <- incl_raw %in% c("是", "YES", "Y", "1", "TRUE", "T")
candidate_vars <- unique(safe_trim(meta[[mc$name_col]][incl_yes]))
candidate_vars <- candidate_vars[!is.na(candidate_vars) & nzchar(candidate_vars)]
candidate_vars <- setdiff(candidate_vars, OUTCOME)

if (length(candidate_vars) != EXPECTED_CANDIDATES) {
  stop("Candidate count = ", length(candidate_vars), "; expected 44.")
}
if (length(setdiff(candidate_vars, names(dat))) > 0L) {
  stop("Some candidate variables are absent from Sheet1.")
}

for (v in candidate_vars) dat[[v]] <- safe_numeric(dat[[v]], v)

# Validate split columns.
needed_outer <- c("row_id", "outcome", "repeat_id", "outer_fold")
needed_inner <- c("row_id", "outcome", "repeat_id", "outer_fold", "inner_fold")
if (length(setdiff(needed_outer, names(outer_splits))) > 0L) {
  stop("Frozen outer split file has missing required columns.")
}
if (length(setdiff(needed_inner, names(inner_splits))) > 0L) {
  stop("Frozen inner split file has missing required columns.")
}

for (v in needed_outer) outer_splits[[v]] <- as.integer(outer_splits[[v]])
for (v in needed_inner) inner_splits[[v]] <- as.integer(inner_splits[[v]])

# ------------------------------------------------------------------
# R5xF5 DESIGN: the frozen split files were originally created for
# 10 repeats. For the prespecified R5xF5 analysis, USE (do not
# regenerate) repeats 1-5 from those frozen files.
# ------------------------------------------------------------------
available_outer_repeats <- sort(unique(outer_splits$repeat_id))
available_inner_repeats <- sort(unique(inner_splits$repeat_id))

if (!all(seq_len(N_REPEATS) %in% available_outer_repeats)) {
  stop(
    "Frozen outer split file does not contain all required repeats 1:",
    N_REPEATS, ". Available repeats: ",
    paste(available_outer_repeats, collapse = ", ")
  )
}
if (!all(seq_len(N_REPEATS) %in% available_inner_repeats)) {
  stop(
    "Frozen inner split file does not contain all required repeats 1:",
    N_REPEATS, ". Available repeats: ",
    paste(available_inner_repeats, collapse = ", ")
  )
}

outer_splits <- outer_splits[
  outer_splits$repeat_id %in% seq_len(N_REPEATS),
  ,
  drop = FALSE
]
inner_splits <- inner_splits[
  inner_splits$repeat_id %in% seq_len(N_REPEATS),
  ,
  drop = FALSE
]

# Hard integrity checks after subsetting the already-frozen R10 files.
if (nrow(outer_splits) != EXPECTED_N * N_REPEATS) {
  stop(
    "Outer split row count mismatch after restricting to repeats 1:",
    N_REPEATS,
    ". Found ", nrow(outer_splits),
    "; expected ", EXPECTED_N * N_REPEATS, "."
  )
}

expected_inner_rows <- sum(
  vapply(
    seq_len(N_REPEATS),
    function(r) {
      sr <- outer_splits[outer_splits$repeat_id == r, , drop = FALSE]
      sum(
        vapply(
          seq_len(N_OUTER_FOLDS),
          function(of) {
            EXPECTED_N - sum(sr$outer_fold == of)
          },
          integer(1)
        )
      )
    },
    integer(1)
  )
)

if (nrow(inner_splits) != expected_inner_rows) {
  stop(
    "Inner split row count mismatch after restricting to repeats 1:",
    N_REPEATS,
    ". Found ", nrow(inner_splits),
    "; expected ", expected_inner_rows, "."
  )
}

if (!all(outer_splits$outcome == y[outer_splits$row_id])) {
  stop("Frozen outer split outcome does not match raw outcome.")
}
if (!all(inner_splits$outcome == y[inner_splits$row_id])) {
  stop("Frozen inner split outcome does not match raw outcome.")
}

log_msg(
  "Frozen R10 split files restricted to prespecified repeats 1-",
  N_REPEATS,
  " without regenerating any split. Outer rows=",
  nrow(outer_splits),
  "; inner rows=",
  nrow(inner_splits),
  "."
)

# ----------------------------------------------------------------------
# 4. MICE HELPERS
# ----------------------------------------------------------------------
build_mice_spec <- function(imp_dat) {
  method <- mice::make.method(imp_dat)
  pred <- mice::make.predictorMatrix(imp_dat)

  # Outcome is observed for training evaluation, but NEVER used to impute X.
  method[OUTCOME] <- ""
  pred[, OUTCOME] <- 0L
  pred[OUTCOME, ] <- 0L

  for (v in candidate_vars) {
    if (all(!is.na(imp_dat[[v]]))) {
      method[v] <- ""
    } else {
      method[v] <- "pmm"
    }
  }

  # Same protection used in the finalized development pipeline.
  if (all(c("pt", "inr") %in% colnames(pred))) {
    pred["pt", "inr"] <- 0L
    pred["inr", "pt"] <- 0L
  }

  list(method = method, predictorMatrix = pred)
}

make_completed_cache <- function(train_ids, eval_ids, seed, label) {
  combined_ids <- c(train_ids, eval_ids)
  ignore <- c(
    rep(FALSE, length(train_ids)),
    rep(TRUE, length(eval_ids))
  )

  imp_dat <- dat[
    combined_ids,
    unique(c(OUTCOME, candidate_vars)),
    drop = FALSE
  ]

  observed_train <- vapply(
    candidate_vars,
    function(v) sum(!is.na(imp_dat[[v]][!ignore])),
    integer(1)
  )
  if (any(observed_train == 0L)) {
    stop(
      label, ": no observed training values for: ",
      paste(names(observed_train)[observed_train == 0L], collapse = ", ")
    )
  }

  spec <- build_mice_spec(imp_dat)

  t0 <- Sys.time()
  imp <- mice::mice(
    data = imp_dat,
    m = M,
    maxit = MICE_MAXIT,
    method = spec$method,
    predictorMatrix = spec$predictorMatrix,
    ignore = ignore,
    seed = seed,
    printFlag = FALSE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  completed <- vector("list", M)
  for (m in seq_len(M)) {
    cc <- as.data.frame(mice::complete(imp, action = m))
    if (anyNA(cc[candidate_vars])) {
      bad <- candidate_vars[
        vapply(cc[candidate_vars], function(z) anyNA(z), logical(1))
      ]
      stop(label, ": residual NA after MICE in ", paste(bad, collapse = ", "))
    }
    completed[[m]] <- list(
      train = cc[!ignore, , drop = FALSE],
      eval = cc[ignore, , drop = FALSE]
    )
  }

  events <- imp$loggedEvents
  if (is.null(events)) events <- data.frame()

  list(
    label = label,
    seed = seed,
    train_ids = train_ids,
    eval_ids = eval_ids,
    method = spec$method,
    predictorMatrix = spec$predictorMatrix,
    completed = completed,
    loggedEvents = events,
    elapsed_seconds = elapsed
  )
}

inner_cache_path <- function(r, of, infold) {
  file.path(
    CACHE_INNER_DIR,
    sprintf("inner_R%02d_O%02d_I%02d.rds", r, of, infold)
  )
}

outer_cache_path <- function(r, of) {
  file.path(
    CACHE_OUTER_DIR,
    sprintf("outer_R%02d_O%02d.rds", r, of)
  )
}

result_path <- function(r, of) {
  file.path(
    RESULT_DIR,
    sprintf("result_R%02d_O%02d.rds", r, of)
  )
}

tuning_path <- function(r, of, model) {
  file.path(
    TUNING_DIR,
    sprintf("tuning_R%02d_O%02d_%s.csv", r, of, model)
  )
}

# ----------------------------------------------------------------------
# 5. BUILD / REUSE ALL NESTED MICE CACHES
# ----------------------------------------------------------------------
CACHE_START <- Sys.time()
cache_timing_rows <- list()
cache_k <- 0L

log_msg("Stage A: building/reusing leakage-free MICE caches.")

for (r in seq_len(N_REPEATS)) {
  split_r <- outer_splits[outer_splits$repeat_id == r, , drop = FALSE]

  for (of in seq_len(N_OUTER_FOLDS)) {
    outer_test_ids <- split_r$row_id[split_r$outer_fold == of]
    outer_train_ids <- setdiff(seq_len(nrow(dat)), outer_test_ids)

    # ---- Inner caches ----
    inner_this <- inner_splits[
      inner_splits$repeat_id == r & inner_splits$outer_fold == of,
      ,
      drop = FALSE
    ]

    for (infold in seq_len(N_INNER_FOLDS)) {
      p <- inner_cache_path(r, of, infold)

      if (file.exists(p)) {
        obj <- readRDS(p)
        status <- "reused"
        elapsed <- obj$elapsed_seconds
      } else {
        inner_val_ids <- inner_this$row_id[inner_this$inner_fold == infold]
        inner_train_ids <- inner_this$row_id[inner_this$inner_fold != infold]
        seed <- seed_from(SEED_MICE_INNER_BASE, r, of, infold)

        log_msg(
          "MICE inner R", r, " O", of, " I", infold,
          " (train=", length(inner_train_ids),
          ", val=", length(inner_val_ids), ")"
        )

        obj <- make_completed_cache(
          train_ids = inner_train_ids,
          eval_ids = inner_val_ids,
          seed = seed,
          label = sprintf("INNER_R%02d_O%02d_I%02d", r, of, infold)
        )
        saveRDS(obj, p, compress = "xz")
        status <- "built"
        elapsed <- obj$elapsed_seconds
      }

      cache_k <- cache_k + 1L
      cache_timing_rows[[cache_k]] <- data.frame(
        cache_type = "inner",
        repeat_id = r,
        outer_fold = of,
        inner_fold = infold,
        status = status,
        elapsed_seconds = elapsed,
        file = p,
        stringsAsFactors = FALSE
      )
    }

    # ---- Outer cache ----
    p <- outer_cache_path(r, of)

    if (file.exists(p)) {
      obj <- readRDS(p)
      status <- "reused"
      elapsed <- obj$elapsed_seconds
    } else {
      seed <- seed_from(SEED_MICE_OUTER_BASE, r, of)

      log_msg(
        "MICE outer R", r, " O", of,
        " (train=", length(outer_train_ids),
        ", test=", length(outer_test_ids), ")"
      )

      obj <- make_completed_cache(
        train_ids = outer_train_ids,
        eval_ids = outer_test_ids,
        seed = seed,
        label = sprintf("OUTER_R%02d_O%02d", r, of)
      )
      saveRDS(obj, p, compress = "xz")
      status <- "built"
      elapsed <- obj$elapsed_seconds
    }

    cache_k <- cache_k + 1L
    cache_timing_rows[[cache_k]] <- data.frame(
      cache_type = "outer",
      repeat_id = r,
      outer_fold = of,
      inner_fold = NA_integer_,
      status = status,
      elapsed_seconds = elapsed,
      file = p,
      stringsAsFactors = FALSE
    )
  }
}

cache_timing <- do.call(rbind, cache_timing_rows)
utils::write.csv(
  cache_timing,
  file.path(ROOT_DIR, "390ML_06C_01_cache_timing_manifest.csv"),
  row.names = FALSE
)

CACHE_ELAPSED <- as.numeric(difftime(Sys.time(), CACHE_START, units = "secs"))
log_msg(
  "Stage A complete. Wall time = ", fmt_time(CACHE_ELAPSED),
  ". Built/reused ", nrow(cache_timing), " cache objects."
)

# ----------------------------------------------------------------------
# 6. MODEL HELPERS
# ----------------------------------------------------------------------
x_matrix <- function(df) {
  x <- as.matrix(df[, candidate_vars, drop = FALSE])
  storage.mode(x) <- "double"
  x
}

y_vector <- function(df) {
  as.integer(df[[OUTCOME]])
}

predict_glmnet_matrix <- function(fit, newx, lambda_values) {
  p <- stats::predict(
    fit,
    newx = newx,
    s = lambda_values,
    type = "response"
  )
  as.matrix(p)
}

choose_1se_lambda <- function(fold_losses, lambda_grid) {
  mean_loss <- colMeans(fold_losses, na.rm = TRUE)
  se_loss <- apply(fold_losses, 2L, stats::sd, na.rm = TRUE) /
    sqrt(rowSums(is.finite(fold_losses)))

  i_min <- which.min(mean_loss)
  cutoff <- mean_loss[i_min] + se_loss[i_min]

  eligible <- which(mean_loss <= cutoff)
  if (length(eligible) == 0L) eligible <- i_min

  # lambda_grid is decreasing. Choose largest lambda satisfying 1-SE.
  i_1se <- eligible[which.max(lambda_grid[eligible])]

  data.frame(
    lambda = lambda_grid,
    mean_logloss = mean_loss,
    se_logloss = se_loss,
    is_min = seq_along(lambda_grid) == i_min,
    is_1se = seq_along(lambda_grid) == i_1se,
    stringsAsFactors = FALSE
  )
}

tune_primary_lasso_5of5 <- function(r, of) {
  tuning_by_m <- vector("list", M)
  chosen_lambda <- numeric(M)

  for (m in seq_len(M)) {
    fold_losses <- matrix(
      NA_real_,
      nrow = N_INNER_FOLDS,
      ncol = length(LAMBDA_GRID)
    )

    for (infold in seq_len(N_INNER_FOLDS)) {
      cc <- readRDS(inner_cache_path(r, of, infold))$completed[[m]]
      xtr <- x_matrix(cc$train)
      ytr <- y_vector(cc$train)
      xva <- x_matrix(cc$eval)
      yva <- y_vector(cc$eval)

      fit <- glmnet::glmnet(
        x = xtr,
        y = ytr,
        family = "binomial",
        alpha = 1,
        lambda = LAMBDA_GRID,
        standardize = TRUE,
        intercept = TRUE
      )

      pred_mat <- predict_glmnet_matrix(fit, xva, LAMBDA_GRID)

      for (j in seq_along(LAMBDA_GRID)) {
        fold_losses[infold, j] <- logloss(yva, pred_mat[, j])
      }
    }

    tab <- choose_1se_lambda(fold_losses, LAMBDA_GRID)
    tab$imputation <- m
    tuning_by_m[[m]] <- tab
    chosen_lambda[m] <- tab$lambda[tab$is_1se][1L]
  }

  list(
    chosen_lambda = chosen_lambda,
    tuning = do.call(rbind, tuning_by_m)
  )
}

fit_primary_outer <- function(r, of, outer_obj, tuned) {
  selected_by_m <- vector("list", M)
  lasso_coefs <- vector("list", M)

  for (m in seq_len(M)) {
    cc <- outer_obj$completed[[m]]
    xtr <- x_matrix(cc$train)
    ytr <- y_vector(cc$train)

    fit <- glmnet::glmnet(
      x = xtr,
      y = ytr,
      family = "binomial",
      alpha = 1,
      lambda = tuned$chosen_lambda[m],
      standardize = TRUE,
      intercept = TRUE
    )

    cf <- as.matrix(stats::coef(fit, s = tuned$chosen_lambda[m]))
    nz <- rownames(cf)[cf[, 1] != 0]
    nz <- setdiff(nz, "(Intercept)")
    selected_by_m[[m]] <- nz

    lasso_coefs[[m]] <- data.frame(
      imputation = m,
      variable = rownames(cf),
      coefficient = cf[, 1],
      lambda = tuned$chosen_lambda[m],
      stringsAsFactors = FALSE
    )
  }

  final_vars <- Reduce(intersect, selected_by_m)
  pred_by_m <- matrix(
    NA_real_,
    nrow = nrow(outer_obj$completed[[1L]]$eval),
    ncol = M
  )
  refit_coefs <- vector("list", M)

  for (m in seq_len(M)) {
    cc <- outer_obj$completed[[m]]
    train_df <- cc$train
    test_df <- cc$eval
    ytr <- y_vector(train_df)

    if (length(final_vars) == 0L) {
      p0 <- mean(ytr)
      pred_by_m[, m] <- rep(p0, nrow(test_df))
      refit_coefs[[m]] <- data.frame(
        imputation = m,
        variable = "(Intercept)",
        coefficient = qlogis(clip_prob(p0)),
        stringsAsFactors = FALSE
      )
    } else {
      fit_df <- data.frame(
        y = ytr,
        train_df[, final_vars, drop = FALSE],
        check.names = FALSE
      )
      form <- stats::as.formula(
        paste("y ~", paste(sprintf("`%s`", final_vars), collapse = " + "))
      )

      fit <- suppressWarnings(
        stats::glm(form, data = fit_df, family = stats::binomial())
      )

      new_df <- test_df[, final_vars, drop = FALSE]
      pp <- suppressWarnings(
        stats::predict(fit, newdata = new_df, type = "response")
      )

      if (any(!is.finite(pp))) {
        stop(
          "Primary logistic refit produced non-finite predictions at R",
          r, " O", of, " M", m
        )
      }

      pred_by_m[, m] <- pp
      cf <- stats::coef(fit)
      refit_coefs[[m]] <- data.frame(
        imputation = m,
        variable = names(cf),
        coefficient = as.numeric(cf),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    prediction = rowMeans(pred_by_m),
    prediction_by_m = pred_by_m,
    selected_by_m = selected_by_m,
    final_vars = final_vars,
    lasso_coefs = do.call(rbind, lasso_coefs),
    refit_coefs = do.call(rbind, refit_coefs)
  )
}

# ----- Elastic Net -----
tune_elastic_net <- function(r, of) {
  grid <- expand.grid(
    alpha = ENET_ALPHA_GRID,
    lambda = LAMBDA_GRID,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  losses <- matrix(
    NA_real_,
    nrow = N_INNER_FOLDS,
    ncol = nrow(grid)
  )

  for (infold in seq_len(N_INNER_FOLDS)) {
    obj <- readRDS(inner_cache_path(r, of, infold))
    yva <- y_vector(obj$completed[[1L]]$eval)

    # One MI-averaged prediction matrix per alpha, then fill lambda columns.
    col_offset <- 0L

    for (a in ENET_ALPHA_GRID) {
      pred_sum <- matrix(
        0,
        nrow = length(yva),
        ncol = length(LAMBDA_GRID)
      )

      for (m in seq_len(M)) {
        cc <- obj$completed[[m]]
        fit <- glmnet::glmnet(
          x = x_matrix(cc$train),
          y = y_vector(cc$train),
          family = "binomial",
          alpha = a,
          lambda = LAMBDA_GRID,
          standardize = TRUE,
          intercept = TRUE
        )
        pred_sum <- pred_sum +
          predict_glmnet_matrix(fit, x_matrix(cc$eval), LAMBDA_GRID)
      }

      pred_avg <- pred_sum / M

      idx <- which(grid$alpha == a)
      for (jj in seq_along(idx)) {
        losses[infold, idx[jj]] <- logloss(yva, pred_avg[, jj])
      }
    }
  }

  grid$mean_logloss <- colMeans(losses)
  grid$se_logloss <- apply(losses, 2L, stats::sd) / sqrt(N_INNER_FOLDS)
  best <- grid[which.min(grid$mean_logloss), , drop = FALSE]

  list(best = best, tuning = grid)
}

fit_elastic_net_outer <- function(outer_obj, best) {
  pred_by_m <- matrix(
    NA_real_,
    nrow = nrow(outer_obj$completed[[1L]]$eval),
    ncol = M
  )
  coef_list <- vector("list", M)

  for (m in seq_len(M)) {
    cc <- outer_obj$completed[[m]]
    fit <- glmnet::glmnet(
      x = x_matrix(cc$train),
      y = y_vector(cc$train),
      family = "binomial",
      alpha = best$alpha,
      lambda = best$lambda,
      standardize = TRUE,
      intercept = TRUE
    )
    pred_by_m[, m] <- as.numeric(
      stats::predict(
        fit,
        newx = x_matrix(cc$eval),
        s = best$lambda,
        type = "response"
      )
    )
    cf <- as.matrix(stats::coef(fit, s = best$lambda))
    coef_list[[m]] <- data.frame(
      imputation = m,
      variable = rownames(cf),
      coefficient = cf[, 1],
      stringsAsFactors = FALSE
    )
  }

  list(
    prediction = rowMeans(pred_by_m),
    prediction_by_m = pred_by_m,
    coefficients = do.call(rbind, coef_list)
  )
}

# ----- Random forest -----
rf_fit_predict <- function(train_df, eval_df, pars, seed) {
  dtrain <- data.frame(
    y = factor(y_vector(train_df), levels = c(0, 1)),
    train_df[, candidate_vars, drop = FALSE],
    check.names = FALSE
  )

  fit <- ranger::ranger(
    dependent.variable.name = "y",
    data = dtrain,
    probability = TRUE,
    num.trees = RF_TREES,
    mtry = as.integer(pars$mtry),
    min.node.size = as.integer(pars$min_node_size),
    sample.fraction = as.numeric(pars$sample_fraction),
    importance = "impurity",
    seed = seed,
    num.threads = N_THREADS
  )

  pp <- stats::predict(
    fit,
    data = eval_df[, candidate_vars, drop = FALSE],
    num.threads = N_THREADS
  )$predictions

  if (is.matrix(pp)) {
    if ("1" %in% colnames(pp)) {
      p <- pp[, "1"]
    } else {
      p <- pp[, ncol(pp)]
    }
  } else {
    p <- pp
  }

  list(pred = as.numeric(p), importance = fit$variable.importance)
}

tune_rf <- function(r, of) {
  losses <- matrix(
    NA_real_,
    nrow = N_INNER_FOLDS,
    ncol = nrow(RF_GRID)
  )

  for (infold in seq_len(N_INNER_FOLDS)) {
    obj <- readRDS(inner_cache_path(r, of, infold))
    yva <- y_vector(obj$completed[[1L]]$eval)

    for (g in seq_len(nrow(RF_GRID))) {
      pred_sum <- rep(0, length(yva))

      for (m in seq_len(M)) {
        cc <- obj$completed[[m]]
        ans <- rf_fit_predict(
          cc$train,
          cc$eval,
          RF_GRID[g, ],
          seed_from(SEED_MODEL_BASE, r, of, infold, m, g)
        )
        pred_sum <- pred_sum + ans$pred
      }

      losses[infold, g] <- logloss(yva, pred_sum / M)
    }
  }

  tab <- RF_GRID
  tab$mean_logloss <- colMeans(losses)
  tab$se_logloss <- apply(losses, 2L, stats::sd) / sqrt(N_INNER_FOLDS)
  best <- tab[which.min(tab$mean_logloss), , drop = FALSE]

  list(best = best, tuning = tab)
}

fit_rf_outer <- function(r, of, outer_obj, best) {
  pred_by_m <- matrix(
    NA_real_,
    nrow = nrow(outer_obj$completed[[1L]]$eval),
    ncol = M
  )
  imp_list <- vector("list", M)

  for (m in seq_len(M)) {
    cc <- outer_obj$completed[[m]]
    ans <- rf_fit_predict(
      cc$train,
      cc$eval,
      best,
      seed_from(SEED_MODEL_BASE + 1000000L, r, of, 0L, m, 1L)
    )
    pred_by_m[, m] <- ans$pred
    imp_list[[m]] <- data.frame(
      imputation = m,
      variable = names(ans$importance),
      importance = as.numeric(ans$importance),
      stringsAsFactors = FALSE
    )
  }

  list(
    prediction = rowMeans(pred_by_m),
    prediction_by_m = pred_by_m,
    importance = do.call(rbind, imp_list)
  )
}

# ----- XGBoost -----
xgb_fit_predict <- function(train_df, eval_df, pars, seed) {
  dtrain <- xgboost::xgb.DMatrix(
    data = x_matrix(train_df),
    label = y_vector(train_df)
  )
  deval <- xgboost::xgb.DMatrix(data = x_matrix(eval_df))

  fit <- xgboost::xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      max_depth = as.integer(pars$max_depth),
      eta = as.numeric(pars$eta),
      min_child_weight = as.numeric(pars$min_child_weight),
      subsample = as.numeric(pars$subsample),
      colsample_bytree = as.numeric(pars$colsample_bytree),
      nthread = N_THREADS,
      seed = seed
    ),
    data = dtrain,
    nrounds = as.integer(pars$nrounds),
    verbose = 0
  )

  pred <- as.numeric(stats::predict(fit, deval))
  importance <- xgboost::xgb.importance(
    feature_names = candidate_vars,
    model = fit
  )

  list(pred = pred, importance = importance)
}

tune_xgb <- function(r, of) {
  losses <- matrix(
    NA_real_,
    nrow = N_INNER_FOLDS,
    ncol = nrow(XGB_GRID)
  )

  for (infold in seq_len(N_INNER_FOLDS)) {
    obj <- readRDS(inner_cache_path(r, of, infold))
    yva <- y_vector(obj$completed[[1L]]$eval)

    for (g in seq_len(nrow(XGB_GRID))) {
      pred_sum <- rep(0, length(yva))

      for (m in seq_len(M)) {
        cc <- obj$completed[[m]]
        ans <- xgb_fit_predict(
          cc$train,
          cc$eval,
          XGB_GRID[g, ],
          seed_from(SEED_MODEL_BASE + 2000000L, r, of, infold, m, g)
        )
        pred_sum <- pred_sum + ans$pred
      }

      losses[infold, g] <- logloss(yva, pred_sum / M)
    }
  }

  tab <- XGB_GRID
  tab$mean_logloss <- colMeans(losses)
  tab$se_logloss <- apply(losses, 2L, stats::sd) / sqrt(N_INNER_FOLDS)
  best <- tab[which.min(tab$mean_logloss), , drop = FALSE]

  list(best = best, tuning = tab)
}

fit_xgb_outer <- function(r, of, outer_obj, best) {
  pred_by_m <- matrix(
    NA_real_,
    nrow = nrow(outer_obj$completed[[1L]]$eval),
    ncol = M
  )
  imp_list <- vector("list", M)

  for (m in seq_len(M)) {
    cc <- outer_obj$completed[[m]]
    ans <- xgb_fit_predict(
      cc$train,
      cc$eval,
      best,
      seed_from(SEED_MODEL_BASE + 3000000L, r, of, 0L, m, 1L)
    )
    pred_by_m[, m] <- ans$pred

    imp <- ans$importance
    if (is.null(imp) || nrow(imp) == 0L) {
      imp_list[[m]] <- data.frame(
        imputation = m,
        Feature = character(0),
        Gain = numeric(0),
        Cover = numeric(0),
        Frequency = numeric(0)
      )
    } else {
      imp$imputation <- m
      imp_list[[m]] <- imp
    }
  }

  list(
    prediction = rowMeans(pred_by_m),
    prediction_by_m = pred_by_m,
    importance = do.call(rbind, imp_list)
  )
}

# ----- SVM RBF -----
svm_fit_predict <- function(train_df, eval_df, pars) {
  xtr <- x_matrix(train_df)
  xev <- x_matrix(eval_df)
  ytr <- factor(y_vector(train_df), levels = c(0, 1))

  # Training-only scaling.
  mu <- colMeans(xtr)
  sds <- apply(xtr, 2L, stats::sd)
  if (any(!is.finite(sds)) || any(sds <= 0)) {
    stop("SVM encountered zero/non-finite training SD.")
  }

  xtr_s <- sweep(sweep(xtr, 2L, mu, "-"), 2L, sds, "/")
  xev_s <- sweep(sweep(xev, 2L, mu, "-"), 2L, sds, "/")

  fit <- kernlab::ksvm(
    x = xtr_s,
    y = ytr,
    type = "C-svc",
    kernel = "rbfdot",
    kpar = list(sigma = as.numeric(pars$sigma)),
    C = as.numeric(pars$C),
    prob.model = TRUE,
    scaled = FALSE
  )

  # kernlab::ksvm is an S4 model. Calling stats::predict() explicitly
  # forces S3 UseMethod dispatch on some R/kernlab versions and can fail with:
  #   no applicable method for class c("ksvm", "vm")
  # Use kernlab's exported S4 predict generic/method instead.
  pp <- kernlab::predict(fit, xev_s, type = "probabilities")
  pp <- as.matrix(pp)

  if ("1" %in% colnames(pp)) {
    p <- pp[, "1"]
  } else {
    p <- pp[, ncol(pp)]
  }

  as.numeric(p)
}

tune_svm <- function(r, of) {
  losses <- matrix(
    NA_real_,
    nrow = N_INNER_FOLDS,
    ncol = nrow(SVM_GRID)
  )

  for (infold in seq_len(N_INNER_FOLDS)) {
    obj <- readRDS(inner_cache_path(r, of, infold))
    yva <- y_vector(obj$completed[[1L]]$eval)

    for (g in seq_len(nrow(SVM_GRID))) {
      pred_sum <- rep(0, length(yva))

      for (m in seq_len(M)) {
        cc <- obj$completed[[m]]
        set.seed(seed_from(SEED_MODEL_BASE + 4000000L, r, of, infold, m, g))
        p <- svm_fit_predict(cc$train, cc$eval, SVM_GRID[g, ])
        pred_sum <- pred_sum + p
      }

      losses[infold, g] <- logloss(yva, pred_sum / M)
    }
  }

  tab <- SVM_GRID
  tab$mean_logloss <- colMeans(losses)
  tab$se_logloss <- apply(losses, 2L, stats::sd) / sqrt(N_INNER_FOLDS)
  best <- tab[which.min(tab$mean_logloss), , drop = FALSE]

  list(best = best, tuning = tab)
}

fit_svm_outer <- function(r, of, outer_obj, best) {
  pred_by_m <- matrix(
    NA_real_,
    nrow = nrow(outer_obj$completed[[1L]]$eval),
    ncol = M
  )

  for (m in seq_len(M)) {
    cc <- outer_obj$completed[[m]]
    set.seed(seed_from(SEED_MODEL_BASE + 5000000L, r, of, 0L, m, 1L))
    pred_by_m[, m] <- svm_fit_predict(cc$train, cc$eval, best)
  }

  list(
    prediction = rowMeans(pred_by_m),
    prediction_by_m = pred_by_m
  )
}

# ----------------------------------------------------------------------
# 7. TRAIN / TUNE ALL MODELS, OUTER FOLD BY OUTER FOLD
# ----------------------------------------------------------------------
MODEL_START <- Sys.time()
log_msg("Stage B: nested tuning and outer-fold model fitting started.")

outer_counter <- 0L

for (r in seq_len(N_REPEATS)) {
  split_r <- outer_splits[outer_splits$repeat_id == r, , drop = FALSE]

  for (of in seq_len(N_OUTER_FOLDS)) {
    outer_counter <- outer_counter + 1L
    rp <- result_path(r, of)

    if (file.exists(rp)) {
      log_msg("Reusing completed model result R", r, " O", of, ".")
      next
    }

    fold_start <- Sys.time()
    log_msg(
      "===== OUTER ", outer_counter, "/",
      N_REPEATS * N_OUTER_FOLDS,
      " : R", r, " O", of, " ====="
    )

    outer_obj <- readRDS(outer_cache_path(r, of))
    test_ids <- outer_obj$eval_ids
    ytest <- y[test_ids]

    fold_result <- list(
      repeat_id = r,
      outer_fold = of,
      train_ids = outer_obj$train_ids,
      test_ids = test_ids,
      outcome_test = ytest,
      algorithms = list(),
      timing = list(),
      status = "running"
    )

    # ---------------- Primary LASSO 5/5 ----------------
    t0 <- Sys.time()
    log_msg("Primary_LASSO_5of5 tuning R", r, " O", of)
    primary_tuned <- tune_primary_lasso_5of5(r, of)
    utils::write.csv(
      primary_tuned$tuning,
      tuning_path(r, of, "Primary_LASSO_5of5"),
      row.names = FALSE
    )
    primary_fit <- fit_primary_outer(r, of, outer_obj, primary_tuned)
    tsec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    fold_result$algorithms$Primary_LASSO_5of5 <- list(
      prediction = primary_fit$prediction,
      prediction_by_m = primary_fit$prediction_by_m,
      chosen_lambda_by_m = primary_tuned$chosen_lambda,
      selected_by_m = primary_fit$selected_by_m,
      final_vars = primary_fit$final_vars,
      lasso_coefs = primary_fit$lasso_coefs,
      refit_coefs = primary_fit$refit_coefs
    )
    fold_result$timing$Primary_LASSO_5of5 <- tsec
    log_msg(
      "Primary done in ", fmt_time(tsec),
      "; final vars = ",
      ifelse(
        length(primary_fit$final_vars) == 0L,
        "<intercept only>",
        paste(primary_fit$final_vars, collapse = ", ")
      )
    )

    # ---------------- Elastic Net ----------------
    t0 <- Sys.time()
    log_msg("ElasticNet tuning R", r, " O", of)
    enet_tuned <- tune_elastic_net(r, of)
    utils::write.csv(
      enet_tuned$tuning,
      tuning_path(r, of, "ElasticNet"),
      row.names = FALSE
    )
    enet_fit <- fit_elastic_net_outer(outer_obj, enet_tuned$best)
    tsec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    fold_result$algorithms$ElasticNet <- list(
      prediction = enet_fit$prediction,
      prediction_by_m = enet_fit$prediction_by_m,
      best = enet_tuned$best,
      coefficients = enet_fit$coefficients
    )
    fold_result$timing$ElasticNet <- tsec
    log_msg("ElasticNet done in ", fmt_time(tsec))

    # ---------------- Random Forest ----------------
    t0 <- Sys.time()
    log_msg("RandomForest tuning R", r, " O", of)
    rf_tuned <- tune_rf(r, of)
    utils::write.csv(
      rf_tuned$tuning,
      tuning_path(r, of, "RandomForest"),
      row.names = FALSE
    )
    rf_fit <- fit_rf_outer(r, of, outer_obj, rf_tuned$best)
    tsec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    fold_result$algorithms$RandomForest <- list(
      prediction = rf_fit$prediction,
      prediction_by_m = rf_fit$prediction_by_m,
      best = rf_tuned$best,
      importance = rf_fit$importance
    )
    fold_result$timing$RandomForest <- tsec
    log_msg("RandomForest done in ", fmt_time(tsec))

    # ---------------- XGBoost ----------------
    t0 <- Sys.time()
    log_msg("XGBoost tuning R", r, " O", of)
    xgb_tuned <- tune_xgb(r, of)
    utils::write.csv(
      xgb_tuned$tuning,
      tuning_path(r, of, "XGBoost"),
      row.names = FALSE
    )
    xgb_fit <- fit_xgb_outer(r, of, outer_obj, xgb_tuned$best)
    tsec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    fold_result$algorithms$XGBoost <- list(
      prediction = xgb_fit$prediction,
      prediction_by_m = xgb_fit$prediction_by_m,
      best = xgb_tuned$best,
      importance = xgb_fit$importance
    )
    fold_result$timing$XGBoost <- tsec
    log_msg("XGBoost done in ", fmt_time(tsec))

    # ---------------- SVM RBF ----------------
    t0 <- Sys.time()
    log_msg("SVM_RBF tuning R", r, " O", of)
    svm_tuned <- tune_svm(r, of)
    utils::write.csv(
      svm_tuned$tuning,
      tuning_path(r, of, "SVM_RBF"),
      row.names = FALSE
    )
    svm_fit <- fit_svm_outer(r, of, outer_obj, svm_tuned$best)
    tsec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    fold_result$algorithms$SVM_RBF <- list(
      prediction = svm_fit$prediction,
      prediction_by_m = svm_fit$prediction_by_m,
      best = svm_tuned$best
    )
    fold_result$timing$SVM_RBF <- tsec
    log_msg("SVM_RBF done in ", fmt_time(tsec))

    # Hard prediction audit.
    for (alg in ALGORITHMS) {
      pp <- fold_result$algorithms[[alg]]$prediction
      if (length(pp) != length(test_ids)) {
        stop("Prediction length mismatch: ", alg, " R", r, " O", of)
      }
      if (any(!is.finite(pp)) || any(pp < 0 | pp > 1)) {
        stop("Invalid prediction: ", alg, " R", r, " O", of)
      }
    }

    fold_result$status <- "PASS"
    fold_result$fold_elapsed_seconds <- as.numeric(
      difftime(Sys.time(), fold_start, units = "secs")
    )

    # Save immediately so interruption does not lose this outer fold.
    saveRDS(fold_result, rp, compress = "xz")

    elapsed_model <- as.numeric(difftime(Sys.time(), MODEL_START, units = "secs"))
    done_n <- length(list.files(RESULT_DIR, pattern = "^result_R.*\\.rds$"))
    avg_per_done <- elapsed_model / max(done_n, 1L)
    remain_est <- avg_per_done * ((N_REPEATS * N_OUTER_FOLDS) - done_n)

    log_msg(
      "Outer R", r, " O", of, " PASS; fold wall time = ",
      fmt_time(fold_result$fold_elapsed_seconds),
      ". Completed ", done_n, "/", N_REPEATS * N_OUTER_FOLDS, ".",
      " Current model-stage elapsed = ", fmt_time(elapsed_model),
      "; rough remaining estimate = ", fmt_time(remain_est)
    )

    gc()
  }
}

MODEL_ELAPSED <- as.numeric(difftime(Sys.time(), MODEL_START, units = "secs"))
log_msg("Stage B complete. Wall time = ", fmt_time(MODEL_ELAPSED))


# Robust row-bind for heterogeneous hyperparameter tables.
# Different algorithms legitimately have different parameter columns
# (e.g., ElasticNet alpha/lambda; RF mtry/min.node.size; XGBoost params;
# SVM C/sigma; Primary imputation/lambda). Base rbind() cannot combine them.
rbind_fill_base <- function(lst) {
  lst <- Filter(function(z) !is.null(z) && nrow(as.data.frame(z)) > 0L, lst)
  if (length(lst) == 0L) return(data.frame())

  dfs <- lapply(lst, function(z) as.data.frame(z, stringsAsFactors = FALSE))
  all_names <- unique(unlist(lapply(dfs, names), use.names = FALSE))

  dfs2 <- lapply(dfs, function(d) {
    miss <- setdiff(all_names, names(d))
    if (length(miss) > 0L) {
      for (nm in miss) d[[nm]] <- NA
    }
    d <- d[, all_names, drop = FALSE]
    rownames(d) <- NULL
    d
  })

  out <- do.call(rbind, dfs2)
  rownames(out) <- NULL
  out
}

# ----------------------------------------------------------------------
# 8. COMBINE ALL PATIENT-LEVEL OOF PREDICTIONS AND INTERMEDIATES
# ----------------------------------------------------------------------
COMBINE_START <- Sys.time()
log_msg("Stage C: combining all saved outer results.")

pred_rows <- list()
fold_metric_rows <- list()
timing_rows <- list()
best_param_rows <- list()
selection_rows <- list()
primary_coef_rows <- list()
enet_coef_rows <- list()
rf_imp_rows <- list()
xgb_imp_rows <- list()

ip <- imet <- itime <- ibest <- isel <- ipc <- iec <- iri <- ixi <- 0L

for (r in seq_len(N_REPEATS)) {
  for (of in seq_len(N_OUTER_FOLDS)) {
    rp <- result_path(r, of)
    if (!file.exists(rp)) stop("Missing outer result: ", rp)

    obj <- readRDS(rp)
    if (!identical(obj$status, "PASS")) stop("Non-PASS outer result: ", rp)

    test_ids <- obj$test_ids
    yy <- obj$outcome_test

    for (alg in ALGORITHMS) {
      pp <- obj$algorithms[[alg]]$prediction

      ip <- ip + 1L
      pred_rows[[ip]] <- data.frame(
        row_id = test_ids,
        subject_id = if ("subject_id" %in% names(dat)) dat$subject_id[test_ids] else NA,
        outcome = yy,
        repeat_id = r,
        outer_fold = of,
        algorithm = alg,
        prediction = pp,
        stringsAsFactors = FALSE
      )

      mm <- metric_row(yy, pp)
      imet <- imet + 1L
      fold_metric_rows[[imet]] <- data.frame(
        repeat_id = r,
        outer_fold = of,
        algorithm = alg,
        test_n = length(yy),
        events = sum(yy == 1L),
        mm,
        stringsAsFactors = FALSE
      )

      itime <- itime + 1L
      timing_rows[[itime]] <- data.frame(
        repeat_id = r,
        outer_fold = of,
        algorithm = alg,
        elapsed_seconds = as.numeric(obj$timing[[alg]]),
        stringsAsFactors = FALSE
      )
    }

    # Primary selections.
    fv <- obj$algorithms$Primary_LASSO_5of5$final_vars
    if (length(fv) == 0L) {
      isel <- isel + 1L
      selection_rows[[isel]] <- data.frame(
        repeat_id = r,
        outer_fold = of,
        variable = "<INTERCEPT_ONLY>",
        selected_final_5of5 = TRUE,
        stringsAsFactors = FALSE
      )
    } else {
      for (v in fv) {
        isel <- isel + 1L
        selection_rows[[isel]] <- data.frame(
          repeat_id = r,
          outer_fold = of,
          variable = v,
          selected_final_5of5 = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }

    cf <- obj$algorithms$Primary_LASSO_5of5$refit_coefs
    cf$repeat_id <- r
    cf$outer_fold <- of
    ipc <- ipc + 1L
    primary_coef_rows[[ipc]] <- cf

    cf2 <- obj$algorithms$ElasticNet$coefficients
    cf2$repeat_id <- r
    cf2$outer_fold <- of
    iec <- iec + 1L
    enet_coef_rows[[iec]] <- cf2

    rfimp <- obj$algorithms$RandomForest$importance
    rfimp$repeat_id <- r
    rfimp$outer_fold <- of
    iri <- iri + 1L
    rf_imp_rows[[iri]] <- rfimp

    ximp <- obj$algorithms$XGBoost$importance
    if (!is.null(ximp) && nrow(ximp) > 0L) {
      ximp$repeat_id <- r
      ximp$outer_fold <- of
      ixi <- ixi + 1L
      xgb_imp_rows[[ixi]] <- ximp
    }

    # Best parameters.
    bp <- list(
      ElasticNet = obj$algorithms$ElasticNet$best,
      RandomForest = obj$algorithms$RandomForest$best,
      XGBoost = obj$algorithms$XGBoost$best,
      SVM_RBF = obj$algorithms$SVM_RBF$best
    )

    for (nm in names(bp)) {
      tmp <- as.data.frame(bp[[nm]])
      tmp$repeat_id <- r
      tmp$outer_fold <- of
      tmp$algorithm <- nm
      ibest <- ibest + 1L
      best_param_rows[[ibest]] <- tmp
    }

    # Primary lambda by MI.
    lam <- obj$algorithms$Primary_LASSO_5of5$chosen_lambda_by_m
    for (m in seq_along(lam)) {
      ibest <- ibest + 1L
      best_param_rows[[ibest]] <- data.frame(
        imputation = m,
        lambda = lam[m],
        repeat_id = r,
        outer_fold = of,
        algorithm = "Primary_LASSO_5of5",
        stringsAsFactors = FALSE
      )
    }
  }
}

predictions <- do.call(rbind, pred_rows)
fold_metrics <- do.call(rbind, fold_metric_rows)
timing_model <- do.call(rbind, timing_rows)
best_params <- rbind_fill_base(best_param_rows)
primary_selections <- do.call(rbind, selection_rows)
primary_coefs <- do.call(rbind, primary_coef_rows)
enet_coefs <- do.call(rbind, enet_coef_rows)
rf_importance <- do.call(rbind, rf_imp_rows)
xgb_importance <- if (length(xgb_imp_rows) > 0L) do.call(rbind, xgb_imp_rows) else data.frame()

utils::write.csv(
  predictions,
  file.path(ROOT_DIR, "390ML_06C_02_patient_level_OOF_predictions.csv"),
  row.names = FALSE
)
utils::write.csv(
  fold_metrics,
  file.path(ROOT_DIR, "390ML_06C_03_outer_fold_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  timing_model,
  file.path(ROOT_DIR, "390ML_06C_04_model_timing_each_outer_fold.csv"),
  row.names = FALSE
)
utils::write.csv(
  best_params,
  file.path(ROOT_DIR, "390ML_06C_05_best_hyperparameters_all_outer.csv"),
  row.names = FALSE
)
utils::write.csv(
  primary_selections,
  file.path(ROOT_DIR, "390ML_06C_06_primary_selected_variables_each_outer.csv"),
  row.names = FALSE
)
utils::write.csv(
  primary_coefs,
  file.path(ROOT_DIR, "390ML_06C_07_primary_refit_coefficients_all_outer.csv"),
  row.names = FALSE
)
utils::write.csv(
  enet_coefs,
  file.path(ROOT_DIR, "390ML_06C_08_elasticnet_coefficients_all_outer.csv"),
  row.names = FALSE
)
utils::write.csv(
  rf_importance,
  file.path(ROOT_DIR, "390ML_06C_09_randomforest_importance_all_outer.csv"),
  row.names = FALSE
)
utils::write.csv(
  xgb_importance,
  file.path(ROOT_DIR, "390ML_06C_10_xgboost_importance_all_outer.csv"),
  row.names = FALSE
)

# Hard OOF integrity: each algorithm must have N=390 predictions per repeat.
oof_count <- aggregate(
  prediction ~ repeat_id + algorithm,
  data = predictions,
  FUN = length
)
if (any(oof_count$prediction != EXPECTED_N)) {
  stop("OOF integrity failure: not exactly N=390 predictions per repeat/algorithm.")
}

# ----------------------------------------------------------------------
# 9. REPEAT-LEVEL PERFORMANCE
# ----------------------------------------------------------------------
repeat_metric_rows <- list()
ir <- 0L

for (r in seq_len(N_REPEATS)) {
  for (alg in ALGORITHMS) {
    dd <- predictions[
      predictions$repeat_id == r & predictions$algorithm == alg,
      ,
      drop = FALSE
    ]
    dd <- dd[order(dd$row_id), ]

    mm <- metric_row(dd$outcome, dd$prediction)

    ir <- ir + 1L
    repeat_metric_rows[[ir]] <- data.frame(
      repeat_id = r,
      algorithm = alg,
      n = nrow(dd),
      events = sum(dd$outcome == 1L),
      mm,
      stringsAsFactors = FALSE
    )
  }
}

repeat_metrics <- do.call(rbind, repeat_metric_rows)
utils::write.csv(
  repeat_metrics,
  file.path(ROOT_DIR, "390ML_06C_11_repeat_level_OOF_metrics.csv"),
  row.names = FALSE
)

# Summary of repeat-level metrics.
summary_rows <- list()
isum <- 0L
for (alg in ALGORITHMS) {
  dd <- repeat_metrics[repeat_metrics$algorithm == alg, ]

  for (metric in c("AUC", "Brier", "Cal_Intercept", "Cal_Slope")) {
    z <- dd[[metric]]
    isum <- isum + 1L
    summary_rows[[isum]] <- data.frame(
      algorithm = alg,
      metric = metric,
      mean = mean(z, na.rm = TRUE),
      sd = stats::sd(z, na.rm = TRUE),
      median = stats::median(z, na.rm = TRUE),
      q025 = stats::quantile(z, 0.025, na.rm = TRUE, names = FALSE),
      q975 = stats::quantile(z, 0.975, na.rm = TRUE, names = FALSE),
      n_defined = sum(is.finite(z)),
      stringsAsFactors = FALSE
    )
  }
}
performance_summary <- do.call(rbind, summary_rows)

utils::write.csv(
  performance_summary,
  file.path(ROOT_DIR, "390ML_06C_12_performance_summary_repeat_distribution.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 10. PAIRED CORRECTED RESAMPLED TESTS (NADEAU-BENGIO STYLE)
# ----------------------------------------------------------------------
# The formal comparison uses paired OUTER-FOLD differences because all
# algorithms are evaluated on exactly the same frozen test folds.
# The corrected-resampled variance inflation accounts for overlapping
# training samples: correction = 1/n_splits + mean(n_test/n_train).
#
# This is NOT a naive DeLong test on pooled repeated OOF predictions.
# ----------------------------------------------------------------------
corrected_resampled_test <- function(diff_vec, test_train_ratio, metric_name, comparator) {
  ok <- is.finite(diff_vec)
  d <- diff_vec[ok]
  n <- length(d)

  if (n < 3L) {
    return(data.frame(
      metric = metric_name,
      comparison = paste(comparator, "-", "Primary_LASSO_5of5"),
      n_outer_folds = n,
      mean_difference = NA_real_,
      corrected_SE = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      t_statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  s2 <- stats::var(d)
  correction <- (1 / n) + mean(test_train_ratio[ok])
  se_corr <- sqrt(correction * s2)
  mean_d <- mean(d)

  if (!is.finite(se_corr) || se_corr == 0) {
    tstat <- NA_real_
    pval <- NA_real_
    lo <- mean_d
    hi <- mean_d
  } else {
    tstat <- mean_d / se_corr
    pval <- 2 * stats::pt(-abs(tstat), df = n - 1L)
    crit <- stats::qt(0.975, df = n - 1L)
    lo <- mean_d - crit * se_corr
    hi <- mean_d + crit * se_corr
  }

  data.frame(
    metric = metric_name,
    comparison = paste(comparator, "-", "Primary_LASSO_5of5"),
    n_outer_folds = n,
    mean_difference = mean_d,
    corrected_SE = se_corr,
    CI_lower = lo,
    CI_upper = hi,
    t_statistic = tstat,
    df = n - 1L,
    p_value = pval,
    stringsAsFactors = FALSE
  )
}

# Match fold metrics by repeat/fold.
primary_fold <- fold_metrics[
  fold_metrics$algorithm == "Primary_LASSO_5of5",
  ,
  drop = FALSE
]
primary_fold <- primary_fold[
  order(primary_fold$repeat_id, primary_fold$outer_fold),
]

test_train_ratio <- primary_fold$test_n / (EXPECTED_N - primary_fold$test_n)

comparators <- setdiff(ALGORITHMS, "Primary_LASSO_5of5")
auc_tests <- list()
brier_tests <- list()

for (i in seq_along(comparators)) {
  alg <- comparators[i]
  dd <- fold_metrics[fold_metrics$algorithm == alg, , drop = FALSE]
  dd <- dd[order(dd$repeat_id, dd$outer_fold), ]

  if (!all(
    dd$repeat_id == primary_fold$repeat_id &
      dd$outer_fold == primary_fold$outer_fold
  )) {
    stop("Fold matching failure for ", alg)
  }

  auc_tests[[i]] <- corrected_resampled_test(
    dd$AUC - primary_fold$AUC,
    test_train_ratio,
    "AUC",
    alg
  )

  # For Brier, lower is better; difference is comparator - primary.
  brier_tests[[i]] <- corrected_resampled_test(
    dd$Brier - primary_fold$Brier,
    test_train_ratio,
    "Brier",
    alg
  )
}

auc_tests <- do.call(rbind, auc_tests)
brier_tests <- do.call(rbind, brier_tests)

auc_tests$p_Holm <- stats::p.adjust(auc_tests$p_value, method = "holm")
brier_tests$p_Holm <- stats::p.adjust(brier_tests$p_value, method = "holm")

utils::write.csv(
  auc_tests,
  file.path(ROOT_DIR, "390ML_06C_13_paired_AUC_corrected_resampled_tests_Holm.csv"),
  row.names = FALSE
)
utils::write.csv(
  brier_tests,
  file.path(ROOT_DIR, "390ML_06C_14_paired_Brier_corrected_resampled_tests_Holm.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 11. SELECTION / IMPORTANCE SUMMARIES
# ----------------------------------------------------------------------
# Primary selection frequency over the 50 outer training sets.
sel_counts <- table(primary_selections$variable)
selection_summary <- data.frame(
  variable = names(sel_counts),
  selected_outer_folds = as.integer(sel_counts),
  selection_frequency_pct = 100 * as.integer(sel_counts) /
    (N_REPEATS * N_OUTER_FOLDS),
  stringsAsFactors = FALSE
)
selection_summary <- selection_summary[
  order(-selection_summary$selection_frequency_pct, selection_summary$variable),
]

utils::write.csv(
  selection_summary,
  file.path(ROOT_DIR, "390ML_06C_15_primary_selection_frequency_outerCV.csv"),
  row.names = FALSE
)

# Random forest mean importance.
rf_imp_summary <- aggregate(
  importance ~ variable,
  data = rf_importance,
  FUN = mean
)
rf_imp_summary <- rf_imp_summary[order(-rf_imp_summary$importance), ]

utils::write.csv(
  rf_imp_summary,
  file.path(ROOT_DIR, "390ML_06C_16_randomforest_mean_importance.csv"),
  row.names = FALSE
)

# XGBoost mean Gain, treating unreported variables as zero is not done here;
# this table summarizes gain among appearances. Raw fold-level importance was saved.
if (nrow(xgb_importance) > 0L && all(c("Feature", "Gain") %in% names(xgb_importance))) {
  xgb_imp_summary <- aggregate(
    Gain ~ Feature,
    data = xgb_importance,
    FUN = mean
  )
  xgb_imp_summary <- xgb_imp_summary[order(-xgb_imp_summary$Gain), ]
} else {
  xgb_imp_summary <- data.frame(Feature = character(0), Gain = numeric(0))
}

utils::write.csv(
  xgb_imp_summary,
  file.path(ROOT_DIR, "390ML_06C_17_xgboost_mean_gain.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 12. TIMING SUMMARY
# ----------------------------------------------------------------------
COMBINE_ELAPSED <- as.numeric(difftime(Sys.time(), COMBINE_START, units = "secs"))
OVERALL_ELAPSED <- as.numeric(difftime(Sys.time(), OVERALL_START, units = "secs"))

model_timing_summary <- aggregate(
  elapsed_seconds ~ algorithm,
  data = timing_model,
  FUN = function(z) c(
    total = sum(z),
    mean = mean(z),
    median = stats::median(z),
    min = min(z),
    max = max(z)
  )
)

model_timing_summary_out <- data.frame(
  algorithm = model_timing_summary$algorithm,
  total_seconds = model_timing_summary$elapsed_seconds[, "total"],
  mean_seconds_per_outer = model_timing_summary$elapsed_seconds[, "mean"],
  median_seconds_per_outer = model_timing_summary$elapsed_seconds[, "median"],
  min_seconds_per_outer = model_timing_summary$elapsed_seconds[, "min"],
  max_seconds_per_outer = model_timing_summary$elapsed_seconds[, "max"],
  stringsAsFactors = FALSE
)

utils::write.csv(
  model_timing_summary_out,
  file.path(ROOT_DIR, "390ML_06C_18_model_timing_summary.csv"),
  row.names = FALSE
)

stage_timing <- data.frame(
  stage = c(
    "MICE_cache_stage_wall_time",
    "Model_nested_tuning_and_outer_fit_wall_time",
    "Combine_metrics_tables_figures_wall_time",
    "Entire_script_wall_time"
  ),
  seconds = c(
    CACHE_ELAPSED,
    MODEL_ELAPSED,
    COMBINE_ELAPSED,
    OVERALL_ELAPSED
  ),
  formatted_hh_mm_ss = vapply(
    c(CACHE_ELAPSED, MODEL_ELAPSED, COMBINE_ELAPSED, OVERALL_ELAPSED),
    fmt_time,
    character(1)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  stage_timing,
  file.path(ROOT_DIR, "390ML_06C_19_STAGE_AND_TOTAL_TIMING.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 13. FIGURES
# ----------------------------------------------------------------------
# Figure 1: repeat-level OOF AUC distributions.
png(
  filename = file.path(FIG_DIR, "390ML_06C_FIG1_repeat_level_AUC_boxplot.png"),
  width = 2000,
  height = 1400,
  res = 220
)
boxplot(
  AUC ~ algorithm,
  data = repeat_metrics,
  las = 2,
  ylab = "Repeated outer-CV OOF AUC",
  xlab = "",
  main = "Nested repeated-CV discrimination"
)
stripchart(
  AUC ~ algorithm,
  data = repeat_metrics,
  vertical = TRUE,
  method = "jitter",
  add = TRUE,
  pch = 16
)
dev.off()

# Figure 2: corrected paired AUC differences vs primary.
png(
  filename = file.path(FIG_DIR, "390ML_06C_FIG2_paired_deltaAUC_corrected_CI.png"),
  width = 1800,
  height = 1200,
  res = 220
)
op <- par(mar = c(5, 10, 4, 2))
yy <- seq_len(nrow(auc_tests))
xlim <- range(
  c(auc_tests$CI_lower, auc_tests$CI_upper, 0),
  finite = TRUE
)
plot(
  auc_tests$mean_difference,
  yy,
  xlim = xlim,
  ylim = c(0.5, nrow(auc_tests) + 0.5),
  yaxt = "n",
  ylab = "",
  xlab = "ΔAUC (comparator - Primary LASSO 5/5)",
  pch = 19,
  main = "Paired corrected-resampled AUC differences"
)
axis(
  2,
  at = yy,
  labels = auc_tests$comparison,
  las = 2
)
abline(v = 0, lty = 2)
segments(
  x0 = auc_tests$CI_lower,
  y0 = yy,
  x1 = auc_tests$CI_upper,
  y1 = yy
)
points(auc_tests$mean_difference, yy, pch = 19)
par(op)
dev.off()

# ----------------------------------------------------------------------
# 14. RUN STATUS, COMPLETE OBJECT, SESSION INFO, MANIFEST
# ----------------------------------------------------------------------
run_status <- data.frame(
  item = c(
    "N",
    "Events",
    "Event_rate",
    "Candidate_predictors",
    "Outer_repeats",
    "Outer_folds",
    "Total_outer_test_folds",
    "Inner_folds",
    "MICE_m",
    "MICE_maxit",
    "Formal_algorithms",
    "Completed_outer_result_files",
    "Patient_level_prediction_rows",
    "OOF_integrity_expected_rows",
    "Cache_stage_seconds",
    "Model_stage_seconds",
    "Combine_stage_seconds",
    "Total_seconds",
    "Total_formatted",
    "Run_status"
  ),
  value = c(
    EXPECTED_N,
    EXPECTED_EVENTS,
    mean(y),
    EXPECTED_CANDIDATES,
    N_REPEATS,
    N_OUTER_FOLDS,
    N_REPEATS * N_OUTER_FOLDS,
    N_INNER_FOLDS,
    M,
    MICE_MAXIT,
    paste(ALGORITHMS, collapse = " | "),
    length(list.files(RESULT_DIR, pattern = "^result_R.*\\.rds$")),
    nrow(predictions),
    EXPECTED_N * N_REPEATS * length(ALGORITHMS),
    CACHE_ELAPSED,
    MODEL_ELAPSED,
    COMBINE_ELAPSED,
    OVERALL_ELAPSED,
    fmt_time(OVERALL_ELAPSED),
    "PASS"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  run_status,
  file.path(ROOT_DIR, "390ML_06C_20_run_status.csv"),
  row.names = FALSE
)

complete_object <- list(
  settings = list(
    work_dir = WORK_DIR,
    data_file = DATA_FILE,
    outer_file = OUTER_FILE,
    inner_file = INNER_FILE,
    outcome = OUTCOME,
    N = EXPECTED_N,
    events = EXPECTED_EVENTS,
    candidates = candidate_vars,
    algorithms = ALGORITHMS,
    repeats = N_REPEATS,
    outer_folds = N_OUTER_FOLDS,
    inner_folds = N_INNER_FOLDS,
    mice_m = M,
    mice_maxit = MICE_MAXIT,
    lambda_grid = LAMBDA_GRID,
    enet_alpha_grid = ENET_ALPHA_GRID,
    rf_grid = RF_GRID,
    rf_trees = RF_TREES,
    xgb_grid = XGB_GRID,
    svm_grid = SVM_GRID,
    n_threads = N_THREADS,
    seed_mice_inner_base = SEED_MICE_INNER_BASE,
    seed_mice_outer_base = SEED_MICE_OUTER_BASE,
    seed_model_base = SEED_MODEL_BASE
  ),
  frozen_outer_splits = outer_splits,
  frozen_inner_splits = inner_splits,
  cache_timing = cache_timing,
  predictions = predictions,
  fold_metrics = fold_metrics,
  repeat_metrics = repeat_metrics,
  performance_summary = performance_summary,
  auc_corrected_resampled_tests = auc_tests,
  brier_corrected_resampled_tests = brier_tests,
  primary_selection_summary = selection_summary,
  primary_selections_raw = primary_selections,
  primary_refit_coefficients = primary_coefs,
  elasticnet_coefficients = enet_coefs,
  randomforest_importance_raw = rf_importance,
  randomforest_importance_summary = rf_imp_summary,
  xgboost_importance_raw = xgb_importance,
  xgboost_importance_summary = xgb_imp_summary,
  best_hyperparameters = best_params,
  model_timing_each_outer = timing_model,
  model_timing_summary = model_timing_summary_out,
  stage_timing = stage_timing,
  run_status = run_status
)

saveRDS(
  complete_object,
  file.path(ROOT_DIR, "390ML_06C_21_COMPLETE_ANALYSIS_OBJECT.rds"),
  compress = "xz"
)

capture.output(
  sessionInfo(),
  file = file.path(ROOT_DIR, "390ML_06C_22_sessionInfo.txt")
)

output_files <- c(
  "390ML_06C_01_cache_timing_manifest.csv",
  "390ML_06C_02_patient_level_OOF_predictions.csv",
  "390ML_06C_03_outer_fold_metrics.csv",
  "390ML_06C_04_model_timing_each_outer_fold.csv",
  "390ML_06C_05_best_hyperparameters_all_outer.csv",
  "390ML_06C_06_primary_selected_variables_each_outer.csv",
  "390ML_06C_07_primary_refit_coefficients_all_outer.csv",
  "390ML_06C_08_elasticnet_coefficients_all_outer.csv",
  "390ML_06C_09_randomforest_importance_all_outer.csv",
  "390ML_06C_10_xgboost_importance_all_outer.csv",
  "390ML_06C_11_repeat_level_OOF_metrics.csv",
  "390ML_06C_12_performance_summary_repeat_distribution.csv",
  "390ML_06C_13_paired_AUC_corrected_resampled_tests_Holm.csv",
  "390ML_06C_14_paired_Brier_corrected_resampled_tests_Holm.csv",
  "390ML_06C_15_primary_selection_frequency_outerCV.csv",
  "390ML_06C_16_randomforest_mean_importance.csv",
  "390ML_06C_17_xgboost_mean_gain.csv",
  "390ML_06C_18_model_timing_summary.csv",
  "390ML_06C_19_STAGE_AND_TOTAL_TIMING.csv",
  "390ML_06C_20_run_status.csv",
  "390ML_06C_21_COMPLETE_ANALYSIS_OBJECT.rds",
  "390ML_06C_22_sessionInfo.txt",
  "figures/390ML_06C_FIG1_repeat_level_AUC_boxplot.png",
  "figures/390ML_06C_FIG2_paired_deltaAUC_corrected_CI.png",
  "logs/390ML_06C_console_log.txt"
)

manifest <- data.frame(
  file = output_files,
  exists = file.exists(file.path(ROOT_DIR, output_files)),
  stringsAsFactors = FALSE
)

utils::write.csv(
  manifest,
  file.path(ROOT_DIR, "390ML_06C_23_output_manifest.csv"),
  row.names = FALSE
)

log_msg("============================================================")
log_msg("06C ALL-IN-ONE NESTED ML COMPLETE: PASS")
log_msg("MICE cache wall time: ", fmt_time(CACHE_ELAPSED))
log_msg("Model wall time:      ", fmt_time(MODEL_ELAPSED))
log_msg("Combine/output time:   ", fmt_time(COMBINE_ELAPSED))
log_msg("TOTAL wall time:       ", fmt_time(OVERALL_ELAPSED))
log_msg("Root output folder: ", ROOT_DIR)
log_msg("Patient-level OOF predictions SAVED.")
log_msg("All split IDs / tuning / selections / timing / complete RDS SAVED.")
log_msg("============================================================")

cat("\n\nFINAL STATUS: PASS\n")
cat("Output folder:", ROOT_DIR, "\n")
cat("Total elapsed:", fmt_time(OVERALL_ELAPSED), "\n")
cat("See 390ML_06C_20_run_status.csv and 390ML_06C_19_STAGE_AND_TOTAL_TIMING.csv\n")
