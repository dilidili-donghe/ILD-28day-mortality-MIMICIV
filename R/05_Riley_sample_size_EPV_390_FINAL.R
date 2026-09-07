# ============================================================================
# 05_Riley_sample_size_EPV_390_FINAL.R
# Retrospective sample-size adequacy assessment for binary prediction model
# Final cohort: N = 390; outcome = death_28days; initial candidate pool from
# Sheet2 where 是否纳入 == "是".
#
# Purpose:
#   1) Report crude events-per-candidate-parameter (EPP/EPV).
#   2) Apply Riley et al. minimum sample-size criteria for binary outcomes.
#   3) Use a conservative anticipated Cox-Snell R2 equal to 15% of the maximum
#      possible Cox-Snell R2, because no independent anticipated R2 was
#      prespecified before model development.
#   4) Save all inputs, calculations, package output, sessionInfo and manifest.
#
# IMPORTANT:
#   This is a retrospective adequacy assessment; it must not be described as a
#   prospective a priori power calculation.
# ============================================================================

options(stringsAsFactors = FALSE)

WORK_DIR <- "D:/R代码/大修1_1"
DATA_FILE <- file.path(WORK_DIR, "390run.xlsx")
OUT_PREFIX <- file.path(WORK_DIR, "390RILEY_")
OUTCOME <- "death_28days"
TARGET_SHRINKAGE <- 0.90
DELTA_NAG_R2 <- 0.05
INTERCEPT_MARGIN <- 0.05
Z_975 <- qnorm(0.975)

required_pkgs <- c("readxl")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       ". Install them before running this script.")
}

if (!file.exists(DATA_FILE)) stop("Data file not found: ", DATA_FILE)

dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[1/5] Reading cohort and candidate-predictor metadata...\n")

dat <- as.data.frame(readxl::read_excel(
  DATA_FILE,
  sheet = "Sheet1",
  na = c("", "NA", "N/A", "NaN", "NULL", ".")
))
meta <- as.data.frame(readxl::read_excel(DATA_FILE, sheet = "Sheet2"))

required_meta_cols <- c("表格中原始名称", "是否纳入", "变量类型")
if (!all(required_meta_cols %in% names(meta))) {
  stop("Sheet2 is missing required metadata columns.")
}
if (!OUTCOME %in% names(dat)) stop("Outcome not found in Sheet1: ", OUTCOME)

candidate_meta <- meta[trimws(as.character(meta$是否纳入)) == "是", , drop = FALSE]
candidate_vars <- as.character(candidate_meta$表格中原始名称)

if (length(candidate_vars) == 0L) stop("No candidate predictors marked 是否纳入 == 是.")
if (anyDuplicated(candidate_vars)) stop("Duplicated candidate predictor names in Sheet2.")
if (!all(candidate_vars %in% names(dat))) {
  stop("Candidate variable(s) absent from Sheet1: ",
       paste(setdiff(candidate_vars, names(dat)), collapse = ", "))
}

# Check parameter count: binary categorical variables contribute 1 parameter;
# continuous variables are assumed to enter as one linear parameter each.
parameter_rows <- lapply(seq_len(nrow(candidate_meta)), function(i) {
  v <- candidate_vars[i]
  typ <- as.character(candidate_meta$变量类型[i])
  x <- dat[[v]]
  observed_levels <- unique(x[!is.na(x)])

  if (grepl("分类", typ)) {
    n_levels <- length(observed_levels)
    if (n_levels != 2L) {
      stop("Categorical candidate '", v, "' has ", n_levels,
           " observed levels. Riley parameter count must use degrees of freedom, not variable count.")
    }
    n_parameters <- 1L
  } else {
    n_levels <- NA_integer_
    n_parameters <- 1L
  }

  data.frame(
    variable = v,
    variable_type = typ,
    observed_levels = n_levels,
    predictor_parameters = n_parameters,
    stringsAsFactors = FALSE
  )
})
parameter_table <- do.call(rbind, parameter_rows)

N <- nrow(dat)
y <- dat[[OUTCOME]]
if (is.factor(y)) y <- as.character(y)
y <- suppressWarnings(as.numeric(y))
if (anyNA(y)) stop("Outcome contains missing/non-numeric values.")
if (!all(y %in% c(0, 1))) stop("Outcome must be binary 0/1.")
EVENTS <- sum(y == 1)
NON_EVENTS <- sum(y == 0)
PREVALENCE <- mean(y == 1)
P <- sum(parameter_table$predictor_parameters)

if (N != 390L) stop("Expected N=390, found N=", N)
if (EVENTS != 111L) stop("Expected 111 deaths, found ", EVENTS)
if (P != 44L) stop("Expected 44 initial candidate predictor parameters, found ", P)

write.csv(parameter_table,
          paste0(OUT_PREFIX, "01_candidate_predictor_parameter_count.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("[2/5] Calculating crude EPP/EPV and conservative anticipated R2...\n")

CRUDE_EPP <- EVENTS / P

# Maximum Cox-Snell R2 for a binary outcome at the observed outcome prevalence.
# ln(L_null)/N = p*ln(p) + (1-p)*ln(1-p)
loglik_null_per_obs <- PREVALENCE * log(PREVALENCE) +
  (1 - PREVALENCE) * log(1 - PREVALENCE)
MAX_CS_R2 <- 1 - exp(2 * loglik_null_per_obs)

# Riley et al. conservative fallback where no reliable anticipated R2 is known:
# Nagelkerke R2 = 0.15, equivalently Cox-Snell R2 = 15% of maximum CS R2.
ANTICIPATED_CS_R2 <- 0.15 * MAX_CS_R2
ANTICIPATED_NAG_R2 <- ANTICIPATED_CS_R2 / MAX_CS_R2

cat("[3/5] Applying Riley minimum sample-size criteria...\n")

# Criterion: target global shrinkage >= 0.90.
# n = P / ((S - 1) * log(1 - R2_CS / S))
N_SHRINKAGE_RAW <- P / ((TARGET_SHRINKAGE - 1) *
                          log(1 - ANTICIPATED_CS_R2 / TARGET_SHRINKAGE))
N_SHRINKAGE <- ceiling(N_SHRINKAGE_RAW)

# Criterion: difference <=0.05 between apparent and adjusted Nagelkerke R2.
# Required shrinkage corresponding to delta in Nagelkerke R2:
S_OPT <- ANTICIPATED_CS_R2 /
  (ANTICIPATED_CS_R2 + DELTA_NAG_R2 * MAX_CS_R2)
N_R2_OPT_RAW <- P / ((S_OPT - 1) * log(1 - ANTICIPATED_CS_R2 / S_OPT))
N_R2_OPT <- ceiling(N_R2_OPT_RAW)

# Criterion: estimate overall outcome proportion with margin of error <= 0.05.
N_INTERCEPT_RAW <- (Z_975 / INTERCEPT_MARGIN)^2 *
  PREVALENCE * (1 - PREVALENCE)
N_INTERCEPT <- ceiling(N_INTERCEPT_RAW)

N_REQUIRED <- max(N_SHRINKAGE, N_R2_OPT, N_INTERCEPT)
DRIVING_CRITERION <- c(
  "Target global shrinkage >=0.90",
  "Nagelkerke R2 optimism <=0.05",
  "Overall outcome proportion precision (MOE <=0.05)"
)[which.max(c(N_SHRINKAGE, N_R2_OPT, N_INTERCEPT))]

EXPECTED_EVENTS_AT_N_REQUIRED <- N_REQUIRED * PREVALENCE
REQUIRED_EPP <- EXPECTED_EVENTS_AT_N_REQUIRED / P
SAMPLE_SIZE_RATIO <- N / N_REQUIRED
SHORTFALL_N <- N_REQUIRED - N

criteria_table <- data.frame(
  criterion = c(
    "Target global shrinkage >=0.90",
    "Nagelkerke R2 optimism <=0.05",
    "Overall outcome proportion precision (MOE <=0.05)"
  ),
  required_n_raw = c(N_SHRINKAGE_RAW, N_R2_OPT_RAW, N_INTERCEPT_RAW),
  required_n_ceiling = c(N_SHRINKAGE, N_R2_OPT, N_INTERCEPT),
  stringsAsFactors = FALSE
)

summary_table <- data.frame(
  item = c(
    "Observed sample size",
    "Observed events",
    "Observed non-events",
    "Observed outcome prevalence",
    "Initial candidate predictors",
    "Initial candidate predictor parameters",
    "Crude events per candidate parameter",
    "Maximum Cox-Snell R2 at observed prevalence",
    "Assumed Nagelkerke R2 for Riley assessment",
    "Assumed Cox-Snell R2 for Riley assessment",
    "Target shrinkage",
    "Required N: shrinkage criterion",
    "Required N: R2 optimism criterion",
    "Required N: overall-risk precision criterion",
    "Overall minimum required N",
    "Driving criterion",
    "Expected events at minimum required N",
    "Required events per predictor parameter",
    "Observed N / required N",
    "Sample-size shortfall"
  ),
  value = c(
    as.character(N),
    as.character(EVENTS),
    as.character(NON_EVENTS),
    sprintf("%.6f", PREVALENCE),
    as.character(length(candidate_vars)),
    as.character(P),
    sprintf("%.6f", CRUDE_EPP),
    sprintf("%.6f", MAX_CS_R2),
    sprintf("%.6f", ANTICIPATED_NAG_R2),
    sprintf("%.6f", ANTICIPATED_CS_R2),
    sprintf("%.2f", TARGET_SHRINKAGE),
    as.character(N_SHRINKAGE),
    as.character(N_R2_OPT),
    as.character(N_INTERCEPT),
    as.character(N_REQUIRED),
    DRIVING_CRITERION,
    sprintf("%.2f", EXPECTED_EVENTS_AT_N_REQUIRED),
    sprintf("%.6f", REQUIRED_EPP),
    sprintf("%.6f", SAMPLE_SIZE_RATIO),
    as.character(SHORTFALL_N)
  ),
  stringsAsFactors = FALSE
)

write.csv(criteria_table,
          paste0(OUT_PREFIX, "02_Riley_criteria.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary_table,
          paste0(OUT_PREFIX, "03_Riley_EPV_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("[4/5] Cross-checking with pmsampsize when available...\n")

pmsamp_status <- "pmsampsize package not installed; formula-based Riley results saved."
pmsamp_text <- character(0)

if (requireNamespace("pmsampsize", quietly = TRUE)) {
  pmsamp_obj <- tryCatch(
    pmsampsize::pmsampsize(
      type = "b",
      csrsquared = ANTICIPATED_CS_R2,
      parameters = P,
      prevalence = PREVALENCE
    ),
    error = function(e) e
  )

  if (inherits(pmsamp_obj, "error")) {
    pmsamp_status <- paste0("pmsampsize ERROR: ", conditionMessage(pmsamp_obj))
    pmsamp_text <- pmsamp_status
  } else {
    pmsamp_status <- "pmsampsize completed successfully."
    pmsamp_text <- capture.output(print(pmsamp_obj))
    saveRDS(pmsamp_obj, paste0(OUT_PREFIX, "04_pmsampsize_object.rds"), compress = "xz")
  }
} else {
  pmsamp_text <- c(
    pmsamp_status,
    "Optional installation command: install.packages('pmsampsize')",
    "The manuscript-facing calculation does not depend on automatic installation."
  )
}

writeLines(pmsamp_text, paste0(OUT_PREFIX, "04_pmsampsize_console_output.txt"), useBytes = TRUE)

# Explicit audit comparing our deterministic formula result with package output is
# left visible in the console; never silently overwrite either result.
cat("\n================ RILEY / EPP SUMMARY ================\n")
cat("N =", N, "\n")
cat("Events =", EVENTS, "\n")
cat("Outcome prevalence =", sprintf("%.4f", PREVALENCE), "\n")
cat("Initial candidate parameters =", P, "\n")
cat("Crude EPP/EPV =", sprintf("%.2f", CRUDE_EPP), "\n")
cat("Maximum Cox-Snell R2 =", sprintf("%.6f", MAX_CS_R2), "\n")
cat("Conservative anticipated Cox-Snell R2 (15% max) =",
    sprintf("%.6f", ANTICIPATED_CS_R2), "\n")
cat("Required N for shrinkage =", N_SHRINKAGE, "\n")
cat("Required N for R2 optimism =", N_R2_OPT, "\n")
cat("Required N for overall-risk precision =", N_INTERCEPT, "\n")
cat("Overall minimum required N =", N_REQUIRED, "\n")
cat("Driving criterion =", DRIVING_CRITERION, "\n")
cat("Observed N / required N =", sprintf("%.3f", SAMPLE_SIZE_RATIO), "\n")
cat("======================================================\n\n")

cat("[5/5] Saving complete audit object, sessionInfo and manifest...\n")

complete_object <- list(
  settings = list(
    data_file = DATA_FILE,
    outcome = OUTCOME,
    target_shrinkage = TARGET_SHRINKAGE,
    delta_nagelkerke_r2 = DELTA_NAG_R2,
    intercept_margin = INTERCEPT_MARGIN,
    anticipated_r2_strategy = "15% of maximum Cox-Snell R2 (Nagelkerke R2 = 0.15)"
  ),
  cohort = list(
    N = N,
    events = EVENTS,
    non_events = NON_EVENTS,
    prevalence = PREVALENCE
  ),
  candidate_parameter_table = parameter_table,
  candidate_predictors = candidate_vars,
  n_candidate_parameters = P,
  crude_epp = CRUDE_EPP,
  max_cox_snell_r2 = MAX_CS_R2,
  anticipated_cox_snell_r2 = ANTICIPATED_CS_R2,
  anticipated_nagelkerke_r2 = ANTICIPATED_NAG_R2,
  criteria = criteria_table,
  overall_required_n = N_REQUIRED,
  driving_criterion = DRIVING_CRITERION,
  expected_events_at_required_n = EXPECTED_EVENTS_AT_N_REQUIRED,
  required_epp = REQUIRED_EPP,
  observed_to_required_ratio = SAMPLE_SIZE_RATIO,
  sample_size_shortfall = SHORTFALL_N,
  pmsampsize_status = pmsamp_status
)

saveRDS(complete_object,
        paste0(OUT_PREFIX, "05_complete_Riley_EPV_object.rds"),
        compress = "xz")

writeLines(capture.output(sessionInfo()),
           paste0(OUT_PREFIX, "06_sessionInfo.txt"), useBytes = TRUE)

manifest_files <- list.files(
  WORK_DIR,
  pattern = "^390RILEY_",
  full.names = FALSE
)
manifest <- data.frame(
  file = manifest_files,
  stringsAsFactors = FALSE
)
write.csv(manifest,
          paste0(OUT_PREFIX, "07_output_manifest.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("DONE. No bootstrap, MICE, LASSO, or model refitting was performed.\n")
