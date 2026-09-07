# ============================================================
# 09_TABLE1_390_BASELINE_P_MISSINGNESS_FINAL.R
# Purpose:
#   Recalculate the FINAL N=390 Table 1 from the ORIGINAL
#   UNIMPUTED dataset:
#     - Overall / 28-day survivors / non-survivors
#     - Missing n (%) for every displayed variable
#     - Exploratory group-comparison P values
#     - Test method used for each variable
#     - Includes all major severity scores:
#         SAPS II, APACHE II, SOFA, OASIS, APS III
#
# IMPORTANT:
#   - Does NOT run MICE.
#   - Does NOT run LASSO.
#   - Does NOT run bootstrap.
#   - Does NOT alter the final model.
#   - P values are descriptive/exploratory only and are NOT used
#     for predictor selection.
#
# Source:
#   D:/R代码/大修1_1/390run.xlsx
# ============================================================

rm(list = ls())

# ----------------------------- PATHS -----------------------------
WORK_DIR <- "D:/R代码/大修1_1"
DATA_FILE <- file.path(WORK_DIR, "390run.xlsx")
OUT_DIR <- file.path(WORK_DIR, "390_TABLE1_FINAL")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ----------------------------- PACKAGES ---------------------------
required_pkgs <- c("readxl")
miss_pkg <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(miss_pkg) > 0L) {
  stop(
    "Missing package(s): ",
    paste(miss_pkg, collapse = ", "),
    ". Install them before running this script."
  )
}

# ----------------------------- HELPERS ----------------------------
clean_numeric <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "na", "Na", "NULL", "null", ".", "-", "--")] <- NA_character_
  suppressWarnings(as.numeric(z))
}

clean_binary <- function(x) {
  z <- trimws(as.character(x))
  out <- rep(NA_real_, length(z))

  out[z %in% c("0", "No", "NO", "no", "N", "n", "Absent", "absent")] <- 0
  out[z %in% c("1", "Yes", "YES", "yes", "Y", "y", "Present", "present")] <- 1

  num <- suppressWarnings(as.numeric(z))
  use_num <- is.na(out) & is.finite(num) & num %in% c(0, 1)
  out[use_num] <- num[use_num]

  out
}

find_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0L) return(hit[1])

  norm <- function(s) gsub("[^a-z0-9]", "", tolower(s))
  nm <- norm(names(df))
  for (cc in candidates) {
    j <- match(norm(cc), nm)
    if (!is.na(j)) return(names(df)[j])
  }

  if (required) {
    stop(
      "Could not find required column. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  NULL
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_num <- function(x, digits = 1) {
  if (length(x) != 1L || !is.finite(x)) return(NA_character_)
  format(round(x, digits), nsmall = digits, trim = TRUE, scientific = FALSE)
}

fmt_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
  paste0(
    fmt_num(q[2], 2),
    " (",
    fmt_num(q[1], 2),
    "–",
    fmt_num(q[3], 2),
    ")"
  )
}

fmt_n_pct <- function(x, denom) {
  n <- sum(x == 1, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_character_)
  paste0(n, " (", sprintf("%.1f", 100 * n / denom), "%)")
}

# Wilcoxon group comparison, robust to missing values.
wilcox_p <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  if (sum(ok) < 3L || length(unique(g[ok])) < 2L) return(NA_real_)
  tryCatch(
    stats::wilcox.test(x[ok] ~ g[ok], exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

# Chi-square/Fisher for categorical variables.
categorical_p <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  if (sum(ok) < 3L || length(unique(g[ok])) < 2L) return(NA_real_)

  tab <- table(x[ok], g[ok], useNA = "no")
  if (any(dim(tab) < 2L)) return(NA_real_)

  # Use Fisher when expected counts are small.
  use_fisher <- FALSE
  chi <- tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL)

  if (!is.null(chi)) {
    exp_min <- min(chi$expected)
    if (is.finite(exp_min) && exp_min < 5) use_fisher <- TRUE
  } else {
    use_fisher <- TRUE
  }

  if (use_fisher) {
    tryCatch(
      stats::fisher.test(tab)$p.value,
      error = function(e) {
        if (!is.null(chi)) chi$p.value else NA_real_
      }
    )
  } else {
    chi$p.value
  }
}

# ----------------------------- READ DATA ---------------------------
if (!file.exists(DATA_FILE)) {
  stop("Data file not found: ", DATA_FILE)
}

dat <- readxl::read_excel(DATA_FILE, sheet = 1)
dat <- as.data.frame(dat, check.names = FALSE)

cat("Rows read from Sheet1:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n")

# ----------------------------- OUTCOME -----------------------------
outcome_col <- find_col(
  dat,
  c("death_28days", "death28days", "death_28_days", "28day_death"),
  required = TRUE
)

outcome <- clean_binary(dat[[outcome_col]])

if (sum(!is.na(outcome)) != 390L) {
  stop(
    "Expected 390 non-missing outcome observations; found ",
    sum(!is.na(outcome)), "."
  )
}

if (sum(outcome == 1, na.rm = TRUE) != 111L) {
  stop(
    "Expected 111 28-day deaths; found ",
    sum(outcome == 1, na.rm = TRUE), "."
  )
}

if (sum(outcome == 0, na.rm = TRUE) != 279L) {
  stop(
    "Expected 279 survivors; found ",
    sum(outcome == 0, na.rm = TRUE), "."
  )
}

group <- factor(
  outcome,
  levels = c(0, 1),
  labels = c("Survivor", "Non-survivor")
)

n_surv <- sum(group == "Survivor")
n_death <- sum(group == "Non-survivor")

# ----------------------- DEFINE TABLE VARIABLES -------------------
# Display order is deliberate: severity scores first, then demographics,
# physiology, laboratories, treatments/comorbidities and ventilation.
var_specs <- list(
  list("SAPS II", "sapsii", "continuous", c("sapsii", "SAPSII", "SAPS II")),
  list("APACHE II", "apacheii", "continuous", c("apacheii", "APACHEII", "APACHE II")),
  list("SOFA", "sofa", "continuous", c("sofa", "SOFA")),
  list("OASIS", "oasis", "continuous", c("oasis", "OASIS")),
  list("APS III", "apsiii", "continuous", c("apsiii", "APSIII", "APS III")),

  list("Age, years", "age", "continuous", c("age", "Age")),
  list("Sex, male", "male", "binary", c("gender", "sex", "Gender", "Sex")),
  list("Weight, kg", "weight", "continuous", c("weight", "Weight")),

  list("Heart rate, beats/min", "hr", "continuous", c("hr", "HR", "heart_rate")),
  list("Respiratory rate, breaths/min", "rr", "continuous", c("rr", "RR", "respiratory_rate")),
  list("SpO2, %", "spo2", "continuous", c("spo2", "SpO2", "SPO2")),
  list("Temperature, °C", "temperature", "continuous", c("temperature", "Temperature", "temp")),

  list("PaCO2, mmHg", "pco2", "continuous", c("pco2", "PCO2", "PaCO2")),
  list("pH", "ph", "continuous", c("ph", "pH")),
  list("PaO2, mmHg", "po2", "continuous", c("po2", "PO2", "PaO2")),
  list("P/F ratio, mmHg", "pf_ratio", "continuous", c("pf_ratio", "P_F_ratio", "P/F Ratio", "PF_ratio")),

  list("Anion gap", "anion_gap", "continuous", c("anion_gap", "aniongap")),
  list("Bicarbonate", "bicarbonate", "continuous", c("bicarbonate", "HCO3", "hco3")),
  list("Lactate, mmol/L", "lactate", "continuous", c("lactate", "Lactate")),
  list("WBC, ×10^9/L", "wbc", "continuous", c("wbc", "WBC", "white_blood_cell")),
  list("Hemoglobin, g/dL", "hb", "continuous", c("hb", "HB", "hemoglobin")),
  list("RDW, %", "rdw", "continuous", c("rdw", "RDW")),
  list("Hematocrit, %", "hematocrit", "continuous", c("hematocrit", "HCT", "hct")),
  list("RBC, ×10^12/L", "rbc", "continuous", c("rbc", "RBC")),
  list("Platelet, ×10^9/L", "platelet", "continuous", c("platelet", "platelets", "PLT")),
  list("BUN, mg/dL", "bun", "continuous", c("bun", "BUN")),
  list("Creatinine, mg/dL", "creatinine", "continuous", c("creatinine", "Creatinine")),
  list("Sodium, mmol/L", "sodiumm", "continuous", c("sodiumm", "sodium", "Na")),
  list("Potassium, mmol/L", "potassium", "continuous", c("potassium", "K")),
  list("Magnesium, mg/dL", "magnesium", "continuous", c("magnesium", "Mg")),
  list("Calcium, mg/dL", "calcium", "continuous", c("calcium", "Ca")),
  list("Glucose, mg/dL", "glucosemg", "continuous", c("glucosemg", "glucose", "Glucose")),
  list("INR", "inr", "continuous", c("inr", "INR")),
  list("PT, s", "pt", "continuous", c("pt", "PT")),
  list("ALT, U/L", "alt", "continuous", c("alt", "ALT")),
  list("AST, U/L", "ast", "continuous", c("ast", "AST")),
  list("Total bilirubin, mg/dL", "bilirubin", "continuous", c("bilirubin", "bilirubin_total", "total_bilirubin")),

  list("Sedation/analgesia use", "sed_analg", "binary", c("sed_analg", "sedanalg", "sedation_analgesia")),
  list("Vasopressor use", "vasopressor", "binary", c("vasopressor", "vasopressor_use")),
  list("Antibiotic use", "antibiotic", "binary", c("antibiotic", "antibiotic_use")),
  list("Glucocorticoid use", "glucocorticoid", "binary", c("glucocorticoid", "glucocorticoid_use")),
  list("Immunosuppressant use", "immunosuppressant", "binary", c("immunosuppressant", "immunosuppressant_use")),
  list("VAP", "vap", "binary", c("vap", "VAP", "ventilator_associated_pneumonia")),

  list("Ventilation mode: IMV", "vent_mode", "binary", c("vent_mode", "Ventilation Mode", "ventilation_mode"))
)

# ------------------------ EXTRACT VARIABLES -----------------------
rows <- vector("list", length(var_specs))
keep_i <- 0L
missing_variable_log <- list()

for (i in seq_along(var_specs)) {
  sp <- var_specs[[i]]

  label <- sp[[1]]
  short_name <- sp[[2]]
  type <- sp[[3]]
  candidates <- sp[[4]]

  cc <- find_col(dat, candidates, required = FALSE)

  if (is.null(cc)) {
    missing_variable_log[[length(missing_variable_log) + 1L]] <- data.frame(
      requested_variable = short_name,
      label = label,
      reason = "Column not found in Sheet1",
      stringsAsFactors = FALSE
    )
    next
  }

  keep_i <- keep_i + 1L
  x_raw <- dat[[cc]]

  if (type == "continuous") {
    x <- clean_numeric(x_raw)
    overall <- fmt_median_iqr(x)
    surv_txt <- fmt_median_iqr(x[group == "Survivor"])
    death_txt <- fmt_median_iqr(x[group == "Non-survivor"])

    p <- wilcox_p(x, group)
    method <- "Wilcoxon rank-sum"
    n_available <- sum(is.finite(x))
    n_missing <- sum(!is.finite(x))
  } else {
    # Sex/gender needs special handling:
    # if source is a text variable, convert male to 1 and female to 0.
    z <- trimws(as.character(x_raw))
    zlow <- tolower(z)

    if (short_name == "male") {
      x <- rep(NA_real_, length(z))
      x[zlow %in% c("male", "m", "1", "yes", "y")] <- 1
      x[zlow %in% c("female", "f", "0", "no", "n")] <- 0
      n_av0 <- sum(!is.na(x))
      n_male_all <- sum(x == 1, na.rm = TRUE)
      n_male_surv <- sum(x[group == "Survivor"] == 1, na.rm = TRUE)
      n_male_death <- sum(x[group == "Non-survivor"] == 1, na.rm = TRUE)

      overall <- fmt_n_pct(x, n_av0)
      surv_txt <- fmt_n_pct(x[group == "Survivor"], sum(!is.na(x[group == "Survivor"])))
      death_txt <- fmt_n_pct(x[group == "Non-survivor"], sum(!is.na(x[group == "Non-survivor"])))

      p <- categorical_p(x, group)
      method <- ifelse(
        any(
          table(
            x[!is.na(x)],
            group[!is.na(x)]
          ) < 5
        ),
        "Fisher's exact",
        "Pearson χ²"
      )
      n_available <- n_av0
      n_missing <- length(x) - n_av0
    } else {
      x <- clean_binary(x_raw)

      # If this is ventilation mode, report IMV as Yes=1.
      overall <- fmt_n_pct(x, sum(!is.na(x)))
      surv_txt <- fmt_n_pct(x[group == "Survivor"], sum(!is.na(x[group == "Survivor"])))
      death_txt <- fmt_n_pct(x[group == "Non-survivor"], sum(!is.na(x[group == "Non-survivor"])))

      p <- categorical_p(x, group)

      # Determine method from expected counts.
      ok <- !is.na(x) & !is.na(group)
      tab <- if (sum(ok) > 0) table(x[ok], group[ok]) else matrix(nrow = 0, ncol = 0)
      method <- "Pearson χ²"
      if (length(tab) > 0) {
        ch <- tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL)
        if (!is.null(ch) && min(ch$expected) < 5) method <- "Fisher's exact"
      }

      n_available <- sum(!is.na(x))
      n_missing <- sum(is.na(x))
    }
  }

  keep_row <- data.frame(
    Variable = label,
    internal_name = short_name,
    source_column = cc,
    type = type,
    Overall = overall,
    Survivor = surv_txt,
    Non_survivor = death_txt,
    Missing_n = n_missing,
    Missing_pct = 100 * n_missing / nrow(dat),
    P_value_numeric = p,
    P_value = fmt_p(p),
    Test = method,
    stringsAsFactors = FALSE
  )

  rows[[keep_i]] <- keep_row
}

table1 <- do.call(rbind, rows[seq_len(keep_i)])
rownames(table1) <- NULL

if (length(missing_variable_log) > 0L) {
  not_found <- do.call(rbind, missing_variable_log)
} else {
  not_found <- data.frame(
    requested_variable = character(0),
    label = character(0),
    reason = character(0),
    stringsAsFactors = FALSE
  )
}

# ---------------------- ADD 44-CANDIDATE MISSINGNESS ---------------
# Sheet2 contains "否纳入=是" in the finalized candidate dictionary.
candidate_missingness <- NULL

if ("Sheet2" %in% readxl::excel_sheets(DATA_FILE)) {
  dict <- readxl::read_excel(DATA_FILE, sheet = "Sheet2")
  dict <- as.data.frame(dict, check.names = FALSE)

  # Try to locate the inclusion flag and candidate-name column.
  names_dict <- names(dict)

  include_col <- NULL
  for (nm in names_dict) {
    vals <- trimws(as.character(dict[[nm]]))
    if (any(vals %in% c("是", "Yes", "YES", "yes"), na.rm = TRUE)) {
      if (grepl("纳入|include|included|入模|变量", nm, ignore.case = TRUE)) {
        include_col <- nm
        break
      }
    }
  }

  if (is.null(include_col)) {
    # broader fallback: choose a column with a majority of 是/否.
    for (nm in names_dict) {
      vals <- trimws(as.character(dict[[nm]]))
      frac <- mean(vals %in% c("是", "否", "Yes", "No", "YES", "NO"), na.rm = TRUE)
      if (is.finite(frac) && frac >= 0.5) {
        include_col <- nm
        break
      }
    }
  }

  var_col <- NULL
  for (nm in names_dict) {
    if (grepl("变量|variable|predictor|name", nm, ignore.case = TRUE)) {
      var_col <- nm
      break
    }
  }
  if (is.null(var_col) && ncol(dict) >= 1L) var_col <- names_dict[1]

  if (!is.null(include_col) && !is.null(var_col)) {
    use <- trimws(as.character(dict[[include_col]])) %in% c("是", "Yes", "YES", "yes")
    candidate_names <- unique(trimws(as.character(dict[[var_col]][use])))
    candidate_names <- candidate_names[nzchar(candidate_names)]

    # Match candidate names to Sheet1.
    cm <- lapply(candidate_names, function(v) {
      cc <- find_col(dat, c(v), required = FALSE)
      if (is.null(cc)) {
        data.frame(
          candidate = v,
          source_column = NA_character_,
          available_n = NA_integer_,
          missing_n = NA_integer_,
          missing_pct = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        xx <- as.character(dat[[cc]])
        miss <- trimws(xx) %in% c("", "NA", "N/A", "na", "Na", "NULL", "null", ".", "-", "--") |
          is.na(dat[[cc]])

        data.frame(
          candidate = v,
          source_column = cc,
          available_n = sum(!miss),
          missing_n = sum(miss),
          missing_pct = 100 * mean(miss),
          stringsAsFactors = FALSE
        )
      }
    })

    candidate_missingness <- do.call(rbind, cm)

    # Preserve the finalized 44-predictor check if exactly 44 are recovered.
    candidate_missingness$candidate_count <- nrow(candidate_missingness)
  }
}

# ------------------------- FINAL SCORE SUMMARY ---------------------
score_names <- c("SAPS II", "APACHE II", "SOFA", "OASIS", "APS III")

score_summary <- table1[
  table1$Variable %in% score_names,
  c("Variable", "Overall", "Survivor", "Non_survivor",
    "Missing_n", "Missing_pct", "P_value", "Test"),
  drop = FALSE
]

# ----------------------------- SAVE --------------------------------
utils::write.csv(
  table1,
  file.path(OUT_DIR, "390_TABLE1_01_full_baseline_with_P_missingness.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

utils::write.csv(
  table1[
    ,
    c("Variable", "Overall", "Survivor", "Non_survivor",
      "Missing_n", "Missing_pct", "P_value", "Test")
  ],
  file.path(OUT_DIR, "390_TABLE1_02_publication_ready.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

utils::write.csv(
  score_summary,
  file.path(OUT_DIR, "390_TABLE1_03_severity_scores_only.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

utils::write.csv(
  not_found,
  file.path(OUT_DIR, "390_TABLE1_04_requested_but_not_found.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (!is.null(candidate_missingness)) {
  utils::write.csv(
    candidate_missingness,
    file.path(OUT_DIR, "390_TABLE1_05_candidate_pool_missingness.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

# Key audit summary
key_summary <- data.frame(
  Item = c(
    "Final cohort N",
    "28-day survivors",
    "28-day deaths",
    "28-day mortality (%)",
    "Number of Table 1 variables found",
    "Number of requested variables not found"
  ),
  Value = c(
    nrow(dat),
    n_surv,
    n_death,
    100 * n_death / nrow(dat),
    nrow(table1),
    nrow(not_found)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  key_summary,
  file.path(OUT_DIR, "390_TABLE1_06_key_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "390_TABLE1_07_sessionInfo.txt")
)

# ----------------------------- CONSOLE -----------------------------
cat("\n============================================================\n")
cat("FINAL N=390 TABLE 1 RECALCULATION COMPLETE\n")
cat("============================================================\n")
cat("N =", nrow(dat), "\n")
cat("Survivors =", n_surv, "\n")
cat("Deaths =", n_death, "\n")
cat("Mortality =", round(100 * n_death / nrow(dat), 1), "%\n")
cat("Variables found =", nrow(table1), "\n")
cat("Variables not found =", nrow(not_found), "\n\n")

cat("Severity scores found:\n")
print(score_summary[, c("Variable", "Overall", "Survivor", "Non_survivor", "P_value", "Test"),
                    drop = FALSE])

cat("\nIMPORTANT:\n")
cat("P values are exploratory descriptive comparisons only.\n")
cat("They were NOT used for predictor selection or model development.\n")
cat("No imputation, LASSO, bootstrap, or model refitting was performed.\n")

if (!is.null(candidate_missingness)) {
  cat("\nRecovered candidate-pool variables from Sheet2:",
      nrow(candidate_missingness), "\n")
  if (nrow(candidate_missingness) != 44L) {
    cat("WARNING: Sheet2 recovery did not yield exactly 44 candidates.\n")
  }
}

cat("\nOutput directory:\n", OUT_DIR, "\n")
cat("============================================================\n")
