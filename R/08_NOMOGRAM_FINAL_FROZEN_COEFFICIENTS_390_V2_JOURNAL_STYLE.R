# ============================================================
# 08_NOMOGRAM_FINAL_FROZEN_COEFFICIENTS_390_V2_JOURNAL_STYLE.R
# Journal-style nomogram based on the FINAL FROZEN coefficients.
# NO model refitting. NO MICE. NO variable reselection.
#
# Frozen model:
# logit(p) =
#   -4.68758
#   + 0.07246  * sapsii
#   + 0.05227  * rr
#   - 0.003365 * pf_ratio
#   + 1.02720  * immunosuppressant
#
# Main improvements over V1:
#   - wider landscape layout
#   - substantially increased vertical spacing
#   - fewer, cleaner tick labels
#   - dynamic thinning of labels to prevent overlap
#   - separate Total Points and Predicted Probability axes
#   - journal-like black/white nomogram appearance
#   - exact frozen coefficients retained
# ============================================================

rm(list = ls())

# ----------------------------- PATHS ----------------------------------
WORK_DIR <- "D:/R代码/大修1_1"
DATA_FILE <- file.path(WORK_DIR, "390run.xlsx")
OUT_DIR <- file.path(WORK_DIR, "390NOMOGRAM_FINAL_V2")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ----------------------------- PACKAGES -------------------------------
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install with install.packages('readxl').")
}

# ----------------------------- FROZEN MODEL ---------------------------
B0 <- -4.68758
BETA <- c(
  sapsii = 0.07246,
  rr = 0.05227,
  pf_ratio = -0.003365,
  immunosuppressant = 1.02720
)

VAR_LABELS <- c(
  sapsii = "SAPS II",
  rr = "Respiratory rate, breaths/min",
  pf_ratio = "P/F ratio",
  immunosuppressant = "Immunosuppressant use"
)

# ----------------------------- READ DATA ------------------------------
dat <- readxl::read_excel(DATA_FILE, sheet = 1)
dat <- as.data.frame(dat, check.names = FALSE)

required_vars <- names(BETA)
miss <- setdiff(required_vars, names(dat))
if (length(miss) > 0L) {
  stop("Missing variable(s): ", paste(miss, collapse = ", "))
}

clean_numeric <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "na", "Na", "NULL", "null", ".")] <- NA_character_
  suppressWarnings(as.numeric(z))
}

for (v in required_vars) dat[[v]] <- clean_numeric(dat[[v]])

if (any(!is.na(dat$immunosuppressant) &
        !dat$immunosuppressant %in% c(0, 1))) {
  stop("immunosuppressant must contain only 0/1.")
}

# ----------------------------- DISPLAY RANGES -------------------------
ranges <- list(
  sapsii = range(dat$sapsii, na.rm = TRUE),
  rr = range(dat$rr, na.rm = TRUE),
  pf_ratio = range(dat$pf_ratio, na.rm = TRUE),
  immunosuppressant = c(0, 1)
)

# A deliberately conservative number of tick labels.
make_clean_ticks <- function(rng, target_n = 5L) {
  x <- pretty(rng, n = target_n)
  x <- x[x >= rng[1] & x <= rng[2]]
  if (length(x) < 3L) {
    x <- seq(rng[1], rng[2], length.out = 4L)
  }
  unique(x)
}

ticks <- list(
  sapsii = make_clean_ticks(ranges$sapsii, 5L),
  rr = make_clean_ticks(ranges$rr, 5L),
  pf_ratio = make_clean_ticks(ranges$pf_ratio, 5L),
  immunosuppressant = c(0, 1)
)

# ----------------------------- POINT SCALE ----------------------------
contrib_limits <- lapply(names(BETA), function(v) {
  a <- BETA[[v]] * ranges[[v]][1]
  b <- BETA[[v]] * ranges[[v]][2]
  c(min(a, b), max(a, b))
})
names(contrib_limits) <- names(BETA)

spans <- vapply(contrib_limits, diff, numeric(1))
if (max(spans) <= 0 || any(!is.finite(spans))) {
  stop("Invalid predictor ranges / contribution spans.")
}

POINTS_PER_LOGIT <- 100 / max(spans)

points_for_value <- function(v, x) {
  lim <- contrib_limits[[v]]
  (BETA[[v]] * x - lim[1]) * POINTS_PER_LOGIT
}

max_points_var <- spans * POINTS_PER_LOGIT
TOTAL_POINTS_MAX <- sum(max_points_var)

MIN_TOTAL_CONTRIB <- sum(vapply(contrib_limits, function(z) z[1], numeric(1)))
LP_AT_ZERO_POINTS <- B0 + MIN_TOTAL_CONTRIB

risk_from_total_points <- function(tp) {
  stats::plogis(LP_AT_ZERO_POINTS + tp / POINTS_PER_LOGIT)
}

# ----------------------------- LABEL THINNING -------------------------
# Remove tick labels that are too close on the common points scale.
# Endpoints are preferentially retained.
thin_ticks <- function(values, point_positions, min_gap_points = 11) {
  if (length(values) <= 2L) {
    return(list(values = values, points = point_positions))
  }

  ord <- order(point_positions)
  values <- values[ord]
  point_positions <- point_positions[ord]

  keep <- rep(FALSE, length(values))
  keep[1] <- TRUE
  last_pt <- point_positions[1]

  if (length(values) > 2L) {
    for (i in 2:(length(values) - 1L)) {
      if ((point_positions[i] - last_pt) >= min_gap_points &&
          (point_positions[length(values)] - point_positions[i]) >= min_gap_points) {
        keep[i] <- TRUE
        last_pt <- point_positions[i]
      }
    }
  }

  keep[length(values)] <- TRUE

  list(
    values = values[keep],
    points = point_positions[keep]
  )
}

for (v in c("sapsii", "rr", "pf_ratio")) {
  tt <- ticks[[v]]
  pp <- points_for_value(v, tt)
  z <- thin_ticks(tt, pp, min_gap_points = 12)
  ticks[[v]] <- z$values
}

# ----------------------------- AUDIT TABLES ---------------------------
range_rows <- do.call(
  rbind,
  lapply(names(BETA), function(v) {
    tv <- ticks[[v]]
    data.frame(
      variable = v,
      label = VAR_LABELS[[v]],
      beta = BETA[[v]],
      observed_min = ranges[[v]][1],
      observed_max = ranges[[v]][2],
      displayed_tick_value = tv,
      displayed_tick_points = points_for_value(v, tv),
      stringsAsFactors = FALSE
    )
  })
)

utils::write.csv(
  range_rows,
  file.path(OUT_DIR, "390NOMOGRAM_V2_01_display_ticks_and_points.csv"),
  row.names = FALSE
)

tp_grid <- seq(0, TOTAL_POINTS_MAX, length.out = 1001L)
risk_map <- data.frame(
  total_points = tp_grid,
  predicted_28d_mortality = risk_from_total_points(tp_grid)
)

utils::write.csv(
  risk_map,
  file.path(OUT_DIR, "390NOMOGRAM_V2_02_total_points_to_risk.csv"),
  row.names = FALSE
)

# ----------------------------- DRAWING HELPERS ------------------------
fmt_num <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-8,
    format(round(x), trim = TRUE, scientific = FALSE),
    format(round(x, 1), trim = TRUE, scientific = FALSE)
  )
}

draw_axis_line <- function(y, x0, x1, lwd = 1.15) {
  segments(x0, y, x1, y, lwd = lwd)
}

draw_ticks <- function(x, y, tick_h = 0.07, lwd = 0.9) {
  segments(x, y - tick_h, x, y + tick_h, lwd = lwd)
}

# ----------------------------- NOMOGRAM -------------------------------
draw_nomogram <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  # Wide journal-like figure.
  par(
    mar = c(3.8, 15.0, 2.8, 3.0),
    xaxs = "i",
    yaxs = "i",
    family = "sans"
  )

  # All variable axes are expressed in NOMOGRAM POINTS.
  XMAX <- max(TOTAL_POINTS_MAX, 100) * 1.02

  # Large vertical gaps between rows.
  y_points <- 8.55
  y_sapsii <- 7.15
  y_rr <- 5.75
  y_pf <- 4.35
  y_immuno <- 2.95
  y_total <- 1.45
  y_risk <- 0.15

  plot(
    NA,
    xlim = c(0, XMAX),
    ylim = c(-0.35, 9.15),
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )

  label_x <- -0.035 * XMAX

  # ---------------- Top Points ruler ----------------
  draw_axis_line(y_points, 0, 100, lwd = 1.2)
  p_ticks <- seq(0, 100, by = 20)
  draw_ticks(p_ticks, y_points)
  text(
    p_ticks,
    y_points + 0.22,
    labels = p_ticks,
    cex = 0.82,
    xpd = NA
  )
  text(
    label_x, y_points,
    "Points",
    adj = 1,
    font = 2,
    cex = 0.95,
    xpd = NA
  )

  # ---------------- Predictor axes ----------------
  var_y <- c(
    sapsii = y_sapsii,
    rr = y_rr,
    pf_ratio = y_pf,
    immunosuppressant = y_immuno
  )

  for (v in names(var_y)) {
    y0 <- var_y[[v]]

    r <- ranges[[v]]
    axis_end <- sort(points_for_value(v, r))

    draw_axis_line(y0, axis_end[1], axis_end[2], lwd = 1.15)

    tv <- ticks[[v]]
    px <- points_for_value(v, tv)

    draw_ticks(px, y0)

    tick_labels <- if (v == "immunosuppressant") {
      c("No", "Yes")
    } else {
      fmt_num(tv)
    }

    # Put tick labels clearly ABOVE each axis.
    text(
      px,
      y0 + 0.24,
      labels = tick_labels,
      cex = 0.82,
      xpd = NA
    )

    text(
      label_x,
      y0,
      VAR_LABELS[[v]],
      adj = 1,
      font = 2,
      cex = 0.92,
      xpd = NA
    )
  }

  # ---------------- Total Points axis ----------------
  draw_axis_line(y_total, 0, TOTAL_POINTS_MAX, lwd = 1.25)

  tp_ticks <- pretty(c(0, TOTAL_POINTS_MAX), n = 7L)
  tp_ticks <- tp_ticks[tp_ticks >= 0 & tp_ticks <= TOTAL_POINTS_MAX]

  # Thin if needed
  ztp <- thin_ticks(tp_ticks, tp_ticks, min_gap_points = 18)
  tp_ticks <- ztp$values

  draw_ticks(tp_ticks, y_total)
  text(
    tp_ticks,
    y_total + 0.23,
    labels = fmt_num(tp_ticks),
    cex = 0.82,
    xpd = NA
  )
  text(
    label_x,
    y_total,
    "Total points",
    adj = 1,
    font = 2,
    cex = 0.95,
    xpd = NA
  )

  # ---------------- Predicted probability axis ----------------
  min_risk <- risk_from_total_points(0)
  max_risk <- risk_from_total_points(TOTAL_POINTS_MAX)

  candidate_risks <- c(
    0.02, 0.05, 0.10, 0.20, 0.30,
    0.50, 0.70, 0.80, 0.90, 0.95
  )

  risk_ticks <- candidate_risks[
    candidate_risks >= min_risk &
      candidate_risks <= max_risk
  ]

  # Map desired risks back onto the total-points scale.
  risk_pts <- (stats::qlogis(risk_ticks) - LP_AT_ZERO_POINTS) *
    POINTS_PER_LOGIT

  # Thin probability labels if needed.
  zr <- thin_ticks(risk_ticks, risk_pts, min_gap_points = 16)
  risk_ticks <- zr$values
  risk_pts <- zr$points

  draw_axis_line(y_risk, 0, TOTAL_POINTS_MAX, lwd = 1.25)
  draw_ticks(risk_pts, y_risk)

  text(
    risk_pts,
    y_risk + 0.24,
    labels = paste0(round(100 * risk_ticks), "%"),
    cex = 0.82,
    xpd = NA
  )

  text(
    label_x,
    y_risk,
    "28-day mortality probability",
    adj = 1,
    font = 2,
    cex = 0.95,
    xpd = NA
  )

  # Minimal journal-style title.
  title(
    main = "Nomogram for prediction of 28-day mortality",
    font.main = 2,
    cex.main = 1.05,
    line = 0.8
  )
}

# ----------------------------- EXPORT --------------------------------
# Wider than V1 to avoid label collisions.
png(
  file.path(OUT_DIR, "390NOMOGRAM_V2_03_journal_style_600dpi.png"),
  width = 11.5,
  height = 8.2,
  units = "in",
  res = 600
)
draw_nomogram()
dev.off()

pdf(
  file.path(OUT_DIR, "390NOMOGRAM_V2_04_journal_style.pdf"),
  width = 11.5,
  height = 8.2,
  useDingbats = FALSE
)
draw_nomogram()
dev.off()

# ----------------------------- METADATA -------------------------------
meta <- c(
  "Journal-style nomogram V2 completed.",
  "",
  "IMPORTANT:",
  "No regression model was re-fitted.",
  "The final frozen coefficients were used exactly as specified.",
  "",
  "Frozen model:",
  "logit(p) = -4.68758 + 0.07246*sapsii + 0.05227*rr - 0.003365*pf_ratio + 1.02720*immunosuppressant",
  "",
  paste0("Source data rows: ", nrow(dat)),
  paste0("Points per logit: ", signif(POINTS_PER_LOGIT, 8)),
  paste0("Maximum total points: ", signif(TOTAL_POINTS_MAX, 8)),
  paste0(
    "Attainable displayed risk range: ",
    signif(risk_from_total_points(0), 6),
    " to ",
    signif(risk_from_total_points(TOTAL_POINTS_MAX), 6)
  ),
  "",
  "Observed display ranges:",
  paste0("SAPS II: ", paste(ranges$sapsii, collapse = " to ")),
  paste0("Respiratory rate: ", paste(ranges$rr, collapse = " to ")),
  paste0("P/F ratio: ", paste(ranges$pf_ratio, collapse = " to ")),
  "Immunosuppressant: 0 to 1",
  "",
  paste0("R version: ", R.version.string),
  paste0("readxl version: ", as.character(utils::packageVersion("readxl")))
)

writeLines(
  meta,
  file.path(OUT_DIR, "390NOMOGRAM_V2_05_run_metadata.txt")
)

cat("\n============================================================\n")
cat("JOURNAL-STYLE NOMOGRAM V2 COMPLETE\n")
cat("Output directory:\n", OUT_DIR, "\n\n")
cat("Main figure:\n")
cat("390NOMOGRAM_V2_03_journal_style_600dpi.png\n")
cat("390NOMOGRAM_V2_04_journal_style.pdf\n\n")
cat("No model was re-fitted.\n")
cat("============================================================\n")
