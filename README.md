# Early 28-day mortality prediction in mechanically ventilated patients with ILD

This repository contains reproducibility materials for the manuscript:

"Early 28-day mortality prediction in mechanically ventilated patients with interstitial lung disease"

## Data source

The study used MIMIC-IV v3.1.

Patient-level MIMIC-IV data are not redistributed in this repository because access is governed by the PhysioNet data use agreement.

Credentialed users can obtain MIMIC-IV directly from PhysioNet.

## Study population

The analytical cohort included adult patients who:

- were aged 18 years or older;
- had a diagnosis of interstitial lung disease;
- were in their first ICU admission;
- received non-invasive ventilation or invasive mechanical ventilation within the first 24 hours after ICU admission.

The prediction landmark was 24 hours after ICU admission.

Patients who died before the 24-hour landmark were excluded.

Final analytical cohort:

- N = 390
- 28-day deaths = 111
- survivors = 279

## Candidate predictors

A total of 44 candidate predictor parameters were considered.

## Primary model-development strategy

1. Multiple imputation by chained equations
2. Five imputed datasets
3. Predictive mean matching
4. Logistic LASSO
5. 10-fold cross-validation
6. lambda.1se
7. Predictor retained only if selected in 5/5 imputed datasets
8. Final ordinary logistic regression
9. Rubin-rule pooling

## Internal validation

Complete-pipeline bootstrap validation was performed using 1,000 bootstrap resamples.

Within each bootstrap resample, the complete development process was repeated, including:

- imputation
- LASSO tuning
- predictor selection
- logistic refitting
- prediction
- performance assessment

## Benchmark analyses

The final model was compared with:

- SAPS II
- APACHE II
- SOFA
- OASIS

Decision-curve analysis and calibration analyses were also performed.

## Machine-learning comparison

A repeated nested cross-validation framework was used.

Algorithms included:

- primary LASSO-logistic strategy
- elastic net
- random forest
- XGBoost
- radial-basis-function support vector machine

All algorithms started from the same 44 candidate predictors and used matched outer resampling partitions.

## Data extraction

Data extraction was performed using the DecisionLinnc 1.0 Software (https://www.statsape.com).

The investigators specified cohort eligibility criteria, variables, and temporal windows through the graphical interface. The platform executes SQL-based queries in the background.Some SQL statements have been included in the cohort_definition.


## Reproducibility

Random seeds, analysis settings, R scripts, and aggregated analysis outputs are provided.

Patient-level data, imputed datasets, patient-level predictions, and other MIMIC-IV-derived individual-level files are not redistributed.

## Software

Analyses were conducted in R.

See the metadata folder for session information and package versions.
