# ============================================================
# DATA GENERATION SCRIPT: Quality of Life (QoL) Dataset
# ============================================================
# PURPOSE:
#   Generate a synthetic dataset (n=200) for linear regression analysis
#   examining the relationship between Quality of Life and multiple
#   predictors including years of work experience, physical activity,
#   and obesity status.
#
# DATASET CHARACTERISTICS:
#   - Sample size: 200 observations
#   - Outcome variable: QoL scores (0-100 scale)
#   - Covariates:
#     * years_working: 0-40 years (continuous, beta-distributed)
#     * phys_activity: 0-10 units (continuous, beta-distributed)
#     * obesity: binary factor (obese/not obese, 35%/65% prevalence)
#   - Interaction term: phys_activity × obesity status
#   - Diagnostic features: 6 influential observations injected for
#     detecting leverage and outlier detection capabilities
#
# OUTPUT: qol_data.csv (comma-separated values file)
# ============================================================

rm(list = ls())
set.seed(42)  # Ensures reproducibility of random number generation

n <- 200  # Sample size for the synthetic dataset

# ============================================================
# COVARIATE GENERATION
# ============================================================
# These predictors simulate realistic distributions with bounded ranges.
# Beta distribution (scaled) provides natural clustering toward extremes,
# reflecting real-world data patterns.

# Years of work experience: uniformly distributed across 0-40 year range
# Beta(2,2) produces symmetric distribution suitable for modeling tenure
years_working  <- rbeta(n, 2, 2) * 40

# Physical activity level: continuous measure on 0-10 scale
# Represents frequency/intensity of exercise behavior
phys_activity  <- rbeta(n, 2, 2) * 10

# Obesity status: categorical variable with realistic prevalence
# 35% prevalence of obesity (obese) and 65% not obese
# Reflects epidemiological patterns in typical populations
obesity        <- sample(c("obese", "not obese"),
                        size = n, replace = TRUE,
                        prob = c(0.35, 0.65))

# Binary coding of obesity: required for interaction term calculation
# (1 = obese, 0 = not obese)
obese_bin <- ifelse(obesity == "obese", 1, 0)

# ============================================================
# LINEAR MODEL WITH INTERACTION: QoL OUTCOME GENERATION
# ============================================================
# Model specification: QoL = β₀ + β₁(years) + β₂(activity) + β₃(obesity) 
#                            + β₄(activity × obesity) + ε
#
# Interpretation of parameters:
#   β₀ = 72      : Intercept (baseline QoL score)
#   β₁ = -0.6    : Years of work → slightly decreases QoL
#   β₂ = 4.2     : Physical activity → substantially improves QoL
#   β₃ = -18.0   : Obesity status → major negative effect on QoL
#   β₄ = -1.8    : Interaction: activity benefit diminished among obese
#   σ = 8        : Residual standard deviation (measurement noise)

beta0    <- 72      # Intercept (baseline QoL)
beta_yr  <- -0.6    # Coefficient for years_working
beta_pa  <-  4.2    # Coefficient for phys_activity
beta_ob  <- -18.0   # Coefficient for obesity status
beta_int <- -1.8    # Interaction coefficient (activity × obesity)
sigma    <-  8      # Standard deviation of random error term

# Generate random error from normal distribution (N(0, σ²))
err <- rnorm(n, mean = 0, sd = sigma)

# Calculate linear predictor plus error
qol <- beta0 +
  beta_yr  * years_working +
  beta_pa  * phys_activity +
  beta_ob  * obese_bin +
  beta_int * phys_activity * obese_bin +
  err

# Bound QoL scores to valid range [0, 100]
# pmin ensures no score exceeds 100; pmax ensures no score falls below 0
qol <- pmin(pmax(qol, 0), 100)

# ============================================================
# DATA FRAME ASSEMBLY AND PREPARATION
# ============================================================
# Assemble generated variables into a single data frame
# Round continuous variables to 1 decimal place for realistic presentation

data <- data.frame(
  qol = round(qol, 1),                  # Quality of Life outcome (0-100)
  years_working = round(years_working, 1),  # Work experience in years
  phys_activity = round(phys_activity, 1),  # Physical activity score
  obesity = obesity,                    # Categorical obesity status
  stringsAsFactors = FALSE
)

# ============================================================
# INJECTION OF INFLUENTIAL OBSERVATIONS
# ============================================================
# PURPOSE:
#   Add 6 strategically-positioned observations with extreme values to
#   serve as test cases for diagnostic regression techniques (leverage,
#   Cook's distance, residual analysis, DFFITS, etc.)
#
# PROPERTIES OF INJECTED OBSERVATIONS:
#   - High leverage: Extreme values on predictor variables
#   - Large residuals: QoL values inconsistent with regression pattern
#   - Combination creates detectability for outlier identification
#   - Simulates real-world data anomalies and measurement errors
#
# OBSERVATIONS:
#   1. High QoL (95) despite obesity & low activity → positive outlier
#   2. Low QoL (15) despite high activity → negative outlier
#   3. Extreme combination: high QoL, extreme tenure, extreme activity
#   4-6. Additional strategic placements for robustness testing

influentials <- data.frame(
  qol = c(95, 15, 90, 10, 5, 98),
  years_working = c(40, 0, 39.5, 38.8, 1.0, 0.5),
  phys_activity = c(0, 10, 9.8, 0.2, 9.9, 0.1),
  obesity = c("obese", "not obese", "obese", "not obese", "obese", "not obese"),
  stringsAsFactors = FALSE
)

# Replace first 6 observations with influential cases
data[1:nrow(influentials), ] <- influentials

# ============================================================
# RANDOMIZATION AND FINALIZATION
# ============================================================
# Shuffle rows to distribute influential observations randomly throughout
# the dataset rather than clustering them at the beginning. This creates
# a more realistic data arrangement for blind diagnostic analysis.

set.seed(123)  # Set seed for reproducible shuffling
data <- data[sample(seq_len(n)), ]  # Random permutation of row indices
row.names(data) <- NULL  # Reset row names to sequential integers

# ============================================================
# OUTPUT EXPORT
# ============================================================
# Save the complete synthetic dataset to CSV format
# No row names exported; includes header row with variable names
write.csv(data, "qol_data.csv", row.names = FALSE)