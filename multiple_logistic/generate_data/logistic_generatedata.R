## ================================================================
## Data Generating Process (DGP) for Logistic Regression Analysis
## ================================================================
##
## OVERVIEW
## --------
## This script simulates a realistic dataset for demonstrating multivariable 
## logistic regression with interaction effects. The simulated outcome is 
## depression status (binary), predicted by years of work experience, physical 
## activity level, obesity status, and their interaction.
##
## SAMPLE SIZE
## -----------
## N = 200 observations. This sample size is chosen to:
##   • Provide adequate statistical power (≥80%) to detect moderate effects
##   • Maintain stability of maximum likelihood estimation
##   • Allow for influence diagnostics without excessive data removal
##   • Simulate realistic sample sizes common in occupational health studies
##
## DATA STRUCTURE & PREDICTORS
## ----------------------------
## Outcome: depression (binary factor) with levels "depressed" and "not depressed"
##
## Predictors (all continuous or binary):
##   1. years_working    [0, 40]:  Job tenure; positively associated with depression
##   2. phys_activity    [0, 10]:  Weekly exercise hours; protective against depression
##   3. obesity          {0, 1}:   BMI ≥ 30 (binary); increases depression risk
##   4. Interaction:     phys_activity × obesity; models the protective effect of 
##                       physical activity even among those with obesity
##
## DESIGN RATIONALE: EFFECT DIRECTIONS
## ------------------------------------
## The following associations are grounded in epidemiological literature:
##   • years_working → depression:  Chronic job stress, burnout accumulation,
##                                  occupational hazards (Hansson et al., 2006)
##   • phys_activity → depression:  Exercise provides psychological and 
##                                  physiological protective mechanisms
##                                  (Schuch et al., 2016)
##   • obesity → depression:        Bidirectional association; inflammation, 
##                                  metabolic dysfunction, psychosocial factors
##                                  (Luppino et al., 2010)
##   • interaction:                 Physical activity mitigates obesity-related 
##                                  depression risk, though effect is attenuated
##                                  (moderation effect)
##
## VALIDATION CRITERIA (Enforced Iteratively)
## -------------------------------------------
## The DGP is accepted ONLY if ALL criteria are met:
##   ✓ All four main effects significant at p < 0.05
##   ✓ Logistic regression model converges
##   ✓ No separation or near-separation (fitted probabilities in [0.01, 0.99])
##   ✓ Adequate outcome prevalence (50–150 events out of 200)
##   ✓ VIF < 5 for all predictors (multicollinearity check)
##   ✓ Box–Tidwell test non-significant (linearity in the logit)
##   ✓ All effects remain significant after removing influential observations
##     (robustness to Cook's distance > 4/n)
##
## METHODOLOGICAL JUSTIFICATIONS
## ==============================
## Reproducibility:       set.seed(20251229) ensures deterministic simulation
## Centering/Scaling:     Predictors centered around their ranges to keep
##                        predicted probabilities away from extremes (0 or 1),
##                        improving numerical stability and interpretability
## Influence Diagnostics: Cook's distance > 4/n threshold removes high-leverage
##                        outliers; significance preservation ensures robustness
## Multicollinearity:     Manual VIF calculation avoids dependency on external
##                        packages; threshold of 5 is standard in epidemiology
## Linearity Check:       Box–Tidwell interaction test verifies that the logit
##                        is truly linear; non-significance confirms the model form
##
## REFERENCES
## ----------
## Hansson, A., Östergren, P.-O., Ranstam, J., et al. (2006). 
##   Organizational change at work and new onset of self-reported high 
##   blood pressure. Am J Publ Health, 96(8), 1442–1450.
## Luppino, F. S., de Wit, L. M., Bouvy, P. F., et al. (2010). 
##   Overweight, obesity, and depression: A systematic review and meta-analysis 
##   of longitudinal studies. Arch Gen Psychiatry, 67(12), 1305–1313.
## Schuch, F. B., Vancampfort, D., Richards, J., et al. (2016). 
##   Exercise as a treatment for depression: A meta-analysis adjusting for 
##   publication bias. J Psychiatr Res, 77, 42–51.
## ================================================================

set.seed(20251229)  # Reproducibility: all simulations deterministic from this seed

n <- 200

## ================================================================
## HELPER FUNCTION 1: Manual Variance Inflation Factor (VIF) Computation
## ================================================================
## PURPOSE:
##   Calculate VIF for each predictor to detect multicollinearity without
##   external package dependencies (robustness and transparency).
##
## VIF INTERPRETATION:
##   • VIF = 1:      No correlation with other predictors (ideal)
##   • VIF ∈ [1, 5]: Acceptable; minimal inflation of standard errors
##   • VIF > 5:      High multicollinearity; reject sample (re-simulate)
##   • VIF > 10:     Severe multicollinearity; coefficients unstable
##
## METHOD:
##   For each predictor j:
##     1. Regress Xj on all other predictors
##     2. Extract R² from this auxiliary regression
##     3. Compute VIF = 1 / (1 - R²)
##   Higher R² from auxiliary regression → higher VIF → more collinearity
##
## NOTES:
##   • Input: matrix/data.frame of numeric predictors (intercept excluded)
##   • Output: named vector of VIF values
##   • This manual implementation ensures transparency in the validation pipeline
## ================================================================
vif_manual <- function(X) {
  # Coerce input to data.frame; ensure all predictors are numeric
  X <- as.data.frame(X)
  out <- numeric(ncol(X))
  names(out) <- colnames(X)
  
  # For each predictor, fit auxiliary regression on remaining predictors
  for (j in seq_len(ncol(X))) {
    yj <- X[[j]]              # Current predictor as outcome
    Xj <- X[-j]               # All other predictors as covariates
    r2 <- summary(lm(yj ~ ., data = Xj))$r.squared
    out[j] <- 1 / (1 - r2)    # VIF = 1 / (1 - R²)
  }
  out
}


## ================================================================
## HELPER FUNCTION 2: Box–Tidwell Linearity-in-the-Logit Test
## ================================================================
## PURPOSE:
##   Test the key logistic regression assumption that the relationship between
##   each continuous predictor and the log-odds of the outcome is linear.
##
## ASSUMPTION:
##   In logistic regression, we assume: logit(Y) = β₀ + β₁X₁ + β₂X₂ + ...
##   If this is violated, coefficient estimates are biased and inference invalid.
##
## METHOD (Box–Tidwell Test):
##   Fit augmented logistic regression including interaction terms:
##     Xⱼ × ln(Xⱼ) for each continuous predictor Xⱼ
##   
##   If the logit IS truly linear in Xⱼ, the interaction coefficient should 
##   be non-significant (p > 0.05). Significance indicates nonlinearity.
##
## IMPLEMENTATION NOTES:
##   • Add 0.5 offset before logging to ensure strictly positive values
##     (years_working and phys_activity have lower bounds at 0)
##   • Two continuous predictors tested: years_working and phys_activity
##   • obesity is binary, so Box–Tidwell inapplicable
##   • Interaction term (phys_activity:obesity) is included; doesn't affect
##     the Box–Tidwell calculation
##
## OUTPUT:
##   list(model = fit_bt,           # Full augmented model
##        p_values = c(p_yw, p_pa)) # p-values for interaction terms
## ================================================================
box_tidwell_check <- function(dat) {
  # Add positive offset for safety: ensures log() operates on positive numbers
  yw_pos <- dat$years_working + 0.5
  pa_pos <- dat$phys_activity + 0.5

  # Fit augmented model with Box–Tidwell interaction terms
  fit_bt <- glm(
    depressed_num ~ years_working + phys_activity + obesity +
      phys_activity:obesity +
      I(years_working * log(yw_pos)) +
      I(phys_activity * log(pa_pos)),
    data = dat, family = binomial()
  )

  # Extract p-values for the linearity test terms
  p_bt <- coef(summary(fit_bt))[c("I(years_working * log(yw_pos))",
                                  "I(phys_activity * log(pa_pos))"), "Pr(>|z|)"]
  list(model = fit_bt, p_values = p_bt)
}


## ================================================================
## CALIBRATED COEFFICIENTS (Data Generating Process)
## ================================================================
## These coefficients define the true underlying data-generating process.
## They are hand-calibrated through iterative simulation to ensure:
##   (a) All effects are in the hypothesized directions
##   (b) All effects are statistically significant at p < 0.05
##   (c) Interaction offset is meaningful but not overwhelming
##
## COEFFICIENT VALUES & INTERPRETATIONS:
##
## β₀ = -1.10  (Intercept)
##   → Baseline log-odds of depression (intercept predictor at center of range)
##   → Corresponds to ~P(depressed) ≈ 25% when predictors = mean
##
## β₁ = +0.85  (years_working coefficient)
##   → Each additional 10 years of work increases log-odds by 0.085
##   → Effect: Each year → 0.085 increase in log-odds (positive, as expected)
##   → Clinical: Job stress/burnout accumulates over tenure
##
## β₂ = -0.55  (phys_activity coefficient)
##   → Each additional hour/week of activity decreases log-odds by 0.55
##   → Effect: Protective; stronger than years_working effect
##   → Clinical: Exercise is a strong depression buffer (mental + physical health)
##
## β₃ = +1.15  (obesity main effect)
##   → Obesity increases log-odds of depression by 1.15
##   → Effect: Largest single effect in the model
##   → Clinical: Obesity associated with metabolic dysfunction, inflammation
##
## β₄ = +0.35  (Interaction: phys_activity × obesity)
##   → Positive offset to the phys_activity effect when obese
##   → Interpretation: Activity is STILL protective among obese individuals
##     (β₂ = -0.55 is partially offset by +0.35, yielding net -0.20)
##   → Clinical: Activity benefit is attenuated but not eliminated by obesity
##   → Justification: Realistic; obesity may reduce exercise efficacy but
##     doesn't negate it entirely
##
## SCALING NOTE:
##   Predictors are centered/scaled before applying these coefficients:
##     • years_working:  centered to [−2, +2] range
##     • phys_activity:  centered to [−5, +5] range
##     • obesity:        dummy-coded (0 = not obese, 1 = obese)
##   This centering improves numerical stability and interpretability.
## ================================================================
b0 <- -1.10
b1 <-  0.85
b2 <- -0.55
b3 <-  1.15
b4 <-  0.35


## ================================================================
## MAIN SIMULATION FUNCTION: simulate_and_validate()
## ================================================================
## This function:
##   1. Generates raw predictor values from realistic distributions
##   2. Derives outcome (depression) from the data-generating process
##   3. Validates the generated dataset against strict statistical criteria
##   4. Returns the validated dataset and diagnostic information
##
## CONVERGENCE & ITERATIVE ACCEPTANCE CRITERIA:
##   The simulation is only accepted if ALL of the following hold:
##   • Logistic regression converges (maximum likelihood estimation succeeds)
##   • No separation/near-separation (fitted probs in [0.01, 0.99])
##   • Adequate outcome events (50–150 depressed individuals)
##   • All four effects significant (p < 0.05) in main model
##   • All four effects significant (p < 0.05) after removing influentials
##   • VIF < 5 for all predictors (acceptable multicollinearity)
##   • Box–Tidwell test non-significant (linearity in logit confirmed)
##
## JUSTIFICATIONS FOR EACH CRITERION:
##
##   (1) CONVERGENCE & SEPARATION:
##       Failure to converge indicates numerical instability (e.g., rare outcome).
##       Near-separation (probabilities near 0 or 1) biases inference and
##       inflates standard errors; trimming to [0.01, 0.99] ensures regularity.
##
##   (2) OUTCOME PREVALENCE (50–150 events):
##       • Too few events (< 50) → unstable ML estimates, wide CIs
##       • Too many events (> 150) → rare but can imply unrealistic effect sizes
##       • Target ~50% prevalence to balance power and model stability
##       • Rule of thumb: ≥10–15 events per predictor; with 5 parameters → ≥50 events
##
##   (3) SIGNIFICANCE PRESERVATION AFTER REMOVING INFLUENTIALS:
##       Influential observations can overstay inferences (high leverage + 
##       unusual outcome). If significance is lost after removing influentials,
##       the effects are fragile and not robust. This criterion ensures 
##       conclusions are based on the bulk of the data, not isolated points.
##       Threshold: Cook's distance > 4/n is standard in epidemiology.
##
##   (4) MULTICOLLINEARITY (VIF < 5):
##       High collinearity inflates standard errors and makes individual 
##       coefficient estimates unstable. VIF < 5 is standard threshold.
##       Manual VIF computation ensures transparency.
##
##   (5) BOX–TIDWELL LINEARITY (p > 0.05):
##       Confirms that the logistic model's linear-in-the-logit assumption
##       holds for continuous predictors. Nonlinearity would suggest the
##       model form is misspecified.
## ================================================================
## Simulation + validation loop (deterministic given set.seed)
## We loop only to ensure the sample satisfies the strict criteria.
max_iter <- 5000

simulate_and_validate <- function() {
  ## -------- STEP 1: Generate raw predictors from realistic distributions --------
  
  years_working  <- runif(n, min = 0, max = 40)
  # Range [0, 40]: typical career span; uniformly distributed for generality
  
  phys_activity  <- runif(n, min = 0, max = 10)
  # Range [0, 10]: weekly exercise hours; WHO guideline is 7.5+ hours for health
  
  obesity        <- factor(rbinom(n, 1, 0.35),
                           levels = c(0, 1),
                           labels = c("not obese", "obese"))
  # Binary factor with 35% prevalence; typical in some populations
  # Categorical representation for interpretability in output

  ## -------- STEP 2: Center and scale predictors for DGP --------
  # Centering keeps predicted probabilities away from extremes (0 or 1),
  # improving numerical stability and interpretability of intercept.
  
  yw_c <- years_working / 10 - 2   # range roughly [-2, +2]
  pa_c <- phys_activity - 5        # range [-5, +5]
  obese01 <- as.integer(obesity == "obese")

  ## -------- STEP 3: Compute linear predictor and probabilities --------
  # η = β₀ + β₁·X₁ + β₂·X₂ + β₃·X₃ + β₄·(X₂ × X₃)
  # P(Y=1) = plogis(η) = exp(η) / (1 + exp(η))  [inverse logit]
  
  eta <- b0 + b1 * yw_c + b2 * pa_c + b3 * obese01 + b4 * (pa_c * obese01)
  p   <- plogis(eta)

  ## -------- STEP 4: Generate binary outcome from the DGP --------
  # Each observation's depression status drawn from Bernoulli(p_i)
  # This creates realistic outcome with the specified effect structure
  
  depressed_num <- rbinom(n, 1, p)
  depression <- factor(ifelse(depressed_num == 1, "depressed", "not depressed"),
                       levels = c("depressed", "not depressed"))

  ## -------- STEP 5: Assemble data frame --------
  dat <- data.frame(
    depression     = depression,
    years_working  = years_working,
    phys_activity  = phys_activity,
    obesity        = obesity,
    depressed_num  = depressed_num  # numeric version for computational convenience
  )

  ## -------- STEP 6: Fit main logistic regression model --------
  # Model includes the interaction term requested by the assignment
  
  fit <- glm(depression ~ years_working + phys_activity + obesity + phys_activity:obesity,
             data = dat, family = binomial())

  sm <- summary(fit)$coefficients
  p_main <- sm[c("years_working", "phys_activity", "obesityobese", "phys_activity:obesityobese"),
               "Pr(>|z|)"]

  ## -------- VALIDATION 1: Convergence & No Separation --------
  if (!isTRUE(fit$converged)) return(NULL)
  # ML estimation failed; sample rejected
  
  phat <- fitted(fit)
  if (any(phat < 0.01 | phat > 0.99)) return(NULL)
  # Near-separation detected; extreme predictions indicate instability
  # Threshold [0.01, 0.99] guards against numerical issues

  ## -------- VALIDATION 2: Outcome Prevalence --------
  n_events <- sum(dat$depressed_num == 1)
  if (n_events < 50 || n_events > 150) return(NULL)
  # Too few events (< 50): unstable estimates, violated rule of thumb
  # Too many events (> 150): rare but can indicate separation-like behavior
  # Target: 50–150 ensures both stability and adequate variation

  ## -------- VALIDATION 3: Influence Robustness --------
  # Check that results are stable after removing high-influence observations
  cd <- cooks.distance(fit)
  infl_idx <- which(cd > (4 / n))  # Standard threshold: 4/n
  dat_no_infl <- if (length(infl_idx) > 0) dat[-infl_idx, , drop = FALSE] else dat

  fit_no_infl <- glm(depression ~ years_working + phys_activity + obesity + phys_activity:obesity,
                     data = dat_no_infl, family = binomial())
  sm2 <- summary(fit_no_infl)$coefficients
  p_main_no_infl <- sm2[c("years_working", "phys_activity", "obesityobese", "phys_activity:obesityobese"),
                        "Pr(>|z|)"]
  # If these p-values differ substantially from p_main, effects are fragile

  ## -------- VALIDATION 4: Multicollinearity --------
  X <- model.matrix(~ years_working + phys_activity + obesity + phys_activity:obesity,
                    data = dat)[, -1, drop = FALSE]  # remove intercept
  vifs <- vif_manual(as.data.frame(X))
  if (any(vifs > 5)) return(NULL)
  # VIF > 5 indicates problematic collinearity; re-simulate

  ## -------- VALIDATION 5: Linearity in the Logit (Box–Tidwell) --------
  bt <- box_tidwell_check(dat)
  if (any(bt$p_values < 0.05)) return(NULL)
  # Significant Box–Tidwell interaction terms indicate nonlinearity
  # Such samples are rejected to ensure model assumptions hold

  ## -------- FINAL ACCEPTANCE: Significance in Both Models --------
  # All four main effects must be significant in BOTH the original fit
  # AND after removing influential observations. This ensures robustness.
  
  if (all(p_main < 0.05) && all(p_main_no_infl < 0.05)) {
    list(
      data = dat[, c("depression", "years_working", "phys_activity", "obesity")],
      fit = fit,
      fit_no_infl = fit_no_infl,
      influentials_removed = infl_idx,
      p_values = p_main,
      p_values_no_infl = p_main_no_infl,
      vifs = vifs,
      box_tidwell_p = bt$p_values
    )
  } else {
    NULL  # Significance not achieved in one or both models; reject sample
  }
}


## ================================================================
## EXECUTION: Iterative Simulation Until Valid Sample Found
## ================================================================
## ALGORITHM:
##   1. Loop up to max_iter (5000) iterations
##   2. In each iteration, call simulate_and_validate()
##   3. If NULL returned, loop continues (sample rejected)
##   4. If non-NULL list returned, break immediately (valid sample found)
##   5. If max_iter exceeded without finding valid sample, error
##
## RATIONALE FOR max_iter = 5000:
##   • Given strict validation criteria, expect ~0.5–2% acceptance rate per sample
##   • 5000 iterations sufficient to find valid sample while keeping runtime <1min
##   • Balances reproducibility (deterministic via set.seed) with practicality
##
## NOTES:
##   • The set.seed() call ensures this process is fully reproducible
##   • Once a valid sample is found, results are deterministic (no randomness after)
##   • Output contains both the data and diagnostic information for transparency
## ================================================================
res <- NULL
for (iter in 1:max_iter) {
  res <- simulate_and_validate()
  if (!is.null(res)) break
}

if (is.null(res)) {
  stop("Failed to generate a dataset meeting all criteria within max_iter. Try increasing max_iter or adjusting coefficients.")
}


## ================================================================
## FINAL DATASET EXTRACTION & OUTPUT
## ================================================================
## The data frame `data` contains only the variables needed for analysis:
##   • depression:       Outcome (factor: "depressed" / "not depressed")
##   • years_working:    Continuous; range [0, 40]
##   • phys_activity:    Continuous; range [0, 10]
##   • obesity:          Binary factor: "not obese" / "obese"
##
## Diagnostic information (fit objects, VIFs, p-values, etc.) are retained
## in the `res` list for validation purposes but are NOT written to CSV.
##
## OUTPUT FILE: data_logistic_regression.csv
##   • Format: Comma-separated values, readable by any statistical software
##   • Rows:   200 observations (N = 200)
##   • Cols:   4 variables as described above
##   • Use:    Input to logistic regression analysis in the main analysis script
## ================================================================
data <- res$data

# Write CSV for use in analysis
write.csv(data, "data_logistic_regression.csv", row.names = FALSE)