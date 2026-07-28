# R Packages for RSF/SDM with Caribou Movement Data

Since you're working with GPS telemetry data and want to do multi-scale RSF (2nd/3rd order, à la Johnson's hierarchy) rather than SSFs, here's how I'd structure the toolkit by workflow stage:

## 1. Movement data handling & availability sampling
- **`amt`** (Animal Movement Tools) — Even without SSFs, `amt` is excellent for RSF workflows. It handles track creation from telemetry, resampling to regular intervals, and generating "available" points via random sampling within a defined domain (home range or study area). Its `random_points()` function works cleanly for both 2nd-order (study-area-wide) and 3rd-order (within-home-range) designs.
- **`adehabitatHR`** — The standard for estimating home ranges (KDE, MCP, LoCoH, Brownian bridge) which you'll need to define the "available" domain for 3rd-order selection (individual home range vs. landscape).
- **`ctmm`** — More rigorous alternative to `adehabitatHR` for utilization distributions; accounts for autocorrelation in GPS fixes (important with caribou data collected at short intervals), giving better-calibrated home range estimates (AKDE) to define availability.

## 2. Spatial covariates
- **`terra`** (replacing the retired `raster`) — Extracting habitat covariates (vegetation, terrain, disturbance layers) at used/available points.
- **`sf`** — Vector data handling (e.g., cutblocks, roads, seismic lines — all relevant disturbance covariates in caribou RSF studies).

## 3. Model fitting
- **`ResourceSelection`** — Purpose-built for RSF/RSPF; includes the Lele & Keim RSPF estimator when you have true absence/presence rather than used-available data, plus RSF fitting with weighted logistic regression.
- **`glmmTMB`** or **`lme4`** — For mixed-effects RSF, essential with caribou data since you typically have repeated observations per individual/herd. The standard approach (Fieberg et al. 2010; Muff et al. 2020) fits weighted logistic regression with a random intercept per individual and random slopes for key covariates to account for individual variation in selection strength.
- **`INLA`** (R-INLA, install from their repo, not CRAN) — Widely used in the caribou RSF literature (Muff, Nielsen, DeCesare) for Bayesian mixed-effects RSF with spatial random effects, especially when you want to properly propagate uncertainty or handle large random-effect structures more efficiently than `glmmTMB`.

## 4. Model selection & multi-scale RSF
- **`MuMIn`** — AIC-based model selection/averaging across candidate RSF models, useful when comparing 2nd- vs 3rd-order model sets or covariate scales.
- **`AICcmodavg`** — Alternative for model averaging with clean output tables, common in wildlife ecology papers.

## 5. Validation & prediction
- **`ResourceSelection::hoslem.test`** or k-fold cross-validation (Boyce index) — often coded manually per Boyce et al. 2002; some functions available in **`ecospat`** if you want the continuous Boyce index for RSF validation.
- **`terra`** again for predicting the fitted RSF surface back onto the raster stack to produce your final habitat suitability/RSF map.

## Suggested core stack
For a clean, current pipeline: **`amt` + `adehabitatHR`/`ctmm` + `terra`/`sf` + `glmmTMB` (or `INLA`) + `MuMIn`**

This mirrors the workflow in most recent caribou RSF papers (e.g., DeCesare et al. 2012 multi-scale approach, Muff et al. 2020 mixed-effects tutorial), and cleanly separates movement/availability definition from model fitting and selection.

A quick question to help me tailor further: are you fitting used-available (weighted logistic regression) or do you have true absence data, and do you want frequentist (`glmmTMB`) or Bayesian (`INLA`) mixed models for the individual-level random effects?


Good — that's the most common setup for caribou RSF work, and it has a well-established recipe (Fieberg et al. 2010; Muff, Signer & Fieberg 2020). Here's how to implement it in `glmmTMB`.

## Core approach: weighted logistic regression with fixed large weight on availability

The trick that makes used-available data behave like a valid likelihood for RSF is assigning a large fixed weight (typically 1000) to available (0) points, while used (1) points get weight 1. This approximates the inhomogeneous Poisson point process likelihood underlying RSF (Fithian & Hastie 2013; Muff et al. 2020).

```r
library(amt)
library(glmmTMB)

# assume 'dat' has columns: case_ (1=used, 0=available), animal_id, covariates...
dat$weight <- ifelse(dat$case_ == 1, 1, 1000)

# Random intercept + random slopes per individual
m <- glmmTMB(
  case_ ~ elevation + canopy_cover + dist_to_road + tree_cover +
    (1 | animal_id) +
    (0 + elevation | animal_id) +
    (0 + canopy_cover | animal_id) +
    (0 + dist_to_road | animal_id),
  family = binomial(),
  weights = weight,
  data = dat
)
```

## Key implementation details specific to this method

**1. Fix the random intercept variance at a large value.** This is the part people miss. Because the intercept absorbs the used:available sampling ratio (not a real biological effect), Muff et al. recommend fixing its variance to something large (e.g., 10^6) rather than estimating it, so it doesn't compete with or bias the slope estimates:

```r
# Set up variance structure manually via glmmTMB's map argument
# Following Muff et al. 2020 tutorial approach:
m <- glmmTMB(
  case_ ~ elevation + canopy_cover + dist_to_road +
    (1 | animal_id) + (0 + elevation | animal_id) + (0 + canopy_cover | animal_id),
  family = binomial(),
  weights = weight,
  data = dat,
  doFit = FALSE
)
m$parameters$theta[1] <- log(1e3)  # fix intercept SD (log scale)
nvarparm <- length(m$parameters$theta)
m$mapArg <- list(theta = factor(c(NA, 1:(nvarparm - 1))))
m_fit <- glmmTMB:::fitTMB(m)
```

This exact recipe (with reproducible code) is laid out step-by-step in the Muff, Signer & Fieberg (2020) *Journal of Animal Ecology* tutorial — worth pulling up directly since the `map` argument syntax is fiddly and version-sensitive in `glmmTMB`.

**2. Availability sampling ratio matters.** With used:available at roughly 1:10 to 1:1000 (common range), keep the weight consistent with your actual sampling ratio design — some authors set the weight to exactly match rather than always using 1000; check what ratio you're using in `amt::random_points()`.

**3. Correlation structure between random slopes.** Decide whether random slopes should be correlated with each other and the intercept (`(1 + elevation + canopy_cover | animal_id)`) or independent (separate `(0 + var | animal_id)` terms, uncorrelated). Independent terms are more common in caribou RSF papers for interpretability and convergence, but test both.

## Multi-scale (2nd vs 3rd order) structuring
Since you want both orders:
- **2nd order**: available points drawn from the broader study area/landscape (e.g., herd range), used points = all caribou locations
- **3rd order**: available points drawn from within each individual's home range (from `adehabitatHR` or `ctmm` UD)

Fit these as **separate models** with the same `glmmTMB` structure above — don't try to combine orders in one model, as the covariate ranges and ecological interpretation differ. Compare coefficient direction/magnitude across orders to assess scale-dependent selection (a common result section in caribou RSF papers — e.g., avoidance of disturbance emerging more strongly at 3rd order).

## Model comparison across orders/candidate sets
```r
library(MuMIn)
model.sel(m1_2ndorder, m2_2ndorder, m3_2ndorder)  # AICc-based ranking within an order
```
Note: `MuMIn` has some compatibility quirks with `glmmTMB` random-effect variance structures when you've fixed `theta` — test on a small candidate set first.

Would it help if I sketched out a fuller end-to-end script (from `amt` availability sampling through to prediction/mapping), or do you already have the used-available dataset built and just need the model-fitting stage?
