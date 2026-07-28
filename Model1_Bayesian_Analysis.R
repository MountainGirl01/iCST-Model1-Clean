# ==============================================================================
# iCST Model 1: Bayesian Analysis of Treatment Effect on ADAS-Cog
# Outcome: ADAS_20_FU2 (26-week follow-up)
# Predictors: Baseline-centred ADAS-Cog, Randomisation (iCST vs TAU Control)
# ==============================================================================


# Load Libraries --------------------------------------------------------

library(readr)
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)


# Import & Prepare Data --------------------------------------------------

dat <- read_csv("~/Downloads/THESIS COMBINED Complete case 26 week follow up data.csv")

# Inspect
head(dat)
names(dat)
table(dat$Randomisation)
summary(dat$c_BASELINE_ADAScog)

# Set factor levels so TAU Control is the reference group
dat$Randomisation <- factor(dat$Randomisation,
                            levels = c("TAU Control", "iCST"))
table(dat$Randomisation)


# Model 1: Primary Bayesian Model -----------------------------------------
# Weakly-informative priors, primary treatment effect prior centred on -1.92

priors <- c(
  prior(normal(20, 7),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(exponential(0.1), class = sigma)
)

fit <- brm(
  formula = ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation,
  data    = dat,
  family  = gaussian(),
  prior   = priors,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  seed    = 42,
  file    = "icst_primary_fit"
)

summary(fit)