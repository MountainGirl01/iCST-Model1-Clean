# ==============================================================================
# iCST Model 1: Bayesian Analysis of Treatment Effect on ADAS-Cog
# Outcome: ADAS_20_FU2 (26-week follow-up)
# Predictors: Baseline-centred ADAS-Cog, Randomisation (iCST vs TAU Control)
# ==============================================================================

# ============================================================
# Load Libraries 
# ============================================================

library(readr)
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)

# ============================================================
# Import & Prepare Data 
# ============================================================

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

# ============================================================
# Model 1: Primary Bayesian Model 
# ============================================================

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

# ============================================================
# Model 2: Optimistic Bayesian Model 
# ============================================================

priors_optimistic_original <- c(
  prior(normal(20, 7),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-2.9, 2),  class = b, coef = RandomisationiCST),
  prior(exponential(0.1), class = sigma)
)

fit_optimistic_original <- brm(
  formula = ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation,
  data    = dat,
  family  = gaussian(),
  prior   = priors_optimistic_original,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  seed    = 42,
  file    = "icst_optimistic_original"
)

summary(fit_optimistic_original)

# ============================================================
# Model 3: Pessimistic Bayesian Model 
# ============================================================

priors_pessimistic_original <- c(
  prior(normal(20, 7),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-0.5, 2),  class = b, coef = RandomisationiCST),
  prior(exponential(0.1), class = sigma)
)

fit_pessimistic_original <- brm(
    formula = ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation,
    data    = dat,
    family  = gaussian(),
    prior   = priors_pessimistic_original,
    chains  = 4,
    iter    = 2000,
    warmup  = 1000,
    cores   = 4,
    seed    = 42,
    file    = "icst_pessimistic_original"
)

summary(fit_pessimistic_original)

# ============================================================
# Run predictive checks
# ============================================================

p_prior <- pp_check(fit_prior_original, ndraws = 100) +
  ggtitle("Prior Predictive Check — Original Intercept Normal(20, 7)") +
  xlab("ADAS-Cog score") +
  xlim(-20, 80) +
  theme_minimal()

p_posterior <- pp_check(fit_original, ndraws = 100) +
  ggtitle("Posterior Predictive Check") +
  xlab("ADAS-Cog score") +
  xlim(-20, 80) +
  theme_minimal()

p_prior + p_posterior

ggsave("predictive_checks_original.png", width = 12, height = 5, dpi = 300)

# ============================================================
# Overlaid Posterior Distributions with Probability Labels
# ============================================================

# Extract posterior draws for the treatment effect from each model
draws_primary     <- as_draws_df(fit)$b_RandomisationiCST
draws_optimistic  <- as_draws_df(fit_optimistic_original)$b_RandomisationiCST
draws_pessimistic <- as_draws_df(fit_pessimistic_original)$b_RandomisationiCST

# Combine into one data frame for plotting
draws_df <- data.frame(
  value = c(draws_primary, draws_optimistic, draws_pessimistic),
  model = rep(c("Primary (−1.92)", "Optimistic (−2.9)", "Pessimistic (−0.5)"),
              each = length(draws_primary))
)
draws_df$model <- factor(draws_df$model,
                         levels = c("Optimistic (−2.9)", "Primary (−1.92)", "Pessimistic (−0.5)"))

# Posterior probability of benefit (treatment effect < 0) for each model
p_benefit_primary     <- mean(draws_primary < 0)
p_benefit_optimistic  <- mean(draws_optimistic < 0)
p_benefit_pessimistic <- mean(draws_pessimistic < 0)

ggplot(draws_df, aes(x = value, fill = model, colour = model)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.8) +
  annotate("text", x = 0.15, y = 0.58, label = "No effect",
           colour = "red", size = 3.5, hjust = 0) +
  annotate("text", x = -5.5, y = 0.50,
           label = paste0("Optimistic: P(benefit) = ", round(p_benefit_optimistic * 100, 1), "%"),
           colour = "#F8766D", size = 3.2, hjust = 0) +
  annotate("text", x = -5.5, y = 0.45,
           label = paste0("Primary: P(benefit) = ", round(p_benefit_primary * 100, 1), "%"),
           colour = "#619CFF", size = 3.2, hjust = 0) +
  annotate("text", x = -5.5, y = 0.40,
           label = paste0("Pessimistic: P(benefit) = ", round(p_benefit_pessimistic * 100, 1), "%"),
           colour = "#00BA38", size = 3.2, hjust = 0) +
  labs(
    title    = "Posterior distribution of iCST treatment effect on ADAS-Cog",
    subtitle = "By prior specification",
    x        = "Treatment effect (ADAS-Cog points)",
    y        = "Density",
    fill     = "Prior specification",
    colour   = "Prior specification",
    caption  = "Negative values indicate improvement (lower ADAS-Cog = better cognitive function).\nP(benefit) = posterior probability of any beneficial effect of iCST versus TAU Control.\nRed dashed line = no effect (0 points)."
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("posterior_distribution.png", width = 10, height = 6, dpi = 300)

# ============================================================
# Prior and Posterior Predictive Check Plots
# ============================================================

fit_prior_only <- brm(
  formula = ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation,
  data    = dat,
  family  = gaussian(),
  prior   = priors,
  sample_prior = "only",
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  seed    = 42,
  file    = "icst_prior_only",
  file_refit = "on_change"
)
p_prior <- pp_check(fit_prior_only, ndraws = 100) +
  ggtitle("Prior Predictive Check") +
  xlab("ADAS-Cog score") +
  xlim(-20, 80) +
  theme_minimal()
p_posterior <- pp_check(fit, ndraws = 100) +
  ggtitle("Posterior Predictive Check") +
  xlab("ADAS-Cog score") +
  xlim(-20, 80) +
  theme_minimal()
p_prior + p_posterior
ggsave("predictive_checks.png", width = 12, height = 5, dpi = 300)

# ============================================================
# Caterpillar Plot — Treatment Effect Across Prior Specifications (95% CI)
# ============================================================

draws_primary     <- as_draws_df(fit)$b_RandomisationiCST
draws_optimistic  <- as_draws_df(fit_optimistic_original)$b_RandomisationiCST
draws_pessimistic <- as_draws_df(fit_pessimistic_original)$b_RandomisationiCST

caterpillar_df <- data.frame(
  model = c("Optimistic (−2.9)", "Primary (−1.92)", "Pessimistic (−0.5)"),
  median = c(median(draws_optimistic), median(draws_primary), median(draws_pessimistic)),
  ci_lower = c(quantile(draws_optimistic, 0.025), quantile(draws_primary, 0.025), quantile(draws_pessimistic, 0.025)),
  ci_upper = c(quantile(draws_optimistic, 0.975), quantile(draws_primary, 0.975), quantile(draws_pessimistic, 0.975))
)

caterpillar_df$model <- factor(caterpillar_df$model,
                               levels = c("Pessimistic (−0.5)", "Primary (−1.92)", "Optimistic (−2.9)"))

ggplot(caterpillar_df, aes(x = median, y = model, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.8) +
  geom_segment(aes(x = ci_lower, xend = ci_upper, y = model, yend = model),
               linewidth = 1.2) +
  geom_point(size = 3, colour = "black") +
  labs(
    title    = "Treatment effect on ADAS-Cog by prior specification",
    subtitle = "Posterior median with 95% credible intervals",
    x        = "Treatment effect (ADAS-Cog points)",
    y        = NULL,
    caption  = "Negative values indicate improvement on ADAS-Cog."
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("caterpillar_plot.png", width = 10, height = 5, dpi = 300)

# ============================================================
# Complete Summary Table — All Three Models
# ============================================================

library(dplyr)

# Extract posterior draws for the treatment effect from each model
draws_primary     <- as_draws_df(fit)$b_RandomisationiCST
draws_optimistic  <- as_draws_df(fit_optimistic_original)$b_RandomisationiCST
draws_pessimistic <- as_draws_df(fit_pessimistic_original)$b_RandomisationiCST

# Helper function to compute all stats for one set of draws
summarise_draws <- function(draws, model_name, prior_used) {
  data.frame(
    Model               = model_name,
    Prior               = prior_used,
    Mean                = mean(draws),
    Median              = median(draws),
    SD                  = sd(draws),
    CI_lower_95         = quantile(draws, 0.025),
    CI_upper_95         = quantile(draws, 0.975),
    P_benefit           = mean(draws < 0),                     # P(any beneficial effect)
    P_harm              = mean(draws > 0),                     # P(any harmful effect)
    P_clinically_meaningful = mean(draws < -1.5),               # adjust threshold if you have one
    row.names = NULL
  )
}

results_table <- bind_rows(
  summarise_draws(draws_optimistic,  "Optimistic",  "Normal(-2.9, 2)"),
  summarise_draws(draws_primary,     "Primary",     "Normal(-1.92, 2)"),
  summarise_draws(draws_pessimistic, "Pessimistic", "Normal(-0.5, 2)")
)

# Round for readability
results_table_rounded <- results_table %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(results_table_rounded)

# Save as CSV
write_csv(results_table_rounded, "icst_results_summary_table.csv")

# ============================================================
# Add convergence diagnostics (Rhat, Bulk_ESS, Tail_ESS)
# ============================================================

extract_diagnostics <- function(fit, model_name) {
  s <- summary(fit)$fixed["RandomisationiCST", ]
  data.frame(
    Model     = model_name,
    Est_Error = s["Est.Error"],
    Rhat      = s["Rhat"],
    Bulk_ESS  = s["Bulk_ESS"],
    Tail_ESS  = s["Tail_ESS"],
    row.names = NULL
  )
}

diagnostics_table <- bind_rows(
  extract_diagnostics(fit_optimistic_original,  "Optimistic"),
  extract_diagnostics(fit,                       "Primary"),
  extract_diagnostics(fit_pessimistic_original,  "Pessimistic")
)

# Merge with results table
full_results_table <- results_table_rounded %>%
  left_join(diagnostics_table, by = "Model") %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(full_results_table)
write_csv(full_results_table, "icst_full_results_table.csv")

# ============================================================
# Render results table as a plot (in RStudio Plots pane)
# ============================================================

library(gridExtra)
library(grid)

# Uses results_table_rounded from the previous step
table_plot <- tableGrob(results_table_rounded, rows = NULL,
                        theme = ttheme_minimal(
                          core = list(fg_params = list(cex = 0.8)),
                          colhead = list(fg_params = list(cex = 0.85, fontface = "bold"))
                        ))

grid.newpage()
grid.draw(table_plot)

ggsave("results_table.png", table_plot, width = 12, height = 3, dpi = 300)

# ============================================================
# Combined Parameter Summary Table — All Three Models
# ============================================================

library(dplyr)

# Extract full posterior summaries (includes b_*, sigma, Intercept, lprior, lp__)
summary_primary     <- as.data.frame(posterior_summary(fit))
summary_optimistic  <- as.data.frame(posterior_summary(fit_optimistic_original))
summary_pessimistic <- as.data.frame(posterior_summary(fit_pessimistic_original))

# Add a Parameter column (from row names) and a Model label to each
summary_primary$Parameter     <- rownames(summary_primary)
summary_optimistic$Parameter  <- rownames(summary_optimistic)
summary_pessimistic$Parameter <- rownames(summary_pessimistic)

summary_primary$Model     <- "Primary"
summary_optimistic$Model  <- "Optimistic"
summary_pessimistic$Model <- "Pessimistic"

# Combine into one table
all_params <- bind_rows(summary_primary, summary_optimistic, summary_pessimistic)

# Reorder columns: Model, Parameter, Estimate, Est.Error, Q2.5, Q97.5
all_params <- all_params %>%
  select(Model, Parameter, Estimate, Est.Error, Q2.5, Q97.5)

# Round for readability
all_params <- all_params %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

print(all_params)

# Save as CSV
write.csv(all_params, "all_model_parameters.csv", row.names = FALSE)

# ============================================================
# Posterior Densities and Trace Plots
# ============================================================

library(brms)
library(bayesplot)
library(ggplot2)

fit <- readRDS("icst_primary_fit.rds")

draws <- as_draws_array(fit)

combo_plot <- mcmc_combo(
  draws,
  combo = c("dens_overlay", "trace"),
  pars  = c("b_Intercept", "b_c_BASELINE_ADAScog", "b_RandomisationiCST", "sigma")
)

combo_plot

ggsave("icst_primary_fit_combo.png", combo_plot,
       width = 10, height = 8, dpi = 300)

plot(fit, variable = "^b_|^sigma", regex = TRUE)

png("icst_primary_fit_histtrace.png", width = 1800, height = 1200, res = 150)
plot(fit, variable = "^b_|^sigma", regex = TRUE)
dev.off()

# ============================================================
# CI Forest Plot
# ============================================================

draws_primary     <- as_draws_df(fit)$b_RandomisationiCST
draws_optimistic  <- as_draws_df(fit_optimistic_original)$b_RandomisationiCST
draws_pessimistic <- as_draws_df(fit_pessimistic_original)$b_RandomisationiCST

caterpillar_df <- data.frame(
  model = c("Optimistic (−2.9)", "Primary (−1.92)", "Pessimistic (−0.5)"),
  median = c(median(draws_optimistic), median(draws_primary), median(draws_pessimistic)),
  ci_lower = c(quantile(draws_optimistic, 0.025), quantile(draws_primary, 0.025), quantile(draws_pessimistic, 0.025)),
  ci_upper = c(quantile(draws_optimistic, 0.975), quantile(draws_primary, 0.975), quantile(draws_pessimistic, 0.975))
)

caterpillar_df$model <- factor(caterpillar_df$model,
                               levels = c("Pessimistic (−0.5)", "Primary (−1.92)", "Optimistic (−2.9)"))

ggplot(caterpillar_df, aes(x = median, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, linewidth = 0.8, colour = "black") +
  geom_point(size = 3, colour = "black", fill = "black", shape = 21) +
  labs(
    title    = "Treatment Effect on ADAS-Cog by Prior Specification",
    subtitle = "Posterior median with 95% credible intervals",
    x        = "Treatment effect (ADAS-Cog points)",
    y        = "Prior specification",
    caption  = "Note. Negative values indicate improvement on ADAS-Cog. Error bars represent 95% credible intervals."
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0, size = 11),
    plot.subtitle     = element_text(hjust = 0, size = 11, colour = "grey30"),
    plot.caption      = element_text(hjust = 0, size = 11, colour = "grey30", face = "italic"),
    axis.title        = element_text(size = 11, face = "bold"),
    axis.text         = element_text(size = 11, colour = "black"),
    axis.line         = element_line(colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3)
  )

ggsave("caterpillar_plot.png", width = 10, height = 5, dpi = 300)

# ============================================================
# Thesis Ready Posterior Distribution Plot
# ============================================================

# Extract posterior draws for the treatment effect from each model
draws_primary     <- as_draws_df(fit)$b_RandomisationiCST
draws_optimistic  <- as_draws_df(fit_optimistic_original)$b_RandomisationiCST
draws_pessimistic <- as_draws_df(fit_pessimistic_original)$b_RandomisationiCST

# Combine into one data frame for plotting
draws_df <- data.frame(
  value = c(draws_primary, draws_optimistic, draws_pessimistic),
  model = rep(c("Primary (−1.92)", "Optimistic (−2.9)", "Pessimistic (−0.5)"),
              each = length(draws_primary))
)

draws_df$model <- factor(draws_df$model,
                         levels = c("Optimistic (−2.9)", "Primary (−1.92)", "Pessimistic (−0.5)"))

# Posterior probability of benefit (treatment effect < 0) for each model
p_benefit_primary     <- mean(draws_primary < 0)
p_benefit_optimistic  <- mean(draws_optimistic < 0)
p_benefit_pessimistic <- mean(draws_pessimistic < 0)

ggplot(draws_df, aes(x = value, fill = model, colour = model)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
  annotate("text", x = 0.15, y = 0.58, label = "No effect",
           colour = "black", size = 3.2, hjust = 0, family = "Arial") +
  scale_fill_manual(values = c(
    "Optimistic (−2.9)"  = "#00BA38",
    "Primary (−1.92)"    = "#619CFF",
    "Pessimistic (−0.5)" = "#F8766D"
  )) +
  scale_colour_manual(values = c(
    "Optimistic (−2.9)"  = "#00BA38",
    "Primary (−1.92)"    = "#619CFF",
    "Pessimistic (−0.5)" = "#F8766D"
  )) +
  labs(
    title    = "Posterior Distribution of iCST Treatment Effect on ADAS-Cog",
    subtitle = "By prior specification",
    x        = "Treatment effect (ADAS-Cog points)",
    y        = "Density",
    fill     = "Prior specification",
    colour   = "Prior specification"
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    legend.position   = "right",
    legend.direction  = "vertical",
    legend.title      = element_text(size = 11, face = "bold"),
    legend.text       = element_text(size = 11),
    plot.title        = element_text(face = "bold", hjust = 0, size = 11),
    plot.subtitle     = element_text(hjust = 0, size = 11, colour = "grey30"),
    axis.title        = element_text(size = 11, face = "bold"),
    axis.text         = element_text(size = 11, colour = "black"),
    axis.line         = element_line(colour = "black"),
    panel.grid.minor  = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1))

ggsave("posterior_distribution.png", width = 11, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Print P(benefit) values to console so you can copy them into your
# APA 7 figure note underneath the image (see note text below)
# ---------------------------------------------------------------------------
cat("Optimistic: P(benefit) =", round(p_benefit_optimistic * 100, 1), "%\n")
cat("Primary: P(benefit) =", round(p_benefit_primary * 100, 1), "%\n")
cat("Pessimistic: P(benefit) =", round(p_benefit_pessimistic * 100, 1), "%\n")

# ============================================================
# Create Characteristics Table for Model 1
# ============================================================

pwd_table <- dat_table1 |>
  select(
    Randomisation,
    pwd_female,
    cRelAgeP,
    BASELINE_ADAS_20,
    pwd_white,
    relationship_spouse,
    relationship_child
  ) |>
  tbl_summary(
    by = Randomisation,
    
    type = list(
      pwd_female ~ "dichotomous",
      pwd_white ~ "dichotomous",
      relationship_spouse ~ "dichotomous",
      relationship_child ~ "dichotomous"
    ),
    
    value = list(
      pwd_female ~ "Yes",
      pwd_white ~ "Yes",
      relationship_spouse ~ "Yes",
      relationship_child ~ "Yes"
    ),
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_dichotomous() ~ "{n} ({p}%)"
    ),
    
    digits = list(
      all_continuous() ~ 1,
      all_dichotomous() ~ c(0, 1)
    ),
    
    missing = "no",
    
    label = list(
      pwd_female ~ "Female, n (%)",
      cRelAgeP ~ "Age, years, mean (SD)",
      BASELINE_ADAS_20 ~ "Baseline ADAS-Cog, mean (SD)",
      pwd_white ~ "White ethnicity, n (%)",
      relationship_spouse ~ "Caregiver relationship: spouse, n (%)",
      relationship_child ~ "Caregiver relationship: son/daughter, n (%)"
    )
  )

caregiver_table <- dat_table1 |>
  select(
    Randomisation,
    caregiver_female,
    cRelAgeC,
    caregiver_white
  ) |>
  tbl_summary(
    by = Randomisation,
    
    type = list(
      caregiver_female ~ "dichotomous",
      caregiver_white ~ "dichotomous"
    ),
    
    value = list(
      caregiver_female ~ "Yes",
      caregiver_white ~ "Yes"
    ),
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_dichotomous() ~ "{n} ({p}%)"
    ),
    
    digits = list(
      all_continuous() ~ 1,
      all_dichotomous() ~ c(0, 1)
    ),
    
    missing = "no",
    
    label = list(
      caregiver_female ~ "Female, n (%)",
      cRelAgeC ~ "Age, years, mean (SD)",
      caregiver_white ~ "White ethnicity, n (%)"
    )
  )

table1 <- tbl_stack(
  tbls = list(
    pwd_table,
    caregiver_table
  ),
  group_header = c(
    "**Person with dementia**",
    "**Caregiver**"
  )
) |>
  modify_header(
    label ~ "**Characteristic**",
    stat_1 ~ "**Treatment group  \nN = {n}**",
    stat_2 ~ "**Control group  \nN = {n}**"
  ) |>
  modify_caption(
    "**Table 1. Baseline Characteristics of People With Dementia and Their Caregivers**"
  ) |>
  modify_footnote(
    all_stat_cols() ~
      "Values are presented as mean (SD) or n (%). Percentages use available data for each variable."
  ) |>
  bold_labels()

table1

# ============================================================
# Model 1 v 2 - Change Priors to Informative and Skeptical 
# ============================================================

## ============================================================================
## iCST Model 1: Bayesian Analysis of Treatment Effect on ADAS-Cog
## Outcome: ADAS_20_FU2 (26-week follow-up)
## Predictors: Baseline-centred ADAS-Cog, Randomisation (iCST vs TAU Control)
## Two-prior structure: Informative vs Skeptical
## ============================================================================

# ---- LIBRARIES ---------------------------------------------------------
library(readr)
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidybayes)
library(gt)


# ---- TRUNCATED FORMULA ---------------------------------------------------
# Bounds set to 0-60: true floor is 0, ceiling padded above the
# observed max (55.67) to avoid cutting off real/predicted values.
adas_formula_trunc <- bf(ADAS_20_FU2 | trunc(lb = 0, ub = 60) ~
                           c_BASELINE_ADAScog + Randomisation)

# ---- FINAL PRIORS ---------------------------------------------------------
# Sigma: normal(9,3) truncated at 0 - centred near the plausible
# residual SD rather than at 0, avoiding unrealistic near-zero draws.

priors_informative_final <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

priors_skeptical_final <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(0, 2),     class = b, coef = RandomisationiCST),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

# ---- FIT BOTH MODELS -------------------------------------------------------

fit_informative_final <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_informative_final,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model1_informative_final"
)

fit_skeptical_final <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_skeptical_final,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model1_skeptical_final"
)

summary(fit_informative_final)
summary(fit_skeptical_final)

# ---- CONVERGENCE DIAGNOSTICS ------------------------------------------------

cat("Max Rhat (Informative):", round(max(brms::rhat(fit_informative_final)), 4), "\n")
cat("Max Rhat (Skeptical):", round(max(brms::rhat(fit_skeptical_final)), 4), "\n")

div_informative <- sum(subset(brms::nuts_params(fit_informative_final),
                              Parameter == "divergent__")$Value)
div_skeptical <- sum(subset(brms::nuts_params(fit_skeptical_final),
                            Parameter == "divergent__")$Value)

cat("Divergent transitions (Informative):", div_informative, "\n")
cat("Divergent transitions (Skeptical):", div_skeptical, "\n")

# ---- PRIOR PREDICTIVE MODELS (for predictive check plots) ------------------

prior_check_informative_final <- brm(
  formula = adas_formula_trunc, data = dat, family = gaussian(),
  prior = priors_informative_final, sample_prior = "only",
  chains = 4, iter = 1000, seed = 42
)

prior_check_skeptical_final <- brm(
  formula = adas_formula_trunc, data = dat, family = gaussian(),
  prior = priors_skeptical_final, sample_prior = "only",
  chains = 4, iter = 1000, seed = 42
)

## ------------------------------------------------------------------
## Model 1: Prior and Posterior Predictive Checks
## Clean side-by-side format (no subtitles), matching reference style.
## Produces two separate panels: Informative and Skeptical.
## ------------------------------------------------------------------

library(ggplot2)
library(patchwork)


p_prior_inf <- brms::pp_check(prior_check_informative_final, ndraws = 100) +
  labs(title = "Prior Predictive Check",
       x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_inf <- brms::pp_check(fit_informative_final, ndraws = 100) +
  labs(title = "Posterior Predictive Check",
       x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_informative <- p_prior_inf | p_post_inf

combined_informative

ggsave("model1_ppc_informative_clean.png", combined_informative,
       width = 10, height = 5, dpi = 300)

p_prior_skep <- brms::pp_check(prior_check_skeptical_final, ndraws = 100) +
  labs(title = "Prior Predictive Check",
       x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_skep <- brms::pp_check(fit_skeptical_final, ndraws = 100) +
  labs(title = "Posterior Predictive Check",
       x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_skeptical <- p_prior_skep | p_post_skep

combined_skeptical

ggsave("model1_ppc_skeptical_clean.png", combined_skeptical,
       width = 10, height = 5, dpi = 300)

## ------------------------------------------------------------------
## Model 1: Overlaid Posterior Distribution Plot
## ------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(tidybayes)

draws_informative <- as_draws_df(fit_informative_final) %>%
  select(b_RandomisationiCST) %>%
  mutate(Prior = "Informative")

draws_skeptical <- as_draws_df(fit_skeptical_final) %>%
  select(b_RandomisationiCST) %>%
  mutate(Prior = "Skeptical")

all_draws_m1 <- bind_rows(draws_informative, draws_skeptical)
all_draws_m1$Prior <- factor(all_draws_m1$Prior, levels = c("Informative", "Skeptical"))


prob_labels_m1 <- all_draws_m1 %>%
  group_by(Prior) %>%
  summarise(prob_benefit = mean(b_RandomisationiCST < 0), .groups = "drop") %>%
  mutate(label = paste0(Prior, ": P(benefit) = ",
                        scales::percent(prob_benefit, accuracy = 0.1)))
library(dplyr)
library(ggplot2)

draws_inf <- as_draws_df(fit_informative_final)$b_RandomisationiCST
draws_skep <- as_draws_df(fit_skeptical_final)$b_RandomisationiCST

forest_data <- tibble::tibble(
  Model = c("Informative", "Skeptical"),
  Median = c(median(draws_inf), median(draws_skep)),
  CI_lower_95 = c(quantile(draws_inf, 0.025), quantile(draws_skep, 0.025)),
  CI_upper_95 = c(quantile(draws_inf, 0.975), quantile(draws_skep, 0.975))
) %>%
  mutate(Model = factor(Model, levels = c("Skeptical", "Informative")))

forest_plot <- ggplot(forest_data, aes(x = Median, y = Model, color = Model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = CI_lower_95, xmax = CI_upper_95), height = 0, linewidth = 1) +
  geom_point(size = 3, color = "black") +
  scale_color_manual(values = c("Informative" = "#4E9FDA", "Skeptical" = "#4CAF50")) +
  labs(
    title = "Treatment effect on ADAS-Cog by prior specification",
    subtitle = "Posterior median with 95% credible intervals",
    x = "Treatment effect (ADAS-Cog points)",
    y = NULL,
    caption = "Negative values indicate improvement on ADAS-Cog."
  ) +
  theme_minimal() +
  theme(legend.position = "none", plot.caption = element_text(hjust = 1, size = 8, color = "grey40"))

forest_plot

ggsave("model1_forest_plot_revised.png", forest_plot, width = 8, height = 4, dpi = 300)
print(prob_labels_m1)


overlay_plot_m1 <- ggplot(all_draws_m1, aes(x = b_RandomisationiCST,
                                            fill = Prior, color = Prior)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  annotate("text", x = 0.15, y = Inf, label = "No effect",
           color = "black", hjust = 0, vjust = 1.5, size = 3.5) +
  scale_fill_manual(
    values = c("Informative" = "#4E9FDA", "Skeptical" = "#4CAF50"),
    labels = c("Informative" = "Informative (\u22121.92)",
               "Skeptical" = "Skeptical (0)")
  ) +
  scale_color_manual(
    values = c("Informative" = "#4E9FDA", "Skeptical" = "#4CAF50"),
    labels = c("Informative" = "Informative (\u22121.92)",
               "Skeptical" = "Skeptical (0)")
  ) +
  labs(
    x = "Treatment effect (ADAS-Cog points)",
    y = "Density",
    fill = "Prior specification",
    color = "Prior specification"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  coord_cartesian(clip = "off")

overlay_plot_m1

ggsave("model1_posterior_treatment_effect_matched_style.png", overlay_plot_m1,
       width = 9, height = 5, dpi = 300)

# ---- FOREST PLOT -------------------------------------------------------------

draws_informative <- as_draws_df(fit_informative_final)$b_RandomisationiCST
draws_skeptical   <- as_draws_df(fit_skeptical_final)$b_RandomisationiCST

caterpillar_df <- data.frame(
  model = c("Informative (\u22121.92)", "Skeptical (0)"),
  median = c(median(draws_informative), median(draws_skeptical)),
  ci_lower = c(quantile(draws_informative, 0.025), quantile(draws_skeptical, 0.025)),
  ci_upper = c(quantile(draws_informative, 0.975), quantile(draws_skeptical, 0.975))
)

caterpillar_df$model <- factor(caterpillar_df$model,
                               levels = c("Skeptical (0)", "Informative (\u22121.92)"))

ggplot(caterpillar_df, aes(x = median, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, linewidth = 0.8, colour = "black") +
  geom_point(size = 3, colour = "black", fill = "black", shape = 21) +
  labs(
    title    = "Treatment Effect on ADAS-Cog by Prior Specification",
    subtitle = "Posterior median with 95% credible intervals",
    x        = "Treatment effect (ADAS-Cog points)",
    y        = "Prior specification",
    caption  = "Note. Negative values indicate improvement on ADAS-Cog. Error bars represent 95% credible intervals."
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0, size = 11),
    plot.subtitle     = element_text(hjust = 0, size = 11, colour = "grey30"),
    plot.caption      = element_text(hjust = 0, size = 11, colour = "grey30", face = "italic"),
    axis.title        = element_text(size = 11, face = "bold"),
    axis.text         = element_text(size = 11, colour = "black"),
    axis.line         = element_line(colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3)
  )

ggsave("caterpillar_plot_bw_arial.png", width = 10, height = 5, dpi = 300)

# ---- FULL RESULTS TABLE (TREATMENT EFFECT ONLY) -------------------------------

mcid_threshold <- -1.4   # replace with your study's actual MCID if different

summary_table <- tibble::tibble(
  Model = c("Informative", "Skeptical"),
  Prior = c("Normal(-1.92, 2)", "Normal(0, 2)"),
  Mean = c(mean(draws_inf), mean(draws_skep)),
  Median = c(median(draws_inf), median(draws_skep)),
  SD = c(sd(draws_inf), sd(draws_skep)),
  CI_lower_95 = c(quantile(draws_inf, 0.025), quantile(draws_skep, 0.025)),
  CI_upper_95 = c(quantile(draws_inf, 0.975), quantile(draws_skep, 0.975)),
  P_benefit = c(mean(draws_inf < 0), mean(draws_skep < 0)),
  P_harm = c(mean(draws_inf > 0), mean(draws_skep > 0)),
  P_clinically_meaningful = c(mean(draws_inf < mcid_threshold), mean(draws_skep < mcid_threshold)),
  Rhat = c(max(brms::rhat(fit_informative_final)), max(brms::rhat(fit_skeptical_final))),
  Bulk_ESS = c(summary(fit_informative_final)$fixed["RandomisationiCST", "Bulk_ESS"],
               summary(fit_skeptical_final)$fixed["RandomisationiCST", "Bulk_ESS"]),
  Tail_ESS = c(summary(fit_informative_final)$fixed["RandomisationiCST", "Tail_ESS"],
               summary(fit_skeptical_final)$fixed["RandomisationiCST", "Tail_ESS"])
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(summary_table)
write.csv(summary_table, "model1_full_results_table.csv", row.names = FALSE)

results_gt <- summary_table %>%
  gt() %>%
  tab_header(
    title = "Model 1: Treatment Effect on ADAS-Cog by Prior Specification",
    subtitle = "iCST vs TAU Control, 26-week follow-up"
  ) %>%
  fmt_number(columns = c(Mean, Median, SD, CI_lower_95, CI_upper_95,
                         P_benefit, P_harm, P_clinically_meaningful, Rhat), decimals = 3) %>%
  fmt_number(columns = c(Bulk_ESS, Tail_ESS), decimals = 0) %>%
  tab_source_note(source_note = paste0(
    "P(benefit) = P(treatment effect < 0). P(clinically meaningful) = P(treatment effect < ",
    mcid_threshold, "). N = ", nrow(dat), "."))

results_gt
gtsave(results_gt, "model1_full_results_table.png")

# ---- FULL PARAMETER TABLE (ALL COEFFICIENTS, BOTH MODELS) ---------------------

extract_all_params <- function(fit, model_label) {
  fixed_df <- as.data.frame(summary(fit)$fixed)
  fixed_df$Parameter <- rownames(fixed_df)
  sigma_df <- as.data.frame(summary(fit)$spec_pars)
  sigma_df$Parameter <- rownames(sigma_df)
  combined <- bind_rows(fixed_df, sigma_df)
  combined$Model <- model_label
  combined %>% select(Model, Parameter, Estimate, `l-95% CI`, `u-95% CI`, Rhat, Bulk_ESS, Tail_ESS)
}

params_informative <- extract_all_params(fit_informative_final, "Informative")
params_skeptical   <- extract_all_params(fit_skeptical_final, "Skeptical")
all_params_table <- bind_rows(params_informative, params_skeptical)

all_params_table$Parameter <- dplyr::recode(all_params_table$Parameter,
                                            "Intercept" = "Intercept", "c_BASELINE_ADAScog" = "Baseline",
                                            "RandomisationiCST" = "Treatment", "sigma" = "Sigma")

param_order <- c("Intercept", "Baseline", "Sigma", "Treatment")
all_params_table$Parameter <- factor(all_params_table$Parameter, levels = param_order)
all_params_table$Model <- factor(all_params_table$Model, levels = c("Informative", "Skeptical"))

all_params_table <- all_params_table %>%
  arrange(Model, Parameter) %>%
  mutate(across(c(Estimate, `l-95% CI`, `u-95% CI`, Rhat), ~round(.x, 2)),
         across(c(Bulk_ESS, Tail_ESS), ~round(.x, 0)))

print(all_params_table)
write.csv(all_params_table, "model1_full_parameter_table.csv", row.names = FALSE)

params_gt <- all_params_table %>%
  rename(`Lower CrI` = `l-95% CI`, `Upper CrI` = `u-95% CI`,
         `Bulk ESS` = Bulk_ESS, `Tail ESS` = Tail_ESS, `R-hat` = Rhat) %>%
  group_by(Model) %>%
  gt() %>%
  tab_header(title = "Table 2", subtitle = "Posterior Estimates of Treatment Effect on ADAS-Cog") %>%
  fmt_number(columns = c(Estimate, `Lower CrI`, `Upper CrI`), decimals = 2) %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Parameter)

params_gt
gtsave(params_gt, "model1_full_parameter_table.png")

## ------------------------------------------------------------------
## Model 1: v2 Sept Excluded missing cases QCPR for consistent participants across models
## ------------------------------------------------------------------

# ============================================================
# Import & Prepare Data
# ============================================================

dat <- read_csv("~/Downloads/THESIS COMBINED Complete case 26 week follow up data.csv")

dat$Randomisation <- factor(dat$Randomisation,
                            levels = c("TAU Control", "iCST"))

# Exclude participants missing baseline or follow-up ADAS-Cog,
# and missing/zero baseline CQ-CPR (0 treated as missing)
dat <- dat %>%
  filter(
    !is.na(c_BASELINE_ADAScog),
    !is.na(ADAS_20_FU2),
    !is.na(BASELINE_CQCPR_20),
    BASELINE_CQCPR_20 != 0
  )

nrow(dat)
table(dat$Randomisation)

# ---- TRUNCATED FORMULA ---------------------------------------------------

adas_formula_trunc <- bf(ADAS_20_FU2 | trunc(lb = 0, ub = 60) ~
                           c_BASELINE_ADAScog + Randomisation)

# ---- FINAL PRIORS ---------------------------------------------------------

priors_informative_final <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

priors_skeptical_final <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(0, 2),     class = b, coef = RandomisationiCST),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

# ---- FIT BOTH MODELS -------------------------------------------------------

fit_informative_matched <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_informative_final,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model1_informative_matched_N258"
)

fit_skeptical_matched <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_skeptical_final,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model1_skeptical_matched_N258"
)

summary(fit_informative_matched)
summary(fit_skeptical_matched)

## ============================================================================
## iCST Model 1 (Final, N=254): Bayesian Analysis of Treatment Effect on ADAS-Cog
## Outcome: ADAS_20_FU2 (26-week follow-up)
## Predictors: Baseline-centred ADAS-Cog, Randomisation (iCST vs TAU Control)
## Two-prior structure: Informative vs Skeptical
##
## Sample matches Model 2's final analytic sample: the 3 participants
## missing baseline QCPR, PLUS the 4 participants with extreme ADAS-Cog
## change scores (identified via inspection of the change score
## distribution, which showed right-skew), have been excluded.
## ============================================================================

# ---- LIBRARIES ---------------------------------------------------------
library(readr)
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidybayes)
library(gt)

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================

# ---- Import ------------------------------------------------------------

dat <- read_csv("~/Downloads/THESIS COMBINED Complete case 26 week follow up data.csv")

dat$Randomisation <- factor(dat$Randomisation,
                            levels = c("TAU Control", "iCST"))

# ---- Exclude participants missing baseline/follow-up ADAS-Cog, and
# missing/zero baseline QCPR (0 treated as missing) -------------------------

dat <- dat %>%
  filter(
    !is.na(c_BASELINE_ADAScog),
    !is.na(ADAS_20_FU2),
    !is.na(BASELINE_CQCPR_20),
    BASELINE_CQCPR_20 != 0
  )

nrow(dat)   # should be 258
table(dat$Randomisation)

# ---- Exclude the 4 outlier cases (extreme ADAS-Cog change scores) --------
# qPIN 14201, 17018, 13044, 11051

dat <- dat[!(dat$qPIN %in% c(14201, 17018, 13044, 11051)), ]

nrow(dat)   # should now be 254
table(dat$Randomisation)

# ============================================================================
# PART 2: MODEL SPECIFICATION
# ============================================================================

# ---- Truncated formula ----------------------------------------------------

adas_formula_trunc <- bf(ADAS_20_FU2 | trunc(lb = 0, ub = 60) ~
                           c_BASELINE_ADAScog + Randomisation)

# ---- Priors -----------------------------------------------------------------

priors_informative_final <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

priors_skeptical_final <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(0, 2),     class = b, coef = RandomisationiCST),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

# ============================================================================
# PART 3: FIT BOTH MODELS
# ============================================================================

fit_informative_final <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_informative_final,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model1_informative_N254"
)

fit_skeptical_final <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_skeptical_final,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model1_skeptical_N254"
)

summary(fit_informative_final)
summary(fit_skeptical_final)

# ============================================================================
# PART 4: CONVERGENCE DIAGNOSTICS
# ============================================================================

cat("Max Rhat (Informative):", round(max(brms::rhat(fit_informative_final)), 4), "\n")
cat("Max Rhat (Skeptical):", round(max(brms::rhat(fit_skeptical_final)), 4), "\n")

div_informative <- sum(subset(brms::nuts_params(fit_informative_final),
                              Parameter == "divergent__")$Value)
div_skeptical <- sum(subset(brms::nuts_params(fit_skeptical_final),
                            Parameter == "divergent__")$Value)

cat("Divergent transitions (Informative):", div_informative, "\n")
cat("Divergent transitions (Skeptical):", div_skeptical, "\n")

# ============================================================================
# PART 5: PRIOR AND POSTERIOR PREDICTIVE CHECKS
# ============================================================================

prior_check_informative_final <- brm(
  formula = adas_formula_trunc, data = dat, family = gaussian(),
  prior = priors_informative_final, sample_prior = "only",
  chains = 4, iter = 1000, seed = 42
)

prior_check_skeptical_final <- brm(
  formula = adas_formula_trunc, data = dat, family = gaussian(),
  prior = priors_skeptical_final, sample_prior = "only",
  chains = 4, iter = 1000, seed = 42
)

p_prior_inf <- brms::pp_check(prior_check_informative_final, ndraws = 100) +
  labs(title = "Prior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_inf <- brms::pp_check(fit_informative_final, ndraws = 100) +
  labs(title = "Posterior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_informative <- p_prior_inf | p_post_inf
combined_informative
ggsave("model1_ppc_informative_N254.png", combined_informative, width = 10, height = 5, dpi = 300)

p_prior_skep <- brms::pp_check(prior_check_skeptical_final, ndraws = 100) +
  labs(title = "Prior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_skep <- brms::pp_check(fit_skeptical_final, ndraws = 100) +
  labs(title = "Posterior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_skeptical <- p_prior_skep | p_post_skep
combined_skeptical
ggsave("model1_ppc_skeptical_N254.png", combined_skeptical, width = 10, height = 5, dpi = 300)

# ============================================================================
# PART 6: OVERLAID POSTERIOR DISTRIBUTION PLOT
# ============================================================================

draws_informative <- as_draws_df(fit_informative_final) %>%
  select(b_RandomisationiCST) %>%
  mutate(Prior = "Informative")

draws_skeptical <- as_draws_df(fit_skeptical_final) %>%
  select(b_RandomisationiCST) %>%
  mutate(Prior = "Skeptical")

all_draws_m1 <- bind_rows(draws_informative, draws_skeptical)
all_draws_m1$Prior <- factor(all_draws_m1$Prior, levels = c("Informative", "Skeptical"))

prob_labels_m1 <- all_draws_m1 %>%
  group_by(Prior) %>%
  summarise(prob_benefit = mean(b_RandomisationiCST < 0), .groups = "drop") %>%
  mutate(label = paste0(Prior, ": P(benefit) = ",
                        scales::percent(prob_benefit, accuracy = 0.1)))

print(prob_labels_m1)

overlay_plot_m1 <- ggplot(all_draws_m1, aes(x = b_RandomisationiCST,
                                            fill = Prior, color = Prior)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  annotate("text", x = 0.15, y = Inf, label = "No effect",
           color = "black", hjust = 0, vjust = 1.5, size = 3.5) +
  scale_fill_manual(
    values = c("Informative" = "#4E9FDA", "Skeptical" = "#4CAF50"),
    labels = c("Informative" = "Informative (\u22121.92)",
               "Skeptical" = "Skeptical (0)")
  ) +
  scale_color_manual(
    values = c("Informative" = "#4E9FDA", "Skeptical" = "#4CAF50"),
    labels = c("Informative" = "Informative (\u22121.92)",
               "Skeptical" = "Skeptical (0)")
  ) +
  labs(
    x = "Treatment effect (ADAS-Cog points)",
    y = "Density",
    fill = "Prior specification",
    color = "Prior specification"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  coord_cartesian(clip = "off")

overlay_plot_m1
ggsave("model1_posterior_overlay_N254.png", overlay_plot_m1, width = 9, height = 5, dpi = 300)

# ============================================================================
# PART 7: FOREST PLOT (black-and-white, Arial style)
# ============================================================================

draws_inf <- as_draws_df(fit_informative_final)$b_RandomisationiCST
draws_skep <- as_draws_df(fit_skeptical_final)$b_RandomisationiCST

caterpillar_df <- data.frame(
  model = c("Informative (\u22121.92)", "Skeptical (0)"),
  median = c(median(draws_inf), median(draws_skep)),
  ci_lower = c(quantile(draws_inf, 0.025), quantile(draws_skep, 0.025)),
  ci_upper = c(quantile(draws_inf, 0.975), quantile(draws_skep, 0.975))
)

caterpillar_df$model <- factor(caterpillar_df$model,
                               levels = c("Skeptical (0)", "Informative (\u22121.92)"))

ggplot(caterpillar_df, aes(x = median, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, linewidth = 0.8, colour = "black") +
  geom_point(size = 3, colour = "black", fill = "black", shape = 21) +
  labs(
    x = "Treatment effect (ADAS-Cog points)",
    y = "Prior specification"
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    legend.position   = "none",
    axis.title        = element_text(size = 11, face = "bold"),
    axis.text         = element_text(size = 11, colour = "black"),
    axis.line         = element_line(colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3)
  )

ggsave("model1_forest_plot_N254.png", width = 10, height = 5, dpi = 300)

# ============================================================================
# PART 8: FULL RESULTS TABLE (TREATMENT EFFECT ONLY)
# ============================================================================

mcid_threshold <- -1.4   # replace with your study's actual MCID if different

summary_table <- tibble::tibble(
  Model = c("Informative", "Skeptical"),
  Prior = c("Normal(-1.92, 2)", "Normal(0, 2)"),
  Mean = c(mean(draws_inf), mean(draws_skep)),
  Median = c(median(draws_inf), median(draws_skep)),
  SD = c(sd(draws_inf), sd(draws_skep)),
  CI_lower_95 = c(quantile(draws_inf, 0.025), quantile(draws_skep, 0.025)),
  CI_upper_95 = c(quantile(draws_inf, 0.975), quantile(draws_skep, 0.975)),
  P_benefit = c(mean(draws_inf < 0), mean(draws_skep < 0)),
  P_harm = c(mean(draws_inf > 0), mean(draws_skep > 0)),
  P_clinically_meaningful = c(mean(draws_inf < mcid_threshold), mean(draws_skep < mcid_threshold)),
  Rhat = c(max(brms::rhat(fit_informative_final)), max(brms::rhat(fit_skeptical_final))),
  Bulk_ESS = c(summary(fit_informative_final)$fixed["RandomisationiCST", "Bulk_ESS"],
               summary(fit_skeptical_final)$fixed["RandomisationiCST", "Bulk_ESS"]),
  Tail_ESS = c(summary(fit_informative_final)$fixed["RandomisationiCST", "Tail_ESS"],
               summary(fit_skeptical_final)$fixed["RandomisationiCST", "Tail_ESS"])
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(summary_table)
write.csv(summary_table, "model1_full_results_table_N254.csv", row.names = FALSE)

results_gt <- summary_table %>%
  gt() %>%
  tab_header(
    title = "Model 1: Treatment Effect on ADAS-Cog by Prior Specification",
    subtitle = "iCST vs TAU Control, 26-week follow-up (N = 254)"
  ) %>%
  fmt_number(columns = c(Mean, Median, SD, CI_lower_95, CI_upper_95,
                         P_benefit, P_harm, P_clinically_meaningful, Rhat), decimals = 3) %>%
  fmt_number(columns = c(Bulk_ESS, Tail_ESS), decimals = 0) %>%
  tab_source_note(source_note = paste0(
    "P(benefit) = P(treatment effect < 0). P(clinically meaningful) = P(treatment effect < ",
    mcid_threshold, "). N = ", nrow(dat), "."))

results_gt
gtsave(results_gt, "model1_full_results_table_N254.png")

# ============================================================================
# PART 9: FULL PARAMETER TABLE (ALL COEFFICIENTS, BOTH MODELS)
# ============================================================================

extract_all_params <- function(fit, model_label) {
  fixed_df <- as.data.frame(summary(fit)$fixed)
  fixed_df$Parameter <- rownames(fixed_df)
  sigma_df <- as.data.frame(summary(fit)$spec_pars)
  sigma_df$Parameter <- rownames(sigma_df)
  combined <- bind_rows(fixed_df, sigma_df)
  combined$Model <- model_label
  combined %>% select(Model, Parameter, Estimate, `l-95% CI`, `u-95% CI`, Rhat, Bulk_ESS, Tail_ESS)
}

params_informative <- extract_all_params(fit_informative_final, "Informative")
params_skeptical   <- extract_all_params(fit_skeptical_final, "Skeptical")
all_params_table <- bind_rows(params_informative, params_skeptical)

all_params_table$Parameter <- dplyr::recode(all_params_table$Parameter,
                                            "Intercept" = "Intercept", "c_BASELINE_ADAScog" = "Baseline",
                                            "RandomisationiCST" = "Treatment", "sigma" = "Sigma")

param_order <- c("Intercept", "Baseline", "Sigma", "Treatment")
all_params_table$Parameter <- factor(all_params_table$Parameter, levels = param_order)
all_params_table$Model <- factor(all_params_table$Model, levels = c("Informative", "Skeptical"))

all_params_table <- all_params_table %>%
  arrange(Model, Parameter) %>%
  mutate(across(c(Estimate, `l-95% CI`, `u-95% CI`, Rhat), ~round(.x, 2)),
         across(c(Bulk_ESS, Tail_ESS), ~round(.x, 0)))

print(all_params_table)
write.csv(all_params_table, "model1_full_parameter_table_N254.csv", row.names = FALSE)

params_gt <- all_params_table %>%
  rename(`Lower CrI` = `l-95% CI`, `Upper CrI` = `u-95% CI`,
         `Bulk ESS` = Bulk_ESS, `Tail ESS` = Tail_ESS, `R-hat` = Rhat) %>%
  group_by(Model) %>%
  gt() %>%
  tab_header(title = "Table 2", subtitle = "Posterior Estimates of Treatment Effect on ADAS-Cog (N = 254)") %>%
  fmt_number(columns = c(Estimate, `Lower CrI`, `Upper CrI`), decimals = 2) %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Parameter)

params_gt
gtsave(params_gt, "model1_full_parameter_table_N254.png")

> # ------------------------------------------------------------------
> # 4. Between-subjects Welch's t-test: does the AMOUNT of change in
  > # QCPR differ between groups?
  > # ------------------------------------------------------------------

dat$qcpr_change <- dat$CQCPR_20_FU2 - dat$BASELINE_CQCPR_20
> t.test(qcpr_change ~ Randomisation, data = dat)

## ============================================================================
## iCST Model 1 (Strict 3 SD Rule, N=256): Bayesian Analysis of Treatment
## Effect on ADAS-Cog
## Outcome: ADAS_20_FU2 (26-week follow-up)
## Predictors: Baseline-centred ADAS-Cog, Randomisation (iCST vs TAU Control)
## Two-prior structure: Informative vs Skeptical
##
## Sample: the 3 participants missing baseline QCPR remain excluded
## (consistent with Model 2's requirements), PLUS only the 2 participants
## whose ADAS-Cog change scores exceeded a STRICT 3 SD threshold
## (qPIN 14201, z=4.15; qPIN 17018, z=4.33) are excluded. The two
## borderline cases (qPIN 13044, 11051; ~2.7 SD) are RETAINED, as they
## do not meet the pre-specified 3 SD criterion.
## ============================================================================

# ---- LIBRARIES ---------------------------------------------------------
library(readr)
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidybayes)
library(gt)

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================

# ---- Import ------------------------------------------------------------

dat <- read_csv("~/Downloads/THESIS COMBINED Complete case 26 week follow up data.csv")

dat$Randomisation <- factor(dat$Randomisation,
                            levels = c("TAU Control", "iCST"))

# ---- Exclude participants missing baseline/follow-up ADAS-Cog, and
# missing/zero baseline QCPR (0 treated as missing) -------------------------

dat <- dat %>%
  filter(
    !is.na(c_BASELINE_ADAScog),
    !is.na(ADAS_20_FU2),
    !is.na(BASELINE_CQCPR_20),
    BASELINE_CQCPR_20 != 0
  )

nrow(dat)   # should be 258
table(dat$Randomisation)

# ---- Exclude ONLY the 2 cases exceeding a strict 3 SD threshold -----------
# qPIN 14201 (z=4.15) and 17018 (z=4.33). qPIN 13044 and 11051 (~2.7 SD)
# are retained.

dat <- dat[!(dat$qPIN %in% c(14201, 17018)), ]

nrow(dat)   # should now be 256
table(dat$Randomisation)

# ============================================================================
# PART 2: MODEL SPECIFICATION
# ============================================================================

adas_formula_trunc <- bf(ADAS_20_FU2 | trunc(lb = 0, ub = 60) ~
                           c_BASELINE_ADAScog + Randomisation)

priors_informative_final <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

priors_skeptical_final <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(0, 2),     class = b, coef = RandomisationiCST),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

# ============================================================================
# PART 3: FIT BOTH MODELS
# ============================================================================

fit_informative_final <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_informative_final,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model1_informative_N256"
)

fit_skeptical_final <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_skeptical_final,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model1_skeptical_N256"
)

summary(fit_informative_final)
summary(fit_skeptical_final)

# ============================================================================
# PART 4: CONVERGENCE DIAGNOSTICS
# ============================================================================

cat("Max Rhat (Informative):", round(max(brms::rhat(fit_informative_final)), 4), "\n")
cat("Max Rhat (Skeptical):", round(max(brms::rhat(fit_skeptical_final)), 4), "\n")

# ============================================================================
# PART 5: POSTERIOR PROBABILITIES
# ============================================================================

hypothesis(fit_informative_final, "RandomisationiCST < 0")
hypothesis(fit_skeptical_final, "RandomisationiCST < 0")

# ============================================================================
# PART 6: FULL PARAMETER TABLE
# ============================================================================

extract_all_params <- function(fit, model_label) {
  fixed_df <- as.data.frame(summary(fit)$fixed)
  fixed_df$Parameter <- rownames(fixed_df)
  sigma_df <- as.data.frame(summary(fit)$spec_pars)
  sigma_df$Parameter <- rownames(sigma_df)
  combined <- bind_rows(fixed_df, sigma_df)
  combined$Model <- model_label
  combined %>% select(Model, Parameter, Estimate, `l-95% CI`, `u-95% CI`, Rhat, Bulk_ESS, Tail_ESS)
}

params_informative <- extract_all_params(fit_informative_final, "Informative")
params_skeptical   <- extract_all_params(fit_skeptical_final, "Skeptical")
all_params_table <- bind_rows(params_informative, params_skeptical)

all_params_table$Parameter <- dplyr::recode(all_params_table$Parameter,
                                            "Intercept" = "Intercept", "c_BASELINE_ADAScog" = "Baseline",
                                            "RandomisationiCST" = "Treatment", "sigma" = "Sigma")

param_order <- c("Intercept", "Baseline", "Sigma", "Treatment")
all_params_table$Parameter <- factor(all_params_table$Parameter, levels = param_order)
all_params_table$Model <- factor(all_params_table$Model, levels = c("Informative", "Skeptical"))

draws_treat_inf <- as_draws_df(fit_informative_final)$b_RandomisationiCST
draws_treat_skep <- as_draws_df(fit_skeptical_final)$b_RandomisationiCST

prob_lookup <- tibble::tibble(
  Model = c("Informative", "Skeptical"),
  Parameter = "Treatment",
  P_benefit = c(mean(draws_treat_inf < 0), mean(draws_treat_skep < 0))
)

all_params_table <- all_params_table %>%
  left_join(prob_lookup, by = c("Model", "Parameter")) %>%
  arrange(Model, Parameter) %>%
  mutate(across(c(Estimate, `l-95% CI`, `u-95% CI`, Rhat, P_benefit), ~round(.x, 2)),
         across(c(Bulk_ESS, Tail_ESS), ~round(.x, 0)))

print(all_params_table)
write.csv(all_params_table, "model1_full_parameter_table_N256.csv", row.names = FALSE)

params_gt <- all_params_table %>%
  rename(`Lower CrI` = `l-95% CI`, `Upper CrI` = `u-95% CI`,
         `Bulk ESS` = Bulk_ESS, `Tail ESS` = Tail_ESS, `R-hat` = Rhat,
         `Posterior Probability` = P_benefit) %>%
  group_by(Model) %>%
  gt() %>%
  tab_header(title = "Table 2", subtitle = "Posterior Estimates of Treatment Effect on ADAS-Cog (N = 256)") %>%
  fmt_number(columns = c(Estimate, `Lower CrI`, `Upper CrI`, `Posterior Probability`), decimals = 2) %>%
  sub_missing(columns = `Posterior Probability`, missing_text = "\u2014") %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Parameter) %>%
  tab_source_note(
    source_note = "Posterior Probability = P(Treatment < 0). Not applicable (\u2014) for Intercept, Baseline, and Sigma."
  )

params_gt
gtsave(params_gt, "model1_full_parameter_table_N256.png")

## ------------------------------------------------------------------
## Model 1 (Strict 3 SD Rule, N=256): Prior and Posterior
## Predictive Checks, both prior specifications.
## ------------------------------------------------------------------

library(brms)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------
# 1. Fit prior-only models (sample_prior = "only")
# ------------------------------------------------------------------

prior_check_informative_final <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_informative_final,
  sample_prior = "only",
  chains  = 4,
  iter    = 1000,
  seed    = 42
)

prior_check_skeptical_final <- brm(
  formula = adas_formula_trunc,
  data    = dat,
  family  = gaussian(),
  prior   = priors_skeptical_final,
  sample_prior = "only",
  chains  = 4,
  iter    = 1000,
  seed    = 42
)

# ------------------------------------------------------------------
# 2. Prior predictive check INFORMATIVE - side by side
# ------------------------------------------------------------------

p_prior_inf <- brms::pp_check(prior_check_informative_final, ndraws = 100) +
  labs(title = "Prior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_inf <- brms::pp_check(fit_informative_final, ndraws = 100) +
  labs(title = "Posterior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_informative <- p_prior_inf | p_post_inf
combined_informative

ggsave("model1_ppc_informative_N256.png", combined_informative,
       width = 10, height = 5, dpi = 300)

# ------------------------------------------------------------------
# 3. Prior predictive check SKEPTICAL - side by side
# ------------------------------------------------------------------

p_prior_skep <- brms::pp_check(prior_check_skeptical_final, ndraws = 100) +
  labs(title = "Prior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_skep <- brms::pp_check(fit_skeptical_final, ndraws = 100) +
  labs(title = "Posterior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_skeptical <- p_prior_skep | p_post_skep
combined_skeptical

ggsave("model1_ppc_skeptical_N256.png", combined_skeptical,
       width = 10, height = 5, dpi = 300)

# ------------------------------------------------------------------
# 4. Traceplots
# ------------------------------------------------------------------

library(bayesplot)
library(ggplot2)

trace_plot_inf <- mcmc_trace(fit_informative_final, pars = c("b_Intercept", "b_c_BASELINE_ADAScog",
                                                             "b_RandomisationiCST", "sigma")) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 11, colour = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

trace_plot_inf

ggsave("model1_traceplots_informative_N256.png", trace_plot_inf, width = 10, height = 8, dpi = 300)


library(bayesplot)
library(ggplot2)

trace_plot_inf <- mcmc_trace(fit_informative_final, pars = c("b_Intercept", "b_c_BASELINE_ADAScog",
                                                             "b_RandomisationiCST", "sigma")) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 11, colour = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

trace_plot_inf

ggsave("model1_traceplots_informative_N256.png", trace_plot_inf, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------------
# 5. Caterpillar Plot
# ------------------------------------------------------------------

library(dplyr)
library(ggplot2)

draws_inf <- as_draws_df(fit_informative_final)$b_RandomisationiCST
draws_skep <- as_draws_df(fit_skeptical_final)$b_RandomisationiCST

caterpillar_df <- data.frame(
  model = c("Informative (\u22121.92)", "Skeptical (0)"),
  median = c(median(draws_inf), median(draws_skep)),
  ci_lower = c(quantile(draws_inf, 0.025), quantile(draws_skep, 0.025)),
  ci_upper = c(quantile(draws_inf, 0.975), quantile(draws_skep, 0.975))
)

caterpillar_df$model <- factor(caterpillar_df$model,
                               levels = c("Skeptical (0)", "Informative (\u22121.92)"))

ggplot(caterpillar_df, aes(x = median, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, linewidth = 0.8, colour = "black") +
  geom_point(size = 3, colour = "black", fill = "black", shape = 21) +
  labs(
    x = "Treatment effect (ADAS-Cog points)",
    y = "Prior specification"
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    legend.position   = "none",
    axis.title        = element_text(size = 11, face = "bold"),
    axis.text         = element_text(size = 11, colour = "black"),
    axis.line         = element_line(colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3)
  )

ggsave("model1_forest_plot_N256.png", width = 10, height = 5, dpi = 300)

# ------------------------------------------------------------------
# 6. Overlay
# ------------------------------------------------------------------

library(dplyr)
library(ggplot2)

draws_informative <- as_draws_df(fit_informative_final) %>%
  select(b_RandomisationiCST) %>%
  mutate(Prior = "Informative")

draws_skeptical <- as_draws_df(fit_skeptical_final) %>%
  select(b_RandomisationiCST) %>%
  mutate(Prior = "Skeptical")

all_draws_m1 <- bind_rows(draws_informative, draws_skeptical)
all_draws_m1$Prior <- factor(all_draws_m1$Prior, levels = c("Informative", "Skeptical"))

prob_labels_m1 <- all_draws_m1 %>%
  group_by(Prior) %>%
  summarise(prob_benefit = mean(b_RandomisationiCST < 0), .groups = "drop") %>%
  mutate(label = paste0(Prior, ": P(benefit) = ",
                        scales::percent(prob_benefit, accuracy = 0.1)))

print(prob_labels_m1)

overlay_plot_m1 <- ggplot(all_draws_m1, aes(x = b_RandomisationiCST,
                                            fill = Prior, color = Prior)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  annotate("text", x = 0.15, y = Inf, label = "No effect",
           color = "black", hjust = 0, vjust = 1.5, size = 3.5) +
  scale_fill_manual(
    values = c("Informative" = "#4E9FDA", "Skeptical" = "#4CAF50"),
    labels = c("Informative" = "Informative (\u22121.92)",
               "Skeptical" = "Skeptical (0)")
  ) +
  scale_color_manual(
    values = c("Informative" = "#4E9FDA", "Skeptical" = "#4CAF50"),
    labels = c("Informative" = "Informative (\u22121.92)",
               "Skeptical" = "Skeptical (0)")
  ) +
  labs(
    x = "Treatment effect (ADAS-Cog points)",
    y = "Density",
    fill = "Prior specification",
    color = "Prior specification"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  coord_cartesian(clip = "off")

overlay_plot_m1

ggsave("model1_posterior_overlay_N256.png", overlay_plot_m1, width = 9, height = 5, dpi = 300)

## ------------------------------------------------------------------
## Model 1: Triplot (informative prior), colour, APA 7 style
## Blue = informative prior, red = likelihood, green = posterior
## Requires these objects in your R session:
##   fit_informative_final  (brms fit, N = 256)
##   dat                    (the N = 256 analysis dataset)
## ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(brms)

# 1. Posterior: draws from the informative-prior model (N = 256)
post_draws <- as_draws_df(fit_informative_final)$b_RandomisationiCST

# 2. Informative prior: Normal(-1.92, 2)
prior_mean <- -1.92
prior_sd   <- 2

# 3. Likelihood: normal curve centred on the maximum likelihood estimate
lm_fit  <- lm(ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation, data = dat)
lik_est <- coef(summary(lm_fit))["RandomisationiCST", "Estimate"]
lik_se  <- coef(summary(lm_fit))["RandomisationiCST", "Std. Error"]

# 4. Common grid
x <- seq(-6, 4, length.out = 1000)

dist_labels <- c("Informative prior", "Likelihood", "Posterior")

tri_df <- data.frame(
  x = rep(x, 3),
  density = c(
    dnorm(x, prior_mean, prior_sd),
    dnorm(x, lik_est, lik_se),
    density(post_draws, from = min(x), to = max(x), n = length(x))$y
  ),
  Distribution = factor(rep(dist_labels, each = length(x)), levels = dist_labels)
)

# 5. APA 7 triplot in colour
triplot_apa_colour <- ggplot(tri_df, aes(x = x, y = density, colour = Distribution)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("Informative prior" = "#4E9FDA",
                                 "Likelihood"        = "#D32F2F",
                                 "Posterior"         = "#4CAF50")) +
  labs(x = "Treatment effect (ADAS-Cog points)", y = "Density", colour = NULL) +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(
    text = element_text(family = "Arial", size = 11, colour = "black"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.5, "cm")
  )

triplot_apa_colour

ggsave("model1_triplot_informative_N256_APA_colour.png", triplot_apa_colour,
       width = 6.5, height = 4, dpi = 300)


  ## ------------------------------------------------------------------
> ## Descriptive statistics: BASELINE_ADAS_20
  > ## Requires `dat` (the N = 256 analysis dataset) in your R session
  > ## ------------------------------------------------------------------
> 
  > library(dplyr)
> 
  > describe_adas <- function(d) {
    +     d %>%
      +         summarise(
        +             n       = sum(!is.na(BASELINE_ADAS_20)),
        +             missing = sum(is.na(BASELINE_ADAS_20)),
        +             M       = mean(BASELINE_ADAS_20, na.rm = TRUE),
        +             SD      = sd(BASELINE_ADAS_20, na.rm = TRUE),
        +             Mdn     = median(BASELINE_ADAS_20, na.rm = TRUE),
        +             Q1      = quantile(BASELINE_ADAS_20, 0.25, na.rm = TRUE),
        +             Q3      = quantile(BASELINE_ADAS_20, 0.75, na.rm = TRUE),
        +             Min     = min(BASELINE_ADAS_20, na.rm = TRUE),
        +             Max     = max(BASELINE_ADAS_20, na.rm = TRUE),
        +             .groups = "drop"
        +         ) %>%
      +         mutate(across(c(M, SD, Mdn, Q1, Q3, Min, Max), ~ round(.x, 2)))
    + }
> 
  > # 1. Whole sample
  > overall <- describe_adas(dat) %>% mutate(Group = "Total", .before = 1)
> 
  > # 2. By trial arm
  > by_arm <- dat %>%
  +     group_by(Group = as.character(Randomisation)) %>%
  +     describe_adas()
> 
  > adas_desc <- bind_rows(by_arm, overall)
> print(adas_desc)

## ------------------------------------------------------------------
> ## Descriptive statistics: BASELINE_CQCPR_20
  > ## Requires `dat` (the N = 256 analysis dataset) in your R session
  > ## ------------------------------------------------------------------
> 
  > library(dplyr)
> 
  > # Skewness as reported by SPSS (adjusted Fisher-Pearson, G1)
  > skew_spss <- function(x) {
    +     x <- x[!is.na(x)]; n <- length(x)
    +     g1 <- (sum((x - mean(x))^3) / n) / (sum((x - mean(x))^2) / n)^(3/2)
    +     g1 * sqrt(n * (n - 1)) / (n - 2)
    + }
> 
  > # Standard error of skewness
  > skew_se <- function(x) {
    +     n <- sum(!is.na(x))
    +     sqrt(6 * n * (n - 1) / ((n - 2) * (n + 1) * (n + 3)))
    + }
> 
  > describe_cqcpr <- function(d) {
    +     d %>%
      +         summarise(
        +             n       = sum(!is.na(BASELINE_CQCPR_20)),
        +             missing = sum(is.na(BASELINE_CQCPR_20)),
        +             M       = mean(BASELINE_CQCPR_20, na.rm = TRUE),
        +             SD      = sd(BASELINE_CQCPR_20, na.rm = TRUE),
        +             Mdn     = median(BASELINE_CQCPR_20, na.rm = TRUE),
        +             Q1      = quantile(BASELINE_CQCPR_20, 0.25, na.rm = TRUE),
        +             Q3      = quantile(BASELINE_CQCPR_20, 0.75, na.rm = TRUE),
        +             Min     = min(BASELINE_CQCPR_20, na.rm = TRUE),
        +             Max     = max(BASELINE_CQCPR_20, na.rm = TRUE),
        +             Skew    = skew_spss(BASELINE_CQCPR_20),
        +             SE_skew = skew_se(BASELINE_CQCPR_20),
        +             .groups = "drop"
        +         ) %>%
      +         mutate(z_skew = Skew / SE_skew) %>%
      +         mutate(across(c(M, SD, Mdn, Q1, Q3, Min, Max, Skew, SE_skew, z_skew), ~ round(.x, 2)))
    + }
> 
  > # 1. Whole sample
  > overall <- describe_cqcpr(dat) %>% mutate(Group = "Total", .before = 1)
> 
  > # 2. By trial arm
  > by_arm <- dat %>%
  +     group_by(Group = as.character(Randomisation)) %>%
  +     describe_cqcpr()
> 
  > cqcpr_desc <- bind_rows(by_arm, overall)
> print(cqcpr_desc)

library(dplyr)

# ---- 0. Check the raw gender codes first ----
table(dat$qRelGenP, useNA = "ifany")
table(dat$qRelGenC, useNA = "ifany")

# ---- 1. Label gender (1 = Male, 2 = Female) ----
dat <- dat %>%
  mutate(
    gender_P = factor(qRelGenP, levels = c(1, 2), labels = c("Male", "Female")),
    gender_C = factor(qRelGenC, levels = c(1, 2), labels = c("Male", "Female"))
  )

# ---- 2. Age descriptives ----
describe_age <- function(d, var) {
  d %>%
    summarise(
      n      = sum(!is.na({{ var }})),
      Mean   = mean({{ var }}, na.rm = TRUE),
      Median = median({{ var }}, na.rm = TRUE),
      SD     = sd({{ var }}, na.rm = TRUE),
      Min    = min({{ var }}, na.rm = TRUE),
      Max    = max({{ var }}, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(c(Mean, Median, SD, Min, Max), ~ round(.x, 2)))
}

age_table <- function(var, who) {
  bind_rows(
    dat %>% group_by(Group = as.character(Randomisation)) %>% describe_age({{ var }}),
    describe_age(dat, {{ var }}) %>% mutate(Group = "Total", .before = 1)
  ) %>%
    mutate(Participant = who, .before = 1)
}

age_desc <- bind_rows(
  age_table(cRelAgeP, "Person with dementia"),
  age_table(cRelAgeC, "Caregiver")
)
print(as.data.frame(age_desc))
View(age_desc, "Age")

# ---- 3. Gender frequencies ----
gender_table <- function(var, who) {
  freq <- function(d) {
    d %>%
      count(Gender = {{ var }}) %>%
      mutate(
        pct_valid = ifelse(is.na(Gender), NA,
                           round(100 * n / sum(n[!is.na(Gender)]), 1)),
        Gender = ifelse(is.na(Gender), "Missing", as.character(Gender))
      )
  }
  bind_rows(
    dat %>% group_by(Group = as.character(Randomisation)) %>% freq(),
    freq(dat) %>% mutate(Group = "Total", .before = 1)
  ) %>%
    ungroup() %>%
    mutate(Participant = who, .before = 1)
}

gender_freq <- bind_rows(
  gender_table(gender_P, "Person with dementia"),
  gender_table(gender_C, "Caregiver")
)
print(as.data.frame(gender_freq))
View(gender_freq, "Gender")

# ---- 4. Save ----
write.csv(age_desc,    "age_descriptives.csv",   row.names = FALSE)
write.csv(gender_freq, "gender_frequencies.csv", row.names = FALSE)

library(dplyr)

# ---- 0. Frequency helper: n and % of valid responses, by arm and total ----
freq_table <- function(var, who) {
  freq <- function(d) {
    d %>%
      count(Category = {{ var }}) %>%
      mutate(
        pct_valid = ifelse(is.na(Category), NA,
                           round(100 * n / sum(n[!is.na(Category)]), 1)),
        Category = ifelse(is.na(Category), "Missing", as.character(Category))
      )
  }
  bind_rows(
    dat %>% group_by(Group = as.character(Randomisation)) %>% freq(),
    freq(dat) %>% mutate(Group = "Total", .before = 1)
  ) %>%
    ungroup() %>%
    mutate(Variable = who, .before = 1)
}

# ---- 1. Raw codes (check against your codebook) ----
raw_freq <- bind_rows(
  freq_table(qRelationship, "Relationship (raw code)"),
  freq_table(cRelEthP,      "Ethnicity PwD (raw code)"),
  freq_table(cRelEthC,      "Ethnicity caregiver (raw code)")
)
print(as.data.frame(raw_freq))

# ---- 2. Labels ----
eth_codes  <- c(1, 2, 3, 4, 5, 6, 13, 15, 22)
eth_labels <- c("White (not specified)", "White British", "White Irish",
                "Other White background", "Black or Black British (not specified)",
                "Caribbean", "Other mixed background", "Indian",
                "Not stated / do not wish to specify")

dat <- dat %>%
  mutate(
    relationship_grp = factor(
      qRelationship,
      levels = c(1, 2, 3, 4, 5, 6, 7, 9),
      labels = c("Wife/husband (spouse)", "Partner", "Son/daughter",
                 "Son/daughter-in-law", "Brother/sister", "Other relative",
                 "Friend", "Other")
    ),
    ethnicity_P = coalesce(eth_labels[match(cRelEthP, eth_codes)],
                           ifelse(is.na(cRelEthP), NA, paste("Code", cRelEthP))),
    ethnicity_C = coalesce(eth_labels[match(cRelEthC, eth_codes)],
                           ifelse(is.na(cRelEthC), NA, paste("Code", cRelEthC)))
  ) %>%
  mutate(across(c(ethnicity_P, ethnicity_C),
                ~ factor(.x, levels = unique(c(eth_labels,
                                               sort(unique(.x[!.x %in% eth_labels])))))))

# ---- 3. Table 1: relationship type ----
relationship_freq <- freq_table(relationship_grp, "Relationship to person with dementia")

# ---- 4. Table 2: ethnicity ----
ethnicity_freq <- bind_rows(
  freq_table(ethnicity_P, "Person with dementia"),
  freq_table(ethnicity_C, "Caregiver")
) %>%
  rename(Participant = Variable)

print(as.data.frame(relationship_freq))
print(as.data.frame(ethnicity_freq))

# Each opens in its own RStudio tab
View(relationship_freq, "Relationship")
View(ethnicity_freq,    "Ethnicity")

# ---- 5. Save: one Excel workbook, a tab for each table ----
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
writexl::write_xlsx(
  list(
    Relationship = relationship_freq,
    Ethnicity    = ethnicity_freq,
    `Raw codes`  = raw_freq
  ),
  "relationship_ethnicity_tables.xlsx"
)

---- 
  library(dplyr)

# ---- 0. Find the baseline MMSE column (ignores capitals) ----
mmse_cols <- grep("mmse", names(dat), ignore.case = TRUE, value = TRUE)
print(mmse_cols)

baseline_cols <- grep("base", mmse_cols, ignore.case = TRUE, value = TRUE)
mmse_var <- if (length(baseline_cols) >= 1) baseline_cols[1] else mmse_cols[1]

if (is.na(mmse_var)) {
  stop("No column with 'mmse' in its name was found in dat. ",
       "Run names(dat) to find the MMSE column, then set mmse_var by hand, ",
       "e.g. mmse_var <- \"YOUR_COLUMN_NAME\"")
}
message("Using column: ", mmse_var)

dat <- dat %>% mutate(mmse = .data[[mmse_var]])

# ---- 1. Skewness (SPSS formula) and its SE ----
skew_spss <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  g1 <- (sum((x - mean(x))^3) / n) / (sum((x - mean(x))^2) / n)^(3/2)
  g1 * sqrt(n * (n - 1)) / (n - 2)
}
skew_se <- function(x) {
  n <- sum(!is.na(x))
  sqrt(6 * n * (n - 1) / ((n - 2) * (n + 1) * (n + 3)))
}

# ---- 2. Descriptives ----
describe_mmse <- function(d) {
  d %>%
    summarise(
      n       = sum(!is.na(mmse)),
      Missing = sum(is.na(mmse)),
      Mean    = mean(mmse, na.rm = TRUE),
      Median  = median(mmse, na.rm = TRUE),
      SD      = sd(mmse, na.rm = TRUE),
      Min     = min(mmse, na.rm = TRUE),
      Max     = max(mmse, na.rm = TRUE),
      Skew    = skew_spss(mmse),
      SE_skew = skew_se(mmse),
      .groups = "drop"
    ) %>%
    mutate(across(c(Mean, Median, SD, Min, Max, Skew, SE_skew), ~ round(.x, 2)))
}

mmse_desc <- bind_rows(
  dat %>% group_by(Group = as.character(Randomisation)) %>% describe_mmse(),
  describe_mmse(dat) %>% mutate(Group = "Total", .before = 1)
)

# ---- 3. Frequencies ----
freq_table <- function(var) {
  freq <- function(d) {
    d %>%
      count(Category = {{ var }}) %>%
      mutate(
        pct_valid = ifelse(is.na(Category), NA,
                           round(100 * n / sum(n[!is.na(Category)]), 1)),
        Category = ifelse(is.na(Category), "Missing", as.character(Category))
      )
  }
  bind_rows(
    dat %>% group_by(Group = as.character(Randomisation)) %>% freq(),
    freq(dat) %>% mutate(Group = "Total", .before = 1)
  ) %>%
    ungroup()
}

# 3a. Every MMSE score
mmse_score_freq <- freq_table(mmse) %>% rename(MMSE_score = Category)

# 3b. MMSE bands: 17 or less, 18-23, 24 and above
dat <- dat %>%
  mutate(mmse_band = cut(mmse,
                         breaks = c(-Inf, 17, 23, Inf),
                         labels = c("17 or less", "18-23", "24 and above")))
mmse_band_freq <- freq_table(mmse_band) %>% rename(Severity = Category)

# ---- 4. Show ----
print(as.data.frame(mmse_desc))
print(as.data.frame(mmse_band_freq))
print(as.data.frame(mmse_score_freq))

View(mmse_desc,       "MMSE descriptives")
View(mmse_band_freq,  "MMSE bands")
View(mmse_score_freq, "MMSE scores")

# ---- 5. Save: one workbook, one tab per table ----
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
writexl::write_xlsx(
  list(
    Descriptives = mmse_desc,
    Bands        = mmse_band_freq,
    Scores       = mmse_score_freq
  ),
  "baseline_MMSE_tables.xlsx"
)