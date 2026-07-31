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
# Create Demographics Table for Model 1
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

