# =============================================================================
# Set-up
# =============================================================================

# -- Load packages ---------------------------------------------------------
library(ggtext)
library(parameters)
library(emmeans)
library(patchwork)
library(ggdist)
library(purrr)
library(tidyverse)

# -- Set common theme ------------------------------------------------------
common_theme <- theme(
  panel.border = element_rect(colour = "grey", fill = NA, linewidth = 1),
  panel.background = element_blank(),
  axis.ticks.y = element_line(colour = "grey", linewidth = 0.5),
  axis.ticks.length.y = unit(-0.1, "cm"),
  axis.ticks.x = element_blank(),
  legend.position = "none",
  axis.text = element_markdown(size = 8),
  axis.title = element_markdown(size = 10),
  plot.title = element_markdown(size = 9, hjust = 0.5, vjust = 0),
  plot.tag = element_markdown(size = 12),
  plot.tag.location = "plot",
  plot.tag.position = c(0.01, 0.97),
  plot.margin = margin(1, 1, 1, 1),
  aspect.ratio = 1
  )

# -- Shared colour and fill scales --------------------------------------------
colour_scale <- scale_colour_manual(values = c(Baseline = "grey80", Spring = "#1b9e77", Summer = "#d95f02"))
fill_scale   <- scale_fill_manual(  values = c(Baseline = "grey80", Spring = "#1b9e77", Summer = "#d95f02"))

# -- Bootstrap function -------------------------------------------------------
# Stratified case resampling: rows are resampled independently within each
# treatment group (preserving group sizes). The model is refit on each
# resample and emmeans extracted on the response scale. Effect size is
# defined as treatment mean - baseline mean (grams).

bootstrap_effects <- function(data_sp, log_transform = FALSE, n_boot = 10000) {
  map_dfr(seq_len(n_boot), function(i) {

    boot_data <- data_sp %>%
      group_by(treatment) %>%
      slice_sample(prop = 1, replace = TRUE) %>%
      ungroup()

    if (log_transform) {
      fit   <- lm(log(dry_uso) ~ treatment, data = boot_data)
      means <- emmeans(fit, ~ treatment, type = "response") %>%
        as.data.frame() %>%
        select(treatment, mean = response)
    } else {
      fit   <- lm(dry_uso ~ treatment, data = boot_data)
      means <- emmeans(fit, ~ treatment) %>%
        as.data.frame() %>%
        select(treatment, mean = emmean)
    }

    baseline_val <- means %>% filter(treatment == "Baseline") %>% pull(mean)

    means %>%
      mutate(
        effect        = mean - baseline_val,
        baseline_mean = baseline_val,
        boot_id       = i
      ) %>%
      select(treatment, effect, baseline_mean, boot_id)
  })
}

# -- Effects plot function ----------------------------------------------------

make_boot_plot <- function(sp_name, sp_label, show_y = TRUE, tag = NULL) {

  # Map treatment to numeric x positions to match panel a's discrete scale
  # (ggplot2 discrete: Baseline = 1, Spring = 2, Summer = 3)
  boot_sp <- boot_all %>%
    filter(species == sp_name, treatment != "Baseline") %>%
    mutate(x_pos = as.numeric(treatment))

  summ_sp <- boot_summary %>%
    filter(species == sp_name) %>%
    mutate(x_pos = as.numeric(treatment))

  ribbon_sp <- baseline_ribbon %>% filter(species == sp_name)

  p <- ggplot() +
    # Baseline ribbon: 95% CI of baseline mean centred at zero
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = ribbon_sp$ribbon_low, ymax = ribbon_sp$ribbon_high,
      fill = "grey80", alpha = 0.5
    ) +
    # Zero reference line
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40", linewidth = 0.6) +
    # Density slab of all 10,000 bootstrap estimates (ggdist)
    stat_slab(
      data = boot_sp,
      aes(x = x_pos, y = effect, fill = treatment),
      alpha     = 0.5,
      scale     = 0.6,
      normalize = "groups",
      adjust    = 3,
      colour    = NA
    ) +
    # 95% empirical CI bar
    geom_errorbar(
      data = summ_sp,
      aes(x = x_pos, ymin = ci_low, ymax = ci_high),
      width     = 0,
      linewidth = 1
    ) +
    # Mean bootstrap estimate
    geom_point(
      data = summ_sp,
      aes(x = x_pos, y = mean_effect, fill = treatment),
      shape = 21, size = 2, stroke = 1
    ) +
    scale_x_continuous(
      breaks = 1:3,
      labels = c("Baseline", "Spring", "Summer")
    ) +
    scale_y_continuous(limits = boot_limits) +
    coord_cartesian(xlim = c(0.5, 3.5)) +
    colour_scale +
    fill_scale +
    labs(
      x     = NULL,
      y     = if (show_y) "Effect size" else NULL,
      title = sp_label,
      tag   = tag
    ) +
    common_theme +
    theme(
      plot.title = element_markdown(colour = "white")
    )

  if (!show_y) {
    p <- p + theme(
      axis.text.y  = element_blank()
    )
  }

  p
}

# -- Read data -------------------------------------------------------------
data <- data.table::fread("data/uso_response.csv") %>%
  filter(
    # Remove NA treatments
    !is.na(treatment)
    ) %>%
  mutate(
    treatment = factor(treatment, levels = c("Baseline", "Spring", "Summer"))
  ) %>%
  glimpse()

# =============================================================================
# Panel (a)
# =============================================================================

# y-axis limits for panel (a)
limits <- c(0, 8)

# -- Arthropodium strictum  -----------------------------------------------------

# Create df for arthropodium strictum 
data_arthropodium <- data %>%
  filter(name == "Arthropodium strictum") %>%
  select(treatment, dry_uso) %>%
  na.omit()

# Fit anova explaining variance in dry weight as a funtion of treatment
model_arthropodium <- lm(log(dry_uso) ~ treatment, data = data_arthropodium)

# Check model assumptions
DHARMa::simulateResiduals(model_arthropodium, plot = TRUE)

# Results
summary_arthropodium <- parameters(model_arthropodium) %>%
  as.data.frame() %>%
  mutate(
    `95% CI` = paste0("[", round(CI_low, 2), ", ", round(CI_high, 2), "]")
  )

# Estimated means the results 
mean_arthropodium <- emmeans(model_arthropodium, ~ treatment, type = "response") %>%
  as.data.frame()

# Plot
plot_arthropodium <- ggplot() +
  geom_point(
    data = data_arthropodium,
    aes(x = treatment, y = dry_uso, colour = treatment),
    shape = 20,
    size = 2,
    alpha = 0.5,
    position = position_jitter(width = 0.3, seed = 1969)
  ) +
  # Plot the fitted means and confidence intervals (model performance)
  geom_errorbar(
    data = mean_arthropodium,
    aes(x = treatment, ymin = lower.CL, ymax = upper.CL),
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = mean_arthropodium,
    aes(x = treatment, y = response, fill = treatment),
    shape = 21,
    size = 2,
    stroke = 1
  ) +
  scale_y_continuous(limits = limits) +
  colour_scale +
  fill_scale +
  labs(
    x = NULL,
    y = "Dry weight (g)",
    title = "***Arthropodium strictum***",
    tag   = "(**a**)"
  ) +
  common_theme +
  theme(
    axis.text.x = element_blank()
  )

# Display the plot
print(plot_arthropodium)

# -- Bulbine bulbosa -----------------------------------------------------------

# Bulbine

# Create df for Bulbine bulbosa 
data_bulbine <- data %>%
  filter(name == "Bulbine bulbosa") %>%
  select(treatment, dry_uso) %>%
  na.omit()

# Fit anova explaining variance in dry weight as a funtion of treatment
model_bulbine <- lm(dry_uso ~ treatment, data = data_bulbine)

# Check model assumptions
DHARMa::simulateResiduals(model_bulbine, plot = TRUE)

# Results
summary_bulbine <- parameters(model_bulbine) %>%
  as.data.frame() %>%
  mutate(
   `95% CI` = paste0("[", round(CI_low, 2), ", ", round(CI_high, 2), "]")
  )

# Summarise the results
mean_bulbine <- emmeans(model_bulbine, ~ treatment) %>%
  as.data.frame()

# Plot
plot_bulbine <- ggplot() +
  geom_point(
    data = data_bulbine,
    aes(x = treatment, y = dry_uso, colour = treatment),
    shape = 20,
    size = 2,
    alpha = 0.5,
    position = position_jitter(width = 0.3, seed = 1969)
) +
  # Plot the fitted means and confidence intervals (model performance)
  geom_errorbar(
    data = mean_bulbine,
    aes(x = treatment, ymin = lower.CL, ymax = upper.CL),
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = mean_bulbine,
    aes(x = treatment, y = emmean, fill = treatment),
    shape = 21,
    size = 2,
    stroke = 1
  ) +
  scale_y_continuous(limits = limits) +
  colour_scale +
  fill_scale +
  labs(
    x = NULL,
    y = NULL,
    title = "***Bulbine bulbosa***"
  ) +
  common_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )

# Display the plot
print(plot_bulbine)

# -- Burchardia umberllata ----------------------------------------------------

# Create df for Burchardia umberllata
data_burchardia <- data %>%
  filter(name == "Burchardia umberllata") %>%
  select(treatment, dry_uso) %>%
  na.omit()

# Fit anova explaining variance in dry weight as a function of treatment
model_burchardia <- lm(dry_uso ~ treatment, data = data_burchardia)

# Check model assumptions
DHARMa::simulateResiduals(model_burchardia, plot = TRUE)

# Results
summary_burchardia <- parameters(model_burchardia) %>%
  as.data.frame() %>%
  mutate(
    `95% CI` = paste0("[", round(CI_low, 2), ", ", round(CI_high, 2), "]")
  )

# Estimated means
mean_burchardia <- emmeans(model_burchardia, ~ treatment) %>%
  as.data.frame()

# Plot
plot_burchardia <- ggplot() +
  geom_point(
    data = data_burchardia,
    aes(x = treatment, y = dry_uso, colour = treatment),
    shape = 20,
    size = 2,
    alpha = 0.5,
    position = position_jitter(width = 0.3, seed = 1969)
  ) +
  # Plot the fitted means and confidence intervals (model performance)
  geom_errorbar(
    data = mean_burchardia,
    aes(x = treatment, ymin = lower.CL, ymax = upper.CL),
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = mean_burchardia,
    aes(x = treatment, y = emmean, fill = treatment),
    shape = 21,
    size = 2,
    stroke = 1
  ) +
  scale_y_continuous(limits = limits) +
  colour_scale +
  fill_scale +
  labs(
    x = NULL,
    y = NULL,
    title = "***Burchardia umbellata***"
  ) +
  common_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )

# Display the plot
print(plot_burchardia)

# -- Eryngium ovinum ----------------------------------------------------------

# Create df for Eryngium ovinum
data_eryngium <- data %>%
  filter(name == "Eryngium ovinum") %>%
  select(treatment, dry_uso) %>%
  na.omit()

# Fit anova explaining variance in dry weight as a function of treatment
model_eryngium <- lm(log(dry_uso) ~ treatment, data = data_eryngium)

# Check model assumptions
DHARMa::simulateResiduals(model_eryngium, plot = TRUE)

# Results
summary_eryngium <- parameters(model_eryngium) %>%
  as.data.frame() %>%
  mutate(
    `95% CI` = paste0("[", round(CI_low, 2), ", ", round(CI_high, 2), "]")
  )

# Estimated means on response scale
mean_eryngium <- emmeans(model_eryngium, ~ treatment, type = "response") %>%
  as.data.frame()

# Plot
plot_eryngium <- ggplot() +
  geom_point(
    data = data_eryngium,
    aes(x = treatment, y = dry_uso, colour = treatment),
    shape = 20,
    size = 2,
    alpha = 0.5,
    position = position_jitter(width = 0.3, seed = 1969)
  ) +
  # Plot the fitted means and confidence intervals (model performance)
  geom_errorbar(
    data = mean_eryngium,
    aes(x = treatment, ymin = lower.CL, ymax = upper.CL),
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = mean_eryngium,
    aes(x = treatment, y = response, fill = treatment),
    shape = 21,
    size = 2,
    stroke = 1
  ) +
  scale_y_continuous(limits = limits) +
  colour_scale +
  fill_scale +
  labs(
    x = NULL,
    y = NULL,
    title = "***Eryngium ovinum***"
  ) +
  common_theme +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
  )

# Display the plot
print(plot_eryngium)

# -- Combine the plots -------------------------------------------------------

panel_a <- wrap_plots(
  plot_arthropodium,
  plot_bulbine,
  plot_burchardia,
  plot_eryngium,
  nrow = 1
)
print(panel_a)

# ============================================================================
# Bootstrap effect sizes (panel b)
# ============================================================================

# -- Run bootstraps ----------------------------------------------------------

set.seed(1986)

boot_arthropodium <- bootstrap_effects(data_arthropodium, log_transform = TRUE)  %>% mutate(species = "Arthropodium strictum")
boot_bulbine      <- bootstrap_effects(data_bulbine,      log_transform = FALSE) %>% mutate(species = "Bulbine bulbosa")
boot_burchardia   <- bootstrap_effects(data_burchardia,   log_transform = FALSE) %>% mutate(species = "Burchardia umbellata")
boot_eryngium     <- bootstrap_effects(data_eryngium,     log_transform = TRUE)  %>% mutate(species = "Eryngium ovinum")

boot_all <- bind_rows(
  boot_arthropodium,
  boot_bulbine,
  boot_burchardia,
  boot_eryngium
) %>%
  mutate(
    treatment = factor(treatment, levels = c("Baseline", "Spring", "Summer"))
    )

# -- Summary stats ------------------------------------------------------------

# Mean + 95% empirical CI for Spring and Summer effects
boot_summary <- boot_all %>%
  filter(treatment != "Baseline") %>%
  group_by(species, treatment) %>%
  summarise(
    mean_effect = mean(effect),
    ci_low      = quantile(effect, 0.025),
    ci_high     = quantile(effect, 0.975),
    .groups     = "drop"
  ) %>%
  mutate(
    treatment = factor(treatment, levels = c("Baseline", "Spring", "Summer"))
    )

# Baseline ribbon: 95% CI of bootstrapped baseline mean, centred at 0.
# Represents estimation uncertainty around the zero reference — effects
# within this ribbon are within baseline sampling noise.
baseline_ribbon <- boot_all %>%
  filter(treatment == "Baseline") %>%
  group_by(species) %>%
  summarise(
    ribbon_low  = quantile(baseline_mean, 0.025) - mean(baseline_mean),
    ribbon_high = quantile(baseline_mean, 0.975) - mean(baseline_mean),
    .groups     = "drop"
  )

# Shared y-axis limits across all panel b plots
boot_limits <- extendrange(
  boot_all$effect[boot_all$treatment != "Baseline"],
  f = 0.01
)

# -- Build panel --------------------------------------------------------------

plot_boot_arthropodium <- make_boot_plot("Arthropodium strictum",  "***Arthropodium strictum***",  show_y = TRUE,  tag = "(**b**)")
plot_boot_bulbine      <- make_boot_plot("Bulbine bulbosa",       "***Bulbine bulbosa***",       show_y = FALSE)
plot_boot_burchardia   <- make_boot_plot("Burchardia umbellata",  "***Burchardia umbellata***",  show_y = FALSE)
plot_boot_eryngium     <- make_boot_plot("Eryngium ovinum",       "***Eryngium ovinum***",       show_y = FALSE)

panel_b <- wrap_plots(
  plot_boot_arthropodium,
  plot_boot_bulbine,
  plot_boot_burchardia,
  plot_boot_eryngium,
  nrow = 1
)

print(panel_b)

# =============================================================================
# Combine and save
# =============================================================================

figure_2 <- panel_a / panel_b

ggsave(
  "output/figure_2.png",
  figure_2,
  width  = 16,
  height = 9.5,
  units  = "cm",
  dpi    = 300
)

# Model coefficients table
lm_summaries <- bind_rows(
  summary_arthropodium %>% mutate(species = "Arthropodium strictum"),
  summary_bulbine      %>% mutate(species = "Bulbine bulbosa"),
  summary_burchardia   %>% mutate(species = "Burchardia umbellata"),
  summary_eryngium     %>% mutate(species = "Eryngium ovinum")
) %>%
  select(
    Model = species, Parameter, Coefficient, `95% CI`,
    df = df_error, `t-value` = t, `p-value` = p
  )

write.csv(lm_summaries, "output/lm_summaries.csv", row.names = FALSE)

# Effect size table: absolute difference and proportional change
# Proportional change = effect / baseline_mean, computed per bootstrap
# iteration so the CI correctly reflects its uncertainty.
effect_table <- boot_all %>%
  filter(treatment != "Baseline") %>%
  mutate(prop_change = effect / baseline_mean) %>%
  group_by(species, treatment) %>%
  summarise(
    effect_mean     = mean(effect),
    effect_ci_low   = quantile(effect, 0.025),
    effect_ci_high  = quantile(effect, 0.975),
    prop_mean       = mean(prop_change),
    prop_ci_low     = quantile(prop_change, 0.025),
    prop_ci_high    = quantile(prop_change, 0.975),
    .groups         = "drop"
  ) %>%
  mutate(
    `Effect size (g)`         = paste0(round(effect_mean, 2),
                                       " [", round(effect_ci_low,  2),
                                       ", ", round(effect_ci_high, 2), "]"),
    `Proportional change (%)`  = paste0(round(prop_mean * 100, 1),
                                        " [", round(prop_ci_low  * 100, 1),
                                        ", ", round(prop_ci_high * 100, 1), "]")
  ) %>%
  select(
    Species   = species,
    Treatment = treatment,
    `Effect size (g)`,
    `Proportional change (%)`
  ) %>%
  arrange(Species, Treatment)

write.csv(effect_table, "output/effect_size_table.csv", row.names = FALSE)
