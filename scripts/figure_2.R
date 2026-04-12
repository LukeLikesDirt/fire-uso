
# Install packages
#library(ggpubr)
library(ggtext)
library(parameters)
library(emmeans)
library(patchwork)
library(ggdist)
library(purrr)
library(tidyverse)

# Set theme
common_theme <- theme(
  panel.border = element_rect(colour = "grey", fill = NA, linewidth = 1),
  panel.background = element_blank(),
  axis.ticks.y = element_line(colour = "grey", linewidth = 0.5),
  axis.ticks.length.y = unit(-0.15, "cm"),
  axis.ticks.x = element_blank(),
  legend.position = "none",
  axis.text = element_markdown(size = 8),
  axis.title = element_markdown(size = 10),
  plot.title = element_markdown(size = 9.5, hjust = 0.5),
  plot.tag = element_markdown(size = 14),
  plot.margin = margin(1, 1, 1, 1),
  aspect.ratio = 1
  )

# y-axis limits
limits <- c(0, 8)

# Read data
data <- data.table::fread("data/uso_response.csv") %>%
  filter(
    # Remove NA treatments
    !is.na(treatment)
    ) %>%
  mutate(
    treatment = factor(treatment, levels = c("Baseline", "Spring", "Summer"))
  ) %>%
  glimpse()

# -- Arthropodium bulbosa -----------------------------------------------------------

# Create df for arthropodium bulbosa 
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
  labs(
    x = NULL,
    y = "Dry weight (g)",
    title = "***Arthropodium bulbosa***"
    #tag = "(**a**)"
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
ggsave(
  "output/figure_2.png",
  panel_a,
  width = 16,
  height = 6,
  units = "cm",
  dpi = 300
  )

# ============================================================================
# Bootstrap effect sizes (panel b)
# ============================================================================

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

# -- Run bootstraps (10,000 per species) -------------------------------------

set.seed(1986)

boot_arthropodium <- bootstrap_effects(data_arthropodium, log_transform = TRUE)  %>% mutate(species = "Arthropodium bulbosa")
boot_bulbine      <- bootstrap_effects(data_bulbine,      log_transform = FALSE) %>% mutate(species = "Bulbine bulbosa")
boot_burchardia   <- bootstrap_effects(data_burchardia,   log_transform = FALSE) %>% mutate(species = "Burchardia umbellata")
boot_eryngium     <- bootstrap_effects(data_eryngium,     log_transform = TRUE)  %>% mutate(species = "Eryngium ovinum")

boot_all <- bind_rows(
  boot_arthropodium,
  boot_bulbine,
  boot_burchardia,
  boot_eryngium
) %>%
  mutate(treatment = factor(treatment, levels = c("Baseline", "Spring", "Summer")))

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
  mutate(treatment = factor(treatment, levels = c("Baseline", "Spring", "Summer")))

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

# -- Plot function ------------------------------------------------------------

make_boot_plot <- function(sp_name, sp_label, show_y = TRUE) {

  boot_sp   <- boot_all        %>% filter(species == sp_name, treatment != "Baseline")
  summ_sp   <- boot_summary    %>% filter(species == sp_name)
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
      aes(x = treatment, y = effect, fill = treatment),
      alpha     = 0.5,
      scale     = 0.7,
      normalize = "groups",
      adjust    = 3,
      colour    = NA
    ) +
    # 95% empirical CI bar
    geom_errorbar(
      data = summ_sp,
      aes(x = treatment, ymin = ci_low, ymax = ci_high),
      width     = 0.08,
      linewidth = 0.8
    ) +
    # Mean bootstrap estimate
    geom_point(
      data = summ_sp,
      aes(x = treatment, y = mean_effect, fill = treatment),
      shape = 21, size = 4, stroke = 1
    ) +
    # Baseline reference point at zero
    geom_point(
      data = data.frame(treatment = factor("Baseline", levels = c("Baseline", "Spring", "Summer")), y = 0),
      aes(x = treatment, y = y),
      shape = 21, size = 4, stroke = 1, fill = "grey60"
    ) +
    scale_x_discrete(limits = c("Baseline", "Spring", "Summer")) +
    labs(
      x     = NULL,
      y     = if (show_y) "Effect size (g)" else NULL,
      title = sp_label
    ) +
    common_theme

  if (!show_y) {
    p <- p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

  p
}

# -- Build panels -------------------------------------------------------------

plot_boot_arthropodium <- make_boot_plot("Arthropodium bulbosa",  "***Arthropodium bulbosa***",  show_y = TRUE)
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

ggsave(
  "output/figure_2b.png",
  panel_b,
  width  = 16,
  height = 6,
  units  = "cm",
  dpi    = 300
)
