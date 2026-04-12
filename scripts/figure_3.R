# Packages
library(tidyverse)
library(ggtext)
library(patchwork)
library(mgcv)

# Common theme (same as figure_2.R)
common_theme <- theme(
  panel.border = element_rect(colour = "grey", fill = NA, linewidth = 0.5),
  panel.background = element_blank(),
  axis.ticks = element_line(colour = "grey", linewidth = 0.25),
  axis.ticks.length = unit(-0.1, "cm"),
  legend.position = "none",
  axis.text = element_markdown(size = 7),
  plot.tag = element_markdown(size = 12),
  plot.tag.location = "plot",
  plot.margin = margin(1, 1, 1, 1),
  strip.background = element_blank(),
  strip.text = element_markdown(size = 8),
  aspect.ratio = 1
)

# Data
data <- read.csv("data/phenology.csv") %>%
  # Remove "Microseris.walterii"
  filter(species != "Microseris.walterii") %>%
  select(treatment:X11.7.2024) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "date",
    values_to = "observation",
    names_prefix = "X"
  ) %>%
  mutate(
    state = case_when(
      observation == 0 ~ NA_character_,
      observation == 6 ~ "Dormant",      # ordinal value = 1
      observation == 7 ~ "Resprouting",  # ordinal value = 2
      observation == 1 ~ "Vegetative",   # ordinal value = 3
      observation == 2 ~ "Budding",      # ordinal value = 4
      observation == 3 ~ "Flowering",    # ordinal value = 5
      observation == 4 ~ "Seed Head",    # ordinal value = 6
      observation == 5 ~ "Seed Shed",    # ordinal value = 7
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(state)) %>%
  mutate(date = as.Date(date, format = "%d.%m.%Y")) %>%
  mutate(state = factor(
    state,
    levels = c("Dormant", "Resprouting", "Vegetative", "Budding", "Flowering", "Seed Head", "Seed Shed")
  )) %>%
  # Replace dots in species names with spaces
  mutate(species = paste0("***", str_replace_all(species, "\\.", " "), "***")) %>%
  # Shift the first observation date a week earlier for the spring treatment group
  group_by(species, treatment, code) %>%
  mutate(date = if_else(treatment == "spring" & date == min(date), date - 7, date)) %>%
  ungroup() %>%
  glimpse()

ordinal_levels <- levels(data$state)
n_levels      <- length(ordinal_levels)

# -- Ordinal GAM (ocat) predictions -------------------------------------------
# mgcv's ocat family models ordered categorical responses directly via a
# proportional-odds / cumulative-link formulation. Predictions are converted
# to an expected category value: E[Y] = sum_k k * P(Y = k), so the smooth
# curve lives on the same 1–7 numeric scale used by the y-axis.
#
# 95% CIs are obtained by posterior simulation of the joint (beta, theta)
# parameter vector, which propagates both smooth and cut-point uncertainty.

ocat_predictions <- function(df, n_levels) {
  df_fit <- df %>%
    mutate(
      date_num  = as.numeric(date),
      state_int = as.integer(state)
    )

  m <- gam(
    state_int ~ s(date_num, bs = "cr"),
    family = ocat(R = n_levels),
    data   = df_fit,
    method = "REML"
  )

  pred_dates <- seq(min(df_fit$date), max(df_fit$date), by = "day")
  nd <- data.frame(date_num = as.numeric(pred_dates))

  # Expected category from posterior mean
  probs   <- predict(m, newdata = nd, type = "response")   # n x R matrix
  ev_mean <- as.numeric(probs %*% seq_len(n_levels))

  # 95% CI via posterior simulation (unconditional = accounts for smoothing
  # parameter uncertainty in addition to coefficient uncertainty)
  set.seed(123)
  B    <- 1000
  sims <- rmvn(B, coef(m), vcov(m, unconditional = TRUE))

  ev_sims <- vapply(seq_len(B), function(b) {
    m_b              <- m
    m_b$coefficients <- sims[b, ]
    p_b              <- predict(m_b, newdata = nd, type = "response")
    as.numeric(p_b %*% seq_len(n_levels))
  }, numeric(nrow(nd)))

  data.frame(
    date   = pred_dates,
    ev     = ev_mean,
    ev_lwr = apply(ev_sims, 1, quantile, 0.025),
    ev_upr = apply(ev_sims, 1, quantile, 0.975)
  )
}

# Fit per species x treatment (group_modify passes the subgroup data frame)
gam_preds <- data %>%
  group_by(treatment, species) %>%
  group_modify(~ ocat_predictions(.x, n_levels)) %>%
  ungroup()

# -- Model summary tables (one CSV per species) -------------------------------
ocat_summary <- function(df, n_levels) {
  df_fit <- df %>%
    mutate(
      date_num  = as.numeric(date),
      state_int = as.integer(state)
    )

  m <- gam(
    state_int ~ s(date_num, bs = "cr"),
    family = ocat(R = n_levels),
    data   = df_fit,
    method = "REML"
  )

  s <- summary(m)
  data.frame(
    term    = rownames(s$s.table),
    edf     = round(s$s.table[, "edf"],     2),
    ref_df  = round(s$s.table[, "Ref.df"],  2),
    chi_sq  = round(s$s.table[, "Chi.sq"],  3),
    p_value = round(s$s.table[, "p-value"], 4)
  )
}

gam_summaries <- data %>%
  group_by(treatment, species) %>%
  group_modify(~ ocat_summary(.x, n_levels)) %>%
  ungroup() %>%
  mutate(
    treatment = str_to_title(treatment),
    species   = str_remove_all(species, "\\*")   # strip markdown italics
  )

# Single CSV with formatted headers
gam_summaries %>%
  mutate(`Model/Species` = paste0(species, " (", treatment, ")")) %>%
  select(
    `Model/Species`,
    `Term`    = term,
    `EDF`     = edf,
    `Ref. df` = ref_df,
    `Chi-sq`  = chi_sq,
    `p-value` = p_value
  ) %>%
  write.csv("output/gam_summaries.csv", row.names = FALSE)

# -- Spring panel (a) ---------------------------------------------------------

plot_spring <- ggplot(
  data %>% filter(treatment == "spring"),
  aes(x = date, y = as.numeric(state), group = code)
) +
  geom_line(
    linewidth = 0.2, alpha = 0.5, colour = "darkgray"
    #position = position_jitter(width = 0, height = 0.05)
    ) +
  geom_ribbon(
    data = gam_preds %>% filter(treatment == "spring"),
    aes(x = date, ymin = ev_lwr, ymax = ev_upr),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = gam_preds %>% filter(treatment == "spring"),
    aes(x = date, y = ev),
    inherit.aes = FALSE,
    linewidth = 0.8, colour = "black"
  ) +
  geom_vline(xintercept = as.Date("2023-11-08"),
             linetype = "dotted", colour = "red") +
  scale_y_continuous(breaks = seq_len(n_levels), labels = ordinal_levels) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(x = NULL, y = NULL, tag = "(**a**)") +
  facet_wrap(~species, nrow = 1) +
  common_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    )

# -- Summer panel (b) ---------------------------------------------------------

plot_summer <- ggplot(
  data %>% filter(treatment == "summer"),
  aes(x = date, y = as.numeric(state), group = code)
) +
  geom_line(
    linewidth = 0.2, alpha = 0.5, colour = "darkgrey"
    #position = position_jitter(width = 0, height = 0.05)
    ) +
  geom_ribbon(
    data = gam_preds %>% filter(treatment == "summer"),
    aes(x = date, ymin = ev_lwr, ymax = ev_upr),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = gam_preds %>% filter(treatment == "summer"),
    aes(x = date, y = ev),
    inherit.aes = FALSE,
    linewidth = 0.8, colour = "black"
  ) +
  geom_vline(xintercept = as.Date("2024-01-16"),
             linetype = "dotted", colour = "red") +
  scale_y_continuous(breaks = seq_len(n_levels), labels = ordinal_levels) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(x = NULL, y = NULL, tag = "(**b**)") +
  facet_wrap(~species, nrow = 1) +
  common_theme +
  theme(
    strip.text = element_markdown(colour = "white"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)
    )

# -- Combine and save ---------------------------------------------------------

figure_3 <- plot_spring / plot_summer
print(figure_3)

ggsave(
  "output/figure_3.png",
  figure_3,
  width = 16,
  height = 9,
  units = "cm",
  dpi = 600
)
