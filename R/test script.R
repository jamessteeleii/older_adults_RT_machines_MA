library(tidyverse)
library(patchwork)
library(metafor)
library(brms)
library(rstan)
library(tidybayes)
library(marginaleffects)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

##### Read csv as data frame into environment 
data <- read_csv("data/data.csv")

# Calculate pre-post SDs from SEs
data$sd_pre <- ifelse(is.na(data$se_pre), data$sd_pre, data$se_pre * sqrt(data$n_post))
data$sd_post <- ifelse(is.na(data$se_post), data$sd_post, data$se_post * sqrt(data$n_post))

# Convert p to t (Change scores)
data$change_t_value <- replmiss(data$change_t_value, with(data, qt(change_p_value/2, df=n_post-1, lower.tail=FALSE)))

# Convert t to SE (Change scores)
data$se_change <- replmiss(data$se_change, with(data, ifelse(is.na(m_change), 
                                                                 (m_post - m_pre)/change_t_value, m_change/change_t_value)))

# Make SE positive (Change scores)
data$se_change <- ifelse(data$se_change < 0, data$se_change * -1, data$se_change)

# Convert CI to SE (Change scores)
data$se_change <- replmiss(data$se_change, with(data, (change_CI_upper - change_CI_lower)/3.92))

# Convert SE to SD (Change scores)
data$sd_change <- replmiss(data$sd_change, with(data, se_change * sqrt(n_post)))

# Calculate pre-post correlation coefficient for those with pre, post, and delta SDs
data$ri <- (data$sd_pre^2 + data$sd_post^2 - data$sd_change^2)/(2 * data$sd_pre * data$sd_post)

# Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
data$ri <- ifelse(between(data$ri,-1,1) == FALSE, NA, data$ri)

# Then we'll convert using Fishers r to z, calculate a meta-analytic point estimate, and impute that across the studies with missing correlations
data <- escalc(measure = "ZCOR", ri = ri, ni = n_post, data = data)

Meta_training_ri <- rma.mv(yi, V=vi, data=filter(data, training_control == "training"),
                     slab=paste(study),
                     random = list(~ 1 | study_code, ~ 1 | group_code, ~ 1 | es_code), method="REML", test="t",
                     control=list(optimizer="optim", optmethod="Nelder-Mead"))

RobuEst_Meta_training_ri <- robust(Meta_training_ri, filter(data, training_control == "training")$study_code)

z2r_training <- psych::fisherz2r(RobuEst_Meta_training_ri$b[1])

Meta_control_ri <- rma.mv(yi, V=vi, data=filter(data, training_control == "control"),
                          slab=paste(study),
                          random = list(~ 1 | study_code, ~ 1 | group_code, ~ 1 | es_code), method="REML", test="t",
                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

RobuEst_Meta_control_ri <- robust(Meta_control_ri, filter(data, training_control == "control")$study)

z2r_control <- psych::fisherz2r(RobuEst_Meta_control_ri$b[1])

data <- data |>
  mutate(ri = if_else(is.na(data$ri) & training_control == "training", z2r_training, 
                      if_else(is.na(data$ri) & training_control == "control", z2r_control, ri))
  ) |>
    # remove yi and vi
    select(-yi,-vi)

# Estimate change score difference SD where only pre-post data available
data$sd_change <- replmiss(data$sd_change, with(data, sqrt(sd_pre^2 + sd_post^2 - (2*ri*sd_pre*sd_post))))


### Standardised mean difference effect size calculations
data_increase <- escalc(measure="SMCR", m1i=m_post, 
                  m2i=m_pre, sd1i=sd_pre, ni=n_post, ri=ri, data = filter(data, increase_decrease == "increase"))

data_decrease <- escalc(measure="SMCR", m1i=m_pre, 
                        m2i=m_post, sd1i=sd_pre, ni=n_post, ri=ri, data = filter(data, increase_decrease == "decrease"))

data_increase_change <- escalc(measure="SMCR", m1i=m_change, 
                        m2i=m_pre, sd1i=sd_pre, ni=n_post, ri=ri, data = filter(data, increase_decrease == "increase" & is.na(m_post)) |>
                          mutate(m_pre = 0)
                        )

data_decrease_change <- escalc(measure="SMCR", m1i=m_pre, 
                        m2i=m_change, sd1i=sd_pre, ni=n_post, ri=ri, data = filter(data, increase_decrease == "decrease" & is.na(m_post)) |>
                          mutate(m_pre = 0)
                        )

data <- rbind(data_increase, data_decrease, data_increase_change, data_decrease_change)

# Functional outcomes
data_function <- filter(data, outcome_group == "function") |>
  filter(!is.na(yi) &
           !is.na(vi))

# meta_function <- rma.mv(yi, V=vi, data=data_function,
#                                        slab=paste(study), mods = ~ 0 + training_control,
#                                        random = list(~ 1 | study_code, ~ 1 | group_code, ~ 1 | es_code), method="REML", test="t",
#                                        control=list(optimizer="optim", optmethod="Nelder-Mead"))
# 
# 
# ### Calculate robust estimate from multi-level model
# RobuEst_meta_function <- robust(meta_function, data_function$study)
# 
# data_function <- cbind(data_function, predict(RobuEst_meta_function))

# data_function |>
#   rowid_to_column() |>
#   ggplot(aes(y=reorder(study, yi))) +
#   geom_vline(xintercept = 0, lty = "dashed") +
#   annotate("rect", ymin = -Inf, ymax = Inf,
#                   xmin = RobuEst_meta_function$ci.lb[1],
#                   xmax = RobuEst_meta_function$ci.ub[1],
#            alpha = 0.25, color = NA, fill = "#D55E00") +
#   geom_vline(xintercept = RobuEst_meta_function$beta[1], color = "#D55E00") +
#   annotate("rect", ymin = -Inf, ymax = Inf,
#            xmin = RobuEst_meta_function$ci.lb[2],
#            xmax = RobuEst_meta_function$ci.ub[2],
#            alpha = 0.25, color = NA, fill  = "#009E73") +
#   geom_vline(xintercept = RobuEst_meta_function$beta[2], color = "#009E73") +
#   geom_pointrange(aes(x=yi, color = training_control,
#                       xmin = yi - (sqrt(vi) * 1.96),
#                       xmax = yi + (sqrt(vi) * 1.96)
#   ), 
#   position = position_jitter(h=0.15, w=0, seed = 1988), alpha = 0.5) +
#   scale_color_manual(values = c("#D55E00","#009E73")) +
#   annotate("text", x = 2.25, y = 3,
#            label = glue::glue("Training Estimate = {round(RobuEst_meta_function$beta[2],2)} [95%CI: {round(RobuEst_meta_function$ci.lb[2],2)} to {round(RobuEst_meta_function$ci.ub[2],2)}]"),
#            size = 3) +
#   annotate("text", x = 2.25, y = 2,
#            label = glue::glue("Control Estimate = {round(RobuEst_meta_function$beta[1],2)} [95%CI: {round(RobuEst_meta_function$ci.lb[1],2)} to {round(RobuEst_meta_function$ci.ub[1],2)}]"),
#            size = 3) +
#   labs(
#     x = "Standardised Mean Change",
#     y = "Study",
#     color = "Group",
#     title = "Functional Outcomes",
#     subtitle = "Vertical lines and bands are the meta-analytic point and 95% confidence interval estimates\nIndividual study level effects shown with point and confidence interval estimates"
#   ) +
#   theme_bw()
# 
# ggsave("function_plot.png", height = 5, width = 10, dpi = 300)

meta_function <- brm(yi | se(sqrt(vi)) ~ training_control + 
                                  (training_control | study_code) + 
                                  (1 | group_code) +
                                  (1 | es_code),
                                data = data_function,
                                chains = 4,
                                cores = 4,
                                seed = 1988,
                                warmup = 2000,
                                iter = 6000,
                                control = list(adapt_delta = 0.99, max_treedepth = 11)
)

plot(meta_function)

pp_check(meta_function)


meta_pred <- predictions(
  meta_function,
  re_formula = NA
) |>
  posterior_draws()

meta_pred_labels <- meta_pred |>
  group_by(training_control) |>
  mean_qi(draw)

study_pred <- predictions(
  meta_function,
  re_formula = NULL,
) |>
  posterior_draws()

study_labels <- data_function |>
  select(1,2) |>
  group_by(study_code) |>
  slice(1)

study_pred <- left_join(study_pred, study_labels, by = "study_code") |>
  group_by(study) |>
  mutate(mean_draw = mean(draw))

study_pred_labels <- study_pred |>
  group_by(study, training_control) |>
  mean_qi(draw)

meta_slopes <- avg_slopes(
  meta_function,
  re_formula = NA,
  variables = "training_control"
) |>
  posterior_draws()

meta_slopes_labels <- meta_slopes |>
  group_by(contrast) |>
  mean_qi(draw)


meta_pred_plot <- ggplot(meta_pred, aes(x = draw, fill = training_control)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  stat_halfeye(slab_alpha = .5, .width = 0.95) +
  facet_grid(training_control~.) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(meta_pred_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      x = draw, y = 0.1
    ),
    size = 3
  ) +
  labs(
    x = "Standardised Mean Change",
    fill = "Condition",
    title = "Global Grand Mean Estimates for Condition"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))



study_pred_plot <- ggplot(study_pred, aes(x = draw, 
                                          y = reorder(study, mean_draw), 
                                          fill = training_control)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  # Add individual study data
  geom_point(
    data = data_function,
    aes(x = yi, y = study, color = training_control),
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1),
    alpha = 0.5
  ) +
  stat_halfeye(slab_alpha = .5, position = position_dodge(width = 0.5), size = 0.25) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(study_pred_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{training_control}: {draw} [{.lower}, {.upper}]"),
      x = 3, y =reorder(study, draw), group = training_control
    ),
    size = 2, position = position_dodge(width = 0.75),
    hjust = "inward"
  ) +
  labs(
    x = "Standardised Mean Change",
    title = "Conditional Estimates for Condition by Study"
  ) +
  guides(
    fill = "none",
    color = "none"
  ) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(title = element_text(size=8))


contrast_plot <- ggplot(meta_slopes, aes(x = draw)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  stat_halfeye(slab_alpha = .5, fill = "black") +
  geom_text(
    data = mutate_if(meta_slopes_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      x = draw, y = 0.1
    ),
    size = 3
  ) +
  labs(
    x = "Standardised Mean Change",
    title = "Contrasts Between Conditions (Training - Control)"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))

meta_plots <- (study_pred_plot / (meta_pred_plot + contrast_plot)) + 
  plot_annotation(title = "Meta-Analysis of Prior Studies Examining the Effects of Machine Based Resistance Training on Function",
                  caption = "Point estimates and 95% quantile intervals reported") +
  plot_layout(guides = "collect", axis_titles = "collect",
              widths = c(2,1,1))  &
  theme(axis.title.x = element_text(size=10),
        legend.position = "bottom")

meta_plots

ggsave("meta_plots_function.tiff", dpi = 300, w=10, h=10)





# Strength outcomes
data_strength <- filter(data, outcome_group == "upper_body_strength" |
                                     outcome_group == "lower_body_strength") |>
  filter(!is.na(yi) &
           !is.na(vi))

# MultiLevelModel_strength <- rma.mv(yi, V=vi, data=data_strength,
#                                    slab=paste(study), mods = ~ 0 + training_control,
#                                    random = list(~ 1 | study_code, ~ 1 | group_code, ~ 1 | es_code), method="REML", test="t",
#                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))
# 
# 
# ### Calculate robust estimate from multi-level model
# RobuEst_MultiLevelModel_strength <- robust(MultiLevelModel_strength, data_strength$study)
# 
# data_strength <- cbind(data_strength, predict(RobuEst_MultiLevelModel_strength))
# 
# data_strength |>
#   rowid_to_column() |>
#   ggplot(aes(y=reorder(study, yi))) +
#   geom_vline(xintercept = 0, lty = "dashed") +
#   annotate("rect", ymin = -Inf, ymax = Inf,
#            xmin = RobuEst_MultiLevelModel_strength$ci.lb[1],
#            xmax = RobuEst_MultiLevelModel_strength$ci.ub[1],
#            alpha = 0.25, color = NA, fill = "#D55E00") +
#   geom_vline(xintercept = RobuEst_MultiLevelModel_strength$beta[1], color = "#D55E00") +
#   annotate("rect", ymin = -Inf, ymax = Inf,
#            xmin = RobuEst_MultiLevelModel_strength$ci.lb[2],
#            xmax = RobuEst_MultiLevelModel_strength$ci.ub[2],
#            alpha = 0.25, color = NA, fill  = "#009E73") +
#   geom_vline(xintercept = RobuEst_MultiLevelModel_strength$beta[2], color = "#009E73") +
#   geom_pointrange(aes(x=yi, color = training_control,
#                       xmin = yi - (sqrt(vi) * 1.96),
#                       xmax = yi + (sqrt(vi) * 1.96)
#   ), 
#   position = position_jitter(h=0.15, w=0, seed = 1988), alpha = 0.5) +
#   scale_color_manual(values = c("#D55E00","#009E73")) +
#   annotate("text", x = 2.25, y = 3,
#            label = glue::glue("Training Estimate = {round(RobuEst_MultiLevelModel_strength$beta[2],2)} [95%CI: {round(RobuEst_MultiLevelModel_strength$ci.lb[2],2)} to {round(RobuEst_MultiLevelModel_strength$ci.ub[2],2)}]"),
#            size = 3) +
#   annotate("text", x = 2.25, y = 2,
#            label = glue::glue("Control Estimate = {round(RobuEst_MultiLevelModel_strength$beta[1],2)} [95%CI: {round(RobuEst_MultiLevelModel_strength$ci.lb[1],2)} to {round(RobuEst_MultiLevelModel_strength$ci.ub[1],2)}]"),
#            size = 3) +
#   labs(
#     x = "Standardised Mean Change",
#     y = "Study",
#     color = "Group",
#     title = "Strength Outcomes",
#     subtitle = "Vertical lines and bands are the meta-analytic point and 95% confidence interval estimates\nIndividual study level effects shown with point and confidence interval estimates"
#   ) +
#   theme_bw()
# 
# 
# ggsave("strength_plot.png", height = 5, width = 10, dpi = 300)

meta_strength <- brm(yi | se(sqrt(vi)) ~ training_control + 
                       (training_control | study_code) + 
                       (1 | group_code) +
                       (1 | es_code),
                     data = data_strength,
                     chains = 4,
                     cores = 4,
                     seed = 1988,
                     warmup = 2000,
                     iter = 6000,
                     control = list(adapt_delta = 0.99, max_treedepth = 11)
)

plot(meta_strength)

pp_check(meta_strength)


meta_pred <- predictions(
  meta_strength,
  re_formula = NA
) |>
  posterior_draws()

meta_pred_labels <- meta_pred |>
  group_by(training_control) |>
  mean_qi(draw)

study_pred <- predictions(
  meta_strength,
  re_formula = NULL,
) |>
  posterior_draws()

study_labels <- data_strength |>
  select(1,2) |>
  group_by(study_code) |>
  slice(1)

study_pred <- left_join(study_pred, study_labels, by = "study_code") |>
  group_by(study) |>
  mutate(mean_draw = mean(draw))

study_pred_labels <- study_pred |>
  group_by(study, training_control) |>
  mean_qi(draw)

meta_slopes <- avg_slopes(
  meta_strength,
  re_formula = NA,
  variables = "training_control"
) |>
  posterior_draws()

meta_slopes_labels <- meta_slopes |>
  group_by(contrast) |>
  mean_qi(draw)


meta_pred_plot <- ggplot(meta_pred, aes(x = draw, fill = training_control)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  stat_halfeye(slab_alpha = .5, .width = 0.95) +
  facet_grid(training_control~.) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(meta_pred_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      x = draw, y = 0.1
    ),
    size = 3
  ) +
  labs(
    x = "Standardised Mean Change",
    fill = "Condition",
    title = "Global Grand Mean Estimates for Condition"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))



study_pred_plot <- ggplot(study_pred, aes(x = draw, 
                                          y = reorder(study, mean_draw), 
                                          fill = training_control)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  # Add individual study data
  geom_point(
    data = data_strength,
    aes(x = yi, y = study, color = training_control),
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1),
    alpha = 0.5
  ) +
  stat_halfeye(slab_alpha = .5, position = position_dodge(width = 0.5), size = 0.25) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(study_pred_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{training_control}: {draw} [{.lower}, {.upper}]"),
      x = 3, y =reorder(study, draw), group = training_control
    ),
    size = 2, position = position_dodge(width = 0.75),
    hjust = "inward"
  ) +
  labs(
    x = "Standardised Mean Change",
    title = "Conditional Estimates for Condition by Study"
  ) +
  guides(
    fill = "none",
    color = "none"
  ) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(title = element_text(size=8))


contrast_plot <- ggplot(meta_slopes, aes(x = draw)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  stat_halfeye(slab_alpha = .5, fill = "black") +
  geom_text(
    data = mutate_if(meta_slopes_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{draw} [{.lower}, {.upper}]"),
      x = draw, y = 0.1
    ),
    size = 3
  ) +
  labs(
    x = "Standardised Mean Change",
    title = "Contrasts Between Conditions (Training - Control)"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))

meta_plots <- (study_pred_plot / (meta_pred_plot + contrast_plot)) + 
  plot_annotation(title = "Meta-Analysis of Prior Studies Examining the Effects of Machine Based Resistance Training on Strength",
                  caption = "Point estimates and 95% quantile intervals reported") +
  plot_layout(guides = "collect", axis_titles = "collect",
              widths = c(2,1,1))  &
  theme(axis.title.x = element_text(size=10),
        legend.position = "bottom")

meta_plots

ggsave("meta_plots_strength.tiff", dpi = 300, w=10, h=10)


