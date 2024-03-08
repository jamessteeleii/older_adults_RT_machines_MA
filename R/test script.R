library(tidyverse)
library(metafor)


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

# Make positive
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

data <- rbind(data_increase, data_decrease)

# Functional outcomes
data_function <- filter(data, outcome_group == "function") |>
  filter(!is.na(yi) &
           !is.na(vi))

MultiLevelModel_function <- rma.mv(yi, V=vi, data=data_function,
                                       slab=paste(study), mods = ~ 0 + training_control,
                                       random = list(~ 1 | study_code, ~ 1 | group_code, ~ 1 | es_code), method="REML", test="t",
                                       control=list(optimizer="optim", optmethod="Nelder-Mead"))


### Calculate robust estimate from multi-level model
RobuEst_MultiLevelModel_function <- robust(MultiLevelModel_function, data_function$study)

data_function <- cbind(data_function, predict(RobuEst_MultiLevelModel_function))

RobuEst_MultiLevelModel_function$ci.lb[1]

data_function |>
  rowid_to_column() |>
  ggplot(aes(y=reorder(study, yi))) +
  annotate("rect", ymin = -Inf, ymax = Inf,
                  xmin = RobuEst_MultiLevelModel_function$ci.lb[1],
                  xmax = RobuEst_MultiLevelModel_function$ci.ub[1],
           alpha = 0.25, color = NA, fill = "#D55E00") +
  geom_vline(xintercept = RobuEst_MultiLevelModel_function$beta[1], color = "#D55E00") +
  annotate("rect", ymin = -Inf, ymax = Inf,
           xmin = RobuEst_MultiLevelModel_function$ci.lb[2],
           xmax = RobuEst_MultiLevelModel_function$ci.ub[2],
           alpha = 0.25, color = NA, fill  = "#009E73") +
  geom_vline(xintercept = RobuEst_MultiLevelModel_function$beta[2], color = "#009E73") +
  geom_pointrange(aes(x=yi, color = training_control,
                      xmin = yi - (sqrt(vi) * 1.96),
                      xmax = yi + (sqrt(vi) * 1.96)
  ), 
  position = position_jitter(h=0.15), alpha = 0.5) +
  scale_color_manual(values = c("#D55E00","#009E73")) +
  annotate("text", x = 2.25, y = 3,
           label = glue::glue("Training Estimate = {round(RobuEst_MultiLevelModel_function$beta[2],2)} [95%CI: {round(RobuEst_MultiLevelModel_function$ci.lb[2],2)} to {round(RobuEst_MultiLevelModel_function$ci.ub[2],2)}]"),
           size = 3) +
  annotate("text", x = 2.25, y = 2,
           label = glue::glue("Control Estimate = {round(RobuEst_MultiLevelModel_function$beta[1],2)} [95%CI: {round(RobuEst_MultiLevelModel_function$ci.lb[1],2)} to {round(RobuEst_MultiLevelModel_function$ci.ub[1],2)}]"),
           size = 3) +
  labs(
    x = "Standardised Mean Change",
    y = "Study",
    color = "Group",
    title = "Functional Outcomes",
    subtitle = "Vertical lines and bands are the meta-analytic point and 95% confidence interval estimates\nIndividual study level effects shown with point and confidence interval estimates"
  ) +
  theme_bw()

ggsave("function_plot.png", height = 5, width = 10, dpi = 300)


# Upper Body Strength outcomes
data_upper_body_strength <- filter(data, outcome_group == "upper_body_strength") |>
  filter(!is.na(yi) &
           !is.na(vi))

MultiLevelModel_upper_body_strength <- rma.mv(yi, V=vi, data=data_upper_body_strength,
                                   slab=paste(study), mods = ~ 0 + training_control,
                                   random = list(~ 1 | study_code, ~ 1 | group_code, ~ 1 | es_code), method="REML", test="t",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))


### Calculate robust estimate from multi-level model
RobuEst_MultiLevelModel_upper_body_strength <- robust(MultiLevelModel_upper_body_strength, data_upper_body_strength$study)

data_upper_body_strength <- cbind(data_upper_body_strength, predict(RobuEst_MultiLevelModel_upper_body_strength))

RobuEst_MultiLevelModel_upper_body_strength$ci.lb[1]

data_upper_body_strength |>
  rowid_to_column() |>
  ggplot(aes(y=reorder(study, yi))) +
  annotate("rect", ymin = -Inf, ymax = Inf,
           xmin = RobuEst_MultiLevelModel_upper_body_strength$ci.lb[1],
           xmax = RobuEst_MultiLevelModel_upper_body_strength$ci.ub[1],
           alpha = 0.25, color = NA, fill = "#D55E00") +
  geom_vline(xintercept = RobuEst_MultiLevelModel_upper_body_strength$beta[1], color = "#D55E00") +
  annotate("rect", ymin = -Inf, ymax = Inf,
           xmin = RobuEst_MultiLevelModel_upper_body_strength$ci.lb[2],
           xmax = RobuEst_MultiLevelModel_upper_body_strength$ci.ub[2],
           alpha = 0.25, color = NA, fill  = "#009E73") +
  geom_vline(xintercept = RobuEst_MultiLevelModel_upper_body_strength$beta[2], color = "#009E73") +
  geom_pointrange(aes(x=yi, color = training_control,
                      xmin = yi - (sqrt(vi) * 1.96),
                      xmax = yi + (sqrt(vi) * 1.96)
  ), 
  position = position_jitter(h=0.15), alpha = 0.5) +
  scale_color_manual(values = c("#D55E00","#009E73")) +
  annotate("text", x = 2.25, y = 3,
           label = glue::glue("Training Estimate = {round(RobuEst_MultiLevelModel_upper_body_strength$beta[2],2)} [95%CI: {round(RobuEst_MultiLevelModel_upper_body_strength$ci.lb[2],2)} to {round(RobuEst_MultiLevelModel_upper_body_strength$ci.ub[2],2)}]"),
           size = 3) +
  annotate("text", x = 2.25, y = 2,
           label = glue::glue("Control Estimate = {round(RobuEst_MultiLevelModel_upper_body_strength$beta[1],2)} [95%CI: {round(RobuEst_MultiLevelModel_upper_body_strength$ci.lb[1],2)} to {round(RobuEst_MultiLevelModel_upper_body_strength$ci.ub[1],2)}]"),
           size = 3) +
  labs(
    x = "Standardised Mean Change",
    y = "Study",
    color = "Group",
    title = "Upper Body Strength Outcomes",
    subtitle = "Vertical lines and bands are the meta-analytic point and 95% confidence interval estimates\nIndividual study level effects shown with point and confidence interval estimates"
  ) +
  theme_bw()


# Lower Body Strength outcomes
data_lower_body_strength <- filter(data, outcome_group == "lower_body_strength") |>
  filter(!is.na(yi) &
           !is.na(vi))

MultiLevelModel_lower_body_strength <- rma.mv(yi, V=vi, data=data_lower_body_strength,
                                              slab=paste(study), mods = ~ 0 + training_control,
                                              random = list(~ 1 | study_code, ~ 1 | group_code, ~ 1 | es_code), method="REML", test="t",
                                              control=list(optimizer="optim", optmethod="Nelder-Mead"))


### Calculate robust estimate from multi-level model
RobuEst_MultiLevelModel_lower_body_strength <- robust(MultiLevelModel_lower_body_strength, data_lower_body_strength$study)

data_lower_body_strength <- cbind(data_lower_body_strength, predict(RobuEst_MultiLevelModel_lower_body_strength))

RobuEst_MultiLevelModel_lower_body_strength$ci.lb[1]

data_lower_body_strength |>
  rowid_to_column() |>
  ggplot(aes(y=reorder(study, yi))) +
  annotate("rect", ymin = -Inf, ymax = Inf,
           xmin = RobuEst_MultiLevelModel_lower_body_strength$ci.lb[1],
           xmax = RobuEst_MultiLevelModel_lower_body_strength$ci.ub[1],
           alpha = 0.25, color = NA, fill = "#D55E00") +
  geom_vline(xintercept = RobuEst_MultiLevelModel_lower_body_strength$beta[1], color = "#D55E00") +
  annotate("rect", ymin = -Inf, ymax = Inf,
           xmin = RobuEst_MultiLevelModel_lower_body_strength$ci.lb[2],
           xmax = RobuEst_MultiLevelModel_lower_body_strength$ci.ub[2],
           alpha = 0.25, color = NA, fill  = "#009E73") +
  geom_vline(xintercept = RobuEst_MultiLevelModel_lower_body_strength$beta[2], color = "#009E73") +
  geom_pointrange(aes(x=yi, color = training_control,
                      xmin = yi - (sqrt(vi) * 1.96),
                      xmax = yi + (sqrt(vi) * 1.96)
  ), 
  position = position_jitter(h=0.15), alpha = 0.5) +
  scale_color_manual(values = c("#D55E00","#009E73")) +
  annotate("text", x = 2.25, y = 3,
           label = glue::glue("Training Estimate = {round(RobuEst_MultiLevelModel_lower_body_strength$beta[2],2)} [95%CI: {round(RobuEst_MultiLevelModel_lower_body_strength$ci.lb[2],2)} to {round(RobuEst_MultiLevelModel_lower_body_strength$ci.ub[2],2)}]"),
           size = 3) +
  annotate("text", x = 2.25, y = 2,
           label = glue::glue("Control Estimate = {round(RobuEst_MultiLevelModel_lower_body_strength$beta[1],2)} [95%CI: {round(RobuEst_MultiLevelModel_lower_body_strength$ci.lb[1],2)} to {round(RobuEst_MultiLevelModel_lower_body_strength$ci.ub[1],2)}]"),
           size = 3) +
  labs(
    x = "Standardised Mean Change",
    y = "Study",
    color = "Group",
    title = "Lower Body Strength Outcomes",
    subtitle = "Vertical lines and bands are the meta-analytic point and 95% confidence interval estimates\nIndividual study level effects shown with point and confidence interval estimates"
  ) +
  theme_bw()

ggsave("lower_body_strength_plot.png", height = 5, width = 10, dpi = 300)




