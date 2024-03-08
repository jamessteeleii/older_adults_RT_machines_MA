# _targets.R file
library(targets)
library(tarchetypes)
source("R/functions.R")
tar_option_set(packages = c("here", "readxl", "janitor", "tidyverse", "base", "scales", "ggtext", "ggh4x", "zoo",
                            "performance", "see", "rstan", "brms", "bayesplot", "marginaleffects", "broom.mixed",
                            "patchwork", "kableExtra", "knitr", "quarto", "officer", "officedown",
                            "lactater"))

list(
  # Load in data
  tar_target(file, here("data","data.csv"), format = "file"),
  tar_target(data, read_prepare_data(file)),
  
  # Fit model
  tar_target(model, fit_model(data)),
  tar_target(tidy_model, get_tidy_model(model)),
  
  # Model checks
  tar_target(model_checks, make_model_checks_tiff(model)),
  
  # # Diagnostic plots
  tar_target(rhat_model, make_rhat_plot(model)),
  tar_target(trace_model, make_trace_plots(model)),
  tar_target(pp_check_model, make_pp_check(model)),
  
  # Calculate thresholds and their agreement
  tar_target(thresholds, calculate_thresholds(data)),
  tar_target(thresholds_agree, calculate_thresholds_agree(thresholds)),
  
  # Make and save plots
  tar_target(individual_data_plot, plot_individual_data(data)),
  tar_target(individual_data_plot_tiff, make_individual_data_plot_tiff(individual_data_plot)),
  
  tar_target(model_plot, plot_model(data, model)),
  tar_target(model_plot_tiff, make_model_plot_tiff(model_plot)),
  
  tar_target(individual_preds_plot, plot_individual_preds(data, model)),
  tar_target(individual_preds_plot_tiff, make_individual_preds_plot_tiff(individual_preds_plot)),
  
  tar_target(main_plot, combine_plots(individual_data_plot, individual_preds_plot, model_plot, thresholds_agree_plot)),
  
  tar_target(thresholds_agree_plot, plot_thresholds_agree(thresholds, thresholds_agree)),
  tar_target(thresholds_agree_plot_tiff, make_thresholds_agree_plot_tiff(thresholds_agree_plot))
  
  # Render the report
  # tar_quarto(report, "report.qmd")
  
  # # Render the supplementary material
  # tar_quarto(diagnostics_plots, "diagnostics_plots.qmd")
  
  
)
