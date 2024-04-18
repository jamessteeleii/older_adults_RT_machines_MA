# _targets.R file
library(targets)
library(tarchetypes)
library(crew)
source("R/functions.R")
tar_option_set(
  controller = crew_controller_local(workers = 2),
  packages = c(
    "here",
    "tidyverse",
    "base",
    "metafor",
    # "scales",
    # "ggtext",
    # "ggh4x",
    # "zoo",
    # "performance",
    # "see",
    "rstan",
    "brms",
    "tidybayes",
    "bayesplot",
    "marginaleffects",
    "broom.mixed",
    "patchwork",
    # "kableExtra",
    # "knitr",
    "quarto"
    # "officer",
    # "officedown",
  )
)

list(
  # Load in and prepare data
  tar_target(file, here("data","data.csv"), format = "file"),
  tar_target(data, read_prepare_data(file)),
  tar_target(data_function, filter_function_data(data)),
  tar_target(data_strength, filter_strength_data(data)),
  
  # Fit models
  tar_target(model_function, fit_model(data_function)),
  tar_target(tidy_model_function, get_tidy_model(model_function)),
  tar_target(model_strength, fit_model(data_strength)),
  tar_target(tidy_model_strength, get_tidy_model(model_strength)),
  
  # Get posterior draws
  tar_target(posterior_draws_meta_function, get_posterior_draws_meta(model_function)),
  tar_target(posterior_draws_study_function, get_posterior_draws_study(data_function, model_function)),
  tar_target(posterior_draws_contrast_function, get_posterior_draws_contrast(model_function)),
  tar_target(posterior_draws_meta_strength, get_posterior_draws_meta(model_strength)),
  tar_target(posterior_draws_study_strength, get_posterior_draws_study(data_strength, model_strength)),
  tar_target(posterior_draws_contrast_strength, get_posterior_draws_contrast(model_strength)),
  
  # Make plots
  tar_target(meta_function_plot, plot_meta(posterior_draws_meta_function)),
  tar_target(study_function_plot, plot_study(data_function, posterior_draws_study_function)),
  tar_target(contrast_function_plot, plot_contrast(posterior_draws_contrast_function)),
  tar_target(meta_strength_plot, plot_meta(posterior_draws_meta_strength)),
  tar_target(study_strength_plot, plot_study(data_strength, posterior_draws_study_strength)),
  tar_target(contrast_strength_plot, plot_contrast(posterior_draws_contrast_strength)),
  
  tar_target(combined_function_plot, combine_plots(study_function_plot, meta_function_plot, contrast_function_plot,
                                                   "Meta-Analysis of Prior Studies Examining the Effects of Machine Based Resistance Training on Function")),
  
  tar_target(combined_strength_plot, combine_plots(study_strength_plot, meta_strength_plot, contrast_strength_plot,
                                                   "Meta-Analysis of Prior Studies Examining the Effects of Machine Based Resistance Training on Strength")),
  
  tar_target(combined_function_plot_tiff, ggsave(combined_function_plot, filename = "plots/meta_plots_function.tiff", dpi = 300, w=10, h=10)),
  tar_target(combined_strength_plot_tiff, ggsave(combined_strength_plot, filename = "plots/meta_plots_strength.tiff", dpi = 300, w=10, h=10)),

  # # Diagnostic plots
  tar_target(rhat_model_function, make_rhat_plot(model_function)),
  tar_target(trace_model_function, make_trace_plots(model_function)),
  tar_target(pp_check_model_function, make_pp_check(model_function)),
  tar_target(rhat_model_strength, make_rhat_plot(model_strength)),
  tar_target(trace_model_strength, make_trace_plots(model_strength)),
  tar_target(pp_check_model_strength, make_pp_check(model_strength))
  
  # # Render the supplementary material
  # tar_quarto(diagnostics_plots, "diagnostics_plots.qmd")
  
  # Render the report
  # tar_quarto(report, "report.qmd")
  

  
)
