#' get range of development patterns from fit_development_pattern_ml function
#'
#' @param ml_results the list object returned by running the fit_development_pattern_ml
#' @param num_fits the number of fits to include in the range. e.g. if set to 25 it will take the best 25 fitting patterns from the ml_results list object.
#' @param metric which goodness of fit metric to use when selecting the best fits. ave_score, cdr_score, huber_loss, rmse, mae or neg_ave_score.
#' @return a list with a plot showing the range of percentage developed at each development month and a tibble with the monthly pattern for each of the fits.
#'
#' @examples
#' \donttest{
#'
#' # set up a small parameters file for illustration
#' cl_parameters <- tidyr::expand_grid(smooth_from = 999,
#'                                     exclude_last_diag = c(TRUE, FALSE),
#'                                     exclude_high = FALSE,
#'                                     exclude_low = FALSE,
#'                                     selected_curve = "weibull",
#'                                     num_periods = c(3, 5),
#'                                     future_dev_periods = 0)
#'
#' # run function (num_cores = 1 for sequential processing)
#' ml_result <- fit_development_pattern_ml(uw_year,
#'                                         dev_year,
#'                                         claim_number,
#'                                         triangle_data, 1, 12,
#'                                         cl_parameters,
#'                                         num_cores = 1)
#'
#' # get range of top 2 best fitting patterns
#' ml_pattern_range(ml_result, 2)
#'
#' }
#'
#' @export

ml_pattern_range <- function(ml_results, num_fits = 25, metric = ave_score) {

  # validations

  allowed_vars <- c("ave_score", "cdr_score", "huber_loss", "rmse", "mae", "neg_ave_score")

  if (!as.character(rlang::ensym(metric)) %in% allowed_vars) {
    stop(paste0("error: 'metric' should be one of ", paste(allowed_vars, collapse = ", "), "."))
  }

  if (is.character(substitute(metric))) {
    stop("error: 'metric' should be provided without quotation marks.")
  }

  if(num_fits > nrow(ml_results$results)) stop("num_fits must be less than or equal to the number of rows in the results dataframe")


  # arrange by metric
  ml_sel <- ml_results[["results"]] |>
    dplyr::mutate(neg_ave_score = -neg_ave_score) |> # as currently this metric has the best fit as the largest number
    dplyr::arrange({{ metric }})

  # get exposure base
  exp_base = ml_results[["exposure_base"]]

  # rename column in exposure base
  if(!is.null(exp_base)) {

    exp_base <- exp_base |>
      dplyr::rename(premium = exposure)

  }

  # check if method column exists in the results and insert if not
  if(!("method" %in% colnames(ml_sel))) {

    ml_sel <- ml_sel |>
      dplyr::mutate(method = "cl")

  }

  # get the data frame stored in ml_inputs
  data <- ml_results$ml_inputs$data


  # fit patterns
  pattern_list <- purrr::map(
    .x = 1:num_fits,
    ~ {
      # Extract the method type from ml_sel for the current iteration
      method <- ml_sel[[.x, "method"]]

      # Fit the pattern based on the method type
      if (method == "cl") {
        fit <- actuary::fit_development_pattern(
          cohort_var = !!rlang::sym(ml_results$ml_inputs$cohort_var),
          dev_var = !!rlang::sym(ml_results$ml_inputs$dev_var),
          weighting_var = !!rlang::sym(ml_results$ml_inputs$weighting_var),
          data = data,
          dev_period_length = ml_results$ml_inputs$dev_period_length,
          dev_period_units = ml_results$ml_inputs$dev_period_units,
          exclude_last_diag = ml_sel[[.x, "exclude_last_diag"]],
          smooth_from = ml_sel[[.x, "smooth_from"]],
          exclude_high = ml_sel[[.x, "exclude_high"]],
          exclude_low = ml_sel[[.x, "exclude_low"]],
          selected_curve = ml_sel[[.x, "selected_curve"]],
          num_periods = ml_sel[[.x, "num_periods"]],
          future_dev_periods = ml_sel[[.x, "future_dev_periods"]]
        )
      } else if (method == "bf") {
        fit <- actuary::fit_development_pattern(
          cohort_var = !!rlang::sym(ml_results$ml_inputs$cohort_var),
          dev_var = !!rlang::sym(ml_results$ml_inputs$dev_var),
          weighting_var = !!rlang::sym(ml_results$ml_inputs$weighting_var),
          data = data,
          dev_period_length = ml_results$ml_inputs$dev_period_length,
          dev_period_units = ml_results$ml_inputs$dev_period_units,
          exclude_last_diag = ml_sel[[.x, "exclude_last_diag"]],
          smooth_from = ml_sel[[.x, "smooth_from"]],
          exclude_high = ml_sel[[.x, "exclude_high"]],
          exclude_low = ml_sel[[.x, "exclude_low"]],
          selected_curve = ml_sel[[.x, "selected_curve"]],
          num_periods = ml_sel[[.x, "num_periods"]],
          future_dev_periods = ml_sel[[.x, "future_dev_periods"]],
          premium = exp_base,
          bf_prior = ml_sel[[.x, "bf_priors"]] # Add bf_prior argument for bf method
        )
      } else if (method == "cc") {
        fit <- actuary::fit_development_pattern(
          cohort_var = !!rlang::sym(ml_results$ml_inputs$cohort_var),
          dev_var = !!rlang::sym(ml_results$ml_inputs$dev_var),
          weighting_var = !!rlang::sym(ml_results$ml_inputs$weighting_var),
          data = data,
          dev_period_length = ml_results$ml_inputs$dev_period_length,
          dev_period_units = ml_results$ml_inputs$dev_period_units,
          exclude_last_diag = ml_sel[[.x, "exclude_last_diag"]],
          smooth_from = ml_sel[[.x, "smooth_from"]],
          exclude_high = ml_sel[[.x, "exclude_high"]],
          exclude_low = ml_sel[[.x, "exclude_low"]],
          selected_curve = ml_sel[[.x, "selected_curve"]],
          num_periods = ml_sel[[.x, "num_periods"]],
          future_dev_periods = ml_sel[[.x, "future_dev_periods"]],
          premium = exp_base,
          cc_decay_factor = ml_sel[[.x, "decay_factors"]] # Add decay_factor argument for cc method
        )
      }

      # Return the fitted pattern with the fit number
      fit$monthly_pattern |> dplyr::mutate(fit_number = .x)
    }
  ) |>
    dplyr::bind_rows()

  # add method back on to patterns list

  pattern_list <- pattern_list |>
    dplyr::left_join(
      ml_sel |>
        dplyr::mutate(fit_number = dplyr::row_number()) |>
        dplyr::select(method, fit_number),
      dplyr::join_by(fit_number)
    )

  # produce plot showing range of development

  summary_data <- pattern_list |>
    dplyr::group_by(dev_month) |>
    dplyr::summarise(
      min_pct_dev = min(pct_dev),
      max_pct_dev = max(pct_dev)
    )

  plot <- ggplot2::ggplot(summary_data, ggplot2::aes(x = dev_month)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = min_pct_dev, ymax = max_pct_dev), fill = "skyblue", alpha = 0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      x = "development month",
      y = "percentage developed"
    )

  # add dev pattern for best fit to the plot
  # add title

  plot <- plot +
    ggplot2::geom_line(
      data = ml_results$best_fit$monthly_pattern,
      ggplot2::aes(y = pct_dev),
      colour = "red3", linetype = "dashed", linewidth = 1
    ) +
    ggplot2::labs(title = "<span style = 'color : skyblue;'>range of development patterns from ML fits</span>  \n<span style = 'color : red3;'>pattern with best score</span>") +
    ggplot2::theme(plot.title = ggtext::element_markdown())

  return(list(plot = plot, patterns = pattern_list))

}


