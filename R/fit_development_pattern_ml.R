#' fit multiple development patterns and assess goodness of fit
#'
#' \enumerate{
#'    \item{creates a holdout sample as the last n diagonals in the triangle data provided where n is specified via the holdout_size parameter.}
#'    \item{fits development patterns to the rest of the data, one for each row in the params_data dataframe.}
#'    \item{the parameters for each of the development patterns are as specified in the params_data dataframe.}
#'    \item{calculates the expected incremental claim amounts in the next diagonal from the holdout sample based on the fitted development pattern.}
#'    \item{calculates the goodness of fit of the fitted development pattern by comparing the expected incremental amounts to the actual incremental amounts.}
#'    \item{six goodness of fit metrics are calculated : AvE, negative AvE, CDR (see reference document or details section below for details on both), \link[yardstick]{huber_loss}, \link[yardstick]{rmse}, \link[yardstick]{mae}}
#'    \item{add the first diagonal from the holdout sample to the training data and repeat steps 2 to 6.}
#'    \item{repeat until all diagonals from the holdout sample are exhausted and then calculate the average goodness of the goodness of fit metrics across all iterations.}
#'    \item{the total reserve implied by the pattern is calculated as the difference between the ultimate and the claim amount at the latest diagonal in the holdout data.}
#' }
#'
#' requires tidy data. one row per cell in the unprojected triangle. \cr
#' missing rows within the triangle will be filled in with the value in the previous period. \cr
#' no rows for development periods beyond the last evaluation point. \cr
#' data provided must be cumulative. \cr
#' AvE score is the weighted average of the AvE in each value of cohort_var where the weighting is the incremental value of weighting_var. \cr
#' Negative AvE score is the negative unweighted average of the AvE. if using this as the eval_metric note that unlike for all other metrics, higher is better. \cr
#' CDR score is the weighted average of the CDR in each value of cohort_var where the weighting is the incremental value of weighting_var. \cr
#' CDR is the AvE plus the change in the IBNR estimate. equivalent to the change in projected ultimate claim value over the time period.
#'
#' @param cohort_var the variable in the input data representing the cohort.
#' @param dev_var the variable in the input data representing the development period. needs to be consistent with dev_period_length and dev_period_units. e.g. if have half annual dev periods and data is in months (6, 12, 18 etc) then dev_period_length = 6 and dev_period_units = 1. if data in years (0.5, 1, 1.5 etc) then dev_period_length = 0.5 and dev_period_units = 12.
#' @param weighting_var the variable in the input data representing the amounts in the triangle e.g. claim counts or claim amounts etc.
#' @param data the cumulative tidy dataset with the triangle data. one row per cohort year and development period. no rows for development periods beyond the last evaluation point.
#' @param dev_period_length the difference between subsequent development periods in the units specified in dev_period_units.
#' @param dev_period_units the units in which dev_var is provided e.g. 1 for months, 12 for years etc.
#' @param params_data a dataframe with all combinations of parameters required for fitting the development patterns. must contain columns with headings : exclude_last_diag, smooth_from, exclude_high, exclude_low, selected_curve, num_periods. set smooth_from to 999 if you don't want to use curves at all. See fit_development_pattern documentation for more detail on the meaning of these parameters.
#' @param holdout_size the number of diagonals to holdout when fitting the development patterns. a
#' @param bf_priors if provided the algorithm will also assess the goodness of fit using the BF method to calculate reserves. can provide a single value or a vector of possible priors.
#' @param decay_factors if provided the algorithm will additionally assess the goodness of fit using the generalised Cape Cod method to set the reserves. can provide a single value or a vector of possible priors.
#' @param exposure_base must be provided if bf_priors or decay_factors is provided. a dataframe with a column called exposure (probably the premium) and another column which must have the same name as cohort_var.
#' @param eval_metric which goodness of fit metric to order the results dataframe and to select the best parameters for the best_fit object.
#' @param num_cores number of cores to use for parallel processing. defaults to the number available - 1.
#'
#' @examples
#' \dontrun{
#'
#' # set up parameters file
#' cl_parameters <- tidyr::expand_grid(smooth_from = c(1, 2),
#'                                     exclude_last_diag = c(TRUE, FALSE),
#'                                     exclude_high = c(TRUE, FALSE),
#'                                     exclude_low = c(TRUE, FALSE),
#'                                     selected_curve = c("weibull", "inverse_power", "exponential_decay"),
#'                                     num_periods = c(1:5),
#'                                     future_dev_periods = c(0, 25))
#'
#' # run function
#' ml_result <- fit_development_pattern_ml(uw_year, dev_year, claim_number, triangle_data, 1, 12, cl_parameters)
#'
#' # view results
#' ml_result$results
#'
#' # plot fitted pattern. note that it doesn't like a good fit because we only supplied two smooth_from values in the params_data.
#' plot_development_pattern(ml_result$best_fit)
#'
#' # view the plot of reserve against the goodness of fit metric
#' ml_result$results_plot
#'
#' # amend the plot so that selected curve is used for the colour of the points
#' ml_result$results_plot + ggplot2::geom_point(ggplot2::aes(colour = selected_curve))
#'
#' }

#'
#' @return
#'
#' \enumerate{
#'    \item{a dataframe with the five goodness of fit metrics, the reserve and the parameter values used for the fit. ordered by descending value of eval_metric.}
#'    \item{a fit_development_pattern object with the best fit according to the eval_metric.}
#'    \item{a plot of reserve against eval_metric.}
#'    \item{the exposure base if provided.}
#'    \item{a tibble with the values of some of the inputs.}
#' }
#'
#'
#' @references Richman, R and Balona, C (2020). The Actuary and IBNR Techniques: A Machine Learning Approach. https://ssrn.com/abstract=3697256
#'
#' @export
#' @seealso \code{\link{fit_development_pattern}}

fit_development_pattern_ml <- function(cohort_var, dev_var, weighting_var, data, dev_period_length, dev_period_units,
                                       params_data, holdout_size = 2, bf_priors = NULL, decay_factors = NULL,
                                       exposure_base = NULL, eval_metric = "ave_score", num_cores = parallel::detectCores() - 1) {

  # quoted strings
  cohort_var_name <- rlang::enquo(cohort_var)
  dev_var_name <- rlang::enquo(dev_var)
  weighting_var_name <- rlang::enquo(weighting_var)

  # if try_bf = TRUE then need to have an exposure_base and bf_priors
  # and exposure_base needs to contain exposure for all values of cohort_var
  # and columns need to be as specified in documentation
  if(!is.null(bf_priors) | !is.null(decay_factors)) {
    if(is.null(exposure_base)) stop("need to supply an exposure base if bf_priors are provided")
    cohort_var_string <- rlang::as_name(cohort_var_name)
    if(!(cohort_var_string %in% colnames(exposure_base))) stop("exposure_base needs to contain the same cohort name as in the triangle data provided")
    data_cohort_values <- data %>% dplyr::pull(!!cohort_var_name) %>% unique()
    exposure_base_cohort_values <- exposure_base %>% dplyr::pull(!!cohort_var_name) %>% unique()

    missing_cohorts <- setdiff(data_cohort_values, exposure_base_cohort_values)

    if(length(missing_cohorts) > 0) {
      stop(paste("The following cohort values are present in the triangle data but missing from exposure_base:",
                 paste(missing_cohorts, collapse = ", ")))
    }

    if(!("exposure" %in% colnames(exposure_base))) {
      stop("exposure_base must contain a column named 'exposure'")
    }

    # Check for additional columns in exposure_base
    expected_columns <- c(cohort_var_string, "exposure")
    extra_columns <- setdiff(colnames(exposure_base), expected_columns)

    if(length(extra_columns) > 0) {
      stop(paste("The exposure_base dataframe contains unexpected columns:", paste(extra_columns, collapse = ", ")))
    }

  }

  # checks on params data

  dev_var_values <- dplyr::pull(data, {{ dev_var }})
  valid_smooth_from_values <- c(999, unique(dev_var_values))

  if(!all(c("exclude_last_diag", "smooth_from", "exclude_high", "exclude_low", "selected_curve", "num_periods", "future_dev_periods") %in% colnames(params_data))) stop("params_data must contain the following columns: exclude_last_diag, smooth_from, exclude_high, exclude_low, selected_curve, num_periods, future_dev_periods")
  if(!all(unique(params_data$selected_curve) %in% c("exponential_decay", "weibull", "inverse_power"))) stop("the selected_curve column in params_data must contain values exponential_decay, inverse_power or weibull")
  if(!all(unique(params_data$exclude_last_diag) %in% c(TRUE, FALSE))) stop("the exclude_last_diag column in params_data must contain values TRUE or FALSE")
  if(!all(unique(params_data$exclude_high) %in% c(TRUE, FALSE))) stop("the exclude_high column in params_data must contain values TRUE or FALSE")
  if(!all(unique(params_data$exclude_low) %in% c(TRUE, FALSE))) stop("the exclude_low column in params_data must contain values TRUE or FALSE")
  if(!all(params_data$smooth_from %in% valid_smooth_from_values)) {
    stop("The smooth_from column in params_data must contain values that correspond to the values of dev_var in the data provided or 999. Set to 999 if you don't want to smooth.")
  }

  # check on eval_metric
  if(!(eval_metric %in% c("ave_score", "cdr_score", "huber_loss", "rmse", "mae", "neg_ave_score"))) stop("eval_metric needs to be one of rmse, huber_loss, mae, cdr_score, ave_score or neg_ave_score")

  # Create a progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = nrow(params_data),
    clear = FALSE,
    width = 60
  )

  # Set up parallel processing
  future::plan(future::multisession, workers = num_cores)

  # Set up progress handler
  # progressr::handlers("txtprogressbar")
  # p <- progressr::progressor(steps = nrow(params_data))

  # Use future_map instead of a for loop
  data_store <- furrr::future_map(1:nrow(params_data), function(i) {

    # Inform progressr about progress
    #p(message = sprintf("Processing parameter set %d of %d", i, nrow(params_data)))

    # Update progress bar
    pb$tick()

    t_values <- 0:(holdout_size - 1)

    metric_store <- furrr::future_map(seq_along(t_values), function(j) {
      t <- t_values[j]
      training_iteration <- t + 1

      # remove diagonals from data
      data_fit <- data |>
        dplyr::group_by(!!cohort_var_name) |>
        dplyr::arrange(!!dev_var_name, .by_group = TRUE) |>
        dplyr::slice_head(n = -(holdout_size - t))

      # fit cl pattern with parameters specified in params_data
      fit <- actuary::fit_development_pattern(
        cohort_var =  !!cohort_var_name,
        dev_var =  !!dev_var_name,
        weighting_var = !!weighting_var_name,
        data = data_fit,
        dev_period_length = dev_period_length,
        dev_period_units = dev_period_units,
        exclude_last_diag = params_data[[i, "exclude_last_diag"]],
        smooth_from = params_data[[i, "smooth_from"]],
        exclude_high = params_data[[i, "exclude_high"]],
        exclude_low = params_data[[i, "exclude_low"]],
        selected_curve = params_data[[i, "selected_curve"]],
        num_periods = params_data[[i, "num_periods"]],
        future_dev_periods = params_data[[i, "future_dev_periods"]]
      )

      # test on next diagonal

      data_test <- data |>
        dplyr::group_by(!!cohort_var_name) |>
        dplyr::arrange(!!dev_var_name, .by_group = TRUE) |>
        dplyr::mutate(row_number = dplyr::row_number(),
                      total_rows = dplyr::n()) |>
        dplyr::filter(row_number <= total_rows - (holdout_size - 1) + t) |>
        dplyr::select(-row_number, -total_rows)

      holdout_flag <- 1 # as just want to flag the last diagonal

      # merge development_pattern back on to original data
      result <- data_test |>
        dplyr::left_join(fit$development_pattern, dplyr::join_by(!!dev_var_name == dev_period)) |>
        dplyr::mutate(pct_dev = tidyr::replace_na(pct_dev, 1)) |>
        dplyr::mutate(ultimate = !!weighting_var_name / pct_dev) |>
        dplyr::group_by(!!cohort_var_name) |>
        dplyr::arrange(!!dev_var_name, .by_group = TRUE) |>
        dplyr::mutate(incr_expected = dplyr::lag(!!weighting_var_name, default = 0) / dplyr::lag(pct_dev, default = 0) * pct_dev -
                        dplyr::lag(!!weighting_var_name, default = 0),
                      incr_actual = !!weighting_var_name - dplyr::lag(!!weighting_var_name, default = 0),
                      ave = incr_actual - incr_expected,
                      expected_amount = dplyr::lag(!!weighting_var_name, default = 0) / dplyr::lag(pct_dev, default = 1) * pct_dev,
                      expected_ibnr = dplyr::lag(ultimate, default = 0) - expected_amount,
                      actual_ibnr = ultimate - !!weighting_var_name,
                      change_in_ibnr = actual_ibnr - expected_ibnr,
                      cdr = ave + change_in_ibnr,
                      reserve = ultimate - !!weighting_var_name) |>
        dplyr::mutate(diag_num = dplyr::n() - dplyr::row_number() + 1,
                      holdout_flag = dplyr::if_else(diag_num <= holdout_flag, 1, 0)) |>
        dplyr::as_tibble()

      # calculate objective functions
      obj_fn_data <- result |>
        dplyr::group_by(!!cohort_var_name) |>
        dplyr::arrange(!!dev_var_name, .by_group = TRUE) |>
        dplyr::mutate(reserve = dplyr::last(reserve)) |>
        dplyr::ungroup() |>
        dplyr::filter(holdout_flag == 1) |>
        dplyr::select(incr_expected, incr_actual, ave, change_in_ibnr, cdr, reserve) |>
        dplyr::filter(!is.nan(incr_expected)) |>
        dplyr::summarise(rmse = yardstick::rmse_vec(truth = incr_actual, estimate = incr_expected),
                         huber_loss = yardstick::mae_vec(truth = incr_actual, estimate = incr_expected),
                         mae = yardstick::mae_vec(truth = incr_actual, estimate = incr_expected),
                         cdr_score = sqrt(sum(abs(incr_actual) * cdr ^ 2) / sum(abs(incr_actual))),
                         ave_score = sqrt(sum(abs(incr_actual) * ave ^ 2) / sum(abs(incr_actual))),
                         neg_ave_score = -sqrt(mean((incr_actual - incr_expected)^2)),
                         reserve = mean(reserve)) |>
        dplyr::mutate(method = "cl")

      if(!is.null(bf_priors)) {

        obj_fn_data_bf <- result |>
          dplyr::select(!!cohort_var_name, !!dev_var_name, !!weighting_var_name, pct_dev, holdout_flag) |>
          dplyr::left_join(exposure_base, dplyr::join_by(!!cohort_var_name)) |>
          tidyr::expand_grid(bf_priors) |>
          dplyr::group_by(!!cohort_var_name, bf_priors) |>
          dplyr::arrange(!!dev_var_name, .by_group = TRUE) |>
          dplyr::mutate(ultimate = !!weighting_var_name + (1 - pct_dev) * exposure * bf_priors,
                        incr_expected = (pct_dev - dplyr::lag(pct_dev, default = 0)) * exposure * bf_priors,
                        incr_actual = !!weighting_var_name - dplyr::lag(!!weighting_var_name, default = 0),
                        cdr = ultimate - dplyr::lag(ultimate, default = 0),
                        ave = incr_actual - incr_expected,
                        change_in_ibnr = cdr - ave,
                        reserve = ultimate - !!weighting_var_name,
                        reserve = dplyr::last(reserve)) |>
          dplyr::ungroup() |>
          dplyr::filter(holdout_flag == 1) |>
          dplyr::group_by(bf_priors) |>
          dplyr::arrange(!!dev_var_name, .by_group = TRUE) |>
          dplyr::summarise(rmse = yardstick::rmse_vec(truth = incr_actual, estimate = incr_expected),
                    huber_loss = yardstick::mae_vec(truth = incr_actual, estimate = incr_expected),
                    mae = yardstick::mae_vec(truth = incr_actual, estimate = incr_expected),
                    cdr_score = sqrt(sum(abs(incr_actual) * cdr ^ 2) / sum(abs(incr_actual))),
                    ave_score = sqrt(sum(abs(incr_actual) * ave ^ 2) / sum(abs(incr_actual))),
                    neg_ave_score = -sqrt(mean((incr_actual - incr_expected)^2)),
                    reserve = mean(reserve)) |>
          dplyr::mutate(method = "bf")

        obj_fn_data <- obj_fn_data |>
          dplyr::bind_rows(obj_fn_data_bf)

      }

      if(!is.null(decay_factors)) {

        # get columns required
        cc_calc <- result |>
          dplyr::select(!!cohort_var_name, !!dev_var_name, !!weighting_var_name, pct_dev, holdout_flag) |>
          dplyr::left_join(exposure_base, join_by(!!cohort_var_name))

        # only need to calculate ultimates in holdout sample and diagonal prior to the holdout sample

        cc_calc <- cc_calc |>
          dplyr::group_by(!!cohort_var_name) |>
          dplyr::mutate(temp = dplyr::row_number() * (1 - holdout_flag),
                        holdout_flag2 = dplyr::if_else(temp == max(temp), 1, 0),
                        holdout_flag3 = holdout_flag + holdout_flag2) |>
          dplyr::filter(holdout_flag3 == 1) |>
          dplyr::mutate(group_num = dplyr::row_number())

        # function to calculate CC lr

        calculate_cape_cod_lr <- function(data, decay_factors) {

          years <- data |>
            dplyr::pull(row_num)

          purrr::map_dfr(decay_factors, function(decay_factors) {
            data |>
              dplyr::mutate(lambda = purrr::map(row_num,
                                                ~ decay_factors ^ base::abs(years - .x))) |>
              dplyr::mutate(lr_numerator = purrr::map_dbl(lambda, ~ base::sum(!!weighting_var_name * .x)),
                            lr_denominator = purrr::map_dbl(lambda, ~ base::sum(exposure * pct_dev * .x))) |>
              dplyr::mutate(lr = lr_numerator / lr_denominator,
                            decay_factors = decay_factors)
          })
        }

        cc_results <- list()

        # for each group calculate cc lr

        for(k in 1:max(cc_calc$group_num)) {

          cc_calc_part <- cc_calc |>
            dplyr::ungroup() |>
            dplyr::filter(group_num == k) |>
            dplyr::mutate(row_num = dplyr::row_number())

          cc_results[[k]] <- calculate_cape_cod_lr(cc_calc_part, decay_factors)

        }

        cc_results <- cc_results |>
          dplyr::bind_rows() |>
          dplyr::select(!!cohort_var_name, !!dev_var_name, decay_factors, !!weighting_var_name, pct_dev, holdout_flag, exposure, lr)


        obj_fn_data_cc <- cc_results |>
          dplyr::group_by(!!cohort_var_name, decay_factors) |>
          dplyr::arrange(!!dev_var_name, .by_group = TRUE) |>
          dplyr::mutate(ultimate = !!weighting_var_name + (1 - pct_dev) * exposure * lr,
                        # for AvE need to calculate expected using cape cod LR at last diagaonal of training data
                        incr_expected = (pct_dev - dplyr::lag(pct_dev, default = 0)) * exposure * first(lr),
                        incr_actual = !!weighting_var_name - dplyr::lag(!!weighting_var_name, default = 0),
                        cdr = ultimate - dplyr::lag(ultimate, default = 0),
                        ave = incr_actual - incr_expected,
                        change_in_ibnr = cdr - ave,
                        reserve = ultimate - !!weighting_var_name,
                        reserve = dplyr::last(reserve)) |>
          dplyr::ungroup() |>
          dplyr::filter(holdout_flag == 1) |>
          dplyr::group_by(decay_factors) |>
          dplyr::arrange(!!dev_var_name, .by_group = TRUE) |>
          dplyr::summarise(rmse = yardstick::rmse_vec(truth = incr_actual, estimate = incr_expected),
                    huber_loss = yardstick::mae_vec(truth = incr_actual, estimate = incr_expected),
                    mae = yardstick::mae_vec(truth = incr_actual, estimate = incr_expected),
                    cdr_score = sqrt(sum(abs(incr_actual) * cdr ^ 2) / sum(abs(incr_actual))),
                    ave_score = sqrt(sum(abs(incr_actual) * ave ^ 2) / sum(abs(incr_actual))),
                    neg_ave_score = -sqrt(mean((incr_actual - incr_expected)^2)),
                    reserve = mean(reserve)) |>
          dplyr::mutate(method = "cc")

        obj_fn_data <- obj_fn_data |>
          dplyr::bind_rows(obj_fn_data_cc)

      }



      return(obj_fn_data)


    })

    # calculate average metrics

    if(!is.null(bf_priors) | !is.null(decay_factors)) {

      metric_store <- metric_store |>
        dplyr::bind_rows()

      if(!("decay_factors" %in% colnames(metric_store))) metric_store[["decay_factors"]] <- NA_real_
      if(!("bf_priors" %in% colnames(metric_store))) metric_store[["bf_priors"]] <- NA_real_

      metric_store <- metric_store |>
        dplyr::group_by(method, bf_priors, decay_factors) |>
        # take reserve at latest diagonal
        dplyr::mutate(reserve = dplyr::last(reserve)) |>
        # average of metrics across all training sets
        dplyr::summarise(dplyr::across(c(huber_loss, rmse, mae, cdr_score, ave_score, neg_ave_score, reserve), mean), .groups = "drop")

    } else {

      metric_store <- metric_store |>
        dplyr::bind_rows() |>
        # take reserve at latest diagonal
        dplyr::mutate(reserve = dplyr::last(reserve)) |>
        # average of metrics across all training sets
        dplyr::summarise(dplyr::across(c(huber_loss, rmse, mae, cdr_score, ave_score, neg_ave_score, reserve), mean))

    }

    # store results
    metric_store |>
      dplyr::mutate(
        exclude_last_diag = params_data[[i, "exclude_last_diag"]],
        smooth_from = params_data[[i, "smooth_from"]],
        exclude_high = params_data[[i, "exclude_high"]],
        exclude_low = params_data[[i, "exclude_low"]],
        selected_curve = params_data[[i, "selected_curve"]],
        num_periods = params_data[[i, "num_periods"]],
        future_dev_periods = params_data[[i, "future_dev_periods"]]
      )
  }, .progress = TRUE)

  # Clean up parallel processing
  future::plan(future::sequential)

  # results dataframe
  data_store <- dplyr::bind_rows(data_store) |>
    dplyr::mutate(neg_ave_score = - neg_ave_score) |>
    dplyr::arrange(!!rlang::sym(eval_metric)) |>
    dplyr::mutate(neg_ave_score = -neg_ave_score) |>
    dplyr::select(exclude_last_diag:future_dev_periods, ave_score, dplyr::everything())

  # fit with best fit

  if ("method" %in% names(data_store)) {

    if(data_store[[1, "method"]] == "cc") {

      best_fit <- actuary::fit_development_pattern(
        cohort_var =  !!cohort_var_name,
        dev_var =  !!dev_var_name,
        weighting_var = !!weighting_var_name ,
        data = data,
        dev_period_length = {{ dev_period_length }},
        dev_period_units = {{ dev_period_units }},
        exclude_last_diag = data_store[[1, "exclude_last_diag"]],
        smooth_from = data_store[[1, "smooth_from"]],
        exclude_high = data_store[[1, "exclude_high"]],
        exclude_low = data_store[[1, "exclude_low"]],
        selected_curve = data_store[[1, "selected_curve"]],
        num_periods = data_store[[1, "num_periods"]],
        future_dev_periods = data_store[[1, "future_dev_periods"]],
        cc_decay_factor = data_store[[1, "decay_factors"]],
        premium = exposure_base |> dplyr::rename(premium = exposure)
      )

    } else if(data_store[[1, "method"]] == "bf") {

      best_fit <- actuary::fit_development_pattern(
        cohort_var =  !!cohort_var_name,
        dev_var =  !!dev_var_name,
        weighting_var = !!weighting_var_name ,
        data = data,
        dev_period_length = {{ dev_period_length }},
        dev_period_units = {{ dev_period_units }},
        exclude_last_diag = data_store[[1, "exclude_last_diag"]],
        smooth_from = data_store[[1, "smooth_from"]],
        exclude_high = data_store[[1, "exclude_high"]],
        exclude_low = data_store[[1, "exclude_low"]],
        selected_curve = data_store[[1, "selected_curve"]],
        num_periods = data_store[[1, "num_periods"]],
        future_dev_periods = data_store[[1, "future_dev_periods"]],
        bf_prior = data_store[[1, "bf_priors"]],
        premium = exposure_base |> dplyr::rename(premium = exposure)
      )
    } else {

      best_fit <- actuary::fit_development_pattern(
        cohort_var =  !!cohort_var_name,
        dev_var =  !!dev_var_name,
        weighting_var = !!weighting_var_name ,
        data = data,
        dev_period_length = {{ dev_period_length }},
        dev_period_units = {{ dev_period_units }},
        exclude_last_diag = data_store[[1, "exclude_last_diag"]],
        smooth_from = data_store[[1, "smooth_from"]],
        exclude_high = data_store[[1, "exclude_high"]],
        exclude_low = data_store[[1, "exclude_low"]],
        selected_curve = data_store[[1, "selected_curve"]],
        num_periods = data_store[[1, "num_periods"]],
        future_dev_periods = data_store[[1, "future_dev_periods"]]
      )

    }

  } else {

    best_fit <- actuary::fit_development_pattern(
      cohort_var =  !!cohort_var_name,
      dev_var =  !!dev_var_name,
      weighting_var = !!weighting_var_name ,
      data = data,
      dev_period_length = {{ dev_period_length }},
      dev_period_units = {{ dev_period_units }},
      exclude_last_diag = data_store[[1, "exclude_last_diag"]],
      smooth_from = data_store[[1, "smooth_from"]],
      exclude_high = data_store[[1, "exclude_high"]],
      exclude_low = data_store[[1, "exclude_low"]],
      selected_curve = data_store[[1, "selected_curve"]],
      num_periods = data_store[[1, "num_periods"]],
      future_dev_periods = data_store[[1, "future_dev_periods"]]
    )

  }

  # create plot

  results_plot <- ggplot2::ggplot(data_store, ggplot2::aes(.data[[eval_metric]], reserve)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(x = eval_metric, y = "reserve") +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "left",
      panel.grid.minor = ggplot2::element_blank()
    )

  # capture the name of the data frame
  data_name <- deparse(substitute(data))

  # return the inputs as well
  ml_inputs <- dplyr::tibble(
    cohort_var = as.character(rlang::ensym(cohort_var)),
    dev_var = as.character(rlang::ensym(dev_var)),
    weighting_var = as.character(rlang::ensym(weighting_var)),
    data_name = data_name,
    dev_period_length = dev_period_length,
    dev_period_units = dev_period_units
  )

  return(list(results = data_store,
              best_fit = best_fit,
              results_plot = results_plot,
              exposure_base = exposure_base,
              ml_inputs = ml_inputs))

}
