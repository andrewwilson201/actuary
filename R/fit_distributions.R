#' fit distributions to data
#'
#' fit multiple distributions and return goodness of fit exhibits together with the fitted parameters\cr
#' all distributions are fitted using maximum likelihood
#'
#' @param data a vector of values to which the distributions should be fitted
#' @param risk_shift the lower bound of the input data. defaults to zero but if specified otherwise then the input data will be filtered to remove any amounts less than the risk_shift
#' @param distributions the distributions to be fitted to the input data. defaults to lnorm, pareto, weibull, gamma and invgauss
#'
#' @return
#'
#'  \enumerate{
#'    \item{the name of the best fitting distribution according to the AIC}
#'    \item{a table showing the AIC for all distributions}
#'    \item{a table showing the parameters for the fitted distributions}
#'    \item{a plot of the CDF}
#'    \item{a Q-Q plot}
#'    \item{a table with percentiles for all distributions}
#'    \item{a summary statistics table based on simulating 100k values from the fitted distributions}
#'    \item{a plot of the summary statistics}
#'
#' }
#'
#' @examples
#' \dontrun{
#' # fit distributions to the losses data
#' losses_sample <- dplyr::slice_sample(losses, n = 10000)
#' fit_distributions(losses_sample$amount, risk_shift = min(losses_sample$amount))
#' }
#'
#' @export

fit_distributions <- function(data, risk_shift = 0, distributions = c("lnorm", "pareto", "weibull", "gamma", "invgauss")) {


  # warn if there are zero or negative values in vector provided
  if(any(data == 0)) warning("the data provided contains zero values. this may result in unreliable fits. check your input data.")
  if(any(data < 0)) warning("the data provided contains negative values. this may result in unreliable fits. check you input data.")

  # shift the input data to filter out anything below the risk_shift
  # and then to shift all amounts by the risk_shift

  shifted_data <- data[data > risk_shift] - risk_shift

  # define custom gamma functions
  dgamma2 <<- function(x, shape, scale, log = FALSE) {
    dgamma(x, shape = shape, scale = scale, log = log)
  }
  pgamma2 <<- function(q, shape, scale) pgamma(q, shape = shape, scale = scale)
  qgamma2 <<- function(p, shape, scale) qgamma(p, shape = shape, scale = scale)
  rgamma2 <<- function(p, shape, scale) rgamma(p, shape = shape, scale = scale)

  # fit all distributions required

  fits <- list()

  for (dist in distributions) {
    tryCatch({
      if (dist == "pareto") {
        # for pareto supply starting values
        fit <- fitdistrplus::fitdist(shifted_data, distr = dpareto_custom, start = list(shape = 1, scale = min(shifted_data)))
      } else if (dist == "gamma") {
        # for gamma manually parameterise the distribution to use shape and scale rather than the defaults of shape and rate (fits more easily this way)
        # also normalise before fitting and supply sensible starting values
        shift_normalised <- shifted_data / mean(shifted_data)
        fit <- fitdistrplus::fitdist(shift_normalised, dgamma2, start = list(shape = 1, scale = 1),
                                     lower = c(0.01, 0.01),
                                     upper = c(100, 100))

        # denormalise the scale parameter after the fit
        original_scale <- fit$estimate["scale"] * mean(shifted_data)

        # update the fit object with the denormalised scale
        fit$estimate["scale"] <- original_scale

        # recalculate log likelihood
        new_loglik <- sum(dgamma2(shifted_data, shape = fit$estimate["shape"], scale = original_scale, log = TRUE))
        fit$loglik <- new_loglik

        # recalculate AIC
        fit$aic <- -2 * new_loglik + 2 * length(fit$estimate)

      } else if (dist == "invgauss") {

        fit <- fitdistrplus::fitdist(shifted_data, distr = dinvgauss_custom, start = list(mean = mean(shifted_data), shape = 1))

      } else {
        # for lognormal and weibull just fit
        fit <- fitdistrplus::fitdist(shifted_data, dist)
      }

      fits[[dist]] <- fit
    }, error = function(e) {
      message(paste("Failed to fit", dist, "distribution:", e$message))
    })
  }

  # create AIC table

  aic_tibble <- dplyr::tibble(
    distribution = names(fits),
    aic = sapply(fits, function(x) x$aic)
  ) |>
    dplyr::arrange(aic)

  # determine best fit based on AIC

  best_fit <- aic_tibble$distribution[1]

  # extract fit parameters

  fit_parameters <- purrr::map_dfr(names(fits), function(dist) {
    params <- fits[[dist]]$estimate
    param_names <- names(params)
    dplyr::tibble(
      distribution = dist,
      parameter_1 = dplyr::if_else(length(params) >= 1, params[1], NA_real_),
      parameter_2 = dplyr::if_else(length(params) >= 2, params[2], NA_real_),
      param_name1 = dplyr::if_else(length(param_names) >= 1, param_names[1], NA_character_),
      param_name2 = dplyr::if_else(length(param_names) >= 2, param_names[2], NA_character_),
      risk_shift = risk_shift
    )
  })

  fit_parameters <- fit_parameters |>
    # swap shape and scale around for pareto2 to be consistent with create_simulations function
    dplyr::mutate(
      parameter_1_new = dplyr::if_else(distribution == "pareto", parameter_2, parameter_1),
      parameter_2_new = dplyr::if_else(distribution == "pareto", parameter_1, parameter_2),
      param_name1_new = dplyr::if_else(distribution == "pareto", param_name2, param_name1),
      param_name2_new = dplyr::if_else(distribution == "pareto", param_name1, param_name2)
    ) |>
    dplyr::transmute(distribution, parameter_1 = parameter_1_new,
                     parameter_2 = parameter_2_new, param_name1 = param_name1_new,
                     param_name2 = param_name2_new)

  # Create data for fitted CDFs

  data_range <- range(data)
  x_range <- seq(risk_shift, data_range[2], length.out = 1000)

  distribution_functions <- list(
    lnorm = function(x, params) plnorm(x - risk_shift, meanlog = params["meanlog"], sdlog = params["sdlog"]),
    pareto = function(x, params) ppareto(x - risk_shift, shape = params["shape"], scale = params["scale"]),
    weibull = function(x, params) pweibull(x - risk_shift, shape = params["shape"], scale = params["scale"]),
    gamma = function(x, params) pgamma2(x - risk_shift, shape = params["shape"], scale = params["scale"]),
    invgauss = function(x, params) pinvgauss(x - risk_shift, mean = params["mean"], shape = params["shape"])
  )

  cdf_data <- lapply(names(fits), function(dist) {
    params <- fits[[dist]]$estimate
    dplyr::tibble(
      x = x_range,
      y = distribution_functions[[dist]](x_range, params),
      distribution = dist
    )
  }) |>
    dplyr::bind_rows()

  # Create CDF plot

  cdf_plot <- ggplot2::ggplot() +
    ggplot2::stat_ecdf(data = data.frame(x = data), ggplot2::aes(x = x, color = "empirical"), linewidth = 1) +
    ggplot2::geom_line(data = cdf_data, ggplot2::aes(x = x, y = y, color = distribution), linewidth = 1, alpha = 0.6) +
    ggplot2::labs(title = "cumulative distribution function",
                  y = "cumulative probability") +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(xlim = c(risk_shift, NA)) +  # Set x-axis limits
    ggplot2::scale_color_manual(name = "distributions",
                                values = c("empirical" = "black", setNames(rainbow(length(fits)), names(fits)))) +
    ggplot2::theme(legend.position = "top",
                   legend.justification = "left",
                   legend.title = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::scale_y_continuous(labels = scales::percent)

  # create q-q plot

  quantile_functions <- list(
    lnorm = function(p, params) qlnorm(p, meanlog = params["meanlog"], sdlog = params["sdlog"]),
    pareto = function(p, params) actuar::qpareto(p, shape = params["shape"], scale = params["scale"]),
    weibull = function(p, params) qweibull(p, shape = params["shape"], scale = params["scale"]),
    gamma = function(p, params) qgamma2(p, shape = params["shape"], scale = params["scale"]),
    invgauss = function(p, params) qinvgauss(p, mean = params["mean"], shape = params["shape"])
  )

  qq_data <- lapply(names(fits), function(dist) {
    params <- fits[[dist]]$estimate
    theoretical_quantiles <- quantile_functions[[dist]](ppoints(length(shifted_data)), params)
    dplyr::tibble(
      theoretical = theoretical_quantiles,
      empirical = sort(shifted_data + risk_shift),
      distribution = dist
    )
  }) |>
    dplyr::bind_rows()

  qq_plot <- ggplot2::ggplot(qq_data, ggplot2::aes(x = theoretical, y = empirical, color = distribution)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_abline(intercept = risk_shift, slope = 1, linetype = "dashed") +
    ggplot2::labs(title = "Q-Q plot", x = "theoretical quantiles", y = "empirical quantiles") +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = setNames(rainbow(length(fits)), names(fits))) +
    ggplot2::theme(legend.position = "top",
                   legend.justification = "left",
                   legend.title = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::scale_y_continuous(labels = scales::comma)

  # create CDF value table

  percentiles <- c(0.01, 0.1, 0.25, 0.35, 0.50, 0.75, 0.85, 0.90, 0.95, 0.975, 0.99)
  cdf_table <- dplyr::tibble(
    percentile = percentiles,
    empirical = quantile(data, probs = percentiles)
  )

  for (dist in names(fits)) {
    params <- fits[[dist]]$estimate
    cdf_table[[dist]] <- sapply(percentiles, function(p) quantile_functions[[dist]](p, params) + risk_shift)
  }

  # format the CDF value table

  cdf_table <- cdf_table |>
    dplyr::mutate(
      percentile = scales::percent(percentile, accuracy = 0.1),
      dplyr::across(c(dplyr::everything(), -percentile), ~ format(round(.), big.mark = ",", scientific = FALSE))
    )

  # calculate summary statistics table

  summary_stats <- function(data) {
    dplyr::tibble(
      mean = mean(data),
      sd = sd(data),
      min = min(data),
      max = max(data)
    )
  }

  empirical_stats <- summary_stats(data)

  distribution_stats <- purrr::map_dfr(names(fits), ~ {
    dist <- .x
    params <- fits[[dist]]$estimate
    simulated_data <- switch(
      dist,
      lnorm = rlnorm(100000, meanlog = params["meanlog"], sdlog = params["sdlog"]) + risk_shift,
      pareto = actuar::rpareto(100000, shape = params["shape"], scale = params["scale"]) + risk_shift,
      weibull = rweibull(100000, shape = params["shape"], scale = params["scale"]) + risk_shift,
      gamma = rgamma2(100000, shape = params["shape"], scale = params["scale"]) + risk_shift,
      invgauss = rinvgauss(100000, mean = params["mean"], shape = params["shape"]) + risk_shift
    )
    summary_stats(simulated_data) |>  dplyr::mutate(distribution = dist)
  })

  summary_table <- dplyr::bind_rows(
    empirical_stats |> dplyr::mutate(distribution = "empirical"),
    distribution_stats
  ) |>
    tidyr::pivot_longer(cols = -distribution, names_to = "metric", values_to = "value") |>
    tidyr::pivot_wider(names_from = distribution, values_from = value) |>
    dplyr::mutate(dplyr::across(c(dplyr::everything(), -metric), ~ format(round(.), big.mark = ",", scientific = FALSE)))

  # create summary statistics plot

  summary_plot <- summary_table |>
    tidyr::pivot_longer(cols = -metric, names_to = "distribution", values_to = "value") |>
    dplyr::mutate(value = as.numeric(gsub(",", "", value))) |>
    dplyr::mutate(metric = forcats::fct_relevel(metric, "mean", "sd", "min", "max")) |>
    ggplot2::ggplot(ggplot2::aes(distribution, value, fill = distribution)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ metric, scales = "free") +
    ggplot2::scale_fill_manual(values = c("empirical" = "black", purrr::set_names(rainbow(length(fits)), names(fits)))) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "summary statistics", x = "distribution", y = "value") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank()) +
    ggplot2::coord_flip()


  # return results

  return(list(
    best_fit = best_fit,
    aic_table = aic_tibble,
    fit_parameters = fit_parameters,
    cdf_plot = cdf_plot,
    qq_plot = qq_plot,
    percentiles_table = cdf_table,
    summary_statistics = summary_table,
    summary_plot = summary_plot
  ))
}



