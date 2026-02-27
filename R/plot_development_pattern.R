#' plot a develoment pattern fitted using the fit_development_pattern function
#'
#' provide the whole object returned by the fit_development_pattern not just the development_data item
#'
#' @param development_data the object returned by the fit_development_pattern function
#' @param show_legend defaults to FALSE. show the legend for the colour variable in the plot
#' @param include_monthly_fit include a separate graph showing the fitted monthly pattern against the chain ladder pattern. defaults to FALSE.
#' @param include_fitted_curve defaults to FALSE. show the Weibull curve fit on the graph
#' @return ggplot2 plot showing the fitted development pattern. if include_monthly_fit set to TRUE then returns a list containing two plots.
#'
#' @examples
#' \donttest{
#'
#' # fit chain ladder to example triangle data
#' fit <- fit_development_pattern(uw_year,
#'                                dev_year,
#'                                claim_number,
#'                                triangle_data)
#'
#' # view the fitted pattern
#' plot_development_pattern(fit)
#'
#' # view the fitted pattern and include the fitted monthly pattern in a
#' # separate graph
#' plot_development_pattern(fit, include_monthly_fit = TRUE)
#'
#' # view the fitted pattern and include the Weibull curve fit
#' plot_development_pattern(fit, include_fitted_curve = TRUE)
#' }
#'
#' @export


plot_development_pattern <- function(development_data, show_legend = FALSE, include_monthly_fit = FALSE,
                                     include_fitted_curve = FALSE) {

  data <- development_data$data
  development_pattern <- development_data$development_pattern
  ultimate <- development_data$ultimates

  # create plot data
  plot <- data |>
    dplyr::left_join(ultimate |>
                       dplyr::select(cohort, ultimate),
                     dplyr::join_by(cohort)) |>
    dplyr::mutate(pct_dev = weighting / ultimate) |>
    dplyr::mutate(cohort = forcats::as_factor(cohort))

  # expand to include an entry at 0,0
  plot <- plot |>
    dplyr::group_by(cohort) |>
    dplyr::slice(1) |>
    dplyr::mutate(dev_period = 0, weighting = 0, imputed_row = 1,
                  ex_numerator = 0, ex_denominator = 0, pct_dev = 0) |>
    dplyr::ungroup() |>
    dplyr::bind_rows(plot) |>
    dplyr::arrange(cohort, dev_period)

  development_pattern <- development_pattern |>
    dplyr::add_row(dev_period = 0, pct_dev = 0, .before = 1)

  # get number of colours in plot
  colours_required <- length(levels(plot$cohort))

  # set colour palette
  colour_palette <- grDevices::colorRampPalette(c("yellow", "red"))(colours_required)

  # prepare horizontal axis title
  x_title <- if (development_data$dev_period_units == 12) {
    "development year"
  } else if (development_data$dev_period_units == 3) {
    "development quarter"
  } else if (development_data$dev_period_units == 1) {
    "development month"
  } else if (development_data$dev_period_units == 6) {
    "development half-year"
  } else {
    "development period"
  }


  # create initial plot

  plot <- plot |>
    ggplot2::ggplot(ggplot2::aes(dev_period, pct_dev, colour = cohort)) +
    ggplot2::geom_line(data = development_pattern, colour = "black", linewidth = 2, linetype = "dashed") +
    ggplot2::scale_y_continuous(labels = scales::percent, name = "percent developed") +
    ggplot2::scale_x_continuous(name = x_title)


  # add lines

  if(show_legend == FALSE) {

    plot <- plot + ggplot2::geom_line(show.legend = FALSE, linewidth = 1, alpha = 0.4) +
      ggplot2::scale_color_manual(values = colour_palette) +
      ggplot2::labs(title = "fitted pattern  \n<span style = 'color: #FFD700;'>low cohort</span> to <span style = 'color: red;'>high cohort") +
      ggplot2::theme_light() +
      ggplot2::theme(
        plot.title = ggtext::element_markdown(),
        legend.position = "top",
        legend.justification = "left",
        panel.grid.minor = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank()
      )

    if (include_fitted_curve == TRUE) {

      fitted_curves <- development_data$development_pattern |>
        dplyr::mutate(weibull_atu = rev(cumprod(rev(weibull_ata))),
                      weibull_pct_dev = 1 / weibull_atu,
                      exponential_atu = rev(cumprod(rev(exponential_ata))),
                      exponential_pct_dev = 1 / exponential_atu,
                      inv_power_atu = rev(cumprod(rev(inv_power_ata))),
                      inv_power_pct_dev = 1 / inv_power_atu) |>
        # add 0,0
        dplyr::add_row(dev_period = 0, weibull_pct_dev = 0, exponential_pct_dev = 0,
                       inv_power_pct_dev = 0, .before = 1)

      weibull_curve <- fitted_curves |> dplyr::select(dev_period, weibull_pct_dev)
      exponential_curve <- fitted_curves |> dplyr::select(dev_period, exponential_pct_dev)
      inv_power_curve <- fitted_curves |> dplyr::select(dev_period, inv_power_pct_dev)

      plot <- plot +
        ggplot2::geom_line(data = weibull_curve, ggplot2::aes(y = weibull_pct_dev), colour = "cornflowerblue", linetype = "dashed",
                           linewidth = 1.5) +
        ggplot2::geom_line(data = exponential_curve, ggplot2::aes(y = exponential_pct_dev), colour = "mediumspringgreen", linetype = "dashed",
                           linewidth = 1.5) +
        ggplot2::geom_line(data = inv_power_curve, ggplot2::aes(y = inv_power_pct_dev), colour = "deeppink2", linetype = "dashed",
                           linewidth = 1.5) +
        ggplot2::labs(title = "fitted pattern  \n<span style = 'color: #FFD700;'>low cohort</span> to <span style = 'color: red;'>high cohort  \n<span style = 'color: cornflowerblue;'>weibull <span style = 'color: mediumspringgreen;'> exponential <span style = 'color: deeppink2;'> inverse power <span style = 'color: black;'>curve fits")

    }

  } else {

    plot <- plot + ggplot2::geom_line(show.legend = TRUE, linewidth = 1, alpha = 0.4) +
      ggplot2::scale_color_manual(values = colour_palette) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
      ggplot2::theme_light() +
      ggplot2::theme(
        plot.title = ggtext::element_markdown(),
        legend.position = "top",
        legend.justification = "left",
        panel.grid.minor = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::ggtitle("fitted pattern")

    if (include_fitted_curve == TRUE) {

      fitted_curves <- development_data$development_pattern |>
        dplyr::mutate(weibull_atu = rev(cumprod(rev(weibull_ata))),
                      weibull_pct_dev = 1 / weibull_atu,
                      exponential_atu = rev(cumprod(rev(exponential_ata))),
                      exponential_pct_dev = 1 / exponential_atu,
                      inv_power_atu = rev(cumprod(rev(inv_power_ata))),
                      inv_power_pct_dev = 1 / inv_power_atu) |>
        # add 0,0
        dplyr::add_row(dev_period = 0, weibull_pct_dev = 0, exponential_pct_dev = 0,
                       inv_power_pct_dev = 0, .before = 1)

      weibull_curve <- fitted_curves |> dplyr::select(dev_period, weibull_pct_dev)
      exponential_curve <- fitted_curves |> dplyr::select(dev_period, exponential_pct_dev)
      inv_power_curve <- fitted_curves |> dplyr::select(dev_period, inv_power_pct_dev)

      plot <- plot +
        ggplot2::geom_line(data = weibull_curve, ggplot2::aes(y = weibull_pct_dev), colour = "cornflowerblue", linetype = "dashed",
                           linewidth = 1.5) +
        ggplot2::geom_line(data = exponential_curve, ggplot2::aes(y = exponential_pct_dev), colour = "mediumspringgreen", linetype = "dashed",
                           linewidth = 1.5) +
        ggplot2::geom_line(data = inv_power_curve, ggplot2::aes(y = inv_power_pct_dev), colour = "deeppink2", linetype = "dashed",
                           linewidth = 1.5) +
        ggplot2::labs(title = "fitted pattern  \n<span style = 'color: #FFD700;'>low cohort</span> to <span style = 'color: red;'>high cohort  \n<span style = 'color: cornflowerblue;'>weibull <span style = 'color: mediumspringgreen;'> exponential <span style = 'color: deeppink2;'> inverse power <span style = 'color: black;'>curve fits")

    }

  }




  # monthly pattern plot

  plot2 <- development_data$development_pattern |>
    dplyr::mutate(dev_month = dev_period * development_data$dev_period_units) |>
    dplyr::add_row(dev_month = 0, pct_dev = 0, .before = 1) |>
    ggplot2::ggplot(ggplot2::aes(dev_month, pct_dev)) +
    ggplot2::geom_line(colour = "black") +
    ggplot2::geom_point(data = development_data$monthly_pattern,
                        colour = "red") +
    ggplot2::scale_y_continuous(labels = scales::percent, name = "percent developed") +
    ggplot2::scale_x_continuous(name = "development month") +
    ggplot2::labs(title = "<span style = 'color : red;'>interpolated monthly pattern</span>  \nchain ladder fitted pattern") +
    ggplot2::theme_light() +
    ggplot2::theme(plot.title = ggtext::element_markdown())

  if (include_monthly_fit == TRUE) {
    return(list(monthly_pattern = plot2, development_pattern = plot))
  } else {
    return(plot)
  }



}




