#' fit chain ladder development pattern
#'
#' requires tidy data. one row per cell in the unprojected triangle. \cr
#' missing rows within the triangle will be filled in with the value in the previous period. \cr
#' no rows for development periods beyond the last evaluation point. \cr
#' data provided must be cumulative. \cr
#' use the plot_cl_pattern function to view fitted pattern.
#'
#' @param cohort_var the variable in the input data representing the cohort
#' @param dev_var the variable in the input data representing the development period. needs to be consistent with dev_period_length and dev_period_units. e.g. if have half annual dev periods and data is in months (6, 12, 18 etc) then dev_period_length = 6 and dev_period_units = 1. if data in years (0.5, 1, 1.5 etc) then dev_period_length = 0.5 and dev_period_units = 12.
#' @param weighting_var the variable in the input data representing the amounts in the triangle e.g. claim counts or claim amounts etc
#' @param data the cumulative tidy dataset with the triangle data. one row per cohort year and development period. no rows for development periods beyond the last evaluation point.
#' @param dev_period_length the difference between subsequent development periods in the units specified in dev_period_units
#' @param dev_period_units the units in which dev_var is provided e.g. 1 for months, 12 for years etc.
#' @param exclude_last_diag defaults to FALSE. set to TRUE to exclude last diagonal
#' @param num_periods the number of prior periods to include in the fit. defaults to the whole triangle but set to an integer to include the last n periods instead. if last diagonal is excluded it takes the n periods prior to these. if ratios are excluded due to any other parameter settings e.g. exclude_high then the number of periods is not adjusted to allow for this e.g. if 4 most recent rows are excluded via exclude_points and num_periods is set to 5 there will only be 1 row included in the ata factor calculation.
#' @param selected_curve which curve to use. weibull, exponential_decay or inverse_power. defaults to weibull. only used if smooth_from is provided.
#' @param smooth_from optional development period to smooth from. if selected a curve is fitted to the CL ata factors and ata factors are selected from this curve beyond the development period specified. the tail factor will be estimated from this fitted curve.
#' @param future_dev_periods how many future development years to use when calculating the tail factor with the fitted curve. defaults to 5 years but should be set to a value such that the ata factors in the tail don't look to be truncated. only has an impact on the atu if smooth_from is provided. set to zero if you want to use a fitted curve but don't want it to estimate the tail.
#' @param exclude_points optional exclude individual ratios. provide as a list of indices into the ratio_table. column index should start after the cohort year column.
#' @param cohort_start all cohorts greater than or equal to this are included in the fit. defaults to the whole triangle if not specified.
#' @param tail_factor tail factor to apply. selected ata factor is set to this value at the final development period in the original data before future_dev_periods are added. hence be wary of using this in addition to using fitted curves to estimate the tail.
#' @param exclude_high set to 1 to exclude the highest ata ratio in each development period.
#' @param exclude_low set to 1 to exclude the lowest ata ratio in each development period.
#' @param bf_prior optional bf prior to use when calculating ultimate. can't be supplied if cc_decay_factor is provided.
#' @param cc_decay_factor optional cape cod decay factor to use when calculating ultimate. can't be supplied if bf_prior is provided.
#' @param premium required if bf_prior or cc_decay_factor is provided. dataframe with with a column called premium and another column which must have the same name as cohort_var.
#'
#' @return a list with the fitted pattern, the original input data, a triangle of individual ratios, the r_squared for the curve fits, a monthly pattern, the development units and the projected ultimates. monthly pattern is fitted using a piecewise cubic Hermite interpolating polynomial.
#'
#' @examples
#' # fit chain ladder to example triangle data
#' fit_cl_pattern(uw_year, dev_year, claim_number, triangle_data)
#'
#' # exclude uw years prior to 2013, exclude last diagonal, use a weibull fit to smooth from dev period 4 onewards
#' fit <- fit_cl_pattern(uw_year,
#' dev_year,
#' claim_number,
#' triangle_data,
#' exclude_last_diag = TRUE,
#' cohort_start = 2013,
#' smooth_from = 4)
#'
#' # view the fitted pattern
#' plot_cl_pattern(fit)
#'
#' # view the table of individual ratios
#' fit$ratio_table
#'
#' # ratios from dev 1 to dev 2 in 2015 and dev 2 to dev 3 in 2009 look high. these are row 13, column 1 and row 7, column 2 in the ratio_table
#' # note that column indexing starts from dev_1 onwards (not uw_year)
#' # refit excluding these points
#' fit <- fit_cl_pattern(uw_year,
#' dev_year,
#' claim_number,
#' triangle_data,
#' exclude_points = list(c(13, 1), c(7, 2)))
#'
#' @export


fit_cl_pattern <- function(cohort_var, dev_var, weighting_var, data,
                           dev_period_length = 1, dev_period_units = 12, exclude_last_diag = FALSE,
                           num_periods = NULL, selected_curve = "weibull", smooth_from = NULL, future_dev_periods = 5,
                           exclude_points = NULL, cohort_start = NULL, tail_factor = 1, exclude_high = FALSE, exclude_low = FALSE,
                           bf_prior = NULL, cc_decay_factor = NULL, premium = NULL) {

  cohort_char <- deparse(substitute(cohort_var))

  # warn if using curve to select tail factor and also provided a manual tail factor
  if(future_dev_periods > 0 & !is.null(smooth_from) & tail_factor != 1) warning("you've provided a manual tail factor but you're also using a curve fit to estimate the tail. inspect the resulting ata factors carefully and check that you know what you're doing! it's very likely you don't want to supply a tail factor if you're already estimating the tail with a curve.")

  # only provide one of bf_prior or cc_decay_factor
  if(!is.null(bf_prior) & !is.null(cc_decay_factor)) stop("choose one of BF or Cape Cod. not possible to supply both a BF prior and a CC decay factor.")

  # if bf prior or cc decay factor provided then premium needs to be provided
  # and premium dataframe needs to contain premium for all values of cohort_var
  # and columns needs to be as specified in the documentation
  if(!is.null(bf_prior) | !is.null(cc_decay_factor)) {
    if(is.null(premium)) stop("need to supply an premium if BF or CC method is used.")
    cohort_var_string <- rlang::as_name(enquo(cohort_var))
    if(!(cohort_var_string %in% colnames(premium))) stop("premium needs to contain the same cohort name as in the triangle data provided")
    data_cohort_values <- data %>% dplyr::pull(!!enquo(cohort_var)) %>% unique()
    exposure_base_cohort_values <- premium %>% dplyr::pull(!!enquo(cohort_var)) %>% unique()

    missing_cohorts <- setdiff(data_cohort_values, exposure_base_cohort_values)

    if(length(missing_cohorts) > 0) {
      stop(paste("The following cohort values are present in the triangle data but missing from exposure_base:",
                 paste(missing_cohorts, collapse = ", ")))
    }

    if(!("premium" %in% colnames(premium))) {
      stop("premium dataframe must contain a column named 'premium'")
    }

  }

  # check that selected_curve is one of the permitted values
  if(!(selected_curve %in% c("weibull", "exponential_decay", "inverse_power"))) stop("selected_curve must be weibull, exponential_decay or inverse_power")

  # if weighting_var contains lots of NA values warn the user as this may indicate that data contains developent periods beyond the last evaluation point

  propn_na <- data |>
    dplyr::ungroup() |>
    dplyr::summarise(na_prop = mean(is.na({{ weighting_var }}))) |>
    dplyr::pull(na_prop)

  if(propn_na > 0.05) warning(paste0(propn_na, " of the claim value rows in the triangle data contain missing values. check that you haven't supplied development periods beyond the last evaluation point."))

  # check that dataframe is complete for all combinations of cohort var and development var

  possible_errors2 <- data |>
    dplyr::group_by({{ cohort_var }}) |>
    dplyr::arrange({{ dev_var }}, .by_group = TRUE) |>
    dplyr::mutate(error = {{ dev_var }} - dplyr::lag({{ dev_var }})) |>
    dplyr::ungroup() |>
    dplyr::summarise(sd_error = sd(error, na.rm = TRUE)) |>
    dplyr::mutate(sd_error = tidyr::replace_na(sd_error, 0))

  if (possible_errors2$sd_error != 0) warning("triangle data doesn't contain one row for each cohort and development period. missing rows have been filled in and can be identified with the imputed_row column.")

  # check that exclude_points is a list with pairs of integers (if it's provided)

  if (!missing(exclude_points)) {

    if (is.list(exclude_points)) {

      all_pairs <- all(sapply(exclude_points, function(x) is.numeric(x) && length(x) == 2))

      if (all_pairs == FALSE) stop("exclude_points contains elements that are not pairs of integers")

    } else {
      print("exclude_points needs to be a list")
    }


  }

  # calculate number of additional rows required for future_dev_periods

  additional_rows <- future_dev_periods / (dev_period_units / 12)

  # ungroup data

  data <- ungroup(data)

  # fill in missing data items

  min_dev <- data |>
    dplyr::summarise(min_dev = min({{ dev_var }})) |>
    dplyr::pull()

  data <- data |>
    dplyr::mutate(imputed_row = 0) |>
    dplyr::group_by({{ cohort_var }}) |>
    tidyr::complete({{ dev_var }} := seq(min_dev, max({{ dev_var }}), by = dev_period_length)) |>
    tidyr::fill({{ weighting_var }}) |>
    tidyr::replace_na(list(imputed_row = 1)) |>
    dplyr::mutate(across({{ weighting_var }}, \(x) tidyr::replace_na(x, 0))) |>
    dplyr::ungroup()


  # create triangle of ratios

  ratio_table <- data |>
    dplyr::group_by({{ cohort_var }}) |>
    dplyr::arrange({{ dev_var }}, .by_group = TRUE) |>
    dplyr::mutate(ratio = dplyr::lead({{ weighting_var }}) / {{ weighting_var }}) |>
    dplyr::select({{ cohort_var }}, {{ dev_var }}, ratio) |>
    tidyr::pivot_wider(names_from = {{ dev_var }}, values_from = ratio, names_prefix = "dev_")

  # flag points to exclude if provided in exclude_points

  data$ex_numerator <- 0
  data$ex_denominator <- 0

  if(!missing(exclude_points)){

    for (point in exclude_points) {

      # find rows in the triangle data that correspond to the index into the ratio_table provided in exclude_points
      row_index <- point[1]
      col_index <- point[2]
      cohort_match <- ratio_table[[cohort_char]][row_index]
      dev_period_match <- as.numeric(sub("dev_", "", names(ratio_table)[col_index + 1]))

      data <- data |>
        dplyr::mutate(ex_denominator = dplyr::if_else({{ cohort_var }} == cohort_match & {{ dev_var }} == dev_period_match, 1, ex_denominator)) |>
        dplyr::mutate(ex_numerator = dplyr::if_else({{ cohort_var }} == cohort_match & {{ dev_var }} == dev_period_match + dev_period_units, 1, ex_numerator))

    }

  }

  # need the magrittr pipe operator for the in-pipe left join below

  `%>%` <- dplyr::`%>%`


  # fit chain ladder pattern

  cl_pattern <- data |>
    # flag how many rows to include in numerator and denominator at each dev period
    dplyr::group_by({{ dev_var }}) |>
    dplyr::arrange({{ cohort_var }}, .by_group = TRUE) |>
    dplyr::mutate(num = dplyr::n()) |>
    dplyr::ungroup() %>%
    dplyr::left_join(
      {.} |>
        dplyr::distinct({{ dev_var }}, num) |>
        dplyr::mutate(dev_temp = {{ dev_var }} - dev_period_length) |>
        dplyr::rename(num_next = num) |>
        dplyr::select(-{{ dev_var }}),
      dplyr::join_by({{ dev_var }} == dev_temp)
    ) |>
    dplyr::group_by({{ dev_var }}) |>
    dplyr::arrange({{ cohort_var }}, .by_group = TRUE) |>
    # flag for which rows to include in the ratio calculation based on exclude_last_diag flag and indicators created in first step
    dplyr::mutate(denominator_include = dplyr::if_else(dplyr::row_number() <= num_next - exclude_last_diag, 1, 0),
                  numerator_include = dplyr::if_else(dplyr::row_number() <= num - exclude_last_diag, 1, 0)) |>
    # adjust for cohort start
    dplyr::mutate(numerator_include = dplyr::if_else({{ cohort_var }} < dplyr::coalesce({{ cohort_start }}, 0),
                                                     0, numerator_include),
                  denominator_include = dplyr::if_else({{ cohort_var }} < dplyr::coalesce({{ cohort_start }}, 0),
                                                       0, denominator_include)) |>
    # calculate ata ratios and flag high and low in each dev period
    dplyr::group_by({{ cohort_var }}) |>
    dplyr::arrange({{ dev_var }}, .by_group = TRUE) |>
    dplyr::mutate(ata_ratio = dplyr::lead({{ weighting_var }}) / {{ weighting_var }},
                  ata_ratio = dplyr::if_else(is.infinite(ata_ratio) | dplyr::lead({{ weighting_var }}) == 0, NA_real_, ata_ratio)) |>
    # flag high / low in each dev year
    dplyr::group_by({{ dev_var }}) |>
    dplyr::mutate(
      any_valid_ratio = any(!is.na(ata_ratio)),
      max_ratio = dplyr::if_else(any_valid_ratio, suppressWarnings(max(ata_ratio, na.rm = TRUE)), NA_real_),
      min_ratio = dplyr::if_else(any_valid_ratio, suppressWarnings(min(ata_ratio, na.rm = TRUE)), NA_real_),
      high_flag = dplyr::if_else(any_valid_ratio & !is.na(ata_ratio) & ata_ratio == max_ratio, 1L, 0L),
      low_flag = dplyr::if_else(any_valid_ratio & !is.na(ata_ratio) & ata_ratio == min_ratio, 1L, 0L)
    ) |>
    dplyr::select(-max_ratio, -min_ratio, -any_valid_ratio) |>
    # only set high / low flag for dev period columns with 3 or more included factors
    dplyr::mutate(
      valid_group = sum(!is.na(ata_ratio)) >= 3,
      high_flag = dplyr::if_else(valid_group, high_flag, 0L),
      low_flag = dplyr::if_else(valid_group, low_flag, 0L)
    ) |>
    # if more than one high / low flag then take the latest uw_year
    dplyr::group_by({{ dev_var }}) |>
    dplyr::mutate(
      max_high_row = max(dplyr::if_else(high_flag == 1, dplyr::row_number(), 0L)),
      max_low_row = max(dplyr::if_else(low_flag == 1, dplyr::row_number(), 0L)),
      high_flag = dplyr::if_else(dplyr::row_number() == max_high_row & max_high_row > 0, 1L, 0L),
      low_flag = dplyr::if_else(dplyr::row_number() == max_low_row & max_low_row > 0, 1L, 0L)
    ) |>
    dplyr::select(-max_high_row, -max_low_row) |>
    dplyr::ungroup() |>
    dplyr::select(-ata_ratio, -valid_group) |>
    # if exclude_high_low then adjust denominator / numerator include flags
    dplyr::group_by({{ cohort_var }}) |>
    dplyr::arrange({{ dev_var }}, .by_group = TRUE) |>
    dplyr::mutate(denominator_include = case_when(
      exclude_high == TRUE & high_flag == 1 ~ 0,
      exclude_low == TRUE & low_flag == 1 ~ 0,
      TRUE ~ denominator_include
    )) |>
    dplyr::mutate(numerator_include = case_when(
      exclude_high == TRUE & dplyr::lag(high_flag, default = 0) == 1 ~ 0,
      exclude_low == TRUE & dplyr::lag(low_flag, default = 0) == 1 ~ 0,
      TRUE ~ numerator_include
    )) |>
    # adjust flag to allow for exclude_point
    dplyr::mutate(denominator_include = dplyr::if_else(ex_denominator == 1, 0, denominator_include),
                  numerator_include = dplyr::if_else(ex_numerator == 1, 0, numerator_include)) |>
    # if preceding row is zero or missing then set numerator include to zero
    dplyr::group_by({{ cohort_var }}) |>
    dplyr::arrange({{ dev_var }}, .by_group = TRUE) |>
    dplyr::mutate(numerator_include = dplyr::if_else(lag({{ weighting_var }}, default = 0) == 0, 0, numerator_include),
                  denominator_include = dplyr::if_else({{ weighting_var }} == 0, 0, denominator_include)) |>
    # adjust flag to allow for num_periods
    dplyr::group_by({{ dev_var }}) |>
    dplyr::arrange(desc({{ cohort_var }}), .by_group = TRUE) |>
    dplyr::mutate(numerator_count = cumsum(numerator_include),
                  denominator_count = cumsum(denominator_include)) |>
    # if num_periods is NULL then set to a large number
    dplyr::mutate(num_periods = dplyr::coalesce(num_periods, 999)) |>
    dplyr::mutate(denominator_include = dplyr::if_else(denominator_count <= num_periods, denominator_include, 0),
                  numerator_include = dplyr::if_else(numerator_count <= num_periods, numerator_include, 0)) |>
    # sum up for total claims in numerator and denominator
    dplyr::summarise(denominator = sum({{ weighting_var }} * denominator_include, na.rm = TRUE),
                     numerator = sum({{ weighting_var }} * numerator_include, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    # calculate age to age ratio at each dev period
    # if NA then set to 1
    dplyr::mutate(numerator = dplyr::lead(numerator),
                  ata = numerator / denominator,
                  ata = tidyr::replace_na(ata, 1)) |>
    # allow for tail factor
    dplyr::mutate(ata = dplyr::if_else(dplyr::row_number() == dplyr::n(), {{ tail_factor }}, ata)) |>
    # if ata is zero then set to 1
    dplyr::mutate(ata = dplyr::if_else(ata == 0, 1, ata)) |>
    # expand for number of future dev periods
    dplyr::mutate(original_data_flag = 1) %>%
    dplyr::bind_rows(
      dplyr::tibble(!!enquo(dev_var) := {
        max_val <- max(dplyr::pull(., !!enquo(dev_var)))
        seq(
          from = max_val + dev_period_length,
          by = dev_period_length,
          length.out = additional_rows
        )
      },
      ata = 1,
      original_data_flag = 0)
    ) |>
    # calculate age to ultimate and percent developed
    dplyr::mutate(atu = rev(cumprod(rev(ata))),
                  pct_dev = 1 / atu) |>
    dplyr::select({{ dev_var }}, ata, atu, pct_dev, original_data_flag) |>
    dplyr::rename(dev_period = {{ dev_var }})

  # fit weibull curve

  weibull_fit <- cl_pattern |>
    dplyr::filter(original_data_flag == 1) |>
    dplyr::mutate(dev_period = dplyr::row_number()) |>
    dplyr::mutate(ata_fit = dplyr::if_else(ata <= 1, NA_real_, ata),
                  y_t = log(log(ata_fit / (ata_fit - 1))),
                  x_t = if_else(ata <= 1, NA_real_, log(dev_period)),
                  ata_count = dplyr::if_else(is.na(ata_fit), 0, 1)) |>
    dplyr::summarise(x = sum(x_t, na.rm = TRUE),
                     y = sum(y_t, na.rm = TRUE),
                     xx = sum(x_t * x_t, na.rm = TRUE),
                     xy = sum(x_t * y_t, na.rm = TRUE),
                     yy = sum(y_t * y_t, na.rm = TRUE),
                     n = sum(ata_count)) |>
    dplyr::mutate(sxx = xx - ((x ^ 2) / n),
                  sxy = xy - ((x * y) / n),
                  syy = yy - ((y ^ 2) / n),
                  b = sxy / sxx,
                  a = exp((y / n) - (b * (x / n))),
                  weibull_r_squared = (sxy / (sqrt(syy * sxx))) ^ 2) |>
    dplyr::transmute(weibull_a = a, weibull_b = b, weibull_r_squared)

  # fit exponential decay

  exponential_fit <- cl_pattern |>
    dplyr::filter(original_data_flag == 1) |>
    dplyr::mutate(dev_period = dplyr::row_number()) |>
    dplyr::mutate(ata_fit = dplyr::if_else(ata <= 1, NA_real_, ata),
                  y_t = log(ata_fit - 1),
                  x_t = if_else(ata <= 1, NA_real_, (dev_period)),
                  ata_count = dplyr::if_else(is.na(ata_fit), 0, 1)) |>
    dplyr::summarise(x = sum(x_t, na.rm = TRUE),
                     y = sum(y_t, na.rm = TRUE),
                     xx = sum(x_t * x_t, na.rm = TRUE),
                     xy = sum(x_t * y_t, na.rm = TRUE),
                     yy = sum(y_t * y_t, na.rm = TRUE),
                     n = sum(ata_count)) |>
    dplyr::mutate(sxx = xx - ((x ^ 2) / n),
                  sxy = xy - ((x * y) / n),
                  syy = yy - ((y ^ 2) / n),
                  b = sxy / sxx,
                  a = exp((y / n) - (b * (x / n))),
                  exponential_r_squared = (sxy / (sqrt(syy * sxx))) ^ 2) |>
    dplyr::transmute(exponential_a = a, exponential_b = b, exponential_r_squared)


  # fit inverse power
  # fit c values -0.5, 0, 1, 3, 5 and select fit with highest r squared

  inv_power_fit <- cl_pattern |>
    dplyr::filter(original_data_flag == 1) |>
    dplyr::mutate(dev_period = dplyr::row_number()) |>
    dplyr::mutate(ata_fit = dplyr::if_else(ata <= 1, NA_real_, ata),
                  y_t = log(ata_fit - 1),
                  x_t_1 = if_else(ata <= 1, NA_real_, log(dev_period - 0.5)),
                  x_t_2 = if_else(ata <= 1, NA_real_, log(dev_period + 0)),
                  x_t_3 = if_else(ata <= 1, NA_real_, log(dev_period + 1)),
                  x_t_4 = if_else(ata <= 1, NA_real_, log(dev_period + 3)),
                  x_t_5 = if_else(ata <= 1, NA_real_, log(dev_period + 5)),
                  ata_count = dplyr::if_else(is.na(ata_fit), 0, 1)) |>
    dplyr::summarise(x_1 = sum(x_t_1, na.rm = TRUE),
                     x_2 = sum(x_t_2, na.rm = TRUE),
                     x_3 = sum(x_t_3, na.rm = TRUE),
                     x_4 = sum(x_t_4, na.rm = TRUE),
                     x_5 = sum(x_t_5, na.rm = TRUE),
                     y = sum(y_t, na.rm = TRUE),
                     xx_1 = sum(x_t_1 * x_t_1, na.rm = TRUE),
                     xx_2 = sum(x_t_2 * x_t_2, na.rm = TRUE),
                     xx_3 = sum(x_t_3 * x_t_3, na.rm = TRUE),
                     xx_4 = sum(x_t_4 * x_t_4, na.rm = TRUE),
                     xx_5 = sum(x_t_5 * x_t_5, na.rm = TRUE),
                     xy_1 = sum(x_t_1 * y_t, na.rm = TRUE),
                     xy_2 = sum(x_t_2 * y_t, na.rm = TRUE),
                     xy_3 = sum(x_t_3 * y_t, na.rm = TRUE),
                     xy_4 = sum(x_t_4 * y_t, na.rm = TRUE),
                     xy_5 = sum(x_t_5 * y_t, na.rm = TRUE),
                     yy = sum(y_t * y_t, na.rm = TRUE),
                     n = sum(ata_count)) |>
    dplyr::mutate(sxx_1 = xx_1 - ((x_1 ^ 2) / n),
                  sxx_2 = xx_2 - ((x_2 ^ 2) / n),
                  sxx_3 = xx_3 - ((x_3 ^ 2) / n),
                  sxx_4 = xx_4 - ((x_4 ^ 2) / n),
                  sxx_5 = xx_5 - ((x_5 ^ 2) / n),
                  sxy_1 = xy_1 - ((x_1 * y) / n),
                  sxy_2 = xy_2 - ((x_2 * y) / n),
                  sxy_3 = xy_3 - ((x_3 * y) / n),
                  sxy_4 = xy_4 - ((x_4 * y) / n),
                  sxy_5 = xy_5 - ((x_5 * y) / n),
                  syy = yy - ((y ^ 2) / n),
                  r_squared_1 = (sxy_1 / (sqrt(syy * sxx_1))) ^ 2,
                  r_squared_2 = (sxy_2 / (sqrt(syy * sxx_2))) ^ 2,
                  r_squared_3 = (sxy_3 / (sqrt(syy * sxx_3))) ^ 2,
                  r_squared_4 = (sxy_4 / (sqrt(syy * sxx_4))) ^ 2,
                  r_squared_5 = (sxy_5 / (sqrt(syy * sxx_5))) ^ 2,
                  max_r_squared = pmax(r_squared_1, r_squared_2, r_squared_3, r_squared_4, r_squared_5),
                  b_1 = sxy_1 / sxx_1,
                  b_2 = sxy_2 / sxx_2,
                  b_3 = sxy_3 / sxx_3,
                  b_4 = sxy_4 / sxx_4,
                  b_5 = sxy_5 / sxx_5,
                  a_1 = exp((y / n) - (b_1 * (x_1 / n))),
                  a_2 = exp((y / n) - (b_2 * (x_2 / n))),
                  a_3 = exp((y / n) - (b_3 * (x_3 / n))),
                  a_4 = exp((y / n) - (b_4 * (x_4 / n))),
                  a_5 = exp((y / n) - (b_5 * (x_5 / n))),
                  a = dplyr::case_when(
                    r_squared_1 == max_r_squared ~ a_1,
                    r_squared_2 == max_r_squared ~ a_2,
                    r_squared_3 == max_r_squared ~ a_3,
                    r_squared_4 == max_r_squared ~ a_4,
                    r_squared_5 == max_r_squared ~ a_5
                  ),
                  b = dplyr::case_when(
                    r_squared_1 == max_r_squared ~ b_1,
                    r_squared_2 == max_r_squared ~ b_2,
                    r_squared_3 == max_r_squared ~ b_3,
                    r_squared_4 == max_r_squared ~ b_4,
                    r_squared_5 == max_r_squared ~ b_5
                  ),
                  inv_power_r_squared = max_r_squared,
                  selected_c = dplyr::case_when(
                    r_squared_1 == max_r_squared ~ -0.5,
                    r_squared_2 == max_r_squared ~ 0,
                    r_squared_3 == max_r_squared ~ 1,
                    r_squared_4 == max_r_squared ~ 3,
                    r_squared_5 == max_r_squared ~ 5,
                  )) |>
    dplyr::transmute(inv_power_a = a, inv_power_b = b, inv_power_r_squared, selected_c)

  # compile all R squared

  r_squared_value <- dplyr::tibble(
    curve = c("weibull", "exponential", "inverse_power"),
    r_squared = c(weibull_fit$weibull_r_squared,
                  exponential_fit$exponential_r_squared,
                  inv_power_fit$inv_power_r_squared)
  )

  # add curves to cl pattern

  cl_pattern <- cl_pattern |>
    dplyr::mutate(weibull_ata = 1 / (1 - exp(-weibull_fit$weibull_a * (dplyr::row_number() ^ weibull_fit$weibull_b))),
                  exponential_ata = 1 + exponential_fit$exponential_a * exp(exponential_fit$exponential_b * dplyr::row_number()),
                  inv_power_ata = 1 + inv_power_fit$inv_power_a * (inv_power_fit$selected_c + dplyr::row_number()) ^ inv_power_fit$inv_power_b,
                  selected_curve_ata = dplyr::case_when(
                    selected_curve == "weibull" ~ weibull_ata,
                    selected_curve == "exponential_decay" ~ exponential_ata,
                    selected_curve == "inverse_power" ~ inv_power_ata,
                    TRUE ~ weibull_ata
                  ),
                  selected_ata = dplyr::if_else(dev_period >= dplyr::coalesce(smooth_from, 999),
                                                selected_curve_ata, ata)) |>
    # if curve(s) have failed to fit then set selected ata to CL ATA
    dplyr::mutate(selected_ata = if_else(is.nan(selected_ata), ata, selected_ata)) |>
    dplyr::mutate(selected_ata = if_else(is.na(selected_ata), ata, selected_ata)) |>
    # if a tail factor has been provided use this at the last development period
    dplyr::mutate(selected_ata = if_else(tail_factor != 1 & row_number() == n(), ata, selected_ata)) |>
    dplyr::mutate(atu = rev(cumprod(rev(selected_ata))),
                  pct_dev = 1 / atu) |>
    dplyr::select(dev_period, ata, weibull_ata, exponential_ata, inv_power_ata, selected_ata, atu, pct_dev)


  # calculate ultimates

  ultimates <- data |>
    dplyr::group_by({{ cohort_var }}) |>
    dplyr::filter({{ dev_var }} == max({{ dev_var }})) |>
    dplyr::select({{ cohort_var }}, {{ dev_var }}, {{ weighting_var }}) |>
    dplyr::left_join(cl_pattern |>
                       dplyr::select(dev_period, pct_dev), dplyr::join_by({{ dev_var }} == dev_period))

  if(!is.null(bf_prior)) {

    ultimates <- ultimates |>
      dplyr::ungroup() |>
      dplyr::left_join(premium, join_by({{ cohort_var }})) |>
      dplyr::mutate(ultimate = {{ weighting_var }} + (1 - pct_dev) * premium * bf_prior)


  } else if(!is.null(cc_decay_factor)) {

    ultimates <- ultimates |>
      dplyr::ungroup() |>
      dplyr::mutate(row_num = row_number()) |>
      dplyr::left_join(premium, join_by({{ cohort_var }}))

    years <- unique(dplyr::pull(ultimates, row_num))

    ultimates <- purrr::map_dfr(cc_decay_factor, function(cc_decay_factor) {
      ultimates |>
        dplyr::mutate(lambda = purrr::map(row_num,
                                          ~ cc_decay_factor ^ base::abs(years - .x))) |>
        dplyr::mutate(lr_numerator = purrr::map_dbl(lambda, ~ base::sum({{ weighting_var }} * .x)),
                      lr_denominator = purrr::map_dbl(lambda, ~ base::sum(premium * pct_dev * .x))) |>
        dplyr::mutate(cc_lr = lr_numerator / lr_denominator,
                      cc_decay_factor = cc_decay_factor)
    }) |>
      dplyr::select(-lambda, -lr_numerator, -lr_denominator) |>
      dplyr::mutate(ultimate = {{ weighting_var }} + (1 - pct_dev) * premium * cc_lr) |>
      dplyr::select(-row_num)

  } else {

    ultimates <- ultimates |>
      dplyr::mutate(ultimate = {{ weighting_var }} / pct_dev) |>
      dplyr::ungroup()

  }

  ultimates <- ultimates |>
    dplyr::rename(cohort = 1, dev_period = 2)

  # get monthly pattern

  monthly_data <- cl_pattern |>
    dplyr::mutate(dev_month = dev_period * {{ dev_period_units }}) |>
    dplyr::select(dev_month, pct_dev) |>
    dplyr::bind_rows(dplyr::tibble(dev_month = 0, pct_dev = 0)) |>
    dplyr::arrange(dev_month) |>
    dplyr::filter(!is.nan(pct_dev))

  interpolated_values <- pracma::pchip(monthly_data$dev_month, monthly_data$pct_dev, 1:max(monthly_data$dev_month))

  monthly_pattern <- dplyr::tibble(
    dev_month = 1:max(monthly_data$dev_month),
    pct_dev = interpolated_values
  )


  return(list(cl_pattern = cl_pattern, data = data |>
                dplyr::rename(cohort = {{ cohort_var }},
                              dev_period = {{ dev_var }},
                              weighting = {{ weighting_var }}),
              ratio_table = ratio_table,
              monthly_pattern = monthly_pattern,
              r_squared = r_squared_value,
              dev_period_units = dev_period_units,
              ultimates = ultimates))

}

