#' create dataset of monte carlo simulations
#'
#' create monte carlo simulations for a given frequency and severity distribution
#'
#' @param num_sims the number of simulated years required in the output
#' @param mean_freq the mean number of claims in each simulated year. unless var_freq is provided assumes poisson
#' @param sev_dist the severity distribution to use. one of lognorm2, gamma, pareto2, weibull or empirical
#' @param sev_param1 the first parameter for the severity distribution. not used if severity distribution is empirical. defaults to 1. for gamma and weibull this should be shape, for pareto2 it should be scale and for lognorm2 it should be meanlog.
#' @param sev_param2 the second parameter for the severity distribution. not used if severity distribution is empirical. defaults to 1. for gamma and weibull this should be scale, for pareto2 it should be shape and for lognorm2 it should be sdlog.
#' @param risk_shift the amount to shift the simulated losses. defaults to zero
#' @param var_freq the variance of the number of claims in each simulated year. if provided and different to mean freq will assume negative binomial. if provided and zero then the number of losses in each year will equal mean_freq
#' @param label optional text variable to include as a column titled "label" in the output. see examples below.
#' @param empirical if empirical distribution selected then this must be provided. dataset with two columns - probability and loss - representing the empirical distribution. probability must be cumulative probability and loss the associated loss amount.
#'
#' @return dataset of monte carlo simulations
#'
#' @examples
#' \donttest{
#'
#' # create 10k simulated years with poisson mean of 5 and
#' # lognormal parameters of 13 and 0.5
#' create_simulations(num_sims = 10e5,
#'                    mean_freq = 5,
#'                    sev_dist = "lognorm2",
#'                    sev_param1 = 13,
#'                    sev_param2 = 0.5)
#'
#' # create 10k simulated years with sames as above except for
#' # negative binomial distribution with mean 5 and variance 10
#' create_simulations(num_sims = 10e5,
#'                    mean_freq = 5,
#'                    sev_dist = "lognorm2",
#'                    sev_param1 = 13,
#'                    sev_param2 = 0.5,
#'                    var_freq = 10)
#'
#' # create 10k simulated years with two different distriubtions using pmap.
#' # use the label variable to differentiate them in the output
#' # note pmap doesn't work if an empirical severity distribution is being used
#' params <- dplyr::tibble(num_sims = 10e5,
#'                         mean_freq = c(5, 2),
#'                         sev_dist = c("lognorm2", "lognorm2"),
#'                         sev_param1 = c(13, 6),
#'                         sev_param2 = c(0.5, 0.7),
#'                         label = c("typeA", "typeB"))
#' purrr::pmap_dfr(params, create_simulations)
#'
#' # create 10k simulated years with poisson frequency and empirical severity
#' # distribution where the loss amounts are distributed evenly between 5-10m.
#' # in practice you would provide a full CDF stored in a separate dataframe.
#' create_simulations(num_sims = 10e5,
#'                    mean_freq = 5,
#'                    sev_dist = "empirical",
#'                    empirical = dplyr::tibble(probability = c(0, 1),
#'                                              loss = c(5e6, 10e6)))
#'
#' }
#'
#' @export

create_simulations <- function(num_sims, mean_freq, sev_dist, sev_param1 = 1, sev_param2 = 1, risk_shift = 0, var_freq = mean_freq, label = NULL,
                               empirical = NULL) {

  # if severity distribution provided isn't one of the options then throw an error
  if(!({{ sev_dist }} %in% c("lognorm2", "gamma", "pareto2", "weibull", "empirical"))) stop("severity distribution needs to be one of lognorm2, pareto2, gamma")

  # if empirical severity distribution selected but no empirical data then throw an error
  if({{ sev_dist }} == "empirical" & is.null({{ empirical }})) stop("empirical distribution selected - please supply an associated dataset")

  # if empirical dataset provided then column names need to be probability and loss
  if({{ sev_dist }} == "empirical" & !(("loss" %in% colnames(empirical)) & ("probability" %in% colnames(empirical)))) stop("empirical dataset needs to have columns titled loss and probability")

  if(var_freq == mean_freq) {

    # if var not specified then simulate required number of poisson variables
    losses <- dplyr::tibble(num = stats::rpois(num_sims, mean_freq)) |>
      # create year
      dplyr::mutate(year = 1:nrow(dplyr::pick(num))) |>
      # expand dataframe so we have one row for each loss
      dplyr::mutate(year = purrr::map2(year, num, function(x, y) rep(x, y))) |>
      tidyr::unnest(year) |>
      dplyr::select(year)

  } else {

    # if var is specified and is zero then simulate fixed n in each year
    if (var_freq == 0) {

      losses <- dplyr::tibble(num = rep(mean_freq, num_sims)) |>
        dplyr::mutate(year = 1:nrow(dplyr::pick(num))) |>
        dplyr::mutate(year = purrr::map2(year, num, function(x, y) rep(x, y))) |>
        tidyr::unnest(year) |>
        dplyr::select(year)


    } else {

      # if var is specified then simulate negative binomial instead
      negbin_theta <- mean_freq / ((var_freq / mean_freq) - 1)

      losses <- dplyr::tibble(num = MASS::rnegbin(num_sims, mean_freq, negbin_theta)) |>
        dplyr::mutate(year = 1:nrow(dplyr::pick(num))) |>
        dplyr::mutate(year = purrr::map2(year, num, function(x, y) rep(x, y))) |>
        tidyr::unnest(year) |>
        dplyr::select(year)

    }

  }


  # add simulated severities using severity pararmeters provided
  random_number <- stats::runif(nrow(losses), 0, 1)
  if(sev_dist == "lognorm2") sev <- dplyr::tibble(loss = stats::rlnorm(nrow(losses), sev_param1, sev_param2) + risk_shift)
  if(sev_dist == "gamma") sev <- dplyr::tibble(loss = stats::rgamma(nrow(losses), shape = sev_param1, scale = sev_param2) + risk_shift)
  if(sev_dist == "pareto2") sev <- dplyr::tibble(loss = gamlss.dist::rPARETO2o(nrow(losses), sev_param1, sev_param2) + risk_shift)
  if(sev_dist == "empirical") sev <- dplyr::tibble(loss = stats::approx(empirical[["probability"]], empirical[["loss"]], xout = random_number, ties = "ordered")$y + risk_shift)
  if(sev_dist == "weibull") sev <- dplyr::tibble(loss = stats::rweibull(nrow(losses), shape = sev_param1, scale = sev_param2) + risk_shift)

  losses <- losses |>
    dplyr::bind_cols(sev)

  if(!missing(label)) losses <- losses |> dplyr::mutate(label = {{ label }})

  return(losses)

}




