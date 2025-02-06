#' example triangle data
#'
#' tidy example triangle data. one row per development year and cohort. no development periods beyond the latest evaluation.
#'
#' @format
#' A dataframe with 441 rows and 3 columns
#' \describe{
#'       \item{uw_year}{underwriting year}
#'       \item{dev_year}{development year}
#'       \item{claim_number}{number of claims}
#' ...
#' }
#' @source <data_raw/create_triangle_data.R>
"triangle_data"

#' example simulated loss data
#'
#' monte carlo simulations with 20k simulated years. based on a poisson frequency distribution and lognormal severity distribution
#'
#' @format 
#' A data frame with 49,936 rows and 3 columns:
#' \describe{
#'   \item{year}{simulated year}
#'   \item{amount}{ground-up simulated loss amount}
#'   \item{territory_peril}{the territory peril in which the loss happened}
#'   ...
#' }
#' @source <data_raw/create example losses.R>
"losses"