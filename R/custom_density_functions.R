#' Pareto density function wrapper for fitdistrplus
#'
#' @param x Numeric vector.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @return Density evaluated at x.
dpareto_custom <- function(x, shape, scale, log = FALSE) {
  actuar::dpareto(x, shape = shape, scale = scale, log = log)
}

#' Inverse Gaussian density function wrapper for fitdistrplus
#'
#' @param x Numeric vector.
#' @param mean Mean parameter.
#' @param shape Shape parameter.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @return Density evaluated at x.
dinvgauss_custom <- function(x, mean, shape, log = FALSE) {
  actuar::dinvgauss(x, mean = mean, shape = shape, log = log)
}
