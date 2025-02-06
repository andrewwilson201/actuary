#' Pareto density function wrapper for fitdistrplus
#'
#' @param x Numeric vector.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @return Density evaluated at x.
#' @export
dpareto <- function(x, shape, scale, log = FALSE) {
  actuar::dpareto(x, shape = shape, scale = scale, log = log)
}

#' Inverse Gaussian density function wrapper for fitdistrplus
#'
#' @param x Numeric vector.
#' @param mean Mean parameter.
#' @param shape Shape parameter.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @return Density evaluated at x.
#' @export
dinvgauss <- function(x, mean, shape, log = FALSE) {
  actuar::dinvgauss(x, mean = mean, shape = shape, log = log)
}

#' Inverse Gaussian density function wrapper for fitdistrplus
#'
#' @param x Numeric vector.
#' @param mean Mean parameter.
#' @param shape Shape parameter.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @return Density evaluated at x.
#' @export
custom_invgauss <- list(
  d = function(x, mean, shape, log = FALSE) {
    dinvgauss(x, mean = mean, shape = shape, log = log)
  },
  p = function(x, mean, shape) {
    actuar::pinvgauss(x, mean = mean, shape = shape)
  },
  q = function(p, mean, shape) {
    actuar::qinvgauss(p, mean = mean, shape = shape)
  },
  r = function(n, mean, shape) {
    actuar::rinvgauss(n, mean = mean, shape = shape)
  },
  name = "invgauss"
)

#' Pareto density function wrapper for fitdistrplus
#'
#' @param x Numeric vector.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @return Density evaluated at x.
#' @export
custom_pareto <- list(
  d = function(x, shape, scale, log = FALSE) {
    dpareto(x, shape = shape, scale = scale, log = log)
  },
  p = function(x, shape, scale) {
    actuar::ppareto(x, shape = shape, scale = scale)
  },
  q = function(p, shape, scale) {
    actuar::qpareto(p, shape = shape, scale = scale)
  },
  r = function(n, shape, scale) {
    actuar::rpareto(n, shape = shape, scale = scale)
  },
  name = "pareto"
)
