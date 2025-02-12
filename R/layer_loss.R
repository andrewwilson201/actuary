#' calculate occurrence or aggregate layer losses
#'
#' create a column with occurrence or aggregate losses into a layer with specified limit, retention and reinstatements
#'
#' @param year the variable in the input data that represents simulation year
#' @param loss the variable in the input data that represents the loss amounts that will feed into the layer
#' @param retention the retention the layer being priced
#' @param limit the limit of the layer being priced
#' @param reinstatements the number of reinstatements. defaults to 999
#' @param type type of layer. occurence or aggregate. defaults to occurrence.
#'
#' @return losses into the layer
#'
#' @examples
#' \dontrun{
#'
#' # calculate 10m x 10m layer with 1 reinstatement
#' losses |> dplyr::mutate(l1 = layer_loss(year, amount, 10e6, 10e6, 1))
#'
#' # calculate 10m xs 10m agg layer with 0 reinstatements
#' losses |> dplyr::mutate(l1 = layer_loss(year, amount, 10e6, 10e6, 0, type = "aggregate"))
#'
#' }
#'
#' @export


layer_loss <- function(year, loss, retention, limit, reinstatements = 999, type = "occurrence") {

<<<<<<< HEAD
=======
  .datatable.aware <<- TRUE

>>>>>>> d0287f3fcdd5818b6bb577ea40f578385a0432fd
  # warn if used within a case_when
  if (any(grepl("case_when", deparse(sys.calls())))) {
    warning("layer_loss used inside a case_when() statement may return unexpected results")
  }

  # if type is not one of the available options then stop with an error
  if(!(type %in% c("occurrence", "aggregate"))) stop("type must be one of occurrence or aggregate")

  if (type == "occurrence") {

    dplyr::tibble(year = {{ year }}, loss = {{ loss }}) |>
      # use data.table calculations
      dtplyr::lazy_dt() |>
      # group by year and calculate
      # the loss to the layer (unlimited reinstatements)
      # the running total (within the year) of the total unlimited layer loss
      # the running total (within the year) of the ul layer loss lagged by 1 row
      dplyr::group_by(year) |>
      dplyr::mutate(ul_layer_loss = pmin(pmax(loss - retention, 0), limit),
                    cumulative_loss = cumsum(ul_layer_loss),
                    cumulative_loss_lag = dplyr::lag(cumsum(ul_layer_loss), default = 0)) |>
      dplyr::ungroup() |>
      # calculate capped layer loss allowing for reinstatements
      dplyr::mutate(capped_layer_loss = dplyr::if_else(cumulative_loss > limit * (reinstatements + 1),
                                                       limit * (reinstatements + 1) - cumulative_loss_lag,
                                                       ul_layer_loss),
                    capped_layer_loss = pmax(capped_layer_loss, 0)) |>
      dplyr::collect() |>
      # extract required variable
      dplyr::pull(capped_layer_loss)

  } else {

    dplyr::tibble(year = {{ year }}, loss = {{ loss }}) |>
      # use data.table calculations
      dtplyr::lazy_dt() |>
      # group by year and calculate cumulative loss
      dplyr::group_by(year) |>
      dplyr::mutate(cumulative_loss = base::cumsum(loss)) |>
      dplyr::ungroup() |>
      # calculate running sum of the loss to the aggregate layer
      dplyr::mutate(layer_loss = base::pmin(base::pmax(cumulative_loss - retention, 0), limit * (reinstatements + 1))) |>
      dplyr::group_by(year) |>
      # calculate lagged running sum
      dplyr::mutate(layer_loss_lag = dplyr::lag(layer_loss, default = 0)) |>
      dplyr::ungroup() |>
      # calculate contribution to the aggregate loss at the individual loss level
      dplyr::mutate(agg_layer_loss = layer_loss - layer_loss_lag) |>
      dplyr::collect() |>
      # extract required variable
      dplyr::pull(agg_layer_loss)

  }
}

<<<<<<< HEAD
.datatable.aware <- TRUE
=======
>>>>>>> d0287f3fcdd5818b6bb577ea40f578385a0432fd
