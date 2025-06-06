#' Splits subsettable objects according to their columns
#'
#' @param x subsettable object to be split
#' @param f vector or list of vectors indicating the grouping of columns
#' @param drop logical indicating if levels that do not occur should be dropped
#' @param ... optional arguments to [split()]
#'
#' @returns List of objects split by columns
#'
#' @examples
#'
#' X <- matrix(1:60, 10, 6)
#' cond <- rep(1:3, each = 2)
#' splitCols(X, cond)
#'
#' @export
#'
splitCols <- function(x, f, drop = FALSE, ...) {
  if (is.list(f)) {
    stopifnot(all(sapply(f, length) == ncol(x)))
  } else {
    stopifnot(length(f) == ncol(x))
  }
  lapply(
    split(seq_len(ncol(x)), f, drop, ...),
    function(ind) x[, ind, drop = FALSE]
  )
}
