# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Title
#'
#' @param res
#' @param r
#' @param n
#' @param w
#' @param m
#' @param s
#' @param J
#'
#' @return
#' @export
#'
#' @examples
draw_indicators_generic <- function(res, r, n, w, m, s, J) {
    .Call('_dspCP_draw_indicators_generic', PACKAGE = 'dspCP', res, r, n, w, m, s, J)
}

#' Title
#'
#' @param row_ind
#' @param col_ind
#' @param val
#' @param n
#' @param num_entries
#' @param linht
#' @param rd
#' @param D
#'
#' @return
#' @export
#'
#' @examples
sample_mat <- function(row_ind, col_ind, val, n, num_entries, linht, rd, D) {
    .Call('_dspCP_sample_mat', PACKAGE = 'dspCP', row_ind, col_ind, val, n, num_entries, linht, rd, D)
}

