#' @title Generate Bates and Watts D Matrix
#'
#' @description Generate a matrix of function derivative information.
#'
#' @param func A function whose first (vector) argument is a parameter vector;
#' @param x The parameter vector passed to `func`;
#' @param method.args Arguments of the algorithm, used for the approximation algorithm;
#' @param ... Additional arguments passed to `func`.
#'
#' @return A list of...
#'
#' @details The derivatives are calculated numerically using Richardson improvement...
#' @export
#'
#' @examples
#' func <- function(x) c(x[1], x[1], x[2]^2)
#' xxx <- c(2, 2, 5)
#' RcppnumDeriv::genD(func, xxx)
genD <- function(func, x, method.args = list(), ...) {
  # Parse default arguments
  args <- list(eps = 1e-4, d = 0.0001, zero.tol = sqrt(.Machine$double.eps / 7e-7), r = 4, v = 2) # default
  args[names(method.args)] <- method.args
  d <- args$d
  r <- args$r
  v <- args$v
  if (v != 2) stop("The current code assumes v is 2 (the default).")

  # Pass a function that incorporate the dots to C++
  funcin <- function(xx = x) func(xx, ...)

  h0 <- abs(d * x) + args$eps * (abs(x) < args$zero.tol)

  # Call internal function
  D <- genD_cpp(h0 = h0, func = funcin, x = x, d = d, r = r, v = v, eps = args$eps, zero_tol = args$zero.tol)
  class(D) <- "Darray"
  return(D)
}
