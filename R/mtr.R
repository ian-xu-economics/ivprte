#' Create an MTR (Marginal Treatment Response) object as a list
#'
#' @param basis A list containing basis functions 'a' and 'b'.
#' @param theta A matrix of coefficients for the basis functions.
#' @return A list containing the basis functions and their coefficients (theta).
#' @export
MTR = function(basis, theta) {
  list(basis = basis, theta = theta)
}

#' Evaluate the MTR (Marginal Treatment Response) at given points
#'
#' @param mtr An MTR object containing basis functions and theta.
#' @param u A vector containing the evaluation points u.
#' @return A numeric vector containing the evaluated MTR values.
#' @export
evaluate_mtr = function(mtr, u) {

  sapply(u,
        function(x){
          sapply(mtr$basis$b,
                 function(f) f(x))
        }) %>%
    t() %*%
    matrix(mtr$theta)

}

#' Evaluate the MTR (Marginal Treatment Response) tuple at given points
#'
#' @param mtrs A list of MTR objects for each treatment level (0 and 1).
#' @param ev A data frame containing the evaluation points (z and u) and treatment assignment (d).
#' @return A numeric vector containing the evaluated MTR values.
#' @export
evaluate_mtr_tuple = function(mtrs, ev) {

  result = numeric(nrow(ev))

  for (d in 0:1) {
    subset_ev = subset(ev, d == d)
    result[d == ev$d] = evaluate_mtr(mtrs[[d + 1]], subset_ev)
  }

  return(result)
}
