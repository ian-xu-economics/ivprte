#' Get Weights from Weight Data Frame
#'
#' @param weights.df Data frame containing weights.
#' @param D Treatment indicator (0 or 1).
#' @param nonZero Accepts logical indicating whether to output only non-zero weights
#'        or all weights.
#' @return A data frame containing weights.
#'
#' @importFrom dplyr %>% rename mutate select filter
#' @importFrom rlang .data
#'
#' @export
get_weights = function(weights.df,
                       D = c(0, 1),
                       nonZero = TRUE){

  stopifnot(D == 0 | D == 1)
  stopifnot(typeof(nonZero) == "logical")

  # Define a small tolerance value
  tol = 1e-12

  if(D == 1){ # Find
    new_weights.df <- weights.df %>%
      dplyr::rename(avgWeight = .data$avgWeightD1) %>%
      dplyr::select(.data$uStart,
                    .data$uEnd,
                    .data$avgWeight)
  } else{
    new_weights.df <- weights.df %>%
      dplyr::rename(avgWeight = .data$avgWeightD0) %>%
      dplyr::select(.data$uStart,
                    .data$uEnd,
                    .data$avgWeight)
  }

  if(nonZero == TRUE){
    new_weights.df <- new_weights.df %>%
      dplyr::filter(abs(.data$avgWeight) > tol)
  }

  rownames(new_weights.df) = 1:nrow(new_weights.df)

  return(new_weights.df)
}


#' Expand Weights from Weight Data Frame
#'
#' @description
#' The `compute_average_weights()` and `compute_average_weights_ivlike()` functions
#' return dataframes with four columns: `uStart`, `uEnd`, `avgWeightD0`, and `avgWeightD1`.
#' This function accepts a dataframe (weights.df) and vector (u), and returns the average
#' weight when d = 0, 1 for specific u values.
#'
#' @param u Treatment indicator (0 or 1).
#' @param weights.df Data frame containing weights.
#' @return A data frame containing u values and their respective average weights when
#' d = 0, 1.
#'
#' @export
expand_weights_df <- function(u, weights.df){

  sapply(u,
         function(x){
           c(u = x,
             avgWeightD1 = ifelse(
               any(x >= weights.df$uStart & (x < weights.df$uEnd | (weights.df$uEnd == 1 & x <= 1))),
               weights.df$avgWeightD1[which(x >= weights.df$uStart & (x < weights.df$uEnd | (weights.df$uEnd == 1 & x <= 1)))],
               NA_real_
             ),
             avgWeightD0 = ifelse(
               any(x >= weights.df$uStart & (x < weights.df$uEnd | (weights.df$uEnd == 1 & x <= 1))),
               weights.df$avgWeightD0[which(x >= weights.df$uStart & (x < weights.df$uEnd | (weights.df$uEnd == 1 & x <= 1)))],
               NA_real_))
         }) %>%
    t() %>%
    data.frame()

}
