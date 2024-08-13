#' Compute Average Weights for Target Parameters
#'
#' This function computes the average weights for different target parameters.
#'
#' @param tp The target parameter to compute the average weights for.
#' @param dgp An optional DGP object containing the data generating process.
#' @param data An optional data.frame from which the estimated weights are calculated from.
#' @param late.lb Required when `tp` is set to "LATE", specifies the lower bound for integration.
#' @param late.ub Required when `tp` is set to "LATE", specifies the upper bound for integration.
#' @return A data frame containing the average weights.
#'
#' @importFrom dplyr %>% group_by summarise n
#' @importFrom utils tail head
#'
#' @export
compute_average_weights = function(tp = c("AUO", "ATO", "ATE", "ATT", "ATU", "LATE"),
                                   dgp = NULL,
                                   data = NULL,
                                   late.lb = NULL,
                                   late.ub = NULL) {
  if(tp == "AUO"){
    weights <- data.frame(uStart = 0,
                          uEnd = 1,
                          avgWeightD1 = 0,
                          avgWeightD0 = 1)
  } else if(tp == "ATO"){
    weights <- data.frame(uStart = 0,
                          uEnd = 1,
                          avgWeightD1 = 1,
                          avgWeightD0 = 0)
  } else if(tp == "ATE"){
    weights <- data.frame(uStart = 0,
                          uEnd = 1,
                          avgWeightD1 = -1,
                          avgWeightD0 = 1)
  } else if(tp == "ATT"){

    if(!is.null(dgp)){

      suppZ <- dgp$suppZ
      probZ <- dgp$densZ
      pscoreZ <- dgp$pscoreZ

      probD1 <- sum(probZ*pscoreZ)

    } else{

      # Calculate P[Z = z] and P[D=1 | Z=z] from data
      summaryZDF = data %>%
        dplyr::group_by(.data$zvals) %>%
        dplyr::summarise(probZ = dplyr::n()/nrow(data),
                         pscoreZ = sum(.data$dvals)/dplyr::n())

      suppZ = summaryZDF$zvals
      probZ = summaryZDF$probZ
      pscoreZ = summaryZDF$pscoreZ

      probD1 = mean(data$dvals)

    }

    weights = NULL

    uPoints <- c(0, pscoreZ, 1)

    for(i in 1:(length(uPoints)-1)){

      avgWeightD1 = ifelse(test = i == 1,
                           yes = (1/probD1) * sum(probZ),
                           no = (1/probD1) * sum(utils::tail(probZ, -i+1)))

      weights <- weights %>%
        rbind(data.frame(uStart = uPoints[i],
                         uEnd = uPoints[i+1],
                         avgWeightD1 = avgWeightD1,
                         avgWeightD0 = -avgWeightD1))
    }

  } else if(tp == "ATU"){

    if(!is.null(dgp)){

      suppZ <- dgp$suppZ
      probZ <- dgp$densZ
      pscoreZ <- dgp$pscoreZ

      probD0 <- sum(probZ*(1-pscoreZ))

    } else{

      # Calculate P[Z = z] and P[D=1 | Z=z] from data
      summaryZDF = data %>%
        dplyr::group_by(.data$zvals) %>%
        dplyr::summarise(probZ = n()/nrow(data),
                         pscoreZ = sum(.data$dvals)/n())

      suppZ = summaryZDF$zvals
      probZ = summaryZDF$probZ
      pscoreZ = summaryZDF$pscoreZ

      probD0 = mean(1 - data$dvals)

    }

    weights = NULL

    uPoints <- c(0, pscoreZ, 1)

    for(i in 1:(length(uPoints)-1)){

      avgWeightD1 = (1/probD0) * sum(utils::head(probZ, i-1))

      weights <- weights %>%
        rbind(data.frame(uStart = uPoints[i],
                         uEnd = uPoints[i+1],
                         avgWeightD1 = avgWeightD1,
                         avgWeightD0 = -avgWeightD1))
    }

  } else if(tp == "LATE"){
    if(is.null(late.lb) || is.null(late.ub)){
      stop("The LATE target parameter requires entry of the integral lower and upper bounds.")
    }

    if(late.lb < 0 || late.lb > 1){
      stop("The integral lower and upper bounds must be within [0,1].")
    }

    weightLower <- if(late.lb > 0){
     data.frame(uStart = 0,
                uEnd = late.lb,
                avgWeightD1 = 0,
                avgWeightD0 = 0)
    } else{
      NULL
    }

    weightMid <- data.frame(uStart = late.lb,
                            uEnd = late.ub,
                            avgWeightD1 = 1/(late.ub - late.lb),
                            avgWeightD0 = -1/(late.ub - late.lb))

    weightUpper <- if(late.ub < 1){
      data.frame(uStart = late.ub,
                 uEnd = 1,
                 avgWeightD1 = 0,
                 avgWeightD0 = 0)
    } else{
      NULL
    }

    weights <- rbind(weightLower, weightMid, weightUpper)
  } else {
    stop("Weight computation is unsupported for this target parameter.")
  }

  return(weights)
}

#' Compute Average Weights for IV-like Parameters
#'
#' This function computes the average weights for IV-like parameters.
#'
#' @param ivlike A list containing the name of the IV-like parameter and its function 's'.
#' @param dgp A DGP object containing the data generating process.
#' @return A data frame containing the average weights.
#'
#' @importFrom dplyr %>%
#'
#' @export
compute_average_weights_ivlike = function(ivlike, dgp) {

  order = order(dgp$pscoreZ)

  s = ivlike$s

  if(ivlike$name == "Saturated"){
    terms = function(d) {
      s_values = sapply(1:length(dgp$suppZ), function(i) s(d, dgp$suppZ[i]))
      return(s_values * dgp$densZ)
    }
  } else{
    terms = function(d){
      s(d, dgp$suppZ) * dgp$densZ
    }
  }

  # d = 1 weights
  d1terms = terms(1)
  summands = c(0, d1terms[rev(order)])  # order by decreasing pscore
  weightsD1 = rev(cumsum(summands))

  # d = 0 weights
  d0terms = terms(0)
  summands = c(0, d0terms[order])  # order by increasing pscore
  weightsD0 = cumsum(summands)

  weights = NULL

  uPoints <- c(0, dgp$pscoreZ, 1)

  for(i in 1:(length(uPoints)-1)){

    weights <- weights %>%
      rbind(data.frame(uStart = uPoints[i],
                       uEnd = uPoints[i+1],
                       avgWeightD1 = weightsD1[i],
                       avgWeightD0 = weightsD0[i]))
  }

  return(weights)
}
