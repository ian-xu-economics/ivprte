#' Create a Target Parameter Structure
#'
#' @param name Name of the target parameter.
#' @param int_limits Function for integration limits.
#' @param int_constant Function for integration constant.
#' @return A list representing the target parameter.
#' @export
TargetParameter <- function(name, int_limits, int_constant) {
  return(list(name=name, int_limits=int_limits, int_constant=int_constant))
}

#' Evaluate Target Parameter
#'
#' @param target.parameter Target parameter object.
#' @param mtrs List of MTR objects.
#' @param dgp Data Generating Process object.
#' @return A numeric value representing the evaluated target parameter.
#' @export
eval_tp <- function(target.parameter, mtrs, dgp) {

  gamma_star <- compute_gamma_star(target.parameter, list(mtrs[[1]]$basis, mtrs[[2]]$basis), dgp)
  total <- 0
  for (l in seq_along(mtrs)) {
    for (d in 0:1) {
      total <- total + sum(gamma_star[[l]][[d + 1]] * mtrs[[d + 1]]$theta)
    }
  }
  return(total)
}

#' Compute Gamma Star for Each Model
#'
#' @param target.parameter Target parameter object.
#' @param bases List of basis functions for each model.
#' @param dgp Data Generating Process object.
#' @return A list of Gamma_star matrices for each model.
#' @export
compute_gamma_star <- function(target.parameter, bases, dgp) {
  lapply(seq_along(bases), function(l) {
    lapply(0:1, function(d) {
      compute_single_gamma_star(target.parameter, bases[[l]], d, l, dgp)
    })
  })
}

#' Compute Single Gamma Star
#'
#' @param target.parameter Target parameter object.
#' @param basis Basis functions for a specific model.
#' @param d Treatment indicator (0 or 1).
#' @param l Model index.
#' @param dgp Data Generating Process object.
#' @return A Gamma_star matrix for the specific basis.
#' @export
compute_single_gamma_star <- function(target.parameter, basis, d, l, dgp) {
  gamma_star <- matrix(0, nrow=1, ncol=length(basis$b))
  for (z in 1:length(dgp$suppZ)) {
    il <- target.parameter$int_limits(dgp$suppZ[z])
    gamma_star <- gamma_star + sapply(basis$ib, function(ibk) ibk(il[1], il[2])) *
        target.parameter$int_constant(l, d, dgp$suppZ[z])

  }
  return(gamma_star)
}

#' Create LATE Target Parameter
#'
#' @param dgp Data Generating Process object.
#' @param u1 Lower limit for LATE.
#' @param u2 Upper limit for LATE.
#' @param l Model index (default is 1).
#' @return A TargetParameter object for LATE.
#'
#' @importFrom glue glue
#'
#' @export
late <- function(dgp, u1, u2, l = 1){
  stopifnot(u1 <= u2)

  name <- "LATE(u1, u2)"
  int_limits <- function(z){
    return(c(u1,u2))
  }
  int_constant <- function(l, d, z){
    return((l==1)*(2*d-1)*dgp$find_density(z)/(u2-u1))
  }

  return(TargetParameter(name, int_limits, int_constant))
}

#' Create ATT Target Parameter
#'
#' @param dgp Data Generating Process object.
#' @param l Model index (default is 1).
#' @return A TargetParameter object for ATT.
#' @export
att <- function(dgp, l = 1){
  prd1 <- sum(dgp$pscoreZ * dgp$densZ)

  name <- "ATT"

  int_limits <- function(z){
    return(c(0, dgp$find_pscore(z)))
  }

  int_constant <- function(l,d,z){
    return((l==1)*(2*d-1)/prd1*dgp$find_density(z))
  }


  return(TargetParameter(name, int_limits, int_constant))

}

#' Create ATU Target Parameter
#'
#' @param dgp Data Generating Process object.
#' @param l Model index (default is 1).
#' @return A TargetParameter object for ATT.
#' @export
atu <- function(dgp, l = 1){
  prd0 <- 1 - sum(dgp$pscoreZ * dgp$densZ)

  name <- "ATT"

  int_limits <- function(z){
    return(c(dgp$find_pscore(z), 1))
  }

  int_constant <- function(l,d,z){
    return((l==1)*(2*d-1)/prd0*dgp$find_density(z))
  }

  return(TargetParameter(name, int_limits, int_constant))

}

#' Create ATE Target Parameter
#'
#' @param dgp Data Generating Process object.
#' @param l Model index (default is 1).
#' @return A TargetParameter object for ATT.
#' @export
ate <- function(dgp, l = 1){

  name <- "ATE"

  int_limits <- function(z){
    return(c(0, 1))
  }

  int_constant <- function(l,d,z){
    return((l==1)*(2*d-1))
  }

  return(TargetParameter(name, int_limits, int_constant))

}

#' Create AUO Target Parameter
#'
#' @param dgp Data Generating Process object.
#' @param l Model index (default is 1).
#' @return A TargetParameter object for ATT.
#' @export
auo <- function(dgp, l = 1){

  name <- "AUO"

  int_limits <- function(z){
    return(c(0, 1))
  }

  int_constant <- function(l,d,z){
    return((l==1)*(1-d))
  }

  return(TargetParameter(name, int_limits, int_constant))

}

#' Create ATO Target Parameter
#'
#' @param dgp Data Generating Process object.
#' @param l Model index (default is 1).
#' @return A TargetParameter object for ATT.
#' @export
ato <- function(dgp, l = 1){

  name <- "ATO"

  int_limits <- function(z){
    return(c(0, 1))
  }

  int_constant <- function(l,d,z){
    return((l==1)*d)
  }

  return(TargetParameter(name, int_limits, int_constant))

}

#' Compute Target Parameter Population Values
#'
#' @param target.parameter A vector of target parameters.
#' @param dgp Data Generating Process object.
#' @param late.lb Lower limit for LATE.
#' @param late.ub Upper limit for LATE.
#'
#' @return The population value for each target parameter.
#' @export
compute_population_value <- function(target.parameter = c("AUO", "ATO", "ATE", "ATT", "ATU", "LATE"),
                                     dgp,
                                     late.lb = NULL,
                                     late.ub = NULL){

  sapply(target.parameter,
         function(x){
           if(x == "AUO"){
             tp <- auo(dgp)
           } else if(x == "ATO"){
             tp <- ato(dgp)
           } else if(x == "ATE"){
             tp <- ate(dgp)
           } else if(x == "ATT"){
             tp <- att(dgp)
           } else if(x == "ATU"){
             tp <- atu(dgp)
           } else if(x == "LATE"){

             if(is.null(late.lb) || is.null(late.ub)){
               cli::cli_abort("'late.lb' and 'late.ub' must be specified when the target parameter is LATE.")
             }

             tp <- late(dgp, late.lb, late.ub)
           } else{
             cli::cli_abort(paste0(x, " is not a recognized target parameter."))
           }

           result <- eval_tp(tp,
                             mtrs = dgp$mtrs,
                             dgp = dgp)

           return(result)
         })

}

#' Compute Target Parameter Population Bounds
#'
#' @param target.parameter A vector of target parameters.
#' @param dgp Data Generating Process object.
#' @param bases A list of basis functions.
#' @param late.lb Lower limit for LATE.
#' @param late.ub Upper limit for LATE.
#' @param assumptions A list of assumptions to be considered. Default is NULL.
#' @param assumptions.extra A list of supplementary parameters that go along with the assumptions. Default is NULL.
#'
#' @return The population value for each target parameter.
#' @export
compute_population_bounds <- function(target.parameter = c("AUO", "ATO", "ATE", "ATT", "ATU", "LATE"),
                                      dgp,
                                      bases,
                                      late.lb,
                                      late.ub,
                                      assumptions,
                                      assumptions.extra){

  sapply(target.parameter,
         function(x){
           if(x == "AUO"){
             tp <- auo(dgp)
           } else if(x == "ATO"){
             tp <- ato(dgp)
           } else if(x == "ATE"){
             tp <- ate(dgp)
           } else if(x == "ATT"){
             tp <- att(dgp)
           } else if(x == "ATU"){
             tp <- atu(dgp)
           } else if(x == "LATE"){

             if(is.null(late.lb) || is.null(late.ub)){
               cli::cli_abort("'late.lb' and 'late.ub' must be specified when the target parameter is LATE.")
             }

             tp <- late(dgp, late.lb, late.ub)
           } else{
             cli::cli_abort(paste0(x, " is not a recognized target parameter."))
           }

           result <- compute_bounds(tp,
                                    bases,
                                    dgp,
                                    assumptions,
                                    assumptions.extra)

           return(c("lower.bound" = result$lower_bound,
                    "upper.bound" = result$upper_bound))
         }) |>
    t()

}



