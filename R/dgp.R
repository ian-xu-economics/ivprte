#' DGP R6 Class
#'
#' An R6 class for representing a Data Generating Process (DGP).
#'
#' @field suppZ Support of Z.
#' @field densZ Density of Z.
#' @field pscoreZ Propensity score.
#' @field mtrs List of Marginal Treatment Response (MTR) functions.
DGP = R6::R6Class("DGP",
                  public = list(
                  suppZ = NULL,
                  densZ = NULL,
                  pscoreZ = NULL,
                  mtrs = NULL,

                  #' Initialize DGP object
                  #'
                  #' @param suppZ Support of Z.
                  #' @param densZ Density of Z.
                  #' @param pscoreZ Propensity score.
                  #' @param mtrs List of MTR functions.
                  #' @return A new DGP object.
                  initialize = function(suppZ, densZ, pscoreZ, mtrs) {
                    # Check that dimensions of suppZ and lengths of densZ and pscoreZ match
                    if (length(suppZ) != length(densZ)) {
                      stop("Dimensions of suppZ and length of densZ do not match.")
                    }

                    if (length(densZ) != length(pscoreZ)) {
                      stop("Lengths of densZ and pscoreZ do not match.")
                    }

                    # Check that densZ values are non-negative
                    if (any(densZ < 0)) {
                      stop("densZ contains negative values.")
                    }

                    # Check that densZ values sum to 1
                    if (sum(densZ) != 1) {
                      stop("densZ values do not sum to 1.")
                    }

                    self$suppZ = suppZ
                    self$densZ = densZ
                    self$pscoreZ = pscoreZ
                    self$mtrs = mtrs
                  },

                  #' Find Propensity Score for Given Z
                  #'
                  #' @param z Value of Z.
                  #' @return Propensity score corresponding to Z.
                  find_pscore = function(z) {
                    idx = which(self$suppZ == z)
                    if (length(idx) != 1) {
                      stop("Expected to find one and only one match for z.")
                    }
                    return(self$pscoreZ[idx])
                  },

                  #' Find Density for Given Z
                  #'
                  #' @param z Value of Z.
                  #' @return Density corresponding to Z.
                  find_density = function(z) {
                    idx = which(self$suppZ == z)
                    if (length(idx) != 1) {
                      stop("Expected to find one and only one match for z.")
                    }
                    return(self$densZ[idx])
                  }
                )
              )

#' Create DGP Object for MST2018 Model
#'
#' This function creates a DGP object based on the MST2018 model.
#'
#' @param suppZ Support of Z.
#' @param densZ Density of Z.
#' @param pscoreZ Propensity score of Z.
#' @param mtrs List of Marginal Treatment Response (MTR) functions.
#' @param MST2018 A Boolean indicating whether to use the DGP in Mogstad, Santos, and Torgovitsy (2018).
#'
#' @return A new DGP object.
#' @export
#'
#' @examples
#' basis0 = bernstein_basis(2)
#' theta0 = c(0.6, 0.4, 0.3)
#' mtr0 = MTR(basis0, theta0)
#'
#' basis1 = bernstein_basis(2)
#' theta1 = c(0.75, 0.5, 0.25)
#' mtr1 = MTR(basis1, theta1)
#'
#' dgp(suppZ = 0:2,
#'     densZ = c(0.5, 0.4, 0.1),
#'     pscoreZ = c(0.35, 0.6, 0.7),
#'     mtrs = list(mtr0, mtr1))
#'
dgp = function(suppZ = NULL,
               densZ = NULL,
               pscoreZ = NULL,
               mtrs = NULL,
               MST2018 = FALSE) {

  if(MST2018 == TRUE &&
       !is.null(suppZ) &&
       !is.null(densZ) &&
       !is.null(pscoreZ) &&
       !is.null(mtrs)){
    stop("`MST2018` must be FALSE if `suppZ`, `densZ`, `pscoreZ`, and `mtrs` are not NULL.")
  }

  if(MST2018 == TRUE){
    suppZ = 0:2
    densZ = c(0.5, 0.4, 0.1)
    pscoreZ = c(0.35, 0.6, 0.7)

    basis0 = bernstein_basis(2)
    theta0 = c(0.6, 0.4, 0.3)
    mtr0 = MTR(basis0, theta0)

    basis1 = bernstein_basis(2)
    theta1 = c(0.75, 0.5, 0.25)
    mtr1 = MTR(basis1, theta1)
  }

  return(DGP$new(suppZ, densZ, pscoreZ, list(mtr0, mtr1)))
}
