#' Compute Bounds for Target Parameters
#'
#' This function computes the upper and lower bounds for a given target parameter
#' under specified assumptions using linear programming.
#'
#' @param tp A list representing the target parameter.
#' @param bases A list of basis functions.
#' @param dgp A list representing the data generating process.
#' @param assumptions A list of assumptions to be considered. Default is NULL.
#'
#' @return A list containing the upper and lower bounds along with their corresponding solutions.
#'
#' @importFrom stringr str_match str_split str_detect
#' @importFrom rlang .data
#' @importFrom lpSolve lp
#'
#' @export
compute_bounds = function(tp,
                          bases,
                          dgp,
                          olsslope = FALSE,
                          ivslope = FALSE,
                          ivslopeind = FALSE,
                          saturated = FALSE,
                          decreasing.MTR = FALSE,
                          params = NULL, # args should be some kind of list (e.g. list(ivslopeind = c()))
                          moment.eqns = NULL){

  objective.in = compute_gamma_star(tp = tp, bases, dgp = dgp)[[1]] %>%
    unlist()

  lenobj = length(objective.in)

  identity_matrix = diag(lenobj)

  const.mat = identity_matrix %>%
    rbind(identity_matrix)

  const.rhs = c(rep.int(1, lenobj),
                rep.int(0, lenobj))

  const.dir = c(rep("<=", lenobj),
                rep(">=", lenobj))

  if(ivslope){
    # If not one of the ones from above, we handle in a different way.
    # Examples of such assumptions are olsslope and ivslope
    # If user wants to plug in their own thing, we will need to run `computer_gamma_s()` and computer_beta_s`, and add to `const.dir`

    const.mat = const.mat %>%
      rbind(compute_gamma_s(bases = bases,
                            dgp = dgp,
                            slist = "ivslope") %>%
              unlist())

    beta_s_vector = compute_beta_s(dgp, "ivslope")

    const.rhs = const.rhs %>%
      c(beta_s_vector)

    const.dir = const.dir %>%
      c(rep("=", length(beta_s_vector)))
  }

  if(olsslope){
    # If not one of the ones from above, we handle in a different way.
    # Examples of such assumptions are olsslope and ivslope
    # If user wants to plug in their own thing, we will need to run `computer_gamma_s()` and computer_beta_s`, and add to `const.dir`

    const.mat = const.mat %>%
      rbind(compute_gamma_s(bases = bases,
                            dgp = dgp,
                            slist = "olsslope") %>%
              unlist())

    beta_s_vector = compute_beta_s(dgp, "olsslope")

    const.rhs = const.rhs %>%
      c(beta_s_vector)

    const.dir = const.dir %>%
      c(rep("=", length(beta_s_vector)))
  }

  if(ivslopeind){

    if(is.null(params$ivslopeind)){
      stop("If 'ivslopeind' is TRUE, include support in the 'params' argument.")
    }

    support = params$ivslopeind

    gamma_s_vector = compute_gamma_s(bases = bases,
                                     dgp = dgp,
                                     slist = "ivslopeind",
                                     param = support) %>%
      unlist()

    for(j in 1:length(support)){
      const.mat = rbind(const.mat,
                        gamma_s_vector[seq(j, length(gamma_s_vector), by = length(support))])
    }

    beta_s_vector = compute_beta_s(dgp, "ivslopeind", param = support)

    const.rhs = c(const.rhs, beta_s_vector)

    const.dir = c(const.dir, rep("=", times = length(beta_s_vector)))

  }

  if(saturated){

    gamma_s_vector = unlist(compute_gamma_s(bases = bases,
                                            dgp = dgp,
                                            slist = "saturated"))


    for(j in 1:(length(dgp$mtrs)*length(dgp$suppZ))){
      const.mat = rbind(const.mat,
                        gamma_s_vector[seq(j, length(gamma_s_vector), by = length(dgp$mtrs)*length(dgp$suppZ))])
    }

    beta_s_vector = compute_beta_s(dgp, "saturated")

    const.rhs = c(const.rhs, beta_s_vector)

    const.dir = c(const.dir, rep("=", times = length(beta_s_vector)))

  }

  if(decreasing.MTR){

    # Create an nxn 0 matrix
    m0 = matrix(0, nrow = lenobj/2, ncol = lenobj/2)

    # Create an nxn identity matrix
    m = diag(1, nrow = lenobj/2, ncol = lenobj/2)

    # Fill in the elements directly to the right of the diagonal with -1s
    # This is how to require decreasing MTR function
    for (i in 1:(nrow(m) - 1)) {
      m[i, i + 1] = -1
    }

    m = m %>%
      cbind(m0) %>%
      rbind(m0 %>%
              cbind(m))

    const.mat = const.mat %>%
      rbind(m)

    const.rhs = const.rhs %>%
      c(rep(0, times = lenobj))

    const.dir = const.dir %>%
      c(rep(">=", times = lenobj))

  }

  if(is.null(moment.eqns)){
    const.mat = rbind(const.mat,
                      moment.eqns$const.mat)

    const.dir = c(const.dir,
                  moment.eqns$const.dir)

    const.rhs = c(const.rhs,
                  moment.eqns$const.rhs)
  }

  Optimum_upper = lpSolve::lp(direction="max",objective.in,const.mat,const.dir,const.rhs)
  Optimum_lower = lpSolve::lp(direction="min",objective.in,const.mat,const.dir,const.rhs)

  return(list(upper_bound = Optimum_upper$objval,
              upper_solution = Optimum_upper$solution,
              upper_solution_d0 = Optimum_upper$solution[1:(length(Optimum_upper$solution)/2)],
              upper_solution_d1 = Optimum_upper$solution[((length(Optimum_upper$solution)/2+1)):length(Optimum_upper$solution)],
              lower_bound = Optimum_lower$objval,
              lower_solution = Optimum_lower$solution,
              lower_solution_d0 = Optimum_lower$solution[1:(length(Optimum_lower$solution)/2)],
              lower_solution_d1 = Optimum_lower$solution[((length(Optimum_lower$solution)/2+1)):length(Optimum_lower$solution)]))

}
