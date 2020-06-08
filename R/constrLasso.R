#' @title The function constrLassoReg
#'
#' @description This function fits a linearly constrained lasso regression, using a predictor matrix X, a response y and a tuning parameter value lambda. It results in a vector of coefficient estimates.
#' The function corresponds to lsq_constrsparsereg of the SparseReg MATLAB-toolbox by Zhou and Gaines, see \href{http://hua-zhou.github.io/SparseReg/}{the project page}.
#' Only the Quadratic Programming Algorithm for ENET is implemented.
#' For more information, see \href{http://hua-zhou.github.io/media/pdf/GainesKimZhou08CLasso.pdf}{\insertCite{gaines2018algorithms;textual}{constrLasso}}.
#'
#' @param X an nxp matrix of p regressors with n observations.
#' @param y an nx1 response vector with n observations.
#' @param lambda a tuning parameter value for the lasso penalty. Default value is lambda=0.
#' @param Aeq a cxp equality constraint matrix, containing c constraints for p regressors. Default value is Aeq=NULL, no equality constraints.
#' @param beq a cx1 equality constraint vector. Default value is beq=NULL, no equality constraints.
#' @param A a cxp inequality constraint matrix, containing c constraints for p regressors. Default value is A=NULL, no inequality constraints.
#' @param b a cx1 inequality constraint vector. Default value is b=NULL, no inequality constraints.
#' @param penidx a logical px1 vector, indicating which coefficients are to be penalized. Default value is penidx=NULL and imposes penalty on all p coefficients.
#' @param method a character string, the method to be used. Possible values are "QP" (default) for Quadratic Programming and "CD" for Coordinate Descent.
#' "CD" uses the glmnet package and only works without equality and inequality constraints.
#'
#' @section Details:
#' The Constrained Lasso as in \insertCite{gaines2018algorithms;textual}{constrLasso} minimizes
#' \deqn{0.5||y - X \beta ||^2_2 + \lambda||\beta||_1,}
#' subject to \eqn{Aeq \beta = beq} and \eqn{A \beta\le b}.
#'
#' @return betahat a px1 vector of estimated coefficients.
#' @return dualEq equality duals. The function returns an empty vector for no equality constraints.
#' @return dualIneq inequality duals. The function returns an empty vector for no inequality constraints.
#'
#' @examples
#' library(constrLasso)
#' library(MASS)
#' set.seed(1234)
#' n <- 200
#' p <- 50
#' Xmat <- matrix(,n,p)
#' for(i in 1:p){Xmat[,i] <- rnorm(n,runif(1,-3,3),runif(1,1,2))}
#' betas <- runif(p,-2,2)
#' nonzeros <- sample(1:p,20,replace=FALSE)
#' yvec <- Xmat[,nonzeros]%*%betas[nonzeros] + rnorm(n,0,2)
#' classoreg_results <- constrLassoReg(Xmat,yvec,lambda=0)
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#'
#' @export constrLassoReg
constrLassoReg <-
  function(X,
           y,
           lambda,
           Aeq = NULL,
           beq = NULL,
           A = NULL,
           b = NULL,
           penidx = NULL,
           method = "QP") {
    X <- as.matrix(X)
    n <- dim(X)[1]
    p <- dim(X)[2]
    y <- as.matrix(y)
    dim(y) <- c(n, 1)

    if (is.null(penidx)) {
      penidx <- matrix(TRUE, p, 1)
    }
    dim(penidx) <- c(p, 1)

    if (is.null(A)) {
      A <- matrix(0, 0, p)
      b <- rep(0, 0)
    }

    if (is.null(Aeq)) {
      Aeq <- matrix(0, 0, p)
      beq <- rep(0, 0)
    }

    m1 <- dim(Aeq)[1]
    m2 <- dim(A)[1]

    if (is.null(lambda)) {
      stop(cat("Please, enter a value for the penalty parameter lambda."))
    }

    ## no constraints
    if (m1 == 0 &&  m2 == 0) {
      #no penalization
      if (abs(lambda) < 1e-16) {
        wt <- matrix(1, n, 1)
        dim(wt) <- c(n, 1)
        Xwt <- X * as.numeric(sqrt(wt))
        ywt <- as.numeric(sqrt(wt)) * y

        betahat <- matrix(lm(ywt ~ Xwt - 1)$coef, p, 1)
        dualEq <-  rep(0, 0)
        dualIneq <-  rep(0, 0)

      } else{
        # with penalization
        if (method == "CD") {
          glmnet_res <-
            glmnet(
              Xwt,
              ywt,
              alpha = 1,
              lambda = lambda,
              penalty.factor = penidx,
              maxit = 1000,
              intercept = FALSE,
              standardize = TRUE
            )

          betahat <- matrix(glmnet_res$beta, p, 1)
          dualEq <-  rep(0, 0)
          dualIneq <-  rep(0, 0)

        } else if (method == "QP") {
          wt <- matrix(1, n, 1)
          dim(wt) <- c(n, 1)
          Xwt <- X * as.numeric(sqrt(wt))
          ywt <- as.numeric(sqrt(wt)) * y

          # quadratic coefficient
          H <- t(Xwt) %*% Xwt
          H <- rbind(cbind(H, -H), cbind(-H, H))

          # linear coefficient
          f <- -t(Xwt) %*% ywt
          f <- rbind(f, -f) + lambda * rbind(penidx, penidx)

          # optimizer
          x <- OP(Q_objective(H, L = t(f)))
          opt <- ROI_solve(x, solver = "qpoases")
          opt_sol <-  opt$message$primal_solution

          betahat <-
            matrix(opt_sol[1:p] - opt_sol[(p + 1):length(opt_sol)], p, 1)
          dualEq <-  rep(0, 0)
          dualIneq <-  rep(0, 0)

        }

      }

    } else{
      ## with constraints

      if (method == "CD") {
        warning("The CD method does not work with constraints. The solution is generated with QP.")
      }

      wt <- matrix(1, n, 1)
      dim(wt) <- c(n, 1)
      Xwt <- X * as.numeric(sqrt(wt))
      ywt <- as.numeric(sqrt(wt)) * y

      # quadratic coefficient
      H <- t(Xwt) %*% Xwt
      H <- rbind(cbind(H, -H), cbind(-H, H))

      # linear coefficient
      f <- -t(Xwt) %*% ywt
      f <- rbind(f, -f) + lambda * rbind(penidx, penidx)

      # constraints
      Amatrix <- rbind(cbind(Aeq, -Aeq), cbind(A, -A))
      bvector <- c(beq, b)

      # optimizer
      x <-
        OP(Q_objective(H, L = t(f)),
           L_constraint(
             L = Amatrix,
             dir = c(rep("==", m1), rep("<=", m2)),
             rhs = bvector
           ))
      opt <- ROI_solve(x, solver = "qpoases")
      opt_sol <-  opt$message$primal_solution

      betahat <-
        matrix(opt_sol[1:p] - opt_sol[(p + 1):length(opt_sol)], p, 1)
      duals <-  opt$message$dual_solution[-(1:(2 * p))]

      if (m1 != 0) {
        dualEq <- -matrix(duals[1:m1], m1, 1)
      } else{
        dualEq <- rep(0, 0)
      }
      if (m2 != 0) {
        dualIneq <- -matrix(duals[(m1 + 1):(m1 + m2)], m2, 1)
      } else{
        dualIneq <- rep(0, 0)
      }
    }

    return(list(
      "betahat" = betahat,
      "dualEq" = dualEq,
      "dualIneq" = dualIneq
    ))
  }


#' @title The function constrLassoPath
#'
#' @description This function performs a Constrained Lasso Solution Path as in \insertCite{gaines2018algorithms;textual}{constrLasso}. It computes the solution path for the constrained lasso problem, using a predictor matrix X and a response y.
#' The constrained lasso solves the standard lasso \insertCite{tibshirani1996regression;textual}{constrLasso}subject to the linear equality constraints \eqn{Aeq \beta = beq} and linear inequality constraints \eqn{A \beta \le b}.
#' The result lambdaPath contains the values of the tuning parameter along the solution path and betaPath - the estimated regression coefficients for each value of lambdaPath.
#' The function corresponds to lsq_classopath of the SparseReg MATLAB-toolboxof by Zhou and Gaines, see \href{http://hua-zhou.github.io/SparseReg/}{the project page}.
#' For more information, see \href{http://hua-zhou.github.io/media/pdf/GainesKimZhou08CLasso.pdf}{\insertCite{gaines2018algorithms;textual}{constrLasso}}.
#'
#' @param X an nxp matrix with p regressors with n observations.
#' @param y an nx1 response vector with n observations.
#' @param Aeq a cxp equality constraint matrix, containing c constraints for p regressors. Default value is Aeq=NULL, no equality constraints.
#' @param beq a cx1 equality constraint vector. Default value is beq=NULL, no equality constraints.
#' @param A a cxp inequality constraint matrix, containing c constraints for p regressors. Default value is A=NULL, no inequality constraints.
#' @param b a cx1 inequality constraint vector. Default value is b=NULL, no inequality constraints.
#' @param penidx a logical px1 vector, indicating which coefficients are to be penalized. Default value is penidx=NULL and allows all p coefficients to be penalized.
#' @param init_method a character string, the initializing method to be used. Possible values are "QP" (default) for Quadratic Programming and "LP" for Linear Programming.
#' "LP" is recommended only when it's reasonable to assume that all coefficient estimates initialize at zero.
#' @param epsilon a tuning parameter for ridge penalty in case of high-dimensional (n>p) regressors matrix X. Default value is 1e-4.
#' @param stop_lambdatol a tolerance value for the tuning lasso parameter. The algorithm stops when hitting this tolerance. Default value is 1e-7.
#' @param ceiling_tol a tolerance value for the change in subgradients. Default value is 1e-10.
#' @param zeros_tol a tolerance value for the zero equality of coefficients. Default value is 1e-20.
#' @param verbose a logical parameter. TRUE prints along the constraint lasso solution path. Default value is FALSE.
#'
#' @section Details:
#' The Constrained Lasso as in \insertCite{gaines2018algorithms;textual}{constrLasso} minimizes
#' \deqn{0.5||y - X \beta ||^2_2 + \lambda||\beta||_1,}
#' subject to \eqn{Aeq \beta = beq} and \eqn{A \beta\le b}.
#'
#' @return lambdaPath a vector of the tuning parameter values along the solution path.
#' @return betaPath  a matrix with estimated regression coefficients for each value of lambdaPath
#' @return dfPath a vector with degrees of freedom along the solution path
#' @return objValPath a vector with values of the objective function for each value of lambdaPath
#' @examples
#'
#' library(constrLasso)
#' set.seed(1234)
#' n <- 200
#' p <- 50
#' Xmat <- matrix(,n,p)
#' for(i in 1:p){Xmat[,i] <- rnorm(n,runif(1,-3,3),runif(1,1,2))}
#' betas <- runif(p,-2, 2)
#' nonzeros <- sample(1:p,20,replace=FALSE)
#' yvec <- Xmat[,nonzeros]%*%betas[nonzeros] + rnorm(n,0,2)
#' classopath_results <- constrLassoPath(Xmat, yvec)
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export constrLassoPath
constrLassoPath <-
  function(X,
           y,
           Aeq = NULL,
           beq = NULL,
           A = NULL,
           b = NULL,
           penidx = NULL,
           init_method = "QP",
           epsilon = 1e-4,
           stop_lambdatol = 1e-7,
           ceiling_tol = 1e-10,
           zeros_tol = 1e-20,
           verbose = FALSE) {
    if (verbose) {
      printer <- print
    } else{
      printer <- function(x) {
      }
    }

    X <- as.matrix(X)
    n <- dim(X)[1]
    n_orig <- n
    p <- dim(X)[2]

    y <- as.matrix(y)
    dim(y) <- c(n, 1)

    if (is.null(penidx)) {
      penidx <- matrix(TRUE, p, 1)
    }
    dim(penidx) <- c(p, 1)

    if (is.null(A)) {
      A <- matrix(0, 0, p)
      b <- rep(0, 0)
    }

    if (is.null(Aeq)) {
      Aeq <- matrix(0, 0, p)
      beq <- rep(0, 0)
    }

    m1 <- dim(Aeq)[1]
    m2 <- dim(A)[1]

    #check if Ridge Penalty has to be included
    if (n < p) {
      warning('Adding a small ridge penalty (default is 1e-4) since n < p')

      if (epsilon <= 0) {
        warning(
          "The ridge tuning parameter epsilon must be positive, switching to default value (1e-4)."
        )
        epsilon <- 1e-4
      }

      # create augmented data for n<p with ridge penalty
      y <- rbind(y, matrix(rep.int(0, p)))
      X <- rbind(X, sqrt(epsilon) * diag(1, p))
      n_orig <- n

    } else {
      #make sure X is full column rank
      R <- qr(X)$qr[1:p, ]
      rankX <-
        sum(matrix(abs(diag(R))) > abs(R[1]) * max(n, p) * .Machine$double.eps)
      rankX <- qr(X)$rank

      if (rankX != p) {
        warning("Adding a small ridge penalty (default is 1e-4) since X is rank deficient")
        if (epsilon <= 0) {
          warning(
            "The ridge tuning parameter epsilon must be positive, switching to default value (1e-4)."
          )
          epsilon <- 1e-4
        }
        # create augmented data for n<p with ridge penalty
        y <- rbind(y, matrix(rep.int(0, p)))
        X <- rbind(X, sqrt(epsilon) * diag(1, p))
        n_orig <- n
      }
    }

    #### Define Iterations and Path help parameters

    maxiters <- 5 * (p + m2) #max number of path segments to consider
    betaPath <- matrix(0, p, maxiters)
    dualpathEq <- matrix(0, m1, maxiters)
    dualpathIneq <- matrix(0, m2, maxiters)
    lambdaPath <- matrix(0, 1, maxiters)
    dfPath <-  matrix(Inf, 1, maxiters)
    objValPath <- matrix(0, 1, maxiters)
    violationsPath <- matrix(Inf, 1, maxiters)

    #### Initialization
    H <- t(X) %*% X

    if (init_method == "LP") {
      warning('LP is used for initialization, assumes initial solution is unique.')

      obj_fun <- matrix(1, 2 * p, 1)
      lb <- 0
      lb_mat <- diag(1, 2 * p)

      constr_mat <- rbind(cbind(Aeq, -Aeq), cbind(A, -A), lb_mat)
      constr_rhs <- c(beq, b, matrix(rep.int(lb, 2 * p), 2 * p, 1))
      constr_dir <- c(rep("=", m1),  rep("<=", m2), rep(">=", 2 * p))

      lpsol <-
        lp(
          direction = "min",
          obj_fun,
          constr_mat,
          constr_dir,
          constr_rhs,
          compute.sens = TRUE
        )
      lpres <- matrix(lpsol$solution, 2 * p, 1)

      betaPath[, 1] <- lpres[1:p] - lpres[(p + 1):(2 * p)]

      dual_eqlin <- matrix(lpsol$duals[1:m1], m1, 1)
      if (m2 != 0) {
        dual_ineqlin <- matrix(lpsol$duals[(m1 + 1):(m1 + m2)], m2, 1)
      } else{
        dual_ineqlin <- rep(0, 0)
      }

      dualpathEq[, 1] <- -dual_eqlin
      dualpathIneq[, 1] <- dual_ineqlin

    } else if (init_method == "QP") {
      obj_fun <- matrix(1, 2 * p, 1)
      lb <- 0
      lb_mat <- diag(1, 2 * p)

      constr_mat <- rbind(cbind(Aeq, -Aeq), cbind(A, -A), lb_mat)
      constr_rhs <- c(beq, b, matrix(rep.int(lb, 2 * p), 2 * p, 1))
      constr_dir <- c(rep("=", m1),  rep("<=", m2), rep(">=", 2 * p))

      lpsol <-
        lp(
          direction = "min",
          obj_fun,
          constr_mat,
          constr_dir,
          constr_rhs,
          compute.sens = TRUE
        )
      lpres <- matrix(lpsol$solution, 2 * p, 1)

      betaPath[, 1] <- lpres[1:p] - lpres[(p + 1):(2 * p)]

      dual_eqlin <- matrix(lpsol$duals[1:m1], m1, 1)
      if (m2 != 0) {
        dual_ineqlin <- matrix(lpsol$duals[(m1 + 1):(m1 + m2)], m2, 1)
      } else{
        dual_ineqlin <- rep(0, 0)
      }

      dualpathEq[, 1] <- -dual_eqlin
      dualpathIneq[, 1] <- dual_ineqlin

      # initialize sets
      dualpathIneq[which(dualpathIneq[, 1] < 0), 1] <-
        0 #fix negative duals
      setActive <- abs(betaPath[, 1]) > 1e-4 | !penidx
      betaPath[!setActive, 1] <- 0

      # find the maximum lambda and initialize subgradient vector
      resid <- y - X %*% betaPath[, 1]
      subgrad <-
        t(X) %*% resid - t(Aeq) %*% dualpathEq[, 1] - t(A) %*% dualpathIneq[, 1]
      lambda_max <- max(abs(subgrad))

      # Use QP at lambda_max to initialize
      constrLassoReg_sol <-
        constrLassoReg(
          X,
          y,
          lambda = lambda_max,
          Aeq = Aeq,
          beq = beq,
          A = A,
          b = b,
          penidx = penidx,
          method = "QP"
        )

      betaPath[, 1] <- constrLassoReg_sol$betahat
      dualpathEq[, 1] <- constrLassoReg_sol$dualEq
      dualpathIneq[, 1] <- constrLassoReg_sol$dualIneq
    }

    # initialize sets
    dualpathIneq[which(dualpathIneq[, 1] < 0), 1] <- 0
    setActive <- abs(betaPath[, 1]) > 1e-4 | !penidx
    betaPath[!setActive, 1] <- 0

    residIneq <- A %*% betaPath[, 1] - b
    setIneqBorder <- residIneq == 0
    nIneqBorder <- length(which(setIneqBorder))

    # find the maximum lambda and initialize subgradient vector
    resid <- y - X %*% betaPath[, 1]
    subgrad <-
      t(X) %*% resid - t(Aeq) %*% dualpathEq[, 1] - t(A) %*% dualpathIneq[, 1]
    lambdaPath[, 1] <- max(abs(subgrad))
    idx <- which.max(abs(subgrad))
    subgrad[setActive] <- sign(betaPath[setActive, 1])
    subgrad[!setActive] <- subgrad[!setActive] / lambdaPath[, 1]
    setActive[idx] <- TRUE
    nActive <- length(which(setActive))
    nInActive <- length(which(!setActive))

    # calculate value for objective function
    objValPath[, 1] <-
      0.5 * (sum((y - X %*% betaPath[, 1]) ^ 2)) + lambdaPath[, 1] * sum(abs(betaPath[, 1]))

    # calculate degrees of freedom
    rankAeq <- qr(Aeq)$rank
    dfPath[, 1] <- nActive - rankAeq - nIneqBorder

    # set initial violations counter to 0
    violationsPath[1] <- 0

    # direction of the algorithm (sign)
    dirsgn = -1


    ####MAIN LOOP FOR PATH FOLLOWING

    for (k in 2:maxiters) {
      printer(k)
      printer(lambdaPath[, k - 1])

      # threshold near-zero lambdas to zero and stop algorithm
      if (lambdaPath[, k - 1] <= (0 + stop_lambdatol)) {
        lambdaPath[, k - 1] <- 0
        printer(paste("BREAK. Previous Lambda < ", stop_lambdatol, ".", sep =
                        ""))
        break
      }

      M <-
        cbind(H[setActive, setActive], t(matrix(Aeq[, setActive], ncol = nActive)), t(matrix(A[setIneqBorder, setActive], ncol =
                                                                                               nActive)))

      M <-
        rbind(M, cbind(
          rbind(matrix(Aeq[, setActive], ncol = nActive), matrix(A[setIneqBorder, setActive], ncol =
                                                                   nActive)),
          matrix(0, m1 + nIneqBorder, m1 + nIneqBorder)
        ))

      ## calculate derivative
      # try using a regular inverse first, otherwise Moore-Penrose Inverse
      error.hand.inv <- function(e) {
        dir <-
          -(ginv(M) %*% rbind(
            matrix(subgrad[setActive], nActive, 1),
            matrix(0, m1 + nIneqBorder, 1)
          ))
        printer("Moore-Penrose-Inverse used.")
        return(dir)
      }
      dir <-
        tryCatch(
          dirsgn * (solve(M, rbind(
            matrix(subgrad[setActive], nActive, 1), matrix(0, m1 + nIneqBorder, 1)
          ))),
          error = error.hand.inv
        )

      if (nInActive != 0) {
        dirSubgrad <-
          -cbind(matrix(H[!setActive, setActive], ncol = nActive), t(matrix(Aeq[, !setActive], ncol =
                                                                              nInActive)), t(matrix(A[setIneqBorder, !setActive], ncol = nInActive))) %*%
          dir
      } else{
        dirSubgrad <- matrix(0, 0, m1 + nIneqBorder)
      }

      ### check additional events related to potential subgradient violations ##

      ## Inactive coefficients moving too slowly

      # Negative subgradient
      inactSlowNegIdx <-
        which((1 * dirsgn - ceiling_tol) <= subgrad[!setActive] &
                subgrad[!setActive] <= (1 * dirsgn + ceiling_tol) &
                1 * dirsgn < dirSubgrad
        )
      # Positive subgradient
      inactSlowPosIdx <-
        which((-1 * dirsgn - ceiling_tol) <= subgrad[!setActive] &
                subgrad[!setActive] <= (-1 * dirsgn + ceiling_tol) &
                dirSubgrad < -1 * dirsgn
        )

      ## "Active" coeficients estimated as 0 with potential sign mismatch
      # Positive subgradient but negative derivative
      signMismatchPosIdx = which((0 - ceiling_tol) <= subgrad[setActive] &
                                   subgrad[setActive] <= (1 + ceiling_tol) &
                                   dirsgn * dir[1:nActive] <= (0 - ceiling_tol)  &
                                   betaPath[setActive, k - 1] == 0
      )
      # Negative subgradient but positive derivative
      signMismatchNegIdx <-
        which((-1 - ceiling_tol) <= subgrad[setActive] &
                subgrad[setActive] <= (0 + ceiling_tol) &
                (0 + ceiling_tol) <= dirsgn * dir[1:nActive] &
                betaPath[setActive, k - 1] == 0
        )

      # reset violation counter (to avoid infinite loops)
      violateCounter <- 0

      ### Outer while loop for checking all conditions together

      while (length(inactSlowNegIdx) != 0 ||
             length(inactSlowPosIdx) != 0  ||
             length(signMismatchPosIdx) != 0  || length(signMismatchNegIdx) != 0) {
        printer("VIOLATIONS DUE TO SLOW ALGO MOVEMENT OR POS/NEG MISMATCH")

        ## Monitor & fix condition 1 violations
        while (length(inactSlowNegIdx) != 0) {
          printer("Violation inactSlowNegIdx")

          # Identify & move problem coefficient

          inactiveCoeffs <-
            which(!setActive) # indices corresponding to inactive coefficients
          viol_coeff <-
            inactiveCoeffs[inactSlowNegIdx] # identify problem coefficient
          setActive[viol_coeff] <-
            TRUE # put problem coefficient back into active set
          nActive <-
            length(which(setActive)) # determine new number of active/inactive coefficients
          nInActive <- length(which(!setActive))
          nIneqBorder <-
            length(which(setIneqBorder)) # determine number of active/binding inequality constraints

          # Recalculate derivative for coefficients & multipliers

          M <-
            cbind(H[setActive, setActive], t(matrix(Aeq[, setActive], ncol = nActive)), t(matrix(A[setIneqBorder, setActive], ncol =
                                                                                                   nActive)))
          M <-
            rbind(M, cbind(
              rbind(matrix(Aeq[, setActive], ncol = nActive), matrix(A[setIneqBorder, setActive], ncol =
                                                                       nActive)),
              matrix(0, m1 + nIneqBorder, m1 + nIneqBorder)
            ))

          error.hand.inv <- function(e) {
            dir <-
              -(ginv(M) %*% rbind(
                matrix(subgrad[setActive], nActive, 1),
                matrix(0, m1 + nIneqBorder, 1)
              ))
            return(dir)
          }
          dir <-
            tryCatch(
              dirsgn * (solve(M, rbind(
                matrix(subgrad[setActive], nActive, 1),
                matrix(0, m1 + nIneqBorder, 1)
              ))),
              error = error.hand.inv
            )

          # calculate derivative for lambda*subgradient
          if (nInActive != 0) {
            dirSubgrad <-
              -cbind(matrix(H[!setActive, setActive], ncol = nActive), t(matrix(Aeq[, !setActive], ncol =
                                                                                  nInActive)), t(matrix(A[setIneqBorder, !setActive], ncol = nInActive))) %*%
              dir
          } else{
            dirSubgrad <- matrix(0, 0, m1 + nIneqBorder)
          }

          # check for violations again

          # Negative subgradient
          inactSlowNegIdx <-
            which((1 * dirsgn - ceiling_tol) <= subgrad[!setActive] &
                    subgrad[!setActive] <= (1 * dirsgn + ceiling_tol) &
                    1 * dirsgn < dirSubgrad
            )
          # Positive subgradient
          inactSlowPosIdx <-
            which((-1 * dirsgn - ceiling_tol) <= subgrad[!setActive] &
                    subgrad[!setActive] <= (-1 * dirsgn + ceiling_tol) &
                    dirSubgrad < -1 * dirsgn
            )
          # Positive subgrad but negative derivative
          signMismatchPosIdx = which((0 - ceiling_tol) <= subgrad[setActive] &
                                       subgrad[setActive] <= (1 + ceiling_tol) &
                                       dirsgn * dir[1:nActive] <= (0 - ceiling_tol) &
                                       betaPath[setActive, k - 1] == 0
          )
          # Negative subgradient but positive derivative
          signMismatchNegIdx <-
            which((-1 - ceiling_tol) <= subgrad[setActive] &
                    subgrad[setActive] <= (0 + ceiling_tol) &
                    (0 + ceiling_tol) <= dirsgn * dir[1:nActive] &
                    betaPath[setActive, k - 1] == 0
            )

          # update violation counter
          violateCounter <- violateCounter + 1
          if (violateCounter >= maxiters) {
            printer("Too many violations.")
            break
          }
        }

        ## Monitor & fix subgradient condition 2 violations

        while (length(inactSlowPosIdx) != 0) {
          printer("violation inactSlowPosIdx")

          # Identify & move problem coefficient

          inactiveCoeffs <-
            which(!setActive) # indices corresponding to inactive coefficients
          viol_coeff <-
            inactiveCoeffs[inactSlowPosIdx] # identify problem coefficient
          setActive[viol_coeff] <-
            TRUE # put problem coefficient back into active set;
          nActive <-
            length(which(setActive)) # determine new number of active/inactive coefficients
          nInActive <- length(which(!setActive))
          nIneqBorder <-
            length(which(setIneqBorder)) # determine number of active/binding inequality constraints

          # Recalculate derivative for coefficients & multipliers

          M <-
            cbind(H[setActive, setActive], t(matrix(Aeq[, setActive], ncol = nActive)), t(matrix(A[setIneqBorder, setActive], ncol =
                                                                                                   nActive)))
          M <-
            rbind(M, cbind(
              rbind(matrix(Aeq[, setActive], ncol = nActive), matrix(A[setIneqBorder, setActive], ncol =
                                                                       nActive)),
              matrix(0, m1 + nIneqBorder, m1 + nIneqBorder)
            ))

          error.hand.inv <- function(e) {
            dir <-
              -(ginv(M) %*% rbind(
                matrix(subgrad[setActive], nActive, 1),
                matrix(0, m1 + nIneqBorder, 1)
              ))
            printer("Moore-Penrose-Inverse used.")
            return(dir)
          }
          dir <-
            tryCatch(
              dirsgn * (solve(M, rbind(
                matrix(subgrad[setActive], nActive, 1),
                matrix(0, m1 + nIneqBorder, 1)
              ))),
              error = error.hand.inv
            )

          if (nInActive != 0) {
            dirSubgrad <-
              -cbind(matrix(H[!setActive, setActive], ncol = nActive), t(matrix(Aeq[, !setActive], ncol =
                                                                                  nInActive)), t(matrix(A[setIneqBorder, !setActive], ncol = nInActive))) %*%
              dir
          } else{
            dirSubgrad <- matrix(0, 0, m1 + nIneqBorder)
          }

          #check for violations again

          # Positive subgradient
          inactSlowPosIdx <-
            which((-1 * dirsgn - ceiling_tol) <= subgrad[!setActive] &
                    subgrad[!setActive] <= (-1 * dirsgn + ceiling_tol) &
                    dirSubgrad < -1 * dirsgn
            )
          # Positive subgrad but negative derivative
          signMismatchPosIdx = which((0 - ceiling_tol) <= subgrad[setActive] &
                                       subgrad[setActive] <= (1 + ceiling_tol) &
                                       dirsgn * dir[1:nActive] <= (0 - ceiling_tol)  &
                                       betaPath[setActive, k - 1] == 0
          )
          # Negative subgradient but positive derivative
          signMismatchNegIdx <-
            which((-1 - ceiling_tol) <= subgrad[setActive] &
                    subgrad[setActive] <= (0 + ceiling_tol) &
                    (0 + ceiling_tol) <= dirsgn * dir[1:nActive] &
                    betaPath[setActive, k - 1] == 0
            )

          # update violation counter
          violateCounter <- violateCounter + 1

          if (violateCounter >= maxiters) {
            printer("Too many violations.")
            break
          }
        }

        ## Monitor & fix subgradient condition 3 violations

        while (length(signMismatchPosIdx) != 0) {
          printer("violation signMismatchPosIdx")

          # Identify & move problem coefficient

          activeCoeffs <-
            which(setActive) # indices corresponding to inactive coefficients
          viol_coeff <-
            activeCoeffs[signMismatchPosIdx] # identify problem coefficient
          setActive[viol_coeff] <-
            FALSE # put problem coefficient back into active set;
          nActive <-
            length(which(setActive)) # determine new number of active/inactive coefficients
          nInActive <- length(which(!setActive))
          nIneqBorder <-
            length(which(setIneqBorder)) # determine number of active/binding inequality constraints


          # Recalculate derivative for coefficients & multipliers

          M <-
            cbind(H[setActive, setActive], t(matrix(Aeq[, setActive], ncol = nActive)), t(matrix(A[setIneqBorder, setActive], ncol =
                                                                                                   nActive)))
          M <-
            rbind(M, cbind(
              rbind(matrix(Aeq[, setActive], ncol = nActive), matrix(A[setIneqBorder, setActive], ncol =
                                                                       nActive)),
              matrix(0, m1 + nIneqBorder, m1 + nIneqBorder)
            ))

          error.hand.inv <- function(e) {
            dir <-
              -(ginv(M) %*% rbind(
                matrix(subgrad[setActive], nActive, 1),
                matrix(0, m1 + nIneqBorder, 1)
              ))
            printer("Moore-Penrose-Inverse used.")
            return(dir)
          }
          dir <-
            tryCatch(
              dirsgn * (solve(M, rbind(
                matrix(subgrad[setActive], nActive, 1),
                matrix(0, m1 + nIneqBorder, 1)
              ))),
              error = error.hand.inv
            )

          if (nInActive != 0) {
            dirSubgrad <-
              -cbind(matrix(H[!setActive, setActive], ncol = nActive), t(matrix(Aeq[, !setActive], ncol =
                                                                                  nInActive)), t(matrix(A[setIneqBorder, !setActive], ncol = nInActive))) %*%
              dir
          } else{
            dirSubgrad <- matrix(0, 0, m1 + nIneqBorder)
          }

          #check for violations again

          # Positive subgrad but negative derivative
          signMismatchPosIdx = which((0 - ceiling_tol) <= subgrad[setActive] &
                                       subgrad[setActive] <= (1 + ceiling_tol) &
                                       dirsgn * dir[1:nActive] <= (0 - ceiling_tol)  &
                                       betaPath[setActive, k - 1] == 0
          )

          # Negative subgradient but positive derivative
          signMismatchNegIdx <-
            which((-1 - ceiling_tol) <= subgrad[setActive] &
                    subgrad[setActive] <= (0 + ceiling_tol) &
                    (0 + ceiling_tol) <= dirsgn * dir[1:nActive] &
                    betaPath[setActive, k - 1] == 0
            )

          # update violation counter
          violateCounter <- violateCounter + 1

          if (violateCounter >= maxiters) {
            printer("Too many violations.")
            break
          }
        }

        ## Monitor & fix subgradient condition 4 violations

        while (length(signMismatchNegIdx) != 0) {
          printer("violation signMismatchNegIdx")

          # Identify & move problem coefficient

          activeCoeffs <-
            which(setActive) # indices corresponding to inactive coefficients
          viol_coeff <-
            activeCoeffs[signMismatchNegIdx] # identify problem coefficient
          setActive[viol_coeff] <-
            FALSE # put problem coefficient back into active set
          nActive <-
            length(which(setActive)) # determine new number of active/inactive coefficients
          nInActive <- length(which(!setActive))
          nIneqBorder <-
            length(which(setIneqBorder)) # determine number of active/binding inequality constraints


          # Recalculate derivative for coefficients & multipliers

          M <-
            cbind(H[setActive, setActive], t(matrix(Aeq[, setActive], ncol = nActive)), t(matrix(A[setIneqBorder, setActive], ncol =
                                                                                                   nActive)))
          M <-
            rbind(M, cbind(
              rbind(matrix(Aeq[, setActive], ncol = nActive), matrix(A[setIneqBorder, setActive], ncol =
                                                                       nActive)),
              matrix(0, m1 + nIneqBorder, m1 + nIneqBorder)
            ))

          error.hand.inv <- function(e) {
            dir <-
              -(ginv(M) %*% rbind(
                matrix(subgrad[setActive], nActive, 1),
                matrix(0, m1 + nIneqBorder, 1)
              ))
            printer("Moore-Penrose-Inverse used.")
            return(dir)
          }
          dir <-
            tryCatch(
              dirsgn * (solve(M, rbind(
                matrix(subgrad[setActive], nActive, 1),
                matrix(0, m1 + nIneqBorder, 1)
              ))),
              error = error.hand.inv
            )

          if (nInActive != 0) {
            dirSubgrad <-
              -cbind(matrix(H[!setActive, setActive], ncol = nActive), t(matrix(Aeq[, !setActive], ncol =
                                                                                  nInActive)), t(matrix(A[setIneqBorder, !setActive], ncol = nInActive))) %*%
              dir
          } else{
            dirSubgrad <- matrix(0, 0, m1 + nIneqBorder)
          }

          #check for violations again

          # Negative subgradient but positive derivative
          signMismatchNegIdx <-
            which((-1 - ceiling_tol) <= subgrad[setActive] &
                    subgrad[setActive] <= (0 + ceiling_tol) &
                    (0 + ceiling_tol) <= dirsgn * dir[1:nActive] &
                    betaPath[setActive, k - 1] == 0
            )

          # update violation counter
          violateCounter <- violateCounter + 1

          if (violateCounter >= maxiters) {
            printer("Too many violations.")
            break
          }
        }

        # update violation trackers to see if any issues persist
        # Negative subgradient
        inactSlowNegIdx <-
          which((1 * dirsgn - ceiling_tol) <= subgrad[!setActive] &
                  subgrad[!setActive] <= (1 * dirsgn + ceiling_tol) &
                  1 * dirsgn < dirSubgrad
          )
        # Positive subgradient
        inactSlowPosIdx <-
          which((-1 * dirsgn - ceiling_tol) <= subgrad[!setActive] &
                  subgrad[!setActive] <= (-1 * dirsgn + ceiling_tol) &
                  dirSubgrad < -1 * dirsgn
          )
        # Positive subgrad but negative derivative
        signMismatchPosIdx = which((0 - ceiling_tol) <= subgrad[setActive] &
                                     subgrad[setActive] <= (1 + ceiling_tol) &
                                     dirsgn * dir[1:nActive] <= (0 - ceiling_tol)  &
                                     betaPath[setActive, k - 1] == 0
        )
        # Negative subgradient but positive derivative
        signMismatchNegIdx <-
          which((-1 - ceiling_tol) <= subgrad[setActive] &
                  subgrad[setActive] <= (0 + ceiling_tol) &
                  (0 + ceiling_tol) <= dirsgn * dir[1:nActive] &
                  betaPath[setActive, k - 1] == 0
          )

        if (violateCounter >= maxiters) {
          printer("Too many violations.")
          break
        }

      } ### end of violation check while-loop

      # store number of violations
      violationsPath[k] <- violateCounter

      # calculate derivative for residual inequality
      dirResidIneq <-
        matrix(A[!setIneqBorder, setActive], nrow = length(which(!setIneqBorder)), ncol =
                 nActive) %*% dir[1:nActive]

      ### Determine lambda for next event (via delta lambda)

      nextlambdaBeta <- matrix(Inf, p, 1)

      ## Events based on changes in coefficient status

      # Active coefficient going inactive
      nextlambdaBeta[setActive, ] = -dirsgn * betaPath[setActive, k - 1] / dir[1:nActive]

      # Inactive coefficient becoming positive
      t1 <-
        dirsgn * lambdaPath[, k - 1] * (1 - subgrad[!setActive]) / (dirSubgrad - 1)
      # threshold values hitting ceiling
      t1[t1 <= (0 + ceiling_tol)] <- Inf

      # Inactive coefficient becoming negative
      t2 = -dirsgn * lambdaPath[, k - 1] * (1 + subgrad[!setActive]) / (dirSubgrad + 1)
      # threshold values hitting ceiling
      t2[t2 <= (0 + ceiling_tol)] <- Inf

      # choose smaller delta lambda out of t1 and t2
      nextlambdaBeta[!setActive, ] <- pmin(t1, t2)

      # ignore delta lambdas numerically equal to zero
      nextlambdaBeta[nextlambdaBeta <= ceiling_tol | !penidx, ] <- Inf

      ## Events based inequality constraints

      # clear previous values
      nextlambdaIneq <- matrix(Inf, m2, 1)

      # Inactive inequality constraint becoming active
      nextlambdaIneq[!setIneqBorder, ] <-
        as.matrix(-dirsgn * residIneq[!setIneqBorder], length(which(!setIneqBorder)), 1) /
        (as.matrix(dirResidIneq, length(which(!setIneqBorder)), 1))

      # Active inequality constraint becoming deactive #
      nextlambdaIneq[setIneqBorder, ] <-
        as.matrix(-dirsgn * dualpathIneq[setIneqBorder, k - 1]) / as.matrix(dir[-(1:(nActive +
                                                                                       m1))], nIneqBorder, 1)

      # ignore delta lambdas equal to zero
      nextlambdaIneq[nextlambdaIneq <= ceiling_tol, ] <- Inf

      # find smallest lambda
      chglambda <-
        min(rbind(nextlambdaBeta, nextlambdaIneq), na.rm = TRUE)

      # find all indices corresponding to this chglambda
      idx <-
        which((rbind(nextlambdaBeta, nextlambdaIneq) - chglambda) <= ceiling_tol)

      # terminate path following if no new event found
      if (is.infinite(chglambda)) {
        chglambda <- lambdaPath[, k - 1]
      }

      # Update values at new lambda and move to next lambda, make sure isnt negative

      if ((lambdaPath[, k - 1] + dirsgn * chglambda) < 0) {
        chglambda <- lambdaPath[, k - 1]
      }

      printer(chglambda)
      ## calculate new value of lambda
      lambdaPath[, k] = lambdaPath[, k - 1] + dirsgn * chglambda

      ## Update parameter and subgradient values

      # new coefficient estimates
      betaPath[setActive, k] <-
        betaPath[setActive, k - 1] + dirsgn * chglambda * dir[1:nActive]

      # force near-zero coefficients to be zero (helps with numerical issues)
      betaPath[abs(betaPath[, k]) < zeros_tol, k]  <- 0

      # new subgradient estimates
      subgrad[!setActive] <-
        as.matrix(lambdaPath[, k - 1] * subgrad[!setActive, ] + dirsgn * chglambda *
                    dirSubgrad) / lambdaPath[, k]

      ## Update dual variables

      # update duals1 (lagrange multipliers for equality constraints)
      dualpathEq[, k] <-
        dualpathEq[, k - 1] + dirsgn * chglambda * as.matrix(dir[(nActive + 1):(nActive +
                                                                                  m1)], m1, 1)

      # update duals2 (lagrange multipliers for inequality constraints)
      dualpathIneq[setIneqBorder, k] <-
        dualpathIneq[setIneqBorder, k - 1] + dirsgn * chglambda * as.matrix(dir[(nActive +
                                                                                   m1 + 1):nrow(dir)], nIneqBorder, 1)

      # update residual inequality
      residIneq <- A %*% betaPath[, k] - b

      # update sets

      for (j in 1:length(idx)) {
        curidx <- idx[j]
        if (curidx <= p && setActive[curidx]) {
          # an active coefficient hits 0, or
          setActive[curidx] <- FALSE
        } else if (curidx <= p && !setActive[curidx]) {
          # a zero coefficient becomes nonzero
          setActive[curidx] <- TRUE
        } else if (curidx > p) {
          # an ineq on boundary becomes strict, or a strict ineq hits boundary
          setIneqBorder[curidx - p] = !setIneqBorder[curidx - p]
        }
      }

      # determine new number of active coefficients
      nActive <- length(which(setActive))
      nInActive <- length(which(!setActive))

      # determine number of active/binding inequality constraints
      nIneqBorder <- length(which(setIneqBorder))

      ## Calcuate and store values of interest along the path

      #calculate value of objective function
      objValPath[k] = 0.5 * (sum((y - X %*% betaPath[, k]) ^ 2)) + lambdaPath[k] *
        sum(abs(betaPath[, k]))

      # calculate degrees of freedom
      dfPath[k] = nActive - rankAeq - nIneqBorder

      # break algorithm when df are exhausted
      if (dfPath[k - 1] >= n_orig) {
        printer("BREAK. No more degrees of freedom.")
        break
      }

      if (length(which(setActive)) == p) {
        printer("BREAK. All coefficients active. No further Sparsity.")
        break
      }

    } #end overall for loop

    betaPath <- betaPath[, 1:(k - 1)]
    lambdaPath <- lambdaPath[1:(k - 1)]
    objValPath <- objValPath[1:(k - 1)]
    dfPath <- dfPath[1:(k - 1)]
    dfPath[dfPath < 0] <- 0

    return(
      list(
        "betaPath" = betaPath,
        "lambdaPath" = lambdaPath,
        "objValPath" = objValPath,
        "dfPath" = dfPath
      )
    )
  }
