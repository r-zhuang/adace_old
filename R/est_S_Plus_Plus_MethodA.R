#' Estimate the treatment effects for population S_++ using Method A
#'
#' @description
#'  The est_S_Plus_Plus_MethodA function produces estimation of treatment
#'  effects for the population that can adhere to both treatments (S_++).
#'  This method (Method A) is based on the potential outcome under the
#'  hypothetical alternative treatment .
#'
#' @param X Matrix of baseline variables. Each row contains the baseline values
#' for each patient.
#' @param A Matrix of indicator for adherence. Each row of A contains the
#' adherence information for each patient.
#' Each column contains the adherence indicator after each intermediate time
#' point.
#' A = 1 means adherence
#' and A = 0 means non-adherence. Monotone missing is assumed.
#' @param Z List of matrices. Intermediate efficacy and safety outcomes that can
#' affect the probability of adherence. For each matrix, the structure is the
#' same as variable X.
#' @param Y Numeric vector of the final outcome (E.g., primary endpoint).
#' @param T Numeric vector of treatment assignment. T = 0 for the control group
#' and T = 1 for the experimental treatment group.
#'
#' @return A list containing the following components:
#'   \item{trt_diff}{Estimate of treatment difference for S_{++} using Method A}
#'   \item{se}{Estimated standard error}
#'   \item{res1}{Estimated mean for the treatment group}
#'   \item{res0}{Estimated mean for the control group}
#'   \item{se_res1}{Estimated standard error for the treatment group}
#'   \item{se_res0}{Estimated standard error for the control group}
#'
#' @details
#' The average treatment difference can be denoted as
#'
#' \deqn{latex}{\mu_{d,++} = E\{Y(1)-Y(0)|A(0) = 1, A(1) = 1\}}
#'
#' The method A exploits the joint distribution of X, Z, and Y by creating a
#' "virtual twin" of the patient from the assigned treatment and estimate the
#' potential outcome of that patient on the alternative treatment for
#' comparison. The variance estimation for the treatment effect is constructed
#' using the sandwich method. Details can be found in the references.
#'
#' The intermediate post-baseline measurements for each intermediate time point
#' are estimated by regressing Z on X
#' using subjects with experimental treatment or placebo. The covariance matrix
#' is estimated based on the residuals of the regression.
#'
#' The probability of adherence is estimated by
#' regressing A on X, Z by using all data. The logistic regression is used in
#' this function.
#'
#' The expected treatment effect is estimated by
#' regression Y on X, Z using subjects with experimental treatment or placebo.
#'
#' The indicator of adherence prior to the first intermediate time point is not
#' included in this model since this function assumes no intercurrent events
#' prior to the first time point. Thus,  the first element of Z should not have
#' missing values.
#'
#' Each element of Z contains the variables at each intermediate time point,
#' i.e., the first element of Z contains the intermediate variables at time
#' point 1, the second element contains the intermediate variables at time point
#'  2, etc.
#'
#' @references
#' Qu, Yongming, et al. "A general framework for treatment effect estimators
#' considering patient adherence."
#' Statistics in Biopharmaceutical Research 12.1 (2020): 1-18.
#'
#' Zhang, Ying, et al. "Statistical inference on the estimators of the adherer
#' average causal effect."
#' Statistics in Biopharmaceutical Research (2021): 1-4.
#'
#' @examples
#'  library(MASS)
#'  n = 1000
#'  alpha1 = matrix(rep(c(2.3, -0.3, -0.01, 0.02, 0.03, 0.04, -0.4),3),ncol=3)
#'  alpha2 = matrix(rep(c(2.3, -0.3, -0.01, 0.02, 0.03, 0.04, -0.9),4),ncol=4)
#'  alpha3 = matrix(rep(c(2.3, -0.3, -0.01, 0.02, 0.03, 0.04, -0.9),5),ncol=5)
#'  beta = c(0.2, -0.3, -0.01, 0.02, 0.03, 0.04, rep(0.02,3),
#'  rep(0.04,4), rep(0.07,5))
#'  beta_T = -0.2
#'  gamma1 = c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-1/3,3))     #setting 1
#'  gamma2 = c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-2/4,4))     #setting 1
#'  gamma3 = c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-2.5/5,5))   #setting 1
#'  sd_z_x = 0.4
#'  X = mvrnorm(n, mu=c(1,5,6,7,8), Sigma=diag(1,5))
#'  T = rbinom(n, size = 1,  prob = 0.5)
#'  Z0_1 = alpha1[1,]+(X%*%alpha1[2:6,])           + mvrnorm(n, mu = rep(0,3),
#'  Sigma = diag(sd_z_x,3))
#'  Z1_1 = alpha1[1,]+(X%*%alpha1[2:6,])+alpha1[7,] + mvrnorm(n, mu = rep(0,3),
#'  Sigma = diag(sd_z_x,3))
#'  Z_1  = Z1_1 * T+Z0_1 * (1-T)
#'
#'  Z0_2 = alpha2[1,]+(X%*%alpha2[2:6,])           + mvrnorm(n, mu = rep(0,4),
#'  Sigma = diag(sd_z_x,4))
#'  Z1_2 = alpha2[1,]+(X%*%alpha2[2:6,])+alpha2[7,] + mvrnorm(n, mu = rep(0,4),
#'  Sigma = diag(sd_z_x,4))
#'  Z_2  = Z1_2 * T + Z0_2 * (1-T)
#'
#'  Z0_3 = alpha3[1,]+(X%*%alpha3[2:6,])           + mvrnorm(n, mu = rep(0,5),
#'  Sigma = diag(sd_z_x,5))
#'  Z1_3 = alpha3[1,]+(X%*%alpha3[2:6,])+alpha3[7,] + mvrnorm(n, mu = rep(0,5),
#'  Sigma = diag(sd_z_x,5))
#'  Z_3  = Z1_3 * T + Z0_3 * (1-T)
#'  Z = list(Z_1, Z_2, Z_3)
#'  Y0 = (beta[1]+(X %*% beta[2:6]) + Z0_1 %*% matrix(beta[7:9], ncol = 1) +
#'  Z0_2 %*% matrix(beta[10:13], ncol = 1) + Z0_3 %*% beta[14:18] +
#'  rnorm(n, mean = 0, sd = 0.3))[,1]
#'  Y1 = (beta[1] + (X %*% beta[2:6]) + Z1_1 %*% matrix(beta[7:9], ncol = 1) +
#'  Z1_2 %*% matrix(beta[10:13], ncol = 1) + Z1_3 %*% beta[14:18] + beta_T +
#'  rnorm(n, mean = 0, sd = 0.3))[,1]
#'  Y  = Y1 * T + Y0 * (1 - T)
#'
#'  A0_1 = rbinom(n, size = 1, prob = 1 / (1 + exp(-(gamma1[1] +
#'  (X %*% gamma1[2:6]) + Z0_1 %*% matrix(gamma1[7:9], ncol = 1))[,1])))
#'  A1_1 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma1[1] +
#'  (X %*% gamma1[2:6]) + Z1_1 %*% matrix(gamma1[7:9], ncol = 1))[,1])))
#'  A_1  = A1_1 * T + A0_1 * (1 - T)
#'
#'  A0_2 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma2[1] +
#'  (X %*% gamma2[2:6]) + Z0_2 %*% matrix(gamma2[7:10], ncol = 1))[,1]))) * A0_1
#'  A1_2 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma2[1] +
#'  (X %*% gamma2[2:6]) + Z1_2 %*% matrix(gamma2[7:10], ncol = 1))[,1]))) * A1_1
#'  A_2  = A1_2 * T + A0_2 * (1 - T)
#'
#'  A0_3 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma3[1] +
#'  (X %*% gamma3[2:6]) + Z0_3 %*% matrix(gamma3[7:11], ncol = 1))[,1]))) * A0_2
#'  A1_3 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma3[1] +
#'  (X %*% gamma3[2:6]) + Z1_3 %*% matrix(gamma3[7:11], ncol = 1))[,1]))) * A1_2
#'  A_3  = A1_3 * T + A0_3 * (1 - T)
#'  A = cbind(A_1, A_2, A_3)
#'
#'  Z[[2]][A_1 == 0] <- NA
#'  Z[[3]][A_2 == 0] <- NA
#'  Y[A_3 == 0]   <- NA
#'  # estimate the treatment difference
#'  fit <- est_S_Plus_Plus_MethodA(X, A, Z, Y, T)
#'  fit
#'  # Calculate the true values
#'  true1 =  mean(Y1[A1_3==1 & A0_3==1])
#'  true1
#'  true0 =  mean(Y0[A1_3==1 & A0_3==1])
#'  true0
#'  true_d  =  true1 - true0
#'  true_d
#' @export
est_S_Plus_Plus_MethodA <- function(X, A, Z, Y, T) {# nolint
  # res is estimator for S_{++}, method A
  # se  is the plug-in variance estimator for S_{++}, method A
  res1 <-  res0 <-  res <-  NULL
  se1 <-  se0 <-  se <-  NULL
    Y <- as.numeric(Y) # nolint
    T <- as.numeric(T) # nolint
    X <- matrix(X, nrow = length(Y)) # nolint
    A <- matrix(A, nrow = length(Y)) # nolint
    n <- dim(X)[1] # nolint
    nX <- dim(X)[2] # nolint
    n_time_points <- length(Z)
    Z <- lapply(Z, as.matrix) # nolint
    # Standardize X
    #X <- apply(X, 2, function(x) {(x-mean(x))/sd(x)}) # nolint
    # Provide column names for X, Z, A
    x_col_names <- paste("X_", 1 : dim(X)[2], sep = "")
    a_col_names <- paste("A_", 1 : dim(A)[2], sep = "")
    colnames(X) <- x_col_names
    colnames(A) <- a_col_names
    for (i in 1:n_time_points) {
      colnames(Z[[i]]) <- paste("Z_", i, "_", 1 : dim(Z[[i]])[2], sep = "")
    }
    z_cbind <- do.call(cbind, Z) # nolint
    data <- data.frame(X, T, z_cbind, Y, A)
    # Model adherence given X, Z using a logistic model
    # fit parametric model for A
    form1 <- formula(paste("A_1 ~ ",
                           paste(c(x_col_names,
                                   colnames(z_cbind)[grep("Z_1",
                                                          colnames(z_cbind))]),
                                 collapse = " + ")))
    models_A_XZ <- list() # nolint
    models_A_XZ[[1]] <- glm(form1, family = "binomial", data = data,
                            control = list(maxit = 50))
    for (i in 2:n_time_points) {
      Z_id <- paste("Z_", i, sep = "") # nolint
      form <- as.formula(paste(a_col_names[i],
                               paste(c(x_col_names,
                                       colnames(z_cbind)[grep(Z_id,
                                                        colnames(z_cbind))]),
                                     collapse = "+"), sep = " ~ "))
      models_A_XZ[[i]] <- glm(form, family = "binomial",
                              data = data[A[, (i - 1)] == 1, ],
                              control = list(maxit = 50))
    }
    coefs_A_XZ <- list() # nolint
    preds_A_XZ <- list() # nolint
    for (i in 1:n_time_points) {
      coefs_A_XZ[[i]] <- c(models_A_XZ[[i]]$coef)
      preds_A_XZ[[i]] <- pmax(pmin(predict(models_A_XZ[[i]], newdata = data,
                                           type = "response"), 0.99), 0.01)
    }

    # Model Z given X using the linear model
    models_Z_X <- list() # nolint
    cov_Z_X <- list() # nolint

    for (i in 1:n_time_points) {
      Z_id <- paste("Z_", i, sep = "") # nolint
      form <- paste("cbind(", paste(colnames(z_cbind)[grep(Z_id,
                    colnames(z_cbind))], collapse = ","), ")~",
                    paste(c(x_col_names, "T"), collapse = "+"), sep = "")
      models_Z_X[[i]] <- lm(form, data = data.frame(data))
      if (dim(Z[[i]])[2] > 1) {
        cov_Z_X[[i]] <- cov(models_Z_X[[i]]$residuals)
        # More robust against model misspecification, but will still
        #increase the variance of the point estimate
      } else {
        cov_Z_X[[i]] <- matrix(var(models_Z_X[[i]]$residuals), 1, 1)
      }
    }
    # Predict Z using models_Z_X
    Zs1_pred <- list() # nolint
    for (i in 1:n_time_points) {
      Zs1_pred[[i]] <- predict(models_Z_X[[i]],
                              newdata = data.frame(X, T = rep(1, nrow(X))))
    }
    # Estimate alpha of Z_i given T for 1 and Xs
    coefs_Z_1X <- list() # nolint
    for (i in 1:n_time_points) {
      coef_time_raw <- matrix(coef(models_Z_X[[i]]), ncol = dim(Z[[i]])[2])
      coef_intercept <- coef_time_raw[1, ] + tail(coef_time_raw, n = 1)
      coef_Xs <- coef_time_raw[-c(1, dim(coef_time_raw)[1]), , drop = FALSE] # nolint
      coefs_Z_1X[[i]] <- rbind(coef_intercept, coef_Xs)
    }

    # fit parametric model for Y
    model_y1_X_Z1 <- paste("Y~", paste("X_", 1 : dim(X)[2], sep = "", # nolint
                                       collapse = "+"), "+",
                           paste(colnames(z_cbind), sep = "",
                                 collapse = "+"), sep = "")
    fit_y1_X_Z1 <- lm(model_y1_X_Z1, data = data[T == 1 &             # nolint
                                                   A[, n_time_points] == 1, ])
    psi1_X_Z1 <- predict(fit_y1_X_Z1, newdata = data) # nolint
    beta1_hat <- c(fit_y1_X_Z1$coef) # nolint
    # Calculate expectations
    Expect_res <- apply(X, 1, Expect_function1D_MA, # nolint
                        n_time_points = n_time_points,
                        gammas = coefs_A_XZ, alphas = coefs_Z_1X,
                        Sigmas = cov_Z_X)
    names_vec <- c("prob1", paste("expz_1_", 1 : dim(Z[[1]])[2], sep = ""),
                   paste("expz^2_1_", 1 : (dim(Z[[1]])[2])^2, sep = ""),
                   "expa1", paste("expaz_1_", 1 : dim(Z[[1]])[2], sep = ""),
                   paste("expaz^2_1_", 1:(dim(Z[[1]])[2])^2, sep = ""))
    for (i in 2:n_time_points) {
      names_temp <- c(paste("prob", i, sep = ""),
                      paste("expz_", i, "_", 1:dim(Z[[i]])[2], sep = ""),
                      paste("expz^2_", i, "_", 1:(dim(Z[[i]])[2])^2, sep = ""),
                      paste("expa", i, sep = ""),
                      paste("expaz_", i, "_", 1:dim(Z[[i]])[2], sep = ""),
                      paste("expaz^2_", i, "_", 1:(dim(Z[[i]])[2])^2, sep = ""))
      names_vec <- c(names_vec, names_temp)
    }
    rownames(Expect_res) <- names_vec

    # Calculate functions of the expectations
    probs_vec <- paste("prob", 1:n_time_points, sep = "")
    Expect_probs <- Expect_res[probs_vec, , drop = FALSE] # nolint
    Expect_A_X <- 1 # nolint

    for (i in 1:n_time_points) {
      Expect_A_X <- Expect_A_X * Expect_probs[i, ] # nolint
    }
    Expect_AY_X <- 0 # nolint
    Expect_AZ_X <- list() # nolint
    Expect_AYZ_X <- list() # nolint
    Expect_AA_X <- list() # nolint
    Expect_AA_Z_X <- list() # nolint
    for (i in 1:n_time_points) {
      prob_remove <- paste("prob", i, sep = "") # nolint
      Expect_AY_X <- Expect_AY_X + Expect_A_X * beta1_hat[paste("Z_", i, "_", # nolint
                                                              1:dim(Z[[i]])[2],
                                                                sep = "")] %*%
        Expect_res[paste("expz_", i, "_", 1:dim(Z[[i]])[2],
            sep = ""), , drop = FALSE] / Expect_res[prob_remove, , drop = FALSE]
      Expect_AZ_X[[i]] <- t(Expect_A_X * t(Expect_res[paste("expz_", i, "_",
                           1:dim(Z[[i]])[2], sep = ""), , drop = FALSE]) /
                              Expect_res[prob_remove, ])
      Expect_AA_X[[i]] <- t(Expect_A_X * t(Expect_res[paste("expa", i,
                              sep = ""), , drop = FALSE]) /
                              Expect_res[prob_remove, ])
      Expect_AA_Z_X[[i]] <- t(Expect_A_X * t(Expect_res[paste("expaz_", i, "_",
                              1:dim(Z[[i]])[2], sep = ""), , drop = FALSE]) /
                                Expect_res[prob_remove, ])
    }
    Expect_AY_X <- as.numeric(cbind(rep(1, n), X) %*% beta1_hat[1:(nX + 1)]) *  # nolint
      Expect_A_X + as.numeric(Expect_AY_X)
# Calculate Expect_AYZ_X
    Expect_AYZ_X <- expect_AYZ_X_MApp(Z, X, n, nX, n_time_points, Expect_res,   # nolint
                                      beta1_hat, Expect_A_X, Expect_AZ_X)
# Calculate Expect_AA_Y_X
    Expect_AA_Y_X <- expect_AA_Y_X_MApp(Z, X, n, nX, n_time_points, Expect_res, # nolint
                                        beta1_hat, Expect_A_X, Expect_AA_X)
# Calculate Expect_AA_YZ_X
    Expect_AA_YZ_X <- expect_AA_YZ_X_MApp(Z, X, n, nX, n_time_points,           # nolint
                                          Expect_res, beta1_hat,
                                          Expect_A_X, Expect_AA_Z_X)


    #-- End: calculate integration -------------------------------------------

    # Calcuate estimator for first term
    res1 <-  sum((1 - T) * A[, n_time_points] * Expect_AY_X) /
      sum((1 - T) * A[, n_time_points] * Expect_A_X)
    # Prepare for plug-in variance estimator
    # g1,g2,g3, g4-res1*g5 are estimating equations
    g1 <- T * A[, n_time_points] * (Y - psi1_X_Z1) *
      cbind(rep(1, n), X, z_cbind)
    g1[A[, n_time_points] == 0, ] <- 0
    covariate_long <- lapply(models_Z_X, model.matrix.lm)
    g2 <- list()
    for (i in 1:n_time_points) {
      g2[[i]] <- array(apply(models_Z_X[[i]]$model[[1]] -
                               models_Z_X[[i]]$fitted.values,
                             2, function(x) {covariate_long[[i]] * x}), # nolint
                       dim = c(dim(covariate_long[[i]]),
                               dim(models_Z_X[[i]]$model[[1]])[2]))
    }
    preds_A_XZ_clean <- lapply(preds_A_XZ, NA_replace) # nolint
    g3 <- list()
    for (i in 1:n_time_points) {
      if (i == 1) {
        g3[[i]] <- (A[, i] - preds_A_XZ_clean[[i]]) *
          cbind(rep(1, n), X, Z[[i]])
      } else {
        g3[[i]] <- A[, (i - 1)] * (A[, i] - preds_A_XZ_clean[[i]]) *
          cbind(rep(1, n), X, Z[[i]])
        g3[[i]][A[, (i - 1)] == 0, ] <- 0
      }
    }

    g4 <- (1 - T) * A[, n_time_points] * Expect_AY_X
    g5 <- (1 - T) * A[, n_time_points] * Expect_A_X


    # Estimating expection of deriatives
    partial_g4_beta <- c(mean((1 - T) * A[, n_time_points] * Expect_A_X),
                  apply(X, 2, function(x){mean((1 - T) *   # nolint
                  A[, n_time_points] * Expect_A_X * x)}),
                  unlist(sapply(Expect_AZ_X,
                  function(x){colMeans((1 - T) * A[, n_time_points] * t(x))}))) # nolint

    partial_g4_alpha <- list()
    part1_partial <- (1 - T) * A[, n_time_points]
    for (i in 1:n_time_points) {
      partial_vec <- t(part1_partial * t(Expect_AYZ_X[[i]] -
                                           t(Expect_AY_X * Zs1_pred[[i]])))
      partial_g4_alpha_z1 <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(partial_vec), ncol = 1)
      partial_g4_alpha_z1x <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(matrix(mapply(outer, split(partial_vec,
                                                   col(partial_vec)),
                                      split(X, row(X))), ncol = n)),
               dim(Z[[i]])[2], dim(X)[2])
      partial_g4_alpha_z1t <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(partial_vec), ncol = 1)
      partial_g4_alpha[[i]] <- cbind(partial_g4_alpha_z1,
                                     partial_g4_alpha_z1x,
                                     partial_g4_alpha_z1t)
    }

    partial_g4_gamma <- list()
    for (i in 1:n_time_points) {
      constant <- (1 - T) * A[, n_time_points]
      vec1 <- Expect_AA_Y_X[[i]]
      vec2 <- Expect_AA_YZ_X[[i]]
      partial_g4_gamma[[i]] <- c(mean(constant * vec1),
                                 colMeans(X * as.numeric(constant * vec1)),
                                 colMeans(constant * t(vec2)))
    }

    partial_g5_beta  <-  0

    partial_g5_alpha <- list()
    for (i in 1:n_time_points) {
      partial_vec <- t(part1_partial * t(Expect_AZ_X[[i]] -
                                           t(Expect_A_X * Zs1_pred[[i]])))
      partial_g5_alpha_z1 <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(partial_vec), ncol = 1)
      partial_g5_alpha_z1x <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(matrix(mapply(outer, split(partial_vec,
                        col(partial_vec)), split(X, row(X))), ncol = n)),
               dim(Z[[i]])[2], dim(X)[2])
      partial_g5_alpha_z1t <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(partial_vec), ncol = 1)
      partial_g5_alpha[[i]] <- cbind(partial_g5_alpha_z1,
                                     partial_g5_alpha_z1x,
                                     partial_g5_alpha_z1t)
    }

    partial_g5_gamma <- list()
    for (i in 1:n_time_points) {
      constant <- (1 - T) * A[, n_time_points]
      vec1 <- Expect_AA_X[[i]]
      vec2 <- Expect_AA_Z_X[[i]]
      partial_g5_gamma[[i]] <- c(mean(constant * vec1),
                                 colMeans(X * as.numeric(constant * vec1)),
                                 colMeans(constant * t(vec2)))
    }
    covariate <- cbind(rep(1, n), X, z_cbind)
    covariate_T1A1 <- covariate[T == 1 & A[, n_time_points] == 1, ] # nolint
    partial_g1_beta <- -t(covariate_T1A1) %*% covariate_T1A1 / n

    partial_g2_alpha <- list()
    for (i in 1:n_time_points) {
      partial_g2_alpha[[i]] <- -t(covariate_long[[i]]) %*%
        covariate_long[[i]] / n
    }

    partial_g3_gamma <- list()
    for (i in 1:n_time_points) {
      if (i == 1) {
        covariate_i <- cbind(rep(1, n), X, Z[[i]])
        partial_g3_gamma[[i]] <- -t(covariate_i) %*%
          (covariate_i * preds_A_XZ[[i]] * (1 - preds_A_XZ[[i]])) / n
      } else {
        covariate_i <- cbind(rep(1, n), X, Z[[i]])[A[, i - 1] != 0, ]
        A_pred_sub <- preds_A_XZ[[i]][A[, i - 1] != 0] # nolint
        partial_g3_gamma[[i]] <- -t(covariate_i) %*%
          (covariate_i * A_pred_sub * (1 - A_pred_sub)) / n
      }
    }
    # Calculate plug-in variance estimator for S_{++} first term
    tau_est1 <- res1

    se_main <- -g1 %*% inv_svd(partial_g1_beta) %*%
      (partial_g4_beta - tau_est1 * partial_g5_beta) / mean(g5)
    se_g2 <- se_g3 <- 0
    for (i in 1:n_time_points) {
      non_miss <- rowMeans(is.na(Z[[i]])) == 0
      for (j in 1:dim(g2[[i]])[3]) {
        temp2 <- matrix(0, nrow = n, ncol = ncol(g2[[i]]))
        temp2[non_miss, ] <- g2[[i]][, , j]
        se_g2 <- se_g2 - temp2 %*% solve(partial_g2_alpha[[i]]) %*%
          (partial_g4_alpha[[i]][j, ] -
             tau_est1 * partial_g5_alpha[[i]][j, ]) / mean(g5)
      }
      se_g3 <- se_g3 - g3[[i]] %*% solve(partial_g3_gamma[[i]]) %*%
        (partial_g4_gamma[[i]] - tau_est1 * partial_g5_gamma[[i]]) / mean(g5)
    }
    se1 <- as.numeric(se_main + se_g2 + se_g3) +
      as.numeric(1 / mean(g5) * (g4 - tau_est1 * g5))


    ##------------------------------------------------------------------------
    ##   Variance calculation for Second Term E(Y0|A^{(h)}==1)
    ##------------------------------------------------------------------------

    # fit parametric model for Y
    model_y0_X_Z0 <- paste("Y~", paste("X_", 1:dim(X)[2], sep = "",   # nolint
                     collapse = "+"), "+", paste(colnames(z_cbind), sep = "",
                                                 collapse = "+"), sep = "")
    fit_y0_X_Z0 <- lm(model_y0_X_Z0,                                 # nolint
                      data = data[T == 0 & A[, n_time_points] == 1, ])
    psi0_X_Z0 <- predict(fit_y0_X_Z0, newdata = data)    # nolint
    beta0_hat <- c(fit_y0_X_Z0$coef)   # nolint

    # Predict Z using models_Z_X
    Zs0_pred <- list() # nolint
    for (i in 1:n_time_points) {
      Zs0_pred[[i]] <- predict(models_Z_X[[i]],
                              newdata = data.frame(X, T = rep(0, nrow(X))))
    }

    # Estimate alpha of Z_i given T for 0 and Xs
    coefs_Z_0X <- list() # nolint
    for (i in 1:n_time_points) {
      coef_time_raw <- matrix(coef(models_Z_X[[i]]), ncol = dim(Z[[i]])[2])
      coef_intercept <- coef_time_raw[1, ]
      coef_Xs <- coef_time_raw[-c(1, dim(coef_time_raw)[1]), , drop = FALSE] # nolint
      coefs_Z_0X[[i]] <- rbind(coef_intercept, coef_Xs)
    }

    #-- calculate integration w.r.t the Second Term (Y0,Z0,A0, etc)------------
    Expect_res <- apply(X, 1, Expect_function1D_MA,                # nolint
                        n_time_points = n_time_points,
                        gammas = coefs_A_XZ,
                        alphas = coefs_Z_0X,
                        Sigmas = cov_Z_X)
    names_vec <- c("prob1", paste("expz_1_", 1:dim(Z[[1]])[2], sep = ""),
                   paste("expz^2_1_", 1:(dim(Z[[1]])[2])^2, sep = ""),
                   "expa1", paste("expaz_1_", 1:dim(Z[[1]])[2], sep = ""),
                   paste("expaz^2_1_", 1:(dim(Z[[1]])[2])^2, sep = ""))
    for (i in 2:n_time_points) {
      names_temp <- c(paste("prob", i, sep = ""),
                      paste("expz_", i, "_", 1:dim(Z[[i]])[2], sep = ""),
                      paste("expz^2_", i, "_", 1:(dim(Z[[i]])[2])^2, sep = ""),
                      paste("expa", i, sep = ""),
                      paste("expaz_", i, "_", 1:dim(Z[[i]])[2], sep = ""),
                      paste("expaz^2_", i, "_", 1:(dim(Z[[i]])[2])^2, sep = ""))
      names_vec <- c(names_vec, names_temp)
    }
    rownames(Expect_res) <- names_vec
    probs_vec <- paste("prob", 1:n_time_points, sep = "")
    Expect_probs <- Expect_res[probs_vec, ] # nolint

    Expect_A_X <- 1 # nolint

    for (i in 1:n_time_points) {
      Expect_A_X <- Expect_A_X * Expect_probs[i, ] # nolint
    }
    Expect_AY_X <- 0 # nolint
    Expect_AZ_X <- list() # nolint
    Expect_AA_X <- list() # nolint
    Expect_AA_Z_X <- list() # nolint
    for (i in 1:n_time_points) {
      prob_remove <- paste("prob", i, sep = "")
      Expect_AY_X <- Expect_AY_X + Expect_A_X *                       # nolint
        beta0_hat[paste("Z_", i, "_", 1:dim(Z[[i]])[2], sep = "")] %*%
        Expect_res[paste("expz_", i, "_", 1:dim(Z[[i]])[2],
                         sep = ""), , drop = FALSE] / Expect_res[prob_remove, ]
      Expect_AZ_X[[i]] <- t(Expect_A_X * t(Expect_res[paste("expz_", i, "_",
            1:dim(Z[[i]])[2], sep = ""), , drop = FALSE]) /
              Expect_res[prob_remove, ])
      Expect_AA_X[[i]] <- t(Expect_A_X * t(Expect_res[paste("expa",
            i, sep = ""), , drop = FALSE]) / Expect_res[prob_remove, ])
      Expect_AA_Z_X[[i]] <-  t(Expect_A_X * t(Expect_res[paste("expaz_", i, "_",
            1:dim(Z[[i]])[2], sep = ""), , drop = FALSE]) /
              Expect_res[prob_remove, ])
    }
    Expect_AY_X <- as.numeric(cbind(rep(1, n), X) %*% beta0_hat[1:(nX + 1)]) *    # nolint
      Expect_A_X + as.numeric(Expect_AY_X)
    # Calculate Expect_AYZ_X
    Expect_AYZ_X <- expect_AYZ_X_MApp(Z, X, n, nX, n_time_points, Expect_res,     # nolint
                                      beta0_hat, Expect_A_X, Expect_AZ_X)


####\int P(A=1|X,Z1)*(1-P(A=1|X,Z1)*E(Y1|Z,X)*f(Z|X) dZ
    # Calculate Expect_AA_Y_X
    Expect_AA_Y_X <- expect_AA_Y_X_MApp(Z, X, n, nX, n_time_points, Expect_res,  # nolint
                                        beta0_hat, Expect_A_X, Expect_AA_X)

    # Calculate Expect_AA_YZ_X
    Expect_AA_YZ_X <- expect_AA_YZ_X_MApp(Z, X, n, nX, n_time_points,            # nolint
                                          Expect_res, beta0_hat,
                                          Expect_A_X, Expect_AA_Z_X)
    #-- End: calculate integration --------------------------------------------

    # Calcuate estimator for the second term
    res0 <- sum(T * A[, n_time_points] * Expect_AY_X) /
      sum(T * A[, n_time_points] * Expect_A_X)
    # Prepare for plug-in variance estimator
    # g1, g4-res1*g5 are estimating equations
    g1 <- (1 - T) * A[, n_time_points] * (Y - psi0_X_Z0) *
      cbind(rep(1, n), X, z_cbind)
    g1[A[, n_time_points] == 0, ] <- 0
    g4 <- T * A[, n_time_points] * Expect_AY_X
    g5 <- T * A[, n_time_points] * Expect_A_X

    # Estimating expection of deriatives
    partial_g4_beta <- c(mean(T * A[, n_time_points] * Expect_A_X),
                         apply(X, 2, function(x) {mean(T * A[,n_time_points] *  # nolint
                                                       Expect_A_X * x)}),
                         unlist(sapply(Expect_AZ_X, function(x) {               # nolint
                           colMeans(T * A[, n_time_points] * t(x))})))

    partial_g4_alpha <- list()
    part1_partial <- T * A[, n_time_points]
    for (i in 1:n_time_points) {
      partial_vec <- t(part1_partial * t(Expect_AYZ_X[[i]] -
                                           t(Expect_AY_X * Zs0_pred[[i]])))
      partial_g4_alpha_z1 <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(partial_vec), ncol = 1)
      partial_g4_alpha_z1x <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(matrix(mapply(outer, split(partial_vec,
              col(partial_vec)), split(X, row(X))), ncol = n)), dim(Z[[i]])[2],
               dim(X)[2])
      partial_g4_alpha_z1t <- matrix(0, nrow = dim(Z[[i]])[2], ncol = 1)
      partial_g4_alpha[[i]] <- cbind(partial_g4_alpha_z1, partial_g4_alpha_z1x,
                                     partial_g4_alpha_z1t)
    }

    partial_g4_gamma <- list()
    for (i in 1:n_time_points) {
      constant <- T * A[, n_time_points]
      vec1 <- Expect_AA_Y_X[[i]]
      vec2 <- Expect_AA_YZ_X[[i]]
      partial_g4_gamma[[i]] <- c(mean(constant * vec1),
                                 colMeans(X * as.numeric(constant * vec1)),
                                 colMeans(constant * t(vec2)))
    }

    partial_g5_beta <- 0

    partial_g5_alpha <- list()
    for (i in 1:n_time_points) {
      partial_vec <- t(part1_partial * t(Expect_AZ_X[[i]] -
                                           t(Expect_A_X * Zs0_pred[[i]])))
      partial_g5_alpha_z1 <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(partial_vec), ncol = 1)
      partial_g5_alpha_z1x <- solve(cov_Z_X[[i]]) %*%
        matrix(rowMeans(matrix(mapply(outer, split(partial_vec,
              col(partial_vec)), split(X, row(X))), ncol = n)),
              dim(Z[[i]])[2], dim(X)[2])
      partial_g5_alpha_z1t <- matrix(0, nrow = dim(Z[[i]])[2], ncol = 1)
      partial_g5_alpha[[i]] <- cbind(partial_g5_alpha_z1,
                                     partial_g5_alpha_z1x,
                                     partial_g5_alpha_z1t)
    }

    partial_g5_gamma <- list()
    for (i in 1:n_time_points) {
      constant <- T * A[, n_time_points]
      vec1 <- Expect_AA_X[[i]]
      vec2 <- Expect_AA_Z_X[[i]]
      partial_g5_gamma[[i]] <- c(mean(constant * vec1),
                                 colMeans(X * as.numeric(constant * vec1)),
                                 colMeans(constant * t(vec2)))
    }

    covariate <- cbind(rep(1, n), X, z_cbind)
    covariate_T0A0 <- covariate[T == 0 & A[, n_time_points] == 1, ] # nolint
    partial_g1_beta <- -t(covariate_T0A0) %*% covariate_T0A0 / n

    # Calculate plug-in variance estimator for S_{++} second term
    tau_est0 <- res0

    se_main <- -g1 %*% inv_svd(partial_g1_beta)  %*% (partial_g4_beta -
                tau_est0 * partial_g5_beta) / mean(g5)
    se_g2 <- se_g3 <- 0
    for (i in 1:n_time_points) {
      non_miss <- rowMeans(is.na(Z[[i]])) == 0
      for (j in 1:dim(g2[[i]])[3]) {
        temp2 <- matrix(0, nrow = n, ncol = ncol(g2[[i]]))
        temp2[non_miss, ] <- g2[[i]][, , j]
        se_g2 <- se_g2 - temp2 %*% solve(partial_g2_alpha[[i]]) %*%
          (partial_g4_alpha[[i]][j, ] - tau_est0 *
             partial_g5_alpha[[i]][j, ]) / mean(g5)
      }
      se_g3 <- se_g3 - g3[[i]] %*% solve(partial_g3_gamma[[i]]) %*%
        (partial_g4_gamma[[i]] - tau_est0 * partial_g5_gamma[[i]]) / mean(g5)
    }
    se0 <- as.numeric(se_main + se_g2 + se_g3) +
      as.numeric(1 / mean(g5) * (g4 - tau_est0 * g5))


    #--------------------------------------------------------------------------
    res <-  res1 - res0
    se  <-  sd(se1 - se0) / sqrt(n)
    se_res1 <-  sd(se1) / sqrt(n)
    se_res0 <-  sd(se0) / sqrt(n)

  RVAL <- list(trt_diff = res,   # nolint
             se = se,
             res1 = res1,
             res0 = res0,
             se_res1 = se_res1,
             se_res0 = se_res0
             )

  return(RVAL)
}

expect_AYZ_X_MApp <- function(Z, X, n, nX, n_time_points, Expect_res, beta_hat, # nolint
                              Expect_A_X, Expect_AZ_X) {                        # nolint
  rval <- list()
  for (i in 1:n_time_points) {
    rval[[i]] <- 0
    for (j in 1:n_time_points) {
      if (j == i) {
        names_target <- paste("expz^2_", j, "_", 1:dim(Z[[j]])[2]^2, sep = "")
        prob_remove <- paste("prob", j, sep = "")
        Expect_res_removed <- Expect_res[prob_remove, ] # nolint
        target_mat <- Expect_res[names_target, , drop = FALSE]
        target_mat <-  apply(target_mat, 2, function(x) {matrix(x, # nolint
                                            dim(Z[[j]])[2], dim(Z[[j]])[2]) %*%
            matrix(beta_hat[paste("Z_", i, "_", 1:dim(Z[[j]])[2], sep = "")],
                   nrow = dim(Z[[j]])[2])
        })
        if (dim(Z[[j]])[2] == 1) {
          target_mat <- matrix(target_mat, nrow = 1)
        }
      } else {
        target_mat <- matrix(beta_hat[paste("Z_", j, "_", 1:dim(Z[[j]])[2],
                                      sep = "")], ncol = dim(Z[[j]])[2]) %*%
          Expect_res[paste("expz_", j, "_",
                           1:dim(Z[[j]])[2], sep = ""), , drop = FALSE]
        target_mat <- t(t(Expect_res[paste("expz_", i, "_", 1:dim(Z[[i]])[2],
                          sep = ""), , drop = FALSE]) * as.numeric(target_mat))
        prob_remove <- c(paste("prob", i, sep = ""),
                         paste("prob", j, sep = ""))
        Expect_res_removed <- apply(Expect_res[prob_remove, ,   # nolint
                                               drop = FALSE], 2, prod)
      }
      rval[[i]] <- rval[[i]] + t(t(target_mat) * Expect_A_X /
                                   Expect_res_removed)
    }
    rval[[i]] <- rval[[i]] + t(t(Expect_AZ_X[[i]]) *
                                 as.numeric(cbind(rep(1, n), X) %*%
                                              beta_hat[1:(nX + 1)]))
  }
  return(rval)
}


expect_AA_Y_X_MApp <- function(Z, X, n, nX, n_time_points, Expect_res, beta_hat, # nolint
                               Expect_A_X, Expect_AA_X) {                        # nolint
  rval <- list()
  for (i in 1:n_time_points) {
    rval[[i]] <- 0
    for (j in 1:n_time_points) {
      if (j == i) {
        prob_remove <- paste("prob", j, sep = "")
        Expect_res_removed <- Expect_res[prob_remove, ] # nolint
        target_mat <- matrix(beta_hat[paste("Z_", j, "_", 1:dim(Z[[j]])[2],
                                        sep = "")], ncol = dim(Z[[j]])[2]) %*%
          Expect_res[paste("expaz_", j, "_",
                           1:dim(Z[[j]])[2], sep = ""), , drop = FALSE]
      } else {
        target_mat <- Expect_res[paste("expa", i, sep = ""), , drop = FALSE] *
          matrix(beta_hat[paste("Z_", j, "_", 1:dim(Z[[j]])[2], sep = "")],
                 ncol = dim(Z[[j]])[2]) %*%
          Expect_res[paste("expz_", j, "_",
                           1:dim(Z[[j]])[2], sep = ""), , drop = FALSE]
        prob_remove <- c(paste("prob", i, sep = ""),
                         paste("prob", j, sep = ""))
        Expect_res_removed <- apply(Expect_res[prob_remove, ,     # nolint
                                               drop = FALSE], 2, prod)
      }
      rval[[i]] <- rval[[i]] +
        t(t(target_mat) * Expect_A_X / Expect_res_removed)
    }
    rval[[i]] <- rval[[i]] + Expect_AA_X[[i]] *
      as.numeric(cbind(rep(1, n), X) %*% beta_hat[1:(nX + 1)])
  }
  return(rval)
}

expect_AA_YZ_X_MApp <- function(Z, X, n, nX, n_time_points,                     # nolint
                                Expect_res, beta_hat,                           # nolint
                                Expect_A_X, Expect_AA_Z_X) {                    # nolint
  rval <- list()
  for (i in 1:n_time_points) {
    rval[[i]] <- 0
    for (j in 1:n_time_points) {
      if (j == i) {
        names_target <- paste("expaz^2_", j, "_", 1:dim(Z[[j]])[2]^2,
                              sep = "")
        prob_remove <- paste("prob", j, sep = "")
        Expect_res_removed <- Expect_res[prob_remove, ] # nolint
        target_mat <- Expect_res[names_target, , drop = FALSE]
        target_mat <- apply(target_mat, 2, function(x) {matrix(x,  # nolint
                                          dim(Z[[j]])[2], dim(Z[[j]])[2]) %*%
            matrix(beta_hat[paste("Z_", i, "_", 1:dim(Z[[j]])[2], sep = "")],
                   nrow = dim(Z[[j]])[2])
        })
        if (dim(Z[[j]])[2] == 1) {
          target_mat <- matrix(target_mat, nrow = 1)
        }
      } else {
        target_mat <- matrix(beta_hat[paste("Z_", j, "_", 1:dim(Z[[j]])[2],
                                        sep = "")], ncol = dim(Z[[j]])[2]) %*%
          Expect_res[paste("expz_", j, "_",
                           1:dim(Z[[j]])[2], sep = ""), , drop = FALSE]
        target_mat <- t(t(Expect_res[paste("expaz_", i, "_", 1:dim(Z[[i]])[2],
                          sep = ""), , drop = FALSE]) * as.numeric(target_mat))
        prob_remove <- c(paste("prob", i, sep = ""),  # nolint
                         paste("prob", j, sep = ""))
        Expect_res_removed <- apply(Expect_res[prob_remove, ,  # nolint
                                               drop = FALSE], 2, prod)
      }
      rval[[i]] <- rval[[i]] +
        t(t(target_mat) * Expect_A_X / Expect_res_removed)
    }
    rval[[i]] <- rval[[i]] +
      t(t(Expect_AA_Z_X[[i]]) *
          as.numeric(cbind(rep(1, n), X) %*% beta_hat[1:(nX + 1)]))
  }

  return(rval)
}
