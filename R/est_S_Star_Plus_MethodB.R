#' Estimate the treatment effects for population S_*+ using Method B
#'
#' @description
#'  The est_S_Star_Plus_MethodB function produces estimation of treatment
#'  effects for the population that can adhere to the experimental
#'  treatment (S_*+). This method (Method B) is based on the potential outcome
#'  under the hypothetical alternative treatment .
#'
#' @param X Matrix of baseline variables. Each row contains the baseline values
#' for each patient across multiple time points.
#' @param A Matrix of indicator for adherence. Each row  of A contains the
#' adherence information for each patient. Each column contains the adherence
#' indicator after each intermediate time point. A = 1 means adherence and A=0
#' means non-adherence. Monotone missing is assumed.
#' @param Z List of matrices. Intermediate efficacy and safety outcomes that can
#' affect the probability of adherence. For each matrix, the structure is the
#' same as variable X.
#' @param Y Numeric vector of the final outcome (E.g., primary endpoint).
#' @param TRT Numeric vector of treatment assignment. TRT=0 for the control
#' group and TRT =1 for the experimental treatment group.
#'
#' @return A list containing the following components:
#'   \item{trt_diff}{Estimate of treatment difference for S_{*+} using Method B}
#'   \item{se}{Estimated standard error}
#'   \item{res1}{Estimated mean for the treatment group}
#'   \item{res0}{Estimated mean for the control group}
#'   \item{se_res1}{Estimated standard error for the treatment group}
#'   \item{se_res0}{Estimated standard error for the control group}
#'
#'#' @details
#' The average treatment difference can be denoted as
#'
#' \deqn{latex}{\mu_{d,*+} = E\{Y(1)-Y(0)|A(1) = 1\}}
#'
#' The method B exploits the joint distribution of X, Z, and A to estimate the
#' probability that a patient would adhere to the hypothetical alternative
#' treatment, and then use IPW to estimate treatment different for a given
#' population. The variance estimation for the treatment effect is constructed
#' using the sandwich method. Details can be found in the references.
#'
#' The intermediate post-baseline measurements for each intermediate time point
#' are estimated by regressing Z on X using subjects with experimental treatment
#' or placebo. The covariance matrix is estimated based on the residuals of
#' the regression.
#'
#' The probability of adherence is estimated by
#' regressing A on X, Z by using all data. The logistic regression is used in
#' this function.
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
#' considering patient adherence." Statistics in Biopharmaceutical Research
#' 12.1 (2020): 1-18.
#'
#' Zhang, Ying, et al. "Statistical inference on the estimators of the adherer
#' average causal effect."
#' Statistics in Biopharmaceutical Research (2021): 1-4.
#'
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
#'  TRT = rbinom(n, size = 1,  prob = 0.5)
#'  Z0_1 = alpha1[1,]+(X%*%alpha1[2:6,])           + mvrnorm(n, mu = rep(0,3),
#'  Sigma = diag(sd_z_x,3))
#'  Z1_1 = alpha1[1,]+(X%*%alpha1[2:6,])+alpha1[7,] + mvrnorm(n, mu = rep(0,3),
#'  Sigma = diag(sd_z_x,3))
#'  Z_1  = Z1_1 * TRT+Z0_1 * (1-TRT)
#'
#'  Z0_2 = alpha2[1,]+(X%*%alpha2[2:6,])           + mvrnorm(n, mu = rep(0,4),
#'  Sigma = diag(sd_z_x,4))
#'  Z1_2 = alpha2[1,]+(X%*%alpha2[2:6,])+alpha2[7,] + mvrnorm(n, mu = rep(0,4),
#'  Sigma = diag(sd_z_x,4))
#'  Z_2  = Z1_2 * TRT + Z0_2 * (1-TRT)
#'
#'  Z0_3 = alpha3[1,]+(X%*%alpha3[2:6,])           + mvrnorm(n, mu = rep(0,5),
#'  Sigma = diag(sd_z_x,5))
#'  Z1_3 = alpha3[1,]+(X%*%alpha3[2:6,])+alpha3[7,] + mvrnorm(n, mu = rep(0,5),
#'  Sigma = diag(sd_z_x,5))
#'  Z_3  = Z1_3 * TRT + Z0_3 * (1-TRT)
#'  Z = list(Z_1, Z_2, Z_3)
#'  Y0 = (beta[1]+(X %*% beta[2:6]) + Z0_1 %*% matrix(beta[7:9], ncol = 1) +
#'  Z0_2 %*% matrix(beta[10:13], ncol = 1) + Z0_3 %*% beta[14:18] +
#'  rnorm(n, mean = 0, sd = 0.3))[,1]
#'  Y1 = (beta[1] + (X %*% beta[2:6]) + Z1_1 %*% matrix(beta[7:9], ncol = 1) +
#'  Z1_2 %*% matrix(beta[10:13], ncol = 1) + Z1_3 %*% beta[14:18] + beta_T +
#'  rnorm(n, mean = 0, sd = 0.3))[,1]
#'  Y  = Y1 * TRT + Y0 * (1 - TRT)
#'
#'  A0_1 = rbinom(n, size = 1, prob = 1 / (1 + exp(-(gamma1[1] +
#'  (X %*% gamma1[2:6]) + Z0_1 %*% matrix(gamma1[7:9], ncol = 1))[,1])))
#'  A1_1 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma1[1] +
#'  (X %*% gamma1[2:6]) + Z1_1 %*% matrix(gamma1[7:9], ncol = 1))[,1])))
#'  A_1  = A1_1 * TRT + A0_1 * (1 - TRT)
#'
#'  A0_2 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma2[1] +
#'  (X %*% gamma2[2:6]) + Z0_2 %*% matrix(gamma2[7:10], ncol = 1))[,1]))) * A0_1
#'  A1_2 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma2[1] +
#'  (X %*% gamma2[2:6]) + Z1_2 %*% matrix(gamma2[7:10], ncol = 1))[,1]))) * A1_1
#'  A_2  = A1_2 * TRT + A0_2 * (1 - TRT)
#'
#'  A0_3 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma3[1] +
#'  (X %*% gamma3[2:6]) + Z0_3 %*% matrix(gamma3[7:11], ncol = 1))[,1]))) * A0_2
#'  A1_3 = rbinom(n, size = 1, prob = 1/(1 + exp(-(gamma3[1] +
#'  (X %*% gamma3[2:6]) + Z1_3 %*% matrix(gamma3[7:11], ncol = 1))[,1]))) * A1_2
#'  A_3  = A1_3 * TRT + A0_3 * (1 - TRT)
#'  A = cbind(A_1, A_2, A_3)
#'
#'  Z[[2]][A_1 == 0] <- NA
#'  Z[[3]][A_2 == 0] <- NA
#'  Y[A_3 == 0]   <- NA
#'  # estimate the treatment difference
#'  fit <- est_S_Star_Plus_MethodA(X, A, Z, Y, TRT)
#'  fit
#'  # Calculate the true values
#'  true1 =  mean(Y1[A1_3 == 1])
#'  true1
#'  true0 =  mean(Y0[A1_3 == 1])
#'  true0
#'  true_d  =  true1 - true0
#'  true_d
#'
#' @export

est_S_Star_Plus_MethodB <- function(X, A, Z, Y, TRT){                  # nolint
  Y <- as.numeric(Y) # nolint
  TRT <- as.numeric(TRT) # nolint
  Y[is.na(Y)] <- 0
  X <- matrix(X, nrow = length(Y)) # nolint
  A <- matrix(A, nrow = length(Y)) # nolint
  n_time_points <- length(Z)
  n <- nrow(X)
  Z <- lapply(Z, as.matrix) # nolint
  Z_cp <- Z # nolint
  for (i in 2:n_time_points) {
    Z_cp[[i]][A[, (i - 1)] == 0] <- 0
  }
  Z_mat <- matrix(unlist(Z), nrow = dim(X)[1]) # nolint

  # Provide column names for X, Z, A
  X_col_names <- paste("X_", 1:dim(X)[2], sep = "")   # nolint
  A_col_names <- paste("A_", 1:dim(A)[2], sep = "")   # nolint
  colnames(X) <- X_col_names
  colnames(A) <- A_col_names
  Z_col_names <- list()                               # nolint
  for (i in 1:n_time_points) {
    part1 <- paste("Z_", i, sep = "")
    Z_col_names[[i]] <- paste(part1, 1:dim(Z[[i]])[2], sep = "")
  }
  colnames(Z_mat) <- unlist(Z_col_names)
  data <- data.frame(X, TRT, Z_mat, Y, A)

  # Model adherence given X, Z using a logistic model
  form1 <- formula(paste("A_1 ~ ",
                         paste(c(X_col_names, Z_col_names[[1]]),
                               collapse = " + ")))
  models_A_XZ <- list()                               # nolint
  models_A_XZ[[1]] <- glm(form1, family = "binomial", data = data)
  for (i in 2:n_time_points) {
    form <- as.formula(paste(A_col_names[i],
                             paste(c(X_col_names, Z_col_names[[i]]),
                                   collapse = "+"), sep = " ~ "))
    models_A_XZ[[i]] <- glm(form, family = "binomial",
                            data = data[A[, (i - 1)] == 1, ],
                            control = list(maxit = 50))
  }
  coefs_A_XZ <- list()  # nolint
  preds_A_XZ <- list()  # nolint
  for (i in 1:n_time_points) {
    coefs_A_XZ[[i]] <- c(models_A_XZ[[i]]$coef)
    preds_A_XZ[[i]] <- predict(models_A_XZ[[i]],
                               newdata = data, type = "response")
  }

  # Model Z given X and estimate variance-covariance of Z given X
  models_Z_X <- list()                                            # nolint
  sigma_mats_Z <- list()                                          # nolint
  for (i in 1:n_time_points) {
    models_Z_X[[i]] <- lm(Z_mat[, Z_col_names[[i]]] ~ X + TRT)
    if (dim(Z[[i]])[2] > 1) {
      sigma_mats_Z[[i]] <- diag(apply(models_Z_X[[i]]$residuals, 2, var))
    } else {
      sigma_mats_Z[[i]] <- matrix(var(models_Z_X[[i]]$residuals), 1, 1)
    }
  }

  # Predict Z using models.Z.X
  Zs_pred <- list()                                             # nolint
  for (i in 1:n_time_points) {

    Zs_pred[[i]] <- predict(models_Z_X[[i]],
                            newdata = data.frame(X, TRT = rep(1, n)))
  }

  # Estimate alpha of Z_i given T for 1 (intercept) and Xs
  coefs_Z_1X <- list()                                         # nolint
  for (i in 1:n_time_points) {
    coef_raw <- matrix(coef(models_Z_X[[i]]), ncol = dim(Z[[i]])[2])
    coefs_Z_1X[[i]] <- matrix(coef_raw[1:(dim(X)[2]+1), ],      # nolint
                              ncol = dim(Z[[i]])[2])
    coefs_Z_1X[[i]][1, ] <- coefs_Z_1X[[i]][1, ] + coef_raw[nrow(coef_raw), ]
  }
  # Calculate required integrations and assigning names for each column
  Expect_res <- apply(X, 1, Expect_function1D_BU, n_time_points = n_time_points,    # nolint
                      gammas = coefs_A_XZ, alphas = coefs_Z_1X,
                      Sigmas = sigma_mats_Z)
  Expect_res_t <- t(Expect_res)                                                     # nolint
  names_vec <- NULL
  for (i in 1:n_time_points) {
    z_dim <- dim(Z[[i]])[2]
    names_prob <- paste("prob", i, sep = "")
    names_expa <- paste("expa", i, sep = "")
    names_other <- paste(c("expz", "expaz"), i, sep = "")
    names_expz <- c(paste(names_other[1], 1:z_dim, sep = ""))
    names_expaz <- c(paste(names_other[2], 1:z_dim, sep = ""))
    names_vec <- c(names_vec, names_prob, names_expz,
                   names_expa, names_expaz)
  }
  colnames(Expect_res_t) <- names_vec

  # Calculate the probability of adherence at the end of study
  probs_vec <- paste("prob", 1:n_time_points, sep = "")
  Expect_probs <- Expect_res_t[, probs_vec]                         # nolint
  Expect_A_X <- 1                                                   # nolint
  for (i in 1:n_time_points) {
    Expect_A_X <- Expect_A_X * Expect_probs[, i]                    # nolint
  }

  # Calculate several integrals
  Expect_AZ_X <- list()   # nolint
  Expect_AA_X <- list()   # nolint
  Expect_AA_Z_X <- list()   # nolint
  for (i in 1:n_time_points) {
    names_target <- paste(c("expz", "expa", "expaz"), i, sep = "")
    prob_remove <- paste("prob", i, sep = "")
    Expect_AZ_X[[i]] <- Expect_A_X *
      Expect_res_t[, grepl(names_target[1], names_vec)] /
      Expect_res_t[, prob_remove]
    Expect_AA_X[[i]] <- Expect_A_X *
      Expect_res_t[, grepl(names_target[2], names_vec)] /
      Expect_res_t[, prob_remove]
    Expect_AA_Z_X[[i]] <-  Expect_A_X *
      Expect_res_t[, grepl(names_target[3], names_vec)] /
      Expect_res_t[, prob_remove]
  }

  # Calcuate a point estimation for S_tar_plus
  preds_A_XZ_clean <- lapply(preds_A_XZ, NA_replace)        # nolint
  prod_A <- 1                    # nolint
  for (i in 1:n_time_points) {
    prod_A <- prod_A * preds_A_XZ_clean[[i]]                 # nolint
  }
  res_main <- sum(TRT * A[, n_time_points] * Y - (1 - TRT) *
                    A[, n_time_points] *
                    Expect_A_X*Y / prod_A) / sum(TRT * A[, n_time_points])    # nolint
  residual <- (sum((1 - TRT) * A[, n_time_points] * Expect_A_X * Y / prod_A) /
                 sum(TRT * A[, n_time_points])) * (1 - sum(TRT) / sum(1 - TRT))
  estimator <- res_main + residual

  # Prepare for the plug-in variance estimator
  # Calcuate sample values of estimation equations

  # At each time point, for each dimension of Z, we have an estimation equation
  # g2 = X(Z_i - X\alpha), thus there will be layers of for loop.
  covariate_long <- lapply(models_Z_X, model.matrix.lm)
  g2 <- list()
  for (i in 1:n_time_points) {
    g2_main <- (models_Z_X[[i]]$model$`Z_mat[, Z_col_names[[i]]]`
                 - predict(models_Z_X[[i]]))
    if (dim(Z[[i]])[2] > 1) {
      dim_z <- dim(g2_main)[2]
      res <- list()
      for (j in 1:dim_z) {
        res[[j]] <- g2_main[, j] * covariate_long[[i]]
      }
      g2[[i]] <- res
    } else {
      g2[[i]] <- g2_main * covariate_long[[i]]
    }

  }

  # At each time point, g3 = (1, x, z)^T(A - Pr(A=1|x, z, \gamma))
  g3 <- list()
  for (i in 1:n_time_points) {
    if (i == 1) {
      g3[[i]] <- (A[, i] - preds_A_XZ_clean[[i]]) * cbind(rep(1, n), X, Z[[i]])
    } else {
      g3[[i]] <- A[, (i - 1)] * (A[, i] - preds_A_XZ_clean[[i]]) *
        cbind(rep(1, n), X, Z[[i]])
      g3[[i]][A[, (i - 1)] == 0, ] <- 0
    }

  }

  # g4 = (1 -T)h(x, \gamma, \alpha)/g(x, z, \gamma), here
  # we use an equivalent form of the above function
  g4 <- TRT * A[, n_time_points] * Y - (1 - TRT) * A[, n_time_points] *
    Expect_A_X * Y / prod_A

  # function g5 = TA
  g5 <- TRT * A[, n_time_points]

  # Estimating expection of deriatives using sample averages
  partial_g4_alpha <- list()
  part1_partial <- (1 - TRT) * A[, n_time_points] * Y / (prod_A)
  for (i in 1:n_time_points) {
    dim_z <- dim(Z[[i]])[2]
    partial_vec <- matrix(part1_partial * (Expect_AZ_X[[i]] - Expect_A_X *
                                             Zs_pred[[i]]), ncol = dim_z)
    partial_g4_alpha_1 <-  apply(partial_vec, 2, mean)
    partial_g4_alpha_x <- matrix(rep(NA, dim_z * dim(X)[2]), ncol = dim_z)
    for (j in 1:dim_z) {
      partial_g4_alpha_x[, j] <- apply(t(X) %*%
                                 diag(as.vector(partial_vec[, j])), 1, mean)
    }
    partial_g4_alpha[[i]] <- -solve(sigma_mats_Z[[i]]) %*%
      t(rbind(partial_g4_alpha_1, partial_g4_alpha_x, partial_g4_alpha_1))
  }

  partial_g4_gamma <- list()
  for (i in 1:n_time_points) {
    constant <- (1 - TRT) * A[, n_time_points] * Y /
      (prod_A / preds_A_XZ_clean[[i]])
    vec1 <- (Expect_AA_X[[i]] * preds_A_XZ_clean[[i]] -
               Expect_A_X * preds_A_XZ_clean[[i]] *
               (1 - preds_A_XZ_clean[[i]])) / (preds_A_XZ_clean[[i]]^2)
    vec2 <- (Expect_AA_X[[i]] * preds_A_XZ_clean[[i]] -
               Expect_A_X * preds_A_XZ_clean[[i]] *
               (1 - preds_A_XZ_clean[[i]])) / (preds_A_XZ_clean[[i]]^2)
    vec3 <- ((Expect_AA_Z_X[[i]] * preds_A_XZ_clean[[i]]) -
               Expect_A_X * preds_A_XZ_clean[[i]] *
              (1 - preds_A_XZ_clean[[i]]) * Z_cp[[i]]) /
              (preds_A_XZ_clean[[i]]^2)
    partial_g4_gamma[[i]] <- -c(mean(constant * vec1),
                                apply(t(X) %*%
                                diag(as.vector(constant * vec2)), 1, mean),
                                apply(constant * vec3, 2, mean))
  }

  partial_g2_alpha <- list()
  for (i in 1:n_time_points) {
    partial_g2_alpha[[i]] <- -t(covariate_long[[i]]) %*% covariate_long[[i]] / n
  }

  partial_g3_gamma <- list()
  for (i in 1:n_time_points) {
    if (i == 1) {
      design_mat <- cbind(rep(1, n), X, Z[[i]])
      partial_g3_gamma[[i]] <- -t(design_mat) %*% (
        design_mat * preds_A_XZ_clean[[i]] * (1 - preds_A_XZ_clean[[i]])) / n
    } else {
      design_mat <- cbind(rep(1, n), X, Z[[i]])
      design_mat_sub <- design_mat[A[, (i - 1)] != 0, ]
      A_pred_sub <- preds_A_XZ_clean[[i]][A[, (i - 1)] != 0]  # nolint
      partial_g3_gamma[[i]] <- -t(design_mat_sub) %*% (design_mat_sub *
                                A_pred_sub * (1 - A_pred_sub)) / n     # nolint
    }
  }

  # Calculate plug-in variance estimator for S_{*+} without residual term
  tau_est <- estimator
  data_long <- cbind(data, subj = seq_len(nrow(data)))
  data_long <- reshape2::melt(data_long, id.vars = c(X_col_names, "TRT", "Y",
                                                    A_col_names, "subj"),
                             variable.name = "Z")
  temp <- list()
  temp[[1]] <- g2[[1]]
  for (i in 2:n_time_points) {
    dim_z <- dim(Z[[i]])[2]
    res <- list()
    if (dim_z > 1) {
      for (j in 1:dim_z) {
        res[[j]] <- matrix(0, nrow = n, ncol = ncol(g2[[i]][[j]]))
        res[[j]][!is.na(data_long$value[data_long$Z == Z_col_names[[i]][j]]),
                 ] <- g2[[i]][[j]]
      }
      temp[[i]] <- res
    } else {
      temp[[i]] <- matrix(0, nrow = n, ncol = ncol(g2[[i]]))
      temp[[i]][!is.na(data_long$value[data_long$Z == Z_col_names[[i]]]),
                ] <- g2[[i]]
    }

  }

  se_main <- 0
  for (i in 1:n_time_points) {
    dim_z <- dim(Z[[i]])[2]
    if (dim_z > 1) {
      for (j in 1:dim_z) {
        se_main <- se_main - temp[[i]][[j]] %*% solve(partial_g2_alpha[[i]]) %*%
          (partial_g4_alpha[[i]][j, ]) / mean(g5)
      }
    } else {
      se_main <- se_main - temp[[i]] %*% solve(partial_g2_alpha[[i]]) %*%
        (t(partial_g4_alpha[[i]])) / mean(g5)
    }
  }

  for (i in 1:n_time_points) {
    se_main <- se_main - g3[[i]] %*% solve(partial_g3_gamma[[i]]) %*%
      partial_g4_gamma[[i]] / mean(g5)
  }
  se_main1 <- se_main + 1 / mean(g5) * (g4 - tau_est * g5)

  tau0_est <- sum((1 - TRT) * A[, n_time_points] * Expect_A_X * Y /
                    (prod_A)) / sum(TRT * A[, n_time_points]) - residual

  #use delta-method for se_residual
  se_residual <- - tau0_est * 1 / (1 - mean(TRT))^2 * (TRT - mean(TRT))
  se <- sd(se_main1 + se_residual) / sqrt(n)
  res1 <- sum(TRT * A[, n_time_points] * Y) / sum(TRT * A[, n_time_points])
  res0 <- sum((1 - TRT) * A[, n_time_points] * Expect_A_X * Y / prod_A) /
    sum(TRT * A[, n_time_points]) * (sum(TRT) / sum(1 - TRT))
  se1  <-  sd(1 / mean(g5) * (TRT * A[, n_time_points] * Y - res1 * g5)) /
    sqrt(n)

  se_main_se0 <- 0
  for (i in 1:n_time_points) {
    dim_z <- dim(Z[[i]])[2]
    if (dim_z > 1) {
      for (j in 1:dim_z) {
        se_main_se0 <- se_main_se0 + temp[[i]][[j]] %*%
          solve(partial_g2_alpha[[i]]) %*%
          (partial_g4_alpha[[i]][j, ]) / mean(g5)
      }
    } else {
      se_main_se0 <- se_main_se0 + temp[[i]] %*%
        solve(partial_g2_alpha[[i]]) %*% (t(partial_g4_alpha[[i]])) / mean(g5)
    }
  }

  for (i in 1:n_time_points) {
    se_main_se0 <- se_main_se0 + g3[[i]] %*% solve(partial_g3_gamma[[i]]) %*%
      partial_g4_gamma[[i]] / mean(g5)
  }

  se0 <- se_main_se0 + ((1 - TRT) * A[, n_time_points] * Expect_A_X * Y /
                          (prod_A) - res0 * g5 - se_residual) / mean(g5)
  se0 <- sd(se0) / sqrt(n)
  RVAL <- list(trt_diff = estimator,                                  # nolint
             se = se,
             res1 = res1,
             res0 = res0,
             se_res1 = se1,
             se_res0 = se0
  )
  return(RVAL)
}
