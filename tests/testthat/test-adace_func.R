test_that("adace main function", {
  #### main function should return a list of 6 elements ####
  library(MASS)
  library(pracma)
  set.seed(123)
  n <-  200
  alpha1 <- matrix(rep(c(2.3, -0.3, -0.01, 0.02,
                         0.03, 0.04, -0.4), 3), ncol = 3)
  alpha2 <- matrix(rep(c(2.3, -0.3, -0.01, 0.02,
                         0.03, 0.04, -0.9), 4), ncol = 4)
  alpha3 <- matrix(rep(c(2.3, -0.3, -0.01, 0.02,
                         0.03, 0.04, -0.9), 5), ncol = 5)
  beta <-  c(0.2, -0.3, -0.01, 0.02, 0.03, 0.04,
             rep(0.02, 3), rep(0.04, 4), rep(0.07, 5))
  beta_t <- -0.2                                               # nolint
  gamma1 <- c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-1 / 3, 3))     #setting 1
  gamma2 <- c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-2 / 4, 4))     #setting 1
  gamma3 <- c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-2.5 / 5, 5))   #setting 1
  sd_z_x <- 0.4
  X <- mvrnorm(n, mu = c(1, 5, 6, 7, 8), Sigma = diag(1, 5))   # nolint
  TRT <- rbinom(n, size = 1,  prob = 0.5)                        # nolint
  Z0_1 <- alpha1[1, ] + (X %*% alpha1[2:6, ]) +                # nolint
    mvrnorm(n, mu = rep(0, 3), Sigma = diag(sd_z_x, 3))
  Z1_1 <- alpha1[1, ] + (X %*% alpha1[2:6, ]) + alpha1[7, ] +           # nolint
    mvrnorm(n, mu = rep(0, 3), Sigma = diag(sd_z_x, 3))
  Z_1  <- Z1_1 * TRT + Z0_1 * (1 - TRT)                                     # nolint
  Z0_2 <- alpha2[1, ] + (X %*% alpha2[2:6, ]) +                         # nolint
    mvrnorm(n, mu = rep(0, 4), Sigma = diag(sd_z_x, 4))
  Z1_2 <- alpha2[1, ] + (X %*% alpha2[2:6, ]) + alpha2[7, ] +           # nolint
    mvrnorm(n, mu = rep(0, 4), Sigma = diag(sd_z_x, 4))
  Z_2  <- Z1_2 * TRT + Z0_2 * (1 - TRT)                                     # nolint

  Z0_3 <- alpha3[1, ] + (X %*% alpha3[2:6, ]) +                         # nolint
    mvrnorm(n, mu = rep(0, 5), Sigma = diag(sd_z_x, 5))
  Z1_3 <- alpha3[1, ] + (X %*% alpha3[2:6, ]) + alpha3[7, ] +           # nolint
    mvrnorm(n, mu = rep(0, 5), Sigma = diag(sd_z_x, 5))
  Z_3  <- Z1_3 * TRT + Z0_3 * (1 - TRT)                                     # nolint
  Z <- list(Z_1, Z_2, Z_3)                                              # nolint
  Y0 <- (beta[1] + (X %*% beta[2:6]) + Z0_1 %*% matrix(beta[7:9], ncol = 1) +   # nolint
           Z0_2 %*% matrix(beta[10:13], ncol = 1) + Z0_3 %*% beta[14:18] +
           rnorm(n, mean = 0, sd = 0.3))[, 1]
  Y1 <- (beta[1] + (X %*% beta[2:6]) + Z1_1 %*% matrix(beta[7:9], ncol = 1) +   # nolint
           Z1_2 %*% matrix(beta[10:13], ncol = 1) + Z1_3 %*% beta[14:18] +
           beta_t + rnorm(n, mean = 0, sd = 0.3))[, 1]
  Y <- Y1 * TRT + Y0 * (1 - TRT)                                            # nolint

  A0_1 <- rbinom(n, size = 1, prob = 1 / (1 + exp(-(gamma1[1] +                 # nolint
          (X %*% gamma1[2:6]) + Z0_1 %*% matrix(gamma1[7:9], ncol = 1))[, 1])))
  A1_1 <- rbinom(n, size = 1, prob = 1 / (1 + exp(-(gamma1[1] +                 # nolint
          (X %*% gamma1[2:6]) + Z1_1 %*% matrix(gamma1[7:9], ncol = 1))[, 1])))
  A_1 <- A1_1 * TRT + A0_1 * (1 - TRT)                                              # nolint

  A0_2 <- rbinom(n, size = 1, prob = 1 / (1 + exp(-(gamma2[1] +                 # nolint
        (X %*% gamma2[2:6]) + Z0_2 %*%
          matrix(gamma2[7:10], ncol = 1))[, 1]))) * A0_1
  A1_2 <- rbinom(n, size = 1, prob = 1 / (1 + exp(-(gamma2[1] +                 # nolint
        (X %*% gamma2[2:6]) + Z1_2 %*%
          matrix(gamma2[7:10], ncol = 1))[, 1]))) * A1_1
  A_2  <- A1_2 * TRT + A0_2 * (1 - TRT)                                             # nolint

  A0_3 <- rbinom(n, size = 1, prob = 1 / (1 + exp(-(gamma3[1] +                 # nolint
          (X %*% gamma3[2:6]) + Z0_3 %*%
            matrix(gamma3[7:11], ncol = 1))[, 1]))) * A0_2
  A1_3 <- rbinom(n, size = 1, prob = 1 / (1 + exp(-(gamma3[1] +                 # nolint
          (X %*% gamma3[2:6]) + Z1_3 %*%
            matrix(gamma3[7:11], ncol = 1))[, 1]))) * A1_2
  A_3 <- A1_3 * TRT + A0_3 * (1 - TRT)                                              # nolint
  A <- cbind(A_1, A_2, A_3)                                                     # nolint

  Z[[2]][A_1 == 0] <- NA
  Z[[3]][A_2 == 0] <- NA
  Y[A_3 == 0] <- NA
  # estimate the treatment difference
  fit1 <- est_S_Plus_Plus_MethodA(X, A, Z, Y, TRT)
  fit2 <- est_S_Star_Plus_MethodA(X, A, Z, Y, TRT)
  fit3 <- est_S_Plus_Plus_MethodB(X, A, Z, Y, TRT)
  fit4 <- est_S_Star_Plus_MethodB(X, A, Z, Y, TRT)
  # Test class
  expect_equal(class(fit1), "list")
  expect_equal(class(fit2), "list")
  expect_equal(class(fit3), "list")
  expect_equal(class(fit4), "list")
  expect_equal(sum(is.na(fit1)), 0)
  expect_equal(sum(is.na(fit2)), 0)
  expect_equal(sum(is.na(fit3)), 0)
  expect_equal(sum(is.na(fit4)), 0)

  # Test number of elements
  expect_equal(length(fit1), 6)
  expect_equal(length(fit2), 6)
  expect_equal(length(fit3), 6)
  expect_equal(length(fit4), 6)
  # Test warning message
  Y <- Y1 * TRT + Y0 * (1 - TRT)                                                    # nolint
  A_3[101:200] <- 0
  Y[A_3 == 0] <- NA
  message1 <- capture_message(fit1 <- est_S_Plus_Plus_MethodA(X, A, Z, Y, TRT))
  expect_equal(message1$message,
               "prediction from a rank-deficient fit may be misleading")
  message2 <- capture_message(fit2 <- est_S_Star_Plus_MethodA(X, A, Z, Y, TRT))
  expect_equal(message2$message,
               "prediction from a rank-deficient fit may be misleading")
  #### Test utility function ####
  A <- mat_vec(1:5)                                                                # nolint
  expect_equal(Rank(A), 5)
  expect_equal(NA_replace(rep(NA, 3)), c(999, 999, 999))
  expect_equal(sum(round(inv_svd(A) - solve(A), 5) == matrix(0, 5, 5)), 25)
  #### Test utitlity function for dim(Z)==1 ####
  al1 <- matrix(rep(c(2.3, -0.3, -0.01, 0.02, 0.03, 0.04), 1), ncol = 1)
  al2 <- matrix(rep(c(2.3, -0.3, -0.01, 0.02, 0.03, 0.04), 1), ncol = 1)
  al3 <- matrix(rep(c(2.3, -0.3, -0.01, 0.02, 0.03, 0.04), 1), ncol = 1)
  alps <- list(al1, al2, al3)
  be <- c(0.2, -0.3, -0.01, 0.02, 0.03, 0.04, rep(0.02, 1),
         rep(0.04, 1), rep(0.07, 1))
  gam1 <- c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-1 / 3, 1))     #setting 1
  gam2 <- c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-2 / 4, 1))     #setting 1
  gam3 <- c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(-2.5 / 5, 1))   #setting 1
  gams <- list(gam1, gam2, gam3)
  sd_z_x <- 0.4
  sigs <- list(diag(sd_z_x, 1), diag(sd_z_x, 1), diag(sd_z_x, 1))
  exp_res <- Expect_function1D_MA(X[1, , drop = FALSE], 3, gams, alps, sigs)
  expect_equal(class(exp_res), "numeric")
  expect_equal(sum(is.na(exp_res)), 0)

  #### Test main functions for dim(Z)==1 ####
  set.seed(12345)
  p_z <- 1
  n_t <- 3
  alphas <- list()
  gammas <- list()
  z_para <- c(-1 / p_z, -1 / p_z, -1 / p_z, -1 / p_z,
              -1.5 / p_z, -1.5 / p_z, -1.5 / p_z, -1.5 / p_z)
  Z <- list() # nolint
  beta <-  c(0.2, -0.3, -0.01, 0.02, 0.03, 0.04,
           rep(rep(0.02, p_z), n_t))
  beta_t <- -0.2
  sd_z_x <- 0.4
  X <- mvrnorm(n, mu = c(1, 5, 6, 7, 8), Sigma = diag(1, 5)) # nolint
  TRT <- rbinom(n, size = 1,  prob = 0.5)                      # nolint
  Y_constant <- beta[1] + (X %*% beta[2:6])                  # nolint
  Y0 <- 0                                                    # nolint
  Y1 <- 0                                                    # nolint
  A <- A1 <- A0 <- matrix(NA, nrow = n, ncol = n_t)          # nolint
  for (i in 1:n_t) {
    alphas[[i]] <- matrix(rep(c(2.3, -0.3, -0.01, 0.02,
                                0.03, 0.04, -0.4), p_z), ncol = p_z)
    gammas[[i]] <- c(1, -0.1, 0.2, 0.2, 0.2, 0.2, rep(z_para[i], p_z))
    Z0 <- alphas[[i]][1, ] + (X %*% alphas[[i]][2:6, ]) +                       # nolint
      mvrnorm(n, mu = rep(0, p_z), Sigma = diag(sd_z_x, p_z))
    Z1 <- alphas[[i]][1, ] + (X %*% alphas[[i]][2:6, ]) + alphas[[i]][7, ] +    # nolint
      mvrnorm(n, mu = rep(0, p_z), Sigma = diag(sd_z_x, p_z))
    Z[[i]] <- Z1 * TRT + Z0 * (1-TRT)                                               # nolint
    Y0 <- (Y0 + Z0 %*% matrix(beta[(7 + (i - 1) * p_z): (6 + p_z * i)],         # nolint
                              ncol = 1))[, 1]                                  # nolint
    Y1 <- (Y1 + Z1 %*% matrix(beta[(7 + (i - 1) * p_z): (6 + p_z * i)],         # nolint
                              ncol = 1))[, 1]
    if (i == 1) {
      A0[, i] <- rbinom(n, size = 1,
                       prob = 1 / (1 + exp(- (gammas[[i]][1] +
                      (X %*% gammas[[i]][2:6]) + Z0 %*%
                        matrix(gammas[[i]][7: (7 + p_z - 1)], ncol = 1))[, 1])))
      A1[, i] <- rbinom(n, size = 1,
                       prob = 1 / (1 + exp(- (gammas[[i]][1] +
                      (X %*% gammas[[i]][2:6]) + Z1 %*%
                        matrix(gammas[[i]][7: (7 + p_z - 1)], ncol = 1))[, 1])))
    } else {
      A0[, i] <- rbinom(n, size = 1,
                       prob = 1 / (1 + exp(- (gammas[[i]][1] +
                      (X %*% gammas[[i]][2:6]) + Z0 %*%
              matrix(gammas[[i]][7: (7 + p_z - 1)], ncol = 1))[, 1]))) *
        A0[, i - 1]
      A1[, i] <- rbinom(n, size = 1,
                       prob = 1 / (1 + exp(- (gammas[[i]][1] +
                      (X %*% gammas[[i]][2:6]) + Z1 %*%
              matrix(gammas[[i]][7: (7 + p_z - 1)], ncol = 1))[, 1]))) *
        A1[, i - 1]

    }
    A[, i] <- A1[, i] * TRT + A0[, i] * (1 - TRT)
  }
  Y0 <- Y0 + rnorm(n, mean = 0, sd = 0.3) + Y_constant                  # nolint
  Y1 <- Y1 + beta_t  + rnorm(n, mean = 0, sd = 0.3) + Y_constant        # nolint
  true1_i <- mean((Y1)[A1[, n_t] == 1 & A0[, n_t] == 1], na.rm = TRUE)
  true0_i <- mean((Y0)[A1[, n_t] == 1 & A0[, n_t] == 1], na.rm = TRUE)

  Y <- as.vector(Y1 * TRT + Y0 * (1 - TRT))                                 # nolint

  for (i in 2:n_t) {
    Z[[i]][A[, (i - 1)] == 0, ] <- NA
  }

  Y[A[, n_t] == 0] <- NA
  # estimate the treatment difference
  fit1 <- est_S_Plus_Plus_MethodA(X, A, Z, Y, TRT)
  fit2 <- est_S_Star_Plus_MethodA(X, A, Z, Y, TRT)
  fit3 <- est_S_Plus_Plus_MethodB(X, A, Z, Y, TRT)
  fit4 <- est_S_Star_Plus_MethodB(X, A, Z, Y, TRT)
  # Test class
  expect_equal(class(fit1), "list")
  expect_equal(class(fit2), "list")
  expect_equal(class(fit3), "list")
  expect_equal(class(fit4), "list")
  expect_equal(sum(is.na(fit1)), 0)
  expect_equal(sum(is.na(fit2)), 0)
  expect_equal(sum(is.na(fit3)), 0)
  expect_equal(sum(is.na(fit4)), 0)

  # Test number of elements
  expect_equal(length(fit1), 6)
  expect_equal(length(fit2), 6)
  expect_equal(length(fit3), 6)
  expect_equal(length(fit4), 6)

})
