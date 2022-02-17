#' @importFrom stats var coef predict qnorm as.formula cov formula
#' @importFrom stats glm lm model.matrix.lm sd
#' @importFrom utils head tail



#Generate k independent vectors given the input vec
mat_vec <- function(vec) {
  k <- length(vec)
  part1 <- cbind(vec[-1], diag(-vec[1], (k - 1)))
  mat <- rbind(vec, part1)
  return(mat)
}

#-- calculate integration------------------------------------------------------
Expect_function1D_MA <- function(x, n_time_points, gammas, alphas, Sigmas) {    # nolint
  res_1 <- NULL
  x <- matrix(x, nrow = 1)
  for (i in 1:n_time_points) {
    z_dim <- dim(alphas[[i]])[2]
    mu_z <- t(cbind(1, x) %*% alphas[[i]])
    part_1X <- (cbind(1, x) %*% head(gammas[[i]], -z_dim))[1, 1]   # nolint
    part_z <- as.numeric(tail(gammas[[i]], z_dim) %*% mu_z)
    Sig_sqrt <- chol(Sigmas[[i]])                                  # nolint
    b1_norm <- norm(Sig_sqrt %*% tail(gammas[[i]], z_dim), "F")
    vec1 <- Sig_sqrt %*% tail(gammas[[i]], z_dim) / b1_norm
    matvec1 <- mat_vec(as.numeric(vec1))
    mat_B <- pracma::gramSchmidt(matvec1)$Q                        # nolint
    seq_v <- seq(1e-4, 1 - 1e-4, length.out = 200)
    ### prob
    fun_prob <- function(v) {
      1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v)))
      }
    prob <- mean(fun_prob(seq_v))
    ### estimate E(Z)
    fun_expz <- function(v) {
      qnorm(v) / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v)))
      }
    intgl_expz <- mean(fun_expz(seq_v))
    expz <- 1 / b1_norm * Sigmas[[i]] %*% tail(gammas[[i]], z_dim) *
      intgl_expz + mu_z * prob
    ### estimate E(Z^2)
    fun_expz2 <- function(v) {
      qnorm(v)^2 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v)))
      }
    if (z_dim > 1) {
      mat_V <- diag(c(mean(fun_expz2(seq_v)), rep(prob, z_dim - 1)))    # nolint
    } else {
      mat_V <- matrix(mean(fun_expz2(seq_v)), 1, 1)                     # nolint
    }
    Sig_matB <- Sig_sqrt %*% mat_B                             # nolint
    mu_intgl <- mu_z %*% t(vec1 * intgl_expz) %*% Sig_sqrt
    expz2 <- Sig_matB %*% mat_V %*% t(Sig_matB) + mu_intgl +
      t(mu_intgl) + mu_z %*% t(mu_z) * prob
    ### proba
    fun_proba <- function(v) {
      1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))) *
        (1 - 1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))))
      }
    proba <- mean(fun_proba(seq_v))
    ### estimate E(aZ)
    fun_expaz <- function(v) {
      qnorm(v) / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))) *
        (1 - 1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))))
      }
    intgl_expaz <- mean(fun_expaz(seq_v))
    expaz <- 1 / b1_norm * Sigmas[[i]] %*% tail(gammas[[i]], z_dim) *
      intgl_expaz + mu_z * proba
    ### estimate E(aZ^2)
    fun_expaz2 <- function(v) {
      qnorm(v)^2 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))) *
        (1 - 1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))))
      }
    if (z_dim > 1) {
      mat_V2 <- diag(c(mean(fun_expaz2(seq_v)), rep(proba, z_dim - 1)))    # nolint
    } else {
      mat_V2 <- matrix(mean(fun_expaz2(seq_v)), 1, 1)                      # nolint
    }
    Sig_matB <- Sig_sqrt %*% mat_B                                      # nolint
    mu_intgla <- mu_z %*% t(vec1 * intgl_expaz) %*% Sig_sqrt
    expaz2 <- Sig_matB %*% mat_V2 %*% t(Sig_matB) + mu_intgla +
      t(mu_intgla) + mu_z %*% t(mu_z) * proba

    ### result at t = i
    res <- c(prob, expz, as.vector(expz2), proba, expaz, as.vector(expaz2))
    ### Append over time points
    res_1 <- c(res_1, res)
  }
  return(res_1)
}

#### Convert NA to 999 ####
NA_replace <- function(vec) {                                           # nolint
  vec[is.na(vec)] <- 999
  return(vec)
}


#### matrix inverse using SVD ####
inv_svd <- function(x) {
  temp <- svd(x)
  RVAL <- temp$v %*% diag(1 / temp$d) %*% t(temp$u)                     # nolint
  return(RVAL)
}

#### expectation approximation function for method B ####
###############################################################################
#
#                         A collection of required functions
###############################################################################

Expect_function1D_BU <- function(x, n_time_points, gammas, alphas, Sigmas) {    # nolint
  res_1 <- NULL
  x <- matrix(x, nrow = 1)
  for (i in 1:n_time_points) {
    z_dim <- dim(alphas[[i]])[2]
    gamma_z <- tail(gammas[[i]], z_dim)
    mu_z <- t(cbind(1, x) %*% alphas[[i]])
    part_1X <- (cbind(1, x) %*% head(gammas[[i]], - z_dim))[1, 1]      # nolint
    part_z <- as.numeric(gamma_z %*% mu_z)
    Sig_sqrt <- chol(Sigmas[[i]])                                      # nolint
    b1_norm <- norm(Sig_sqrt %*% tail(gammas[[i]], z_dim), "F")
    vec1 <- Sig_sqrt %*% tail(gammas[[i]], z_dim) / b1_norm
    matvec1 <- mat_vec(as.numeric(vec1))
    mat_B <- pracma::gramSchmidt(matvec1)$Q                            # nolint
    seq_v <- seq(1e-4, 1 - 1e-4, length.out = 1000)
    ### prob
    fun_prob <- function(v) {
      1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v)))
      }
    prob <- mean(fun_prob(seq_v))
    ### estimate E(Z)
    fun_expz <- function(v) {
      qnorm(v) / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v)))
      }
    intgl_expz <- mean(fun_expz(seq_v))
    expz <- 1 / b1_norm * Sigmas[[i]] %*% gamma_z * intgl_expz + mu_z * prob
    ### proba
    fun_proba <- function(v) {
      1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))) *
        (1 - 1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))))
      }
    proba <- mean(fun_proba(seq_v))
    ### estimate E(aZ)
    fun_expaz <- function(v) {
      qnorm(v) / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))) *
        (1 - 1 / (1 + exp(-part_1X - part_z - b1_norm * qnorm(v))))
      }
    intgl_expaz <- mean(fun_expaz(seq_v))
    expaz <- 1 / b1_norm * Sigmas[[i]] %*% gamma_z * intgl_expaz + mu_z * proba
    ### result at t = i
    res <- c(prob, expz, proba, expaz)
    ### Append over time points
    res_1 <- c(res_1, res)
  }
  return(res_1)
}
