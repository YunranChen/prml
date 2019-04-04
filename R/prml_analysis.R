# logsum and logmean
#' Calculate the log(sum(exp(x)))
#'
#' @param x A vector.
#' @return A number, the \eqn{ log(sum(exp(x))) }

logsum <- function(x) {
    res <- max(x) + log(sum(exp(x - max(x))))
    return(res)
}

#' Calculate the log(mean(exp(x)))
#'
#' @param x A vector.
#' @return A number, the \eqn{ log(mean(exp(x))) }

logmean <- function(x) {
    # to calculate log(mean(exp(x)))
    res <- max(x) + log(sum(exp(x - max(x)))) - log(length(x))
    return(res)
}

# truncated gamma.

#' Generate random variable from truncated gamma
#'
#' @param size A number, size of the random variable.
#' @param shape_ A number, shape of the truncated gamma.
#' @param rate_ A number, rate of the truncated gamma.
#' @param a lower bound of the range
#' @param b upper bound of the range
#' @return A vector, iid samples from \eqn{ Ga_{[a,b]}(shape_,rate_) }
rtgamma <- function(size, shape_, rate_, a, b) {
    u <- runif(n = size)
    c_inv <- pgamma(q = b, shape = shape_, rate = rate_) -
        pgamma(q = a, shape = shape_, rate = rate_)
    x <- qgamma(p = u * c_inv + pgamma(q = a, shape = shape_, rate = rate_),
                shape = shape_, rate = rate_)
    return(x)
}

#' Calculate pdf of the truncated gamma
#'
#' @param x_ A number, value.
#' @param shape_ A number, shape of the truncated gamma.
#' @param rate_ A number, rate of the truncated gamma.
#' @param a lower bound of the range
#' @param b upper bound of the range
#' @return A vector, pdf of \eqn{ Ga_{[a,b]}(x_; shape_,rate_) }
dtgamma <- function(x_, shape_, rate_, a, b) {
    c_inv <- pgamma(q = b, shape = shape_, rate = rate_) -
        pgamma(q = a, shape = shape_, rate = rate_)
    x <- dgamma(x = x_, shape = shape_, rate = rate_) / c_inv
    return(x)
}

# Pprml for single trial
prml_int <- function(xs_bn, n_gq = 20, n_per = 100, alpha = 0.5) {
    # return log(m_x) permutation
    n <- length(xs_bn)
    ind <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })
    # guass.quad
    out <- gauss.quad(n_gq)
    ## weight seq
    w <- 1 / (1:n + 1)
    ## Recursion
    Xs_bn_p <- sapply(1:n_per, function(i) {
        xs_bn[ind[, i]]
    })  #each col is a permutation
    lq <- quantile(xs_bn, 0.25)
    uq <- quantile(xs_bn, 0.75)
    a <- lq - alpha * (uq - lq)
    b <- uq + alpha * (uq - lq)
    names(a) <- names(b) <- NULL
    xs_gq <- out$nodes * (b - a) / 2 + (a + b) / 2  #change interval
    ## Initial Guess: Uniform
    f_gq <- rep(1 / (b - a), n_gq)
    res_t <- map(1:n_per, function(ix) {
        m_x <- 0
        for (j in 1:n) {
            px_grid <- map_dbl(xs_gq, ~dpois(x = Xs_bn_p[j, ix], lambda = .x))
            mi_x <- ((b - a) / 2 * sum(out$weights * px_grid * f_gq))
            f_gq <- (1 - w[j]) * f_gq + w[j] * px_grid * f_gq / mi_x
            m_x <- log(mi_x) + m_x
        }
        return(list(m_x = m_x))
    })
    res_py <- logmean(map_dbl(res_t, ~.x$m_x))
    return(res_m = res_py)
}
cml <- function(xs_bn, alpha = 0.5) {
    yi <- xs_bn
    alpha <- 0.5
    lq <- quantile(yi, 0.25)
    uq <- quantile(yi, 0.75)
    a <- lq - alpha * (uq - lq)
    b <- uq + alpha * (uq - lq)
    names(a) <- names(b) <- NULL
    res <- -log(b - a) - sum(lgamma(yi + 1)) +
        lgamma(sum(yi) + 1) - (sum(yi) + 1) * log(length(yi)) +
        log(pgamma(q = b, shape = sum(yi) + 1, rate = length(yi)) -
                pgamma(q = a, shape = sum(yi) + 1, rate = length(yi)))
    return(res)
}

#' PRML filter: Bayes factor of Poisson versus Poisson mixtures.
#'
#' @param xs_bn A vector. Spike counts of repeated single-stimulus trial data.
#' @param n_gq A number. 20 by default. Number of grids in Gaussion quadrature.
#' @param n_per A number. 100 by default. Permutation of likihood estimation to obtain the order-invariant estimator.
#' @param alpha 0.5 by default. The range of the spike counts estimator \eqn{ [Y_{0.25}-\alpha IQR,Y_{0.75}+\alpha IQR] }
#' @return A number. Bayes factor of Poisson versus Poisson mixtures by PRML algorithm.
prml_filter <- function(xs_bn, n_gq = 20, n_per = 100, alpha = 0.5) {
    res <- exp(cml(xs_bn, alpha) - prml_int(xs_bn, n_gq, n_per, alpha))
    return(res)
}

# Pprml for multiplexing
prml_outA_lp <- function(xs_bn, xs_a, r_a = 0.5, s_a = 2e-10, mu_l = "min",
                         n_gq = 20, n_per = 100) {
    # default is Jeffery
    r_A <- sum(xs_a) + r_a
    s_A <- length(xs_a) + s_a
    imu_a <- mean(xs_a)
    if (mu_l == "min") {
        mu_l <- max(0, min(min(xs_bn) - 2 * sd(xs_bn),
                           min(xs_a) - 2 * sd(xs_a)))
    }
    n <- length(xs_bn)
    indx <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })
    # guass.quad
    out <- gauss.quad(n_gq)

    allf <- function(phi) {
        a <- 0
        b <- 1
        xs_gq <- out$nodes * (b - a) / 2 + (a + b) / 2  #change interval
        ## Initial Guess:
        d_fgq <- rep(0, n_gq)
        f_gq <- rep(1 / (b - a), n_gq)
        ## weight seq
        w <- 1 / (1:n + 1)

        ## Recursion
        Xs_bn_p <- sapply(1:n_per, function(i) {
            xs_bn[indx[, i]]
        })  #each col is a permutation
        res_t <- map(1:n_per, function(ix) {
            m_x <- 0
            grad <- 0
            for (j in 1:n) {
                # step1:prml
                px_grid <- map_dbl(xs_gq,
                                   ~dpois(x = Xs_bn_p[j, ix],
                                          lambda = .x * (phi - mu_l) + mu_l))
                mi_x <- ((b - a) / 2 * sum(out$weights * px_grid * f_gq))
                # step2:gradient
                d_px_grid <- px_grid *
                    (Xs_bn_p[j, ix] * xs_gq /(mu_l + xs_gq * (phi - mu_l)) - xs_gq)
                G1_gq <- d_px_grid * f_gq + px_grid * d_fgq
                G2_gq <- px_grid * f_gq
                d_logmi_x <- (b - a) / 2 * sum(out$weights * G1_gq) / mi_x
                # step3:update f,gr_f
                f_gq <- (1 - w[j]) * f_gq + w[j] * px_grid * f_gq / mi_x
                d_fgq <- (1 - w[j]) * d_fgq + w[j] * (G1_gq - G2_gq * d_logmi_x) / mi_x
                # store the result
                m_x <- log(mi_x) + m_x
                grad <- grad + d_logmi_x
            }
            return(list(m_x = m_x, grad = grad))
        })
        pm_x <- map_dbl(res_t, ~.x$m_x) %>% logmean(.)
        pgrad <- map_dbl(res_t, function(res_t) {
            res_t$grad
        }) %>% logmean(.)
        return(list(m_x = pm_x, grad = pgrad))
    }
    obj <- function(theta) {
        phi <- exp(theta) + mu_l
        p_mua <- dgamma(x = phi, shape = r_A, rate = s_A, log = TRUE)
        return(-(allf(phi)$m_x + p_mua))
    }
    grad <- function(theta) {
        phi <- exp(theta) + mu_l
        g_mua <- ((r_A - 1) / phi - s_A)
        prml <- allf(phi)$grad
        res <- -(prml + g_mua) * exp(theta)
        return(res)
    }
    res <- optim(par = log(imu_a - mu_l), fn = obj,
                 gr = grad, method = "BFGS", hessian = TRUE)
    # optim is for minimization
    H <- (-res$hessian - grad(res$par)) * exp(-2 * res$par)
    logmx <- 0.5 * log(2 * pi) + 0.5 * log((-H)^(-1)) - res$value
    res_py <- logmx
    return(res_py)
}
prml_outB_lp <- function(xs_bn, xs_b, r_b = 0.5, s_b = 2e-10,
                         mu_u = "max", n_gq = 20, n_per = 100) {

    r_B <- sum(xs_b) + r_b
    s_B <- length(xs_b) + s_b
    imu_b <- mean(xs_b)
    # set the outside B range: max + 2sigma
    if (mu_u == "max") {
        mu_u <- max(max(xs_bn) + 2 * sd(xs_bn), max(xs_b) + 2 * sd(xs_b))
    }
    n <- length(xs_bn)
    indx <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })

    # guass.quad
    out <- gauss.quad(n_gq)

    allf <- function(phi) {
        a <- 0
        b <- 1
        xs_gq <- out$nodes * (b - a) / 2 + (a + b) / 2  #change interval
        ## Initial Guess:
        d_fgq <- rep(0, n_gq)
        f_gq <- rep(1 / (b - a), n_gq)

        ## weight seq
        w <- 1 / (1:n + 1)

        ## Recursion

        Xs_bn_p <- sapply(1:n_per, function(i) {
            xs_bn[indx[, i]]
        })  #each col is a permutation
        res_t <- map(1:n_per, function(ix) {

            m_x <- 0
            grad <- 0
            for (j in 1:n) {
                # step1:prml
                px_grid <- map_dbl(xs_gq, ~dpois(x = Xs_bn_p[j, ix],
                                                 lambda = .x * (mu_u - phi) + phi))
                mi_x <- ((b - a) / 2 * sum(out$weights * px_grid * f_gq))
                # step2:gradient
                d_px_grid <- px_grid * (Xs_bn_p[j, ix] * (1 - xs_gq) /
                                            (phi + xs_gq * (mu_u - phi)) - (1 - xs_gq))
                G1_gq <- d_px_grid * f_gq + px_grid * d_fgq
                G2_gq <- px_grid * f_gq
                d_logmi_x <- (b - a) / 2 * sum(out$weights * G1_gq) / mi_x
                # step3:update f,gr_f
                f_gq <- (1 - w[j]) * f_gq + w[j] * px_grid * f_gq / mi_x
                d_fgq <- (1 - w[j]) * d_fgq + w[j] * (G1_gq - G2_gq * d_logmi_x) / mi_x
                # store the result
                m_x <- log(mi_x) + m_x
                grad <- grad + d_logmi_x
            }
            # return(list(f=f,m_x=m_x,grad=grad))
            return(list(m_x = m_x, grad = grad))
        })
        pm_x <- map_dbl(res_t, ~.x$m_x) %>% logmean(.)
        pgrad <- map_dbl(res_t, function(res_t) {
            res_t$grad
        }) %>% logmean(.)
        return(list(m_x = pm_x, grad = pgrad))
    }

    obj <- function(theta) {
        phi <- exp(theta) * mu_u / (1 + exp(theta))
        p_mua <- dgamma(x = phi, shape = r_B, rate = s_B, log = TRUE)
        return(-(allf(phi)$m_x + p_mua))
    }
    grad <- function(theta) {
        phi <- exp(theta) * mu_u / (1 + exp(theta))
        g_mua <- ((r_B - 1) / phi - s_B)
        prml <- allf(phi)$grad
        res <- -(prml + g_mua) * phi * (1 - phi / mu_u)
        return(res)
    }
    res <- optim(par = log(imu_b / (mu_u - imu_b)), fn = obj, gr = grad,
                 method = "BFGS", hessian = TRUE)
    # optim is for minimization
    phi_r <- exp(res$par) * mu_u / (1 + exp(res$par))
    H <- (res$hessian * mu_u + grad(res$par) * (2 * phi_r - mu_u)) *
        (phi_r * (mu_u - phi_r))^(-2) * mu_u
    logmx <- 0.5 * log(2 * pi) + 0.5 * log((H)^(-1)) - res$value
    res_py <- logmx
    return(res_py)
}
prml_sin_lp <- function(xs_bn, xs_a, r_a = 0.5, s_a = 2e-10, n_gq = 20) {

    r_A <- sum(xs_a) + r_a
    s_A <- length(xs_a) + s_a
    n <- length(xs_bn)
    r <- sum(xs_bn) + r_A
    s <- n + s_A
    logmx <- -sum(lgamma(xs_bn + 1)) + lgamma(r) - r * log(s) + r_A * log(s_A) - lgamma(r_A)
    res_py <- logmx
    return(res_py)
}
prml_int_lp <- function(xs_bn, xs_a, xs_b, e = 0, r_a = 0.5, s_a = 2e-10,
                        r_b = 0.5, s_b = 2e-10, n_gq = 20, n_per = 100) {

    r_A <- sum(xs_a) + r_a
    s_A <- length(xs_a) + s_a
    r_B <- sum(xs_b) + r_b
    s_B <- length(xs_b) + s_b
    imu_a <- mean(xs_a)
    imu_b <- mean(xs_b)

    n <- length(xs_bn)
    indx <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })

    # guass.quad
    out <- gauss.quad(n_gq)

    allf <- function(phi) {
        a <- e
        b <- 1 - e
        xs_gq <- out$nodes * (b - a) / 2 + (a + b) / 2  #change interval

        ## Initial Guess:
        f_gq <- rep(1 / (b - a), n_gq)
        d_fgq <- matrix(0, nrow = 2, ncol = n_gq)
        ## weight seq
        w <- 1 / (1:n + 1)

        ## Recursion
        Xs_bn_p <- sapply(1:n_per, function(i) {
            xs_bn[indx[, i]]
        })  #each col is a permutation
        res_t <- map(1:n_per, function(ix) {
            m_x <- 0
            grad <- c(0, 0)
            for (j in 1:n) {
                # step1:prml
                px_grid <- map_dbl(xs_gq, ~dpois(x = Xs_bn_p[j, ix],
                                                 lambda = phi[1] + .x * (phi[2] - phi[1])))
                mi_x <- ((b - a) / 2 * sum(out$weights * px_grid * f_gq))
                # step2:gradient
                d_px_grid_a <- px_grid * (Xs_bn_p[j, ix] * (1 - xs_gq) /
                                              (phi[1] + xs_gq * (phi[2] - phi[1])) - (1 - xs_gq))
                d_px_grid_b <- px_grid * (Xs_bn_p[j, ix] * xs_gq /
                                              (phi[1] + xs_gq * (phi[2] - phi[1])) - xs_gq)
                G_a_gq <- d_px_grid_a * f_gq + px_grid * d_fgq[1, ]
                G_b_gq <- d_px_grid_b * f_gq + px_grid * d_fgq[2, ]
                G2_gq <- px_grid * f_gq
                d_logmi_x_a <- (b - a) / 2 * sum(out$weights * G_a_gq) / mi_x
                d_logmi_x_b <- (b - a) / 2 * sum(out$weights * G_b_gq) / mi_x
                d_logmi_x <- c(d_logmi_x_a, d_logmi_x_b)
                # step3:update f,gr_f
                d_fgq[1, ] <- (1 - w[j]) * d_fgq[1, ] + w[j] *
                    (G_a_gq - G2_gq * d_logmi_x_a) / mi_x
                d_fgq[2, ] <- (1 - w[j]) * d_fgq[2, ] + w[j] *
                    (G_b_gq - G2_gq * d_logmi_x_b) / mi_x
                f_gq <- (1 - w[j]) * f_gq + w[j] * px_grid * f_gq / mi_x
                # store the result
                m_x <- log(mi_x) + m_x
                grad <- grad + d_logmi_x
            }
            return(list(m_x = m_x, grad = grad))
        })
        pm_x <- map_dbl(res_t, ~.x$m_x) %>% logmean(.)
        pgrad <- sapply(res_t, function(res_t) {
            res_t$grad
        }) %>% apply(., 1, logmean)
        return(list(m_x = pm_x, grad = pgrad))
    }
    obj <- function(theta) {
        phi <- exp(sort(theta))
        p_mua <- dgamma(x = phi[1], shape = r_A, rate = s_A, log = TRUE)
        p_mub <- dgamma(x = phi[2], shape = r_B, rate = s_B, log = TRUE)
        return(-(allf(phi)$m_x + p_mua + p_mub))
    }
    grad <- function(theta) {
        phi <- exp(sort(theta))
        g_mua <- (r_A - 1) / phi[1] - s_A
        g_mub <- (r_B - 1) / phi[2] - s_B
        prml <- allf(phi)$grad
        res <- -(prml + c(g_mua, g_mub)) * phi
        return(res[order(theta)])
    }
    res <- optim(par = c(log(imu_a), log(imu_b)), fn = obj, gr = grad,
                 method = "BFGS", hessian = TRUE)
    # optim is for minimization
    H <- matrix(0, nrow = 2, ncol = 2)
    res_par <- res$par %>% sort(.)
    H[1, 1] <- (-res$hessian[1, 1] - grad(res_par)[1]) * exp(-2 * res_par[1])
    H[2, 2] <- (-res$hessian[2, 2] - grad(res_par)[2]) * exp(-2 * res_par[2])
    H[1, 2] <- -res$hessian[1, 2] * exp(-sum(res_par))
    H[2, 1] <- -res$hessian[2, 1] * exp(-sum(res_par))

    logmx <- log(2 * pi) - 0.5 * log(det(-H)) - res$value
    res_py <- logmx
    return(res_py)
}
prml_mix_lp <- function(xs_bn, xs_a, xs_b, e = 0, r_a = 0.5, s_a = 2e-10,
                        r_b = 0.5, s_b = 2e-10, n_gq = 20, n_per = 100) {

    r_A <- sum(xs_a) + r_a
    s_A <- length(xs_a) + s_a
    r_B <- sum(xs_b) + r_b
    s_B <- length(xs_b) + s_b
    imu_a <- mean(xs_a)
    imu_b <- mean(xs_b)
    del <- 2 * e

    n <- length(xs_bn)
    indx <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })

    allf <- function(phi) {
        a <- 0
        b <- 1

        z <- c(a, b)

        ## Initial Guess:
        f <- c(0.5, 0.5)
        d_f <- matrix(0, nrow = 2, ncol = 2)

        ## weight seq
        w <- 1 / (1:n + 1)

        ## Recursion
        Xs_bn_p <- sapply(1:n_per, function(i) {
            xs_bn[indx[, i]]
        })  #each col is a permutation
        res_t <- map(1:n_per, function(ix) {
            m_x <- 0
            grad <- c(0, 0)
            for (j in 1:n) {
                # step1:prml
                poi_x <- map_dbl(z, ~(1 - del) * (dpois(x = Xs_bn_p[j, ix],
                                                        lambda = phi[1] + .x * (phi[2] - phi[1]))) +
                                     del / 2 *(dpois(x = Xs_bn_p[j, ix], lambda = phi[1]) +
                                                   dpois(x = Xs_bn_p[j, ix], lambda = phi[2])))
                mi_x <- sum(poi_x * f)

                # step2:gradient
                d_poi_x_a <- poi_x * (Xs_bn_p[j, ix] * (1 - z) / (phi[1] + z * (phi[2] - phi[1])) - (1 - z)) * (1 - del) +
                    del / 2 * poi_x * (Xs_bn_p[j, ix] / (phi[1]) - (1))
                d_poi_x_b <- poi_x * (Xs_bn_p[j, ix] * z / (phi[1] + z * (phi[2] - phi[1])) - z) * (1 - del) +
                    del / 2 * poi_x * (Xs_bn_p[j, ix] / (phi[2]) - (1))
                G_a <- d_poi_x_a * f + poi_x * d_f[1, ]
                G_b <- d_poi_x_b * f + poi_x * d_f[2, ]
                G2 <- poi_x * f
                d_logmi_x_a <- sum(G_a) / mi_x
                d_logmi_x_b <- sum(G_b) / mi_x
                d_logmi_x <- c(d_logmi_x_a, d_logmi_x_b)
                # step3:update f,gr_f
                d_f[1, ] <- (1 - w[j]) * d_f[1, ] + w[j] * (G_a - G2 * d_logmi_x_a) / mi_x
                d_f[2, ] <- (1 - w[j]) * d_f[2, ] + w[j] * (G_b - G2 * d_logmi_x_b) / mi_x
                f <- (1 - w[j]) * f + w[j] * poi_x * f / mi_x

                # store the result
                m_x <- log(mi_x) + m_x
                grad <- grad + d_logmi_x
            }
            return(list(m_x = m_x, grad = grad))
        })
        pm_x <- map_dbl(res_t, ~.x$m_x) %>% logmean(.)
        pgrad <- sapply(res_t, function(res_t) {
            res_t$grad
        }) %>% apply(., 1, logmean)
        return(list(m_x = pm_x, grad = pgrad))
    }

    obj <- function(theta) {
        phi <- exp(sort(theta))
        p_mua <- dgamma(x = phi[1], shape = r_A, rate = s_A, log = TRUE)
        p_mub <- dgamma(x = phi[2], shape = r_B, rate = s_B, log = TRUE)
        return(-(allf(phi)$m_x + p_mua + p_mub))
    }
    grad <- function(theta) {
        phi <- exp(sort(theta))
        g_mua <- (r_A - 1) / phi[1] - s_A
        g_mub <- (r_B - 1) / phi[2] - s_B
        prml <- allf(phi)$grad
        res <- -(prml + c(g_mua, g_mub)) * phi
        return(res[order(theta)])
    }
    res <- optim(par = c(log(imu_a), log(imu_b)), fn = obj, gr = grad,
                 method = "BFGS", hessian = TRUE)
    # optim is for minimization
    H <- matrix(0, nrow = 2, ncol = 2)
    res_par <- res$par %>% sort(.)
    H[1, 1] <- (-res$hessian[1, 1] - grad(res_par)[1]) * exp(-2 * res_par[1])
    H[2, 2] <- (-res$hessian[2, 2] - grad(res_par)[2]) * exp(-2 * res_par[2])
    H[1, 2] <- -res$hessian[1, 2] * exp(-sum(res_par))
    H[2, 1] <- -res$hessian[2, 1] * exp(-sum(res_par))

    logmx <- log(2 * pi) - 0.5 * log(det(-H)) - res$value
    res_py <- logmx
    return(res_py)
}

#' PRML classifier: Posterior probability under each Poisson mixtures hypotheses.
#'
#' @param xs_bn A vector. Spike counts of repeated dual-stimuli trial data AB.
#' @param xs_a A vector. Spike counts of repeated single-stimulus trial data A.
#' @param xs_b A vector. Spike counts of repeated single-stimulus trial data B.
#' @param mu_l A number. Lower bound of spike counts. "min" by default. Indicating \eqn{ max(0,\underset{j=A,B,AB}{\min}(\min(Y_j)-2\text{std}(Y_j))) }
#' @param mu_u A number. Upper bound of spike counts. "max" by default. Indicating \eqn{ \underset{j=A,B,AB}{\max}(\max(Y_j)+2\text{std}(Y_j)) }
#' @param e A number. 0 by default. Shringkage on the domain and meansurement of mixing density f under the Intermediate and Mixture hypothese.
#' @param r_a A number. The parameter in gamma prior of spike rate mu_A. rate. Jeffereys' prior by default.
#' @param s_a A number. The parameter in gamma prior of spike rate mu_A. shape. Jeffereys' prior by default.
#' @param r_b A number. The parameter in gamma prior of spike rate mu_B. rate. Jeffereys' prior by default.
#' @param s_b A number. The parameter in gamma prior of spike rate mu_B. shape. Jeffereys' prior by default.
#' @param n_gq A number. 20 by default. Number of grids in Gaussion quadrature.
#' @param n_per A number. 100 by default. Permutation of likihood estimation to obtain the order-invariant estimator.
#' @return A list.
#' \describe{
#'   \item{post.prob}{posterior probabilities under Mixture, Intermediate, Outside, Single hypotheses.}
#'   \item{win.model}{the model has largest post.prob.}
#' }
prml_classifier <- function(xs_bn, xs_a, xs_b, mu_l = "min", mu_u = "max", e = 0,
                     r_a = 0.5, s_a = 2e-10, r_b = 0.5, s_b = 2e-10, n_gq = 20, n_per = 100) {
    if (mean(xs_a) <= mean(xs_b)) {
        outA <- prml_outA_lp(xs_bn, xs_a, r_a, s_a, mu_l, n_gq, n_per)
        outB <- prml_outB_lp(xs_bn, xs_b, r_b, s_b, mu_u, n_gq, n_per)
        sinA <- prml_sin_lp(xs_bn, xs_a, r_a, s_a, n_gq)
        sinB <- prml_sin_lp(xs_bn, xs_b, r_b, s_b, n_gq)
        int <- prml_int_lp(xs_bn, xs_a, xs_b, e, r_a, s_a, r_b, s_b, n_gq, n_per)
        mix <- prml_mix_lp(xs_bn, xs_a, xs_b, e, r_a, s_a, r_b, s_b, n_gq, n_per)
    } else {
        outA <- prml_outA_lp(xs_bn, xs_b, r_b, s_b, mu_l, n_gq, n_per)
        outB <- prml_outB_lp(xs_bn, xs_a, r_a, s_a, mu_u, n_gq, n_per)
        sinA <- prml_sin_lp(xs_bn, xs_b, r_b, s_b, n_gq)
        sinB <- prml_sin_lp(xs_bn, xs_a, r_a, s_a, n_gq)
        int <- prml_int_lp(xs_bn, xs_b, xs_a, e, r_b, s_b, r_a, s_a, n_gq, n_per)
        mix <- prml_mix_lp(xs_bn, xs_b, xs_a, e, r_b, s_b, r_a, s_a, n_gq, n_per)
    }

    mxs <- c(mix, int, logmean(c(outA, outB)), logmean(c(sinA, sinB)))
    post.prob <- exp(mxs - logsum(mxs))
    # models <- c('Mix', 'Int', 'Out', 'Sing')
    win.model <- which.max(post.prob)
    out <- list(post.prob = post.prob, win.model = win.model)
    return(out)
}


log.pm <- function(x, a, b) {
    a.x <- a + sum(x)
    b.x <- b + length(x)
    mu.map <- a.x / b.x
    return(sum(dpois(x, mu.map, log = TRUE)) -
               diff(dgamma(mu.map, c(a, a.x), c(b, b.x), log = TRUE)))
}

#' Result of PRML classifier and PRML filter.
#'
#' @param xA A vector. Spike counts of repeated dual-stimuli trial data AB.
#' @param xB A vector. Spike counts of repeated single-stimulus trial data A.
#' @param xAB A vector. Spike counts of repeated single-stimulus trial data B.
#' @param labels A vector. labels for the trials.
#' @param remove.zeros A logical value. Whether to remove 0s in spike counts.
#' @param mu_l A number. Lower bound of spike counts. "min" by default. Indicating \eqn{ \text{max}(0,\underset{j=A,B,AB}{\text{min}}(\text{min}(Y_j)-2\text{std}(Y_j))) }
#' @param mu_u A number. Upper bound of spike counts. "max" by default. Indicating \eqn{ \underset{ j=A,B,AB }{\text{max}}(\text{max}(Y_j)+2\text{std}(Y_j)) }
#' @param e A number. 0 by default. Shringkage on the domain and meansurement of mixing density f under the Intermediate and Mixture hypothese.
#' @param gamma.pars A length 2 vector. The shape and rate of gamma prior for spike rate mu_A and mu_B. Jeffereys' prior by default.
#' @param n_gq A number. 20 by default. Number of grids in Gaussion quadrature.
#' @param n_per A number. 100 by default. Permutation of likihood estimation to obtain the order-invariant estimator.
#' @param alpha 0.5 by default. (For PRML filter) The range of the spike counts estimator \eqn{ [Y_{0.25}-\alpha \text{IQR},Y_{0.75}+\alpha \text{IQR}] }
#' @return A list.
#' \describe{
#'   \item{separation.logBF}{log Bayes factor for the hypothesis \eqn{ mu_A=mu_B } versus \eqn{ mu_A \neq mu_B }.}
#'   \item{post.prob}{posterior probabilities under Mixture, Intermediate, Outside, Single hypotheses.}
#'   \item{win.model}{the model has largest post.prob.}
#'   \item{prml.filter.bf}{Bayes factor of PRML filter for single-stimulus trial A and B}
#'   \item{samp.sizes}{number of repeated trials under condition A, B, AB}
#' }
prml_tests <- function(xA, xB, xAB, labels = c("A", "B", "AB"), remove.zeros = FALSE,
                       mu_l = "min", mu_u = "max", e = 0, gamma.pars = c(0.5, 2e-10),
                       n_gq = 20, n_per = 100, alpha = 0.5) {

    a <- r_a <- r_b <- gamma.pars[1]
    b <- s_a <- s_b <- gamma.pars[2]

    if (remove.zeros) {
        xA <- xA[xA != 0]
        xB <- xB[xB != 0]
        xAB <- xAB[xAB != 0]
    }

    nA <- length(xA)
    nB <- length(xB)
    nAB <- length(xAB)
    if (nA == 0 | nB == 0)
        stop("not enough data in single sound")

    ## how different are the two pure trials?

    two.poi.ibf <- Vectorize(function(i, j)
        return(log.pm(xA[-i], xA[i] + a, 1 + b) +
                   log.pm(xB[-j], xB[j] + a, 1 + b) -
                   log.pm(c(xA[-i], xB[-j]), xA[i] + xB[j] + a, 2 + b)))
    lbf.pure <- mean(c(outer(1:length(xA), 1:length(xB), two.poi.ibf)))

    pvls <- c(prml_filter(xA, n_gq, n_per, alpha), prml_filter(xB, n_gq, n_per, alpha))
    prmlmix_res <- prml_classifier(xAB, xA, xB, mu_l, mu_u, e, r_a, s_a, r_b, s_b, n_gq, n_per)
    out <- list(separation.logBF = lbf.pure, post.prob = prmlmix_res$post.prob,
                win.model = prmlmix_res$win.model, prml.filter.bf = pvls,
                samp.sizes = c(nA, nB, nAB))
    return(out)
}

## decide the range for x according to rule-of-thumb
x_range <- function(xs_vec, n_mu, l, u) {
    n <- length(xs_vec)
    bw_cut <- sd(xs_vec) * (4 / 3 / n)^(1 / 5) * 3
    ub <- max(xs_vec) + bw_cut
    lb <- min(xs_vec) - bw_cut
    xxs <- seq(max(lb, l), min(ub, u), length.out = n_mu)
    return(xxs)
}


# prml for multiplexing
prml_outA_lp_f <- function(xs_bn, xs_a, r_a = 0.5, s_a = 2e-10, mu_l = "min",
                           n_gq = 20, n_mu = 100, n_per = 100) {
    # default is Jeffery
    r_A <- sum(xs_a) + r_a
    s_A <- length(xs_a) + s_a
    imu_a <- mean(xs_a)
    if (mu_l == "min") {
        mu_l <- max(0, min(min(xs_bn) - 2 * sd(xs_bn), min(xs_a) - 2 * sd(xs_a)))
    }
    n <- length(xs_bn)
    indx <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })

    # guass.quad
    out <- gauss.quad(n_gq)

    allf <- function(phi) {
        a <- 0
        b <- 1
        xs_gq <- out$nodes * (b - a) / 2 + (a + b) / 2  #change interval
        ind <- ceiling((xs_gq - a) / ((b - a) / n_mu))
        z <- seq(a, b, length.out = n_mu)

        ## Initial Guess:

        f <- rep(1 / (b - a), n_mu)  #uniform
        f[n_mu] <- 0
        d_f <- 0
        d_fgq <- rep(0, n_gq)
        f_gq <- rep(1 / (b - a), n_gq)
        ## weight seq
        w <- 1 / (1:n + 1)

        ## Recursion
        Xs_bn_p <- sapply(1:n_per, function(i) {
            xs_bn[indx[, i]]
        })  #each col is a permutation
        res_t <- map(1:n_per, function(ix) {
            m_x <- 0
            grad <- 0
            for (j in 1:n) {
                # step1:prml
                px_grid <- map_dbl(xs_gq,
                                   ~dpois(x = Xs_bn_p[j, ix],
                                          lambda = .x * (phi - mu_l) + mu_l))
                mi_x <- ((b - a) / 2 * sum(out$weights * px_grid * f_gq))
                poi_x <- map_dbl(z,
                                 ~dpois(x = Xs_bn_p[j, ix],
                                        lambda = .x * (phi - mu_l) + mu_l))
                # step2:gradient
                d_poi_x <- poi_x * (Xs_bn_p[j, ix] * z / (mu_l + z * (phi - mu_l)) - z)
                d_px_grid <- px_grid *
                    (Xs_bn_p[j, ix] * xs_gq / (mu_l + xs_gq * (phi - mu_l)) - xs_gq)
                G1 <- d_poi_x * f + poi_x * d_f
                G1_gq <- d_px_grid * f_gq + px_grid * d_fgq
                G2 <- poi_x * f
                G2_gq <- px_grid * f_gq
                d_logmi_x <- (b - a) / 2 * sum(out$weights * G1_gq) / mi_x
                # step3:update f,gr_f
                f <- (1 - w[j]) * f + w[j] * poi_x * f / mi_x
                f_gq <- (1 - w[j]) * f_gq + w[j] * px_grid * f_gq / mi_x
                d_f <- (1 - w[j]) * d_f + w[j] * (G1 - G2 * d_logmi_x) / mi_x
                d_fgq <- (1 - w[j]) * d_fgq + w[j] * (G1_gq - G2_gq * d_logmi_x) / mi_x
                # store the result
                m_x <- log(mi_x) + m_x
                grad <- grad + d_logmi_x
            }
            return(list(f = f, m_x = m_x, grad = grad))
        })
        pf <- sapply(res_t, function(res_t) {
            res_t$f
        }) %>% apply(., 1, mean)
        pm_x <- map_dbl(res_t, ~.x$m_x) %>% logmean(.)
        pgrad <- map_dbl(res_t, function(res_t) {
            res_t$grad
        }) %>% logmean(.)
        return(list(f = pf, m_x = pm_x, grad = pgrad))
    }
    obj <- function(theta) {
        phi <- exp(theta) + mu_l
        p_mua <- dgamma(x = phi, shape = r_A, rate = s_A, log = TRUE)
        return(-(allf(phi)$m_x + p_mua))
    }
    grad <- function(theta) {
        phi <- exp(theta) + mu_l
        g_mua <- ((r_A - 1) / phi - s_A)
        prml <- allf(phi)$grad
        res <- -(prml + g_mua) * exp(theta)
        return(res)
    }
    res <- optim(par = log(imu_a - mu_l), fn = obj, gr = grad,
                 method = "BFGS", hessian = TRUE)
    # optim is for minimization
    H <- (-res$hessian - grad(res$par)) * exp(-2 * res$par)
    logmx <- 0.5 * log(2 * pi) + 0.5 * log((-H)^(-1)) - res$value
    res_py <- logmx
    res_phi <- exp(res$par) + mu_l
    res_pf <- allf(res_phi)$f
    res <- list(res_py = res_py, res_pf = res_pf)
    return(res)
}
prml_outB_lp_f <- function(xs_bn, xs_b, r_b = 0.5, s_b = 2e-10, mu_u = "max",
                           n_gq = 20, n_mu = 100, n_per = 100) {

    r_B <- sum(xs_b) + r_b
    s_B <- length(xs_b) + s_b
    imu_b <- mean(xs_b)
    # set the outside B range: max + 2sigma
    if (mu_u == "max") {
        mu_u <- max(max(xs_bn) + 2 * sd(xs_bn), max(xs_b) + 2 * sd(xs_b))
    }
    n <- length(xs_bn)
    indx <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })

    # guass.quad
    out <- gauss.quad(n_gq)

    allf <- function(phi) {
        a <- 0
        b <- 1
        xs_gq <- out$nodes * (b - a) / 2 + (a + b) / 2  #change interval
        z <- seq(a, b, length.out = n_mu)

        ## Initial Guess:

        f <- rep(1 / (b - a), n_mu)  #uniform
        f[1] <- 0
        d_f <- 0
        d_fgq <- rep(0, n_gq)
        f_gq <- rep(1 / (b - a), n_gq)

        ## weight seq
        w <- 1 / (1:n + 1)

        ## Recursion

        Xs_bn_p <- sapply(1:n_per, function(i) {
            xs_bn[indx[, i]]
        })  #each col is a permutation
        res_t <- map(1:n_per, function(ix) {

            m_x <- 0
            grad <- 0
            for (j in 1:n) {
                # step1:prml
                px_grid <- map_dbl(xs_gq,
                                   ~dpois(x = Xs_bn_p[j, ix], lambda = .x * (mu_u - phi) + phi))
                mi_x <- ((b - a) / 2 * sum(out$weights * px_grid * f_gq))
                poi_x <- map_dbl(z,
                                 ~dpois(x = Xs_bn_p[j, ix], lambda = .x * (mu_u - phi) + phi))
                # step2:gradient
                d_poi_x <- poi_x *
                    (Xs_bn_p[j, ix] * (1 - z) / (phi + z * (mu_u - phi)) - (1 - z))
                d_px_grid <- px_grid *
                    (Xs_bn_p[j, ix] * (1 - xs_gq) / (phi + xs_gq * (mu_u - phi)) - (1 - xs_gq))
                G1 <- d_poi_x * f + poi_x * d_f
                G1_gq <- d_px_grid * f_gq + px_grid * d_fgq
                G2 <- poi_x * f
                G2_gq <- px_grid * f_gq
                d_logmi_x <- (b - a) / 2 * sum(out$weights * G1_gq) / mi_x
                # step3:update f,gr_f
                f <- (1 - w[j]) * f + w[j] * poi_x * f / mi_x
                f_gq <- (1 - w[j]) * f_gq + w[j] * px_grid * f_gq / mi_x
                d_f <- (1 - w[j]) * d_f + w[j] * (G1 - G2 * d_logmi_x) / mi_x
                d_fgq <- (1 - w[j]) * d_fgq + w[j] * (G1_gq - G2_gq * d_logmi_x) / mi_x
                # store the result
                m_x <- log(mi_x) + m_x
                grad <- grad + d_logmi_x
            }
            return(list(f = f, m_x = m_x, grad = grad))
        })
        pf <- sapply(res_t, function(res_t) {
            res_t$f
        }) %>% apply(., 1, mean)
        pm_x <- map_dbl(res_t, ~.x$m_x) %>% logmean(.)
        pgrad <- map_dbl(res_t, function(res_t) {
            res_t$grad
        }) %>% logmean(.)
        return(list(f = pf, m_x = pm_x, grad = pgrad))
    }

    obj <- function(theta) {
        phi <- exp(theta) * mu_u / (1 + exp(theta))
        p_mua <- dgamma(x = phi, shape = r_B, rate = s_B, log = TRUE)
        return(-(allf(phi)$m_x + p_mua))
    }
    grad <- function(theta) {
        phi <- exp(theta) * mu_u / (1 + exp(theta))
        g_mua <- ((r_B - 1) / phi - s_B)
        prml <- allf(phi)$grad
        res <- -(prml + g_mua) * phi * (1 - phi / mu_u)
        return(res)
    }
    res <- optim(par = log(imu_b / (mu_u - imu_b)), fn = obj, gr = grad,
                 method = "BFGS", hessian = TRUE)
    # optim is for minimization
    phi_r <- exp(res$par) * mu_u / (1 + exp(res$par))
    H <- (res$hessian * mu_u + grad(res$par) * (2 * phi_r - mu_u)) *
        (phi_r * (mu_u - phi_r))^(-2) * mu_u
    logmx <- 0.5 * log(2 * pi) + 0.5 * log((H)^(-1)) - res$value
    res_py <- logmx
    res_pf <- allf(phi_r)$f
    res <- list(res_py = res_py, res_pf = res_pf)
    return(res)
}
prml_sin_lp_f <- function(xs_bn, xs_a, r_a = 0.5, s_a = 2e-10, n_gq = 20) {

    r_A <- sum(xs_a) + r_a
    s_A <- length(xs_a) + s_a
    n <- length(xs_bn)
    r <- sum(xs_bn) + r_A
    s <- n + s_A
    logmx <- -sum(lgamma(xs_bn + 1)) + lgamma(r) - r * log(s) + r_A * log(s_A) - lgamma(r_A)
    res_py <- logmx
    return(res_py)
}
prml_int_lp_f <- function(xs_bn, xs_a, xs_b, e = 0, r_a = 0.5, s_a = 2e-10,
                          r_b = 0.5, s_b = 2e-10, n_gq = 20, n_mu = 100, n_per = 100) {

    r_A <- sum(xs_a) + r_a
    s_A <- length(xs_a) + s_a
    r_B <- sum(xs_b) + r_b
    s_B <- length(xs_b) + s_b
    imu_a <- mean(xs_a)
    imu_b <- mean(xs_b)

    n <- length(xs_bn)
    indx <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })

    # guass.quad
    out <- gauss.quad(n_gq)

    allf <- function(phi) {
        a <- e
        b <- 1 - e
        xs_gq <- out$nodes * (b - a) / 2 + (a + b) / 2  #change interval
        z <- seq(a, b, length.out = n_mu)

        ## Initial Guess:

        f <- rep(1 / (b - a), n_mu)
        f_gq <- rep(1 / (b - a), n_gq)

        d_f <- matrix(0, nrow = 2, ncol = n_mu)
        d_fgq <- matrix(0, nrow = 2, ncol = n_gq)
        ## weight seq
        w <- 1 / (1:n + 1)

        ## Recursion
        Xs_bn_p <- sapply(1:n_per, function(i) {
            xs_bn[indx[, i]]
        })  #each col is a permutation
        res_t <- map(1:n_per, function(ix) {
            m_x <- 0
            grad <- c(0, 0)
            for (j in 1:n) {
                # step1:prml
                px_grid <- map_dbl(xs_gq,
                                   ~dpois(x = Xs_bn_p[j, ix],
                                          lambda = phi[1] + .x * (phi[2] - phi[1])))
                mi_x <- ((b - a) / 2 * sum(out$weights * px_grid * f_gq))
                poi_x <- map_dbl(z,
                                 ~dpois(x = Xs_bn_p[j, ix],
                                        lambda = phi[1] + .x * (phi[2] - phi[1])))
                # step2:gradient
                d_poi_x_a <- poi_x *
                    (Xs_bn_p[j, ix] * (1 - z) / (phi[1] + z * (phi[2] - phi[1])) - (1 - z))
                d_poi_x_b <- poi_x *
                    (Xs_bn_p[j, ix] * z / (phi[1] + z * (phi[2] - phi[1])) - z)
                d_px_grid_a <- px_grid *
                    (Xs_bn_p[j, ix] * (1 - xs_gq) / (phi[1] + xs_gq * (phi[2] - phi[1])) -
                         (1 - xs_gq))
                d_px_grid_b <- px_grid *
                    (Xs_bn_p[j, ix] * xs_gq / (phi[1] + xs_gq * (phi[2] - phi[1])) - xs_gq)
                G_a <- d_poi_x_a * f + poi_x * d_f[1, ]
                G_b <- d_poi_x_b * f + poi_x * d_f[2, ]
                G_a_gq <- d_px_grid_a * f_gq + px_grid * d_fgq[1, ]
                G_b_gq <- d_px_grid_b * f_gq + px_grid * d_fgq[2, ]
                G2 <- poi_x * f
                G2_gq <- px_grid * f_gq
                d_logmi_x_a <- (b - a) / 2 * sum(out$weights * G_a_gq) / mi_x
                d_logmi_x_b <- (b - a) / 2 * sum(out$weights * G_b_gq) / mi_x
                d_logmi_x <- c(d_logmi_x_a, d_logmi_x_b)
                # step3:update f,gr_f
                d_f[1, ] <- (1 - w[j]) * d_f[1, ] + w[j] * (G_a - G2 * d_logmi_x_a) / mi_x
                d_f[2, ] <- (1 - w[j]) * d_f[2, ] + w[j] * (G_b - G2 * d_logmi_x_b) / mi_x
                f <- (1 - w[j]) * f + w[j] * poi_x * f / mi_x
                d_fgq[1, ] <- (1 - w[j]) * d_fgq[1, ] +
                    w[j] * (G_a_gq - G2_gq * d_logmi_x_a) / mi_x
                d_fgq[2, ] <- (1 - w[j]) * d_fgq[2, ] +
                    w[j] * (G_b_gq - G2_gq * d_logmi_x_b) / mi_x
                f_gq <- (1 - w[j]) * f_gq + w[j] * px_grid * f_gq / mi_x
                # store the result
                m_x <- log(mi_x) + m_x
                grad <- grad + d_logmi_x
            }
            return(list(f = f, m_x = m_x, grad = grad))
        })
        pf <- sapply(res_t, function(res_t) {
            res_t$f
        }) %>% apply(., 1, mean)
        pm_x <- map_dbl(res_t, ~.x$m_x) %>% logmean(.)
        pgrad <- sapply(res_t, function(res_t) {
            res_t$grad
        }) %>% apply(., 1, logmean)
        return(list(f = pf, m_x = pm_x, grad = pgrad))
    }
    obj <- function(theta) {
        phi <- exp(sort(theta))
        p_mua <- dgamma(x = phi[1], shape = r_A, rate = s_A, log = TRUE)
        p_mub <- dgamma(x = phi[2], shape = r_B, rate = s_B, log = TRUE)
        return(-(allf(phi)$m_x + p_mua + p_mub))
    }
    grad <- function(theta) {
        phi <- exp(sort(theta))
        g_mua <- (r_A - 1) / phi[1] - s_A
        g_mub <- (r_B - 1) / phi[2] - s_B
        prml <- allf(phi)$grad
        res <- -(prml + c(g_mua, g_mub)) * phi
        return(res[order(theta)])
    }
    res <- optim(par = c(log(imu_a), log(imu_b)), fn = obj, gr = grad,
                 method = "BFGS", hessian = TRUE)
    # optim is for minimization
    H <- matrix(0, nrow = 2, ncol = 2)
    res_par <- res$par %>% sort(.)
    H[1, 1] <- (-res$hessian[1, 1] - grad(res_par)[1]) * exp(-2 * res_par[1])
    H[2, 2] <- (-res$hessian[2, 2] - grad(res_par)[2]) * exp(-2 * res_par[2])
    H[1, 2] <- -res$hessian[1, 2] * exp(-sum(res_par))
    H[2, 1] <- -res$hessian[2, 1] * exp(-sum(res_par))

    logmx <- log(2 * pi) - 0.5 * log(det(-H)) - res$value
    res_py <- logmx
    res_phi <- exp(res_par)
    res_pf <- allf(res_phi)$f
    res <- list(res_py = res_py, res_pf = res_pf)
    return(res)
}
prml_mix_lp_f <- function(xs_bn, xs_a, xs_b, e = 0, r_a = 0.5, s_a = 2e-10,
                          r_b = 0.5, s_b = 2e-10, n_gq = 20, n_mu = 100, n_per = 100) {

    r_A <- sum(xs_a) + r_a
    s_A <- length(xs_a) + s_a
    r_B <- sum(xs_b) + r_b
    s_B <- length(xs_b) + s_b
    imu_a <- mean(xs_a)
    imu_b <- mean(xs_b)
    del <- 2 * e

    n <- length(xs_bn)
    indx <- sapply(1:n_per, function(x) {
        sample.int(n = n, size = n, replace = FALSE)
    })

    allf <- function(phi) {
        a <- 0
        b <- 1

        z <- c(a, b)

        ## Initial Guess:
        f <- c(0.5, 0.5)
        d_f <- matrix(0, nrow = 2, ncol = 2)

        ## weight seq
        w <- 1 / (1:n + 1)

        ## Recursion
        Xs_bn_p <- sapply(1:n_per, function(i) {
            xs_bn[indx[, i]]
        })  #each col is a permutation
        res_t <- map(1:n_per, function(ix) {
            m_x <- 0
            grad <- c(0, 0)
            for (j in 1:n) {
                # step1:prml
                poi_x <- map_dbl(z,
                                 ~(1 - del) *
                                     (dpois(x = Xs_bn_p[j, ix],
                                            lambda = phi[1] + .x * (phi[2] - phi[1]))) +
                                     del / 2 * (dpois(x = Xs_bn_p[j, ix], lambda = phi[1]) +
                                                    dpois(x = Xs_bn_p[j, ix],
                                                          lambda = phi[2])))
                mi_x <- sum(poi_x * f)

                # step2:gradient
                d_poi_x_a <- poi_x *
                    (Xs_bn_p[j, ix] * (1 - z) / (phi[1] + z * (phi[2] - phi[1])) - (1 - z)) *
                    (1 - del) + del / 2 * poi_x * (Xs_bn_p[j, ix] / (phi[1]) - (1))
                d_poi_x_b <- poi_x *
                    (Xs_bn_p[j, ix] * z / (phi[1] + z * (phi[2] - phi[1])) - z) * (1 - del) +
                    del / 2 * poi_x * (Xs_bn_p[j, ix] / (phi[2]) - (1))
                G_a <- d_poi_x_a * f + poi_x * d_f[1, ]
                G_b <- d_poi_x_b * f + poi_x * d_f[2, ]
                G2 <- poi_x * f
                d_logmi_x_a <- sum(G_a) / mi_x
                d_logmi_x_b <- sum(G_b) / mi_x
                d_logmi_x <- c(d_logmi_x_a, d_logmi_x_b)
                # step3:update f,gr_f
                d_f[1, ] <- (1 - w[j]) * d_f[1, ] + w[j] * (G_a - G2 * d_logmi_x_a) / mi_x
                d_f[2, ] <- (1 - w[j]) * d_f[2, ] + w[j] * (G_b - G2 * d_logmi_x_b) / mi_x
                f <- (1 - w[j]) * f + w[j] * poi_x * f / mi_x

                # store the result
                m_x <- log(mi_x) + m_x
                grad <- grad + d_logmi_x
            }
            return(list(f = f, m_x = m_x, grad = grad))
        })
        pf <- sapply(res_t, function(res_t) {
            res_t$f
        }) %>% apply(., 1, mean)
        pm_x <- map_dbl(res_t, ~.x$m_x) %>% logmean(.)
        pgrad <- sapply(res_t, function(res_t) {
            res_t$grad
        }) %>% apply(., 1, logmean)
        return(list(f = pf, m_x = pm_x, grad = pgrad))
    }

    obj <- function(theta) {
        phi <- exp(sort(theta))
        p_mua <- dgamma(x = phi[1], shape = r_A, rate = s_A, log = TRUE)
        p_mub <- dgamma(x = phi[2], shape = r_B, rate = s_B, log = TRUE)
        return(-(allf(phi)$m_x + p_mua + p_mub))
    }
    grad <- function(theta) {
        phi <- exp(sort(theta))
        g_mua <- (r_A - 1) / phi[1] - s_A
        g_mub <- (r_B - 1) / phi[2] - s_B
        prml <- allf(phi)$grad
        res <- -(prml + c(g_mua, g_mub)) * phi
        return(res[order(theta)])
    }
    res <- optim(par = c(log(imu_a), log(imu_b)), fn = obj, gr = grad,
                 method = "BFGS", hessian = TRUE)
    # optim is for minimization
    H <- matrix(0, nrow = 2, ncol = 2)
    res_par <- res$par %>% sort(.)
    H[1, 1] <- (-res$hessian[1, 1] - grad(res_par)[1]) * exp(-2 * res_par[1])
    H[2, 2] <- (-res$hessian[2, 2] - grad(res_par)[2]) * exp(-2 * res_par[2])
    H[1, 2] <- -res$hessian[1, 2] * exp(-sum(res_par))
    H[2, 1] <- -res$hessian[2, 1] * exp(-sum(res_par))

    logmx <- log(2 * pi) - 0.5 * log(det(-H)) - res$value

    res_py <- logmx
    res_phi <- exp(res_par)
    res_pf <- allf(res_phi)$f
    res <- list(res_py = res_py, res_pf = res_pf)
    return(res)
}

#' PRML classifier together with density estimation of mixing density.
#' @param xs_bn A vector. Spike counts of repeated dual-stimuli trial data AB.
#' @param xs_a A vector. Spike counts of repeated single-stimulus trial data A.
#' @param xs_b A vector. Spike counts of repeated single-stimulus trial data B.
#' @param mu_l A number. Lower bound of spike counts. "min" by default. Indicating \eqn{ max(0,\underset{j=A,B,AB}{\min}(\min(Y_j)-2\text{std}(Y_j))) }
#' @param mu_u A number. Upper bound of spike counts. "max" by default. Indicating \eqn{ \underset{j=A,B,AB}{\max}(\max(Y_j)+2\text{std}(Y_j)) }
#' @param e A number. 0 by default. Shringkage on the domain and meansurement of mixing density f under the Intermediate and Mixture hypothese.
#' @param r_a A number. The parameter in gamma prior of spike rate mu_A. rate. Jeffereys' prior by default.
#' @param s_a A number. The parameter in gamma prior of spike rate mu_A. shape. Jeffereys' prior by default.
#' @param r_b A number. The parameter in gamma prior of spike rate mu_B. rate. Jeffereys' prior by default.
#' @param s_b A number. The parameter in gamma prior of spike rate mu_B. shape. Jeffereys' prior by default.
#' @param n_gq A number. 20 by default. Number of grids in Gaussion quadrature.
#' @param n_mu A number. 100 by default. The number of grids used to represent the pdf of f.
#' @param n_per A number. 100 by default. Permutation of likihood estimation to obtain the order-invariant estimator.
#' @return A list.
#' \describe{
#'   \item{out1}{Result of \link{prml_tests}}
#'   \item{out2}{density estimation of mixing density f under Mixture, Intermediate, OutsideA, OutsideB hypotheses}
#' }
#' @seealso \code{\link{prml_classifier}}
prml_classifier_f <- function(xs_bn, xs_a, xs_b, mu_l = 0, mu_u = 180, e = 0, r_a = 0.5, s_a = 2e-10,
                       r_b = 0.5, s_b = 2e-10, n_gq = 20, n_mu = 100, n_per = 100) {
    if (mean(xs_a) <= mean(xs_b)) {
        outA <- prml_outA_lp_f(xs_bn, xs_a, r_a, s_a, mu_l, n_gq, n_mu, n_per)
        outB <- prml_outB_lp_f(xs_bn, xs_b, r_b, s_b, mu_u, n_gq, n_mu, n_per)
        sinA <- prml_sin_lp_f(xs_bn, xs_a, r_a, s_a, n_gq)
        sinB <- prml_sin_lp_f(xs_bn, xs_b, r_b, s_b, n_gq)
        int <- prml_int_lp_f(xs_bn, xs_a, xs_b, e, r_a, s_a, r_b, s_b, n_gq, n_mu, n_per)
        mix <- prml_mix_lp_f(xs_bn, xs_a, xs_b, e, r_a, s_a, r_b, s_b, n_gq, n_mu, n_per)
    } else {
        outA <- prml_outA_lp_f(xs_bn, xs_b, r_b, s_b, mu_l, n_gq, n_mu, n_per)
        outB <- prml_outB_lp_f(xs_bn, xs_a, r_a, s_a, mu_u, n_gq, n_mu, n_per)
        sinA <- prml_sin_lp_f(xs_bn, xs_b, r_b, s_b, n_gq)
        sinB <- prml_sin_lp_f(xs_bn, xs_a, r_a, s_a, n_gq)
        int <- prml_int_lp_f(xs_bn, xs_b, xs_a, e, r_b, s_b, r_a, s_a, n_gq, n_mu, n_per)
        mix <- prml_mix_lp_f(xs_bn, xs_b, xs_a, e, r_b, s_b, r_a, s_a, n_gq, n_mu, n_per)
    }

    mxs <- c(mix$res_py, int$res_py, logmean(c(outA$res_py, outB$res_py)),
             logmean(c(sinA, sinB)))
    post.prob <- exp(mxs - logsum(mxs))
    # models <- c('Mix', 'Int', 'Out', 'Sing')
    win.model <- which.max(post.prob)
    out <- list(post.prob = post.prob, win.model = win.model, mix_pf = mix$res_pf,
                int_pf = int$res_pf, outA_pf = outA$res_pf, outB_pf = outB$res_pf)
    return(out)
}

#' PRML filter and PRML classifier together with density estimation of mixing density.
#' @param xA A vector. Spike counts of repeated dual-stimuli trial data AB.
#' @param xB A vector. Spike counts of repeated single-stimulus trial data A.
#' @param xAB A vector. Spike counts of repeated single-stimulus trial data B.
#' @param labels A vector. labels for the trials.
#' @param remove.zeros A logical value. Whether to remove 0s in spike counts.
#' @param mu_l A number. Lower bound of spike counts. "min" by default. Indicating \eqn{ \text{max}(0,\underset{j=A,B,AB}{\text{min}}(\text{min}(Y_j)-2\text{std}(Y_j))) }
#' @param mu_u A number. Upper bound of spike counts. "max" by default. Indicating \eqn{ \underset{ j=A,B,AB }{\text{max}}(\text{max}(Y_j)+2\text{std}(Y_j)) }
#' @param e A number. 0 by default. Shringkage on the domain and meansurement of mixing density f under the Intermediate and Mixture hypothese.
#' @param gamma.pars A length 2 vector. The shape and rate of gamma prior for spike rate mu_A and mu_B. Jeffereys' prior by default.
#' @param n_gq A number. 20 by default. Number of grids in Gaussion quadrature.
#' @param n_mu A number. 100 by default. The number of grids used to represent the pdf of f.
#' @param n_per A number. 100 by default. Permutation of likihood estimation to obtain the order-invariant estimator.
#' @param alpha 0.5 by default. (For PRML filter) The range of the spike counts estimator \eqn{ [Y_{0.25}-\alpha \text{IQR},Y_{0.75}+\alpha \text{IQR}] }
#' @return A list.
#' \describe{
#'   \item{out1}{Result of \link{prml_tests}}
#'   \item{out2}{density estimation of mixing density f under Mixture, Intermediate, OutsideA, OutsideB hypotheses}
#' }
#' @seealso \code{\link{prml_tests}}
prml_tests_f <- function(xA, xB, xAB, labels = c("A", "B", "AB"), remove.zeros = FALSE,
                         mu_l = "min", mu_u = "max", e = 0, gamma.pars = c(0.5, 2e-10),
                         n_gq = 20, n_mu = 100, n_per = 100, alpha = 0.5) {

    a <- r_a <- r_b <- gamma.pars[1]
    b <- s_a <- s_b <- gamma.pars[2]

    if (remove.zeros) {
        xA <- xA[xA != 0]
        xB <- xB[xB != 0]
        xAB <- xAB[xAB != 0]
    }

    nA <- length(xA)
    nB <- length(xB)
    nAB <- length(xAB)
    if (nA == 0 | nB == 0)
        stop("not enough data in single sound")

    ## how different are the two pure trials?

    two.poi.ibf <- Vectorize(function(i, j)
        return(log.pm(xA[-i], xA[i] + a, 1 + b) +
                   log.pm(xB[-j], xB[j] + a, 1 + b) -
                   log.pm(c(xA[-i], xB[-j]), xA[i] + xB[j] + a, 2 + b)))
    lbf.pure <- mean(c(outer(1:length(xA), 1:length(xB), two.poi.ibf)))

    pvls <- c(prml_filter(xA, n_gq, n_per, alpha),
              prml_filter(xB, n_gq, n_per, alpha))
    prmlmix_res <- prml_classifier_f(xAB, xA, xB, mu_l, mu_u, e,
                              r_a, s_a, r_b, s_b, n_gq, n_mu, n_per)
    out1 <- list(separation.logBF = lbf.pure, post.prob = prmlmix_res$post.prob,
                 win.model = prmlmix_res$win.model, prml.filter.bf = pvls,
                 samp.sizes = c(nA, nB, nAB))
    out2 <- list(mix_pf = prmlmix_res$mix_pf, int_pf = prmlmix_res$int_pf,
                 outA_pf = prmlmix_res$outA_pf, outB_pf = prmlmix_res$outB_pf)
    out <- list(out1 = out1, out2 = out2)
    return(out)
}


