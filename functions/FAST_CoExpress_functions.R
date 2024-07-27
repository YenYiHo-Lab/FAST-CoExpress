library(GJRM)
library(gamlss.dist)
library(copula)
library(pracma)
library(mvtnorm)
library(pscl)



# sigmoid <- function(x) {
#   exp(x) / (exp(x) + 1)
# }
# 
# logit <- function(x) {
#   log(x / (1 - x))
# }


copnb.sim.individual <- function(p1, p2, mu1, mu2, sig1, sig2, rho, copulafun){
  
  Cop <- do.call(copulafun, args = list(param=rho)) 
  
  spec <- mvdc(copula = Cop, c("NBI", "NBI"),
               list(list( mu = mu1, sigma = sig1 ),
                    list( mu = mu2, sigma = sig2 )) )
  
  y1y2 <- rMvdc(1, spec) * rbinom(2, size=1, prob = 1-c(p1, p2) )
  
  y1y2
  
}


copnb.sim <- function(p1, p2, mu1, mu2, sig1, sig2, rho, copula){
  
  if (copula == "Gaussian"){
    copulafun <- "normalCopula"
  } else if (copula == "Clayton"){
    copulafun <- "claytonCopula"
  } else if (copula == "Frank"){
    copulafun <- "frankCopula"
  } else if (copula == "Joe"){
    copulafun <- "joeCopula"
  } else if (copula == "Gumbel"){
    copulafun <- "gumbelCopula"
  } else {
    stop("copula needs to be one of \"Gaussian\", \"Clayton\", \"Frank\", \"Joe\", \"Gumbel\"")
  }
  
  param.mat <- cbind(p1, p2, mu1, mu2, sig1, sig2, rho)
  
  y1y2 <- t(apply(param.mat, 1, function(a) copnb.sim.individual(a[1], a[2], a[3], a[4], a[5], a[6], a[7], copulafun)) )
  
  colnames(y1y2) <- c("y1", "y2")
  
  y1y2
  
}






## Modifying function from rCoCoA to give standard errors ## 

get_params_reml <- function (dat_xyz, init_mle = TRUE, rescale_zt = FALSE, optim_type = "optim") 
{
  
  n = nrow(dat_xyz)
  p = ncol(dat_xyz) - 1
  if (init_mle) {
    params_mle = get_params_mle(dat_xyz)$params
    start_vals = params_mle[-(1:2)]
  }
  else {
    start_vals = c(sd(dat_xyz$X), sd(dat_xyz$Y), rep(0.1, 
                                                     p))
  }
  dat_xyz2 <- dat_xyz
  if (rescale_zt) {
    xy_ind = which(colnames(dat_xyz2) %in% c("X", "Y"))
    means_zt = apply(dat_xyz2[, -xy_ind], 2, mean)
    sds_zt = apply(dat_xyz2[, -xy_ind], 2, sd)
    dat_xyz2[, -xy_ind] <- apply(dat_xyz2[, -xy_ind], 2, 
                                 scale2)
  }
  t0 = Sys.time()
  if (optim_type == "optim") {
    opt = optim(par = start_vals, loglik_reml, dat_xyz = dat_xyz2, 
                method = "L-BFGS-B", lower = c(0.001, 0.001, rep(-5, 
                                                                 p)), upper = c(Inf, Inf, rep(5, p)), hessian = TRUE, 
                control = list(fnscale = -1))
    ll = 1 * (opt$value)
  }
  else {
    opt = optimr::optimr(par = start_vals, loglik_reml, 
                         dat_xyz = dat_xyz2, method = "Rvmmin", lower = c(0.001, 
                                                                          0.001, rep(-5, p)), upper = c(Inf, Inf, rep(5, 
                                                                                                                      p)), hessian = TRUE, control = list(maximize = TRUE))
    ll = -1 * (opt$value)
  }
  
  se_use <- sqrt(-diag(solve(opt$hessian)))
  se_use <- se_use[(length(se_use)-(p-1)):length(se_use)]
  
  t1 = Sys.time() - t0
  pars = opt$par
  final_beta_hat = rCoCoA:::reml_get_beta(pars, dat_xyz)
  param_names <- rCoCoA:::grk_letters[1:(length(pars) - 2)]
  var_names <- c("1", colnames(dat_xyz)[-(1:2)])
  if (rescale_zt) {
    alpha0 = pars[3]
    beta_gamma0 = pars[-(1:3)]
    alpha2 = alpha0 - sum((means_zt * beta_gamma0)/sds_zt)
    beta_gamma2 = beta_gamma0/sds_zt
    out2 = data.frame(Coefficient = c("mu_x", "mu_y", "sigma_x", 
                                      "sigma_y", param_names), Variable = c("X", "Y", 
                                                                            "X", "Y", var_names), estimate = c(final_beta_hat, 
                                                                                                               pars[1:2], alpha2, beta_gamma2))
    out2 = dplyr::mutate(out2, san.se = as.numeric(NA), 
                         wald = as.numeric(NA), p = as.numeric(NA), Model = c(rep("mean", 
                                                                                  2), rep("scale", 2), rep("correlation", length(param_names))))
  }
  else {
    out2 = data.frame(Coefficient = c("mu_x", "mu_y", "sigma_x", 
                                      "sigma_y", param_names), Variable = c("X", "Y", 
                                                                            "X", "Y", var_names), estimate = c(final_beta_hat, 
                                                                                                               pars))
    out2 = dplyr::mutate(out2, san.se = as.numeric(NA), 
                         wald = as.numeric(NA), p = as.numeric(NA), Model = c(rep("mean", 
                                                                                  2), rep("scale", 2), rep("correlation", length(param_names))))
    
    rhorows <- (nrow(out2)-(p-1)):nrow(out2)
    
    out2[rhorows, "san.se"] <- se_use
    
    out2[rhorows, "wald"] <- apply(out2[rhorows, 3:4], 1, function(a) sign(a[1] + 1.96*a[2]) == sign(a[1] - 1.96*a[2]))
    
    
  }
  return(list(params = as.matrix(pars), out2 = out2, runtime = t1, 
              loglik = ll))
}





overall.svG <- function (formula, data, ngc, margins, M, vo, gam1, gam2, type = "copR", 
          inde = NULL, c.gam2 = NULL, gam3 = NULL, knots = NULL) 
{
  X3 <- X4 <- X5 <- X6 <- X7 <- X8 <- X9 <- NULL
  X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- X9.d2 <- NULL
  gp3 <- gp4 <- gp5 <- gp6 <- gp7 <- gp8 <- gp9 <- NULL
  gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- gam9 <- NULL
  l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- 0
  sp3 <- sp4 <- sp5 <- sp6 <- sp7 <- sp8 <- sp9 <- NULL
  X3s <- X4s <- X5s <- X6s <- X7s <- X8s <- X9s <- NULL
  Sl.sf2 <- Sl.sf3 <- NULL
  if (type == "ROY") {
    X4 <- X5 <- X6 <- X7 <- X8 <- X9 <- X4s <- X5s <- X6s <- X7s <- X8s <- X9s <- matrix(1, 
                                                                                         vo$n, 1)
    gp3 <- gp4 <- gp5 <- gp6 <- gp7 <- gp8 <- gp9 <- 0
    if (margins[2] %in% c(M$bl, M$m1d) && margins[3] %in% 
        c(M$bl, M$m1d)) {
      formula.eq4 <- formula[[4]]
      nad <- "theta12"
      formula.eq4 <- as.formula(paste(nad, "~", formula.eq4[2], 
                                      sep = ""))
      set.seed(1)
      theta12 <- rnorm(vo$n, vo$i.rho1, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                  subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq5 <- formula[[5]]
      nad <- "theta13"
      formula.eq5 <- as.formula(paste(nad, "~", formula.eq5[2], 
                                      sep = ""))
      set.seed(1)
      theta13 <- rnorm(vo$n, vo$i.rho2, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                  subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      l.sp4 <- length(gam4$sp)
      l.sp5 <- length(gam5$sp)
      if (l.sp4 != 0) {
        ngc <- 2
        while (any(round(summary(gam4)$edf, 1) > 1)) {
          gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                        1, subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp5 != 0) {
        ngc <- 2
        while (any(round(summary(gam5)$edf, 1) > 1)) {
          gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                        1, subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]
      X5 <- model.matrix(gam5)
      X5.d2 <- dim(X5)[2]
      if (l.sp4 != 0) 
        sp4 <- gam4$sp
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf
      if (l.sp5 != 0) 
        sp5 <- gam5$sp
      environment(gam5$formula) <- environment(gam2$formula)
      gp5 <- gam5$nsdf
      X4s <- try(predict.gam(gam4, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X4s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X5s <- try(predict.gam(gam5, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X5s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      start.v <- c(gam1$coefficients, gam2$coefficients, 
                   gam3$coefficients, gam4$coefficients, gam5$coefficients)
    }
    if (margins[2] %in% c(M$m2d, M$m2) && margins[3] %in% 
        c(M$m2d, M$m2)) {
      formula.eq4 <- formula[[4]]
      nad <- "sigma2"
      formula.eq4 <- as.formula(paste(nad, "~", formula.eq4[2], 
                                      sep = ""))
      set.seed(1)
      sigma2 <- rnorm(vo$n, vo$log.sig1, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                  subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq5 <- formula[[5]]
      nad <- "sigma3"
      formula.eq5 <- as.formula(paste(nad, "~", formula.eq5[2], 
                                      sep = ""))
      set.seed(1)
      sigma3 <- rnorm(vo$n, vo$log.sig2, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                  subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq6 <- formula[[6]]
      nad <- "theta12"
      formula.eq6 <- as.formula(paste(nad, "~", formula.eq6[2], 
                                      sep = ""))
      set.seed(1)
      theta12 <- rnorm(vo$n, vo$i.rho1, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                  subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq7 <- formula[[7]]
      nad <- "theta13"
      formula.eq7 <- as.formula(paste(nad, "~", formula.eq7[2], 
                                      sep = ""))
      set.seed(1)
      theta13 <- rnorm(vo$n, vo$i.rho2, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam7 <- gam(formula.eq7, data = data, gamma = ngc, 
                  subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      l.sp4 <- length(gam4$sp)
      l.sp5 <- length(gam5$sp)
      l.sp6 <- length(gam6$sp)
      l.sp7 <- length(gam7$sp)
      if (l.sp4 != 0) {
        ngc <- 2
        while (any(round(summary(gam4)$edf, 1) > 1)) {
          gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                        1, subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp5 != 0) {
        ngc <- 2
        while (any(round(summary(gam5)$edf, 1) > 1)) {
          gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                        1, subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp6 != 0) {
        ngc <- 2
        while (any(round(summary(gam6)$edf, 1) > 1)) {
          gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                        1, subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp7 != 0) {
        ngc <- 2
        while (any(round(summary(gam7)$edf, 1) > 1)) {
          gam7 <- gam(formula.eq7, data = data, gamma = ngc + 
                        1, subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]
      X5 <- model.matrix(gam5)
      X5.d2 <- dim(X5)[2]
      X6 <- model.matrix(gam6)
      X6.d2 <- dim(X6)[2]
      X7 <- model.matrix(gam7)
      X7.d2 <- dim(X7)[2]
      if (l.sp4 != 0) 
        sp4 <- gam4$sp
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf
      if (l.sp5 != 0) 
        sp5 <- gam5$sp
      environment(gam5$formula) <- environment(gam2$formula)
      gp5 <- gam5$nsdf
      if (l.sp6 != 0) 
        sp6 <- gam6$sp
      environment(gam6$formula) <- environment(gam2$formula)
      gp6 <- gam6$nsdf
      if (l.sp7 != 0) 
        sp7 <- gam7$sp
      environment(gam7$formula) <- environment(gam2$formula)
      gp7 <- gam7$nsdf
      X4s <- try(predict.gam(gam4, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X4s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X5s <- try(predict.gam(gam5, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X5s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X6s <- try(predict.gam(gam6, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X6s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X7s <- try(predict.gam(gam7, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X7s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      start.v <- c(gam1$coefficients, gam2$coefficients, 
                   gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                   gam6$coefficients, gam7$coefficients)
    }
    if (margins[2] %in% c(M$m3) && margins[3] %in% c(M$m3)) {
      formula.eq4 <- formula[[4]]
      nad <- "sigma2"
      formula.eq4 <- as.formula(paste(nad, "~", formula.eq4[2], 
                                      sep = ""))
      set.seed(1)
      sigma2 <- rnorm(vo$n, vo$log.sig1, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                  subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq5 <- formula[[5]]
      nad <- "sigma3"
      formula.eq5 <- as.formula(paste(nad, "~", formula.eq5[2], 
                                      sep = ""))
      set.seed(1)
      sigma3 <- rnorm(vo$n, vo$log.sig2, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                  subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq6 <- formula[[6]]
      nad <- "nu2"
      formula.eq6 <- as.formula(paste(nad, "~", formula.eq6[2], 
                                      sep = ""))
      set.seed(1)
      nu2 <- rnorm(vo$n, vo$log.nu1, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                  subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq7 <- formula[[7]]
      nad <- "nu3"
      formula.eq7 <- as.formula(paste(nad, "~", formula.eq7[2], 
                                      sep = ""))
      set.seed(1)
      nu3 <- rnorm(vo$n, vo$log.nu2, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam7 <- gam(formula.eq7, data = data, gamma = ngc, 
                  subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq8 <- formula[[8]]
      nad <- "theta12"
      formula.eq8 <- as.formula(paste(nad, "~", formula.eq8[2], 
                                      sep = ""))
      set.seed(1)
      theta12 <- rnorm(vo$n, vo$i.rho1, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam8 <- gam(formula.eq8, data = data, gamma = ngc, 
                  subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      formula.eq9 <- formula[[9]]
      nad <- "theta13"
      formula.eq9 <- as.formula(paste(nad, "~", formula.eq9[2], 
                                      sep = ""))
      set.seed(1)
      theta13 <- rnorm(vo$n, vo$i.rho2, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam9 <- gam(formula.eq9, data = data, gamma = ngc, 
                  subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      l.sp4 <- length(gam4$sp)
      l.sp5 <- length(gam5$sp)
      l.sp6 <- length(gam6$sp)
      l.sp7 <- length(gam7$sp)
      l.sp8 <- length(gam8$sp)
      l.sp9 <- length(gam9$sp)
      if (l.sp4 != 0) {
        ngc <- 2
        while (any(round(summary(gam4)$edf, 1) > 1)) {
          gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                        1, subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp5 != 0) {
        ngc <- 2
        while (any(round(summary(gam5)$edf, 1) > 1)) {
          gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                        1, subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp6 != 0) {
        ngc <- 2
        while (any(round(summary(gam6)$edf, 1) > 1)) {
          gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                        1, subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp7 != 0) {
        ngc <- 2
        while (any(round(summary(gam7)$edf, 1) > 1)) {
          gam7 <- gam(formula.eq7, data = data, gamma = ngc + 
                        1, subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp8 != 0) {
        ngc <- 2
        while (any(round(summary(gam8)$edf, 1) > 1)) {
          gam8 <- gam(formula.eq8, data = data, gamma = ngc + 
                        1, subset = vo$inde0, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp9 != 0) {
        ngc <- 2
        while (any(round(summary(gam9)$edf, 1) > 1)) {
          gam9 <- gam(formula.eq9, data = data, gamma = ngc + 
                        1, subset = vo$inde1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]
      X5 <- model.matrix(gam5)
      X5.d2 <- dim(X5)[2]
      X6 <- model.matrix(gam6)
      X6.d2 <- dim(X6)[2]
      X7 <- model.matrix(gam7)
      X7.d2 <- dim(X7)[2]
      X8 <- model.matrix(gam8)
      X8.d2 <- dim(X8)[2]
      X9 <- model.matrix(gam9)
      X9.d2 <- dim(X9)[2]
      if (l.sp4 != 0) 
        sp4 <- gam4$sp
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf
      if (l.sp5 != 0) 
        sp5 <- gam5$sp
      environment(gam5$formula) <- environment(gam2$formula)
      gp5 <- gam5$nsdf
      if (l.sp6 != 0) 
        sp6 <- gam6$sp
      environment(gam6$formula) <- environment(gam2$formula)
      gp6 <- gam6$nsdf
      if (l.sp7 != 0) 
        sp7 <- gam7$sp
      environment(gam7$formula) <- environment(gam2$formula)
      gp7 <- gam7$nsdf
      if (l.sp8 != 0) 
        sp8 <- gam8$sp
      environment(gam8$formula) <- environment(gam2$formula)
      gp8 <- gam8$nsdf
      if (l.sp9 != 0) 
        sp9 <- gam9$sp
      environment(gam9$formula) <- environment(gam2$formula)
      gp9 <- gam9$nsdf
      X4s <- try(predict.gam(gam4, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X4s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X5s <- try(predict.gam(gam5, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X5s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X6s <- try(predict.gam(gam6, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X6s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X7s <- try(predict.gam(gam7, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X7s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X8s <- try(predict.gam(gam8, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X8s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      X9s <- try(predict.gam(gam9, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X9s) == "try-error")) 
        stop("Check that the factor variables' levels\nin the selected sample for the third margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
      start.v <- c(gam1$coefficients, gam2$coefficients, 
                   gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                   gam6$coefficients, gam7$coefficients, gam8$coefficients, 
                   gam9$coefficients)
    }
  }
  if (type != "triv") 
    gam3 <- NULL
  if (type == "triv") {
    formula.eq4 <- formula[[4]]
    nad <- "theta12"
    formula.eq4 <- as.formula(paste(nad, "~", formula.eq4[2], 
                                    sep = ""))
    set.seed(1)
    theta12 <- rnorm(vo$n, vo$theta12, 0.001)
    rm(list = ".Random.seed", envir = globalenv())
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, subset = inde, 
                knots = knots, drop.unused.levels = vo$drop.unused.levels)
    formula.eq5 <- formula[[5]]
    nad <- "theta13"
    formula.eq5 <- as.formula(paste(nad, "~", formula.eq5[2], 
                                    sep = ""))
    set.seed(1)
    theta13 <- rnorm(vo$n, vo$theta13, 0.001)
    rm(list = ".Random.seed", envir = globalenv())
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, subset = inde, 
                knots = knots, drop.unused.levels = vo$drop.unused.levels)
    formula.eq6 <- formula[[6]]
    nad <- "theta23"
    formula.eq6 <- as.formula(paste(nad, "~", formula.eq6[2], 
                                    sep = ""))
    set.seed(1)
    theta23 <- rnorm(vo$n, vo$theta23, 0.001)
    rm(list = ".Random.seed", envir = globalenv())
    gam6 <- gam(formula.eq6, data = data, gamma = ngc, subset = inde, 
                knots = knots, drop.unused.levels = vo$drop.unused.levels)
    l.sp5 <- length(gam5$sp)
    l.sp4 <- length(gam4$sp)
    l.sp6 <- length(gam6$sp)
    if (l.sp4 != 0) {
      ngc <- 2
      while (any(round(summary(gam4)$edf, 1) > 1)) {
        gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                      1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
        ngc <- ngc + 1
        if (ngc > 5) 
          break
      }
    }
    if (l.sp5 != 0) {
      ngc <- 2
      while (any(round(summary(gam5)$edf, 1) > 1)) {
        gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                      1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
        ngc <- ngc + 1
        if (ngc > 5) 
          break
      }
    }
    if (l.sp6 != 0) {
      ngc <- 2
      while (any(round(summary(gam6)$edf, 1) > 1)) {
        gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                      1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
        ngc <- ngc + 1
        if (ngc > 5) 
          break
      }
    }
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2]
    if (l.sp4 != 0) 
      sp4 <- gam4$sp
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf
    if (l.sp5 != 0) 
      sp5 <- gam5$sp
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf
    if (l.sp6 != 0) 
      sp6 <- gam6$sp
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf
    start.v <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, 
                 gam4$coefficients, gam5$coefficients, gam6$coefficients)
  }
  if (type == "biv") {
    if (M$l.flist == 3) {
      formula.eq3 <- formula[[3]]
      nad <- "theta"
      formula.eq3 <- as.formula(paste(nad, "~", formula.eq3[2], 
                                      sep = ""))
      set.seed(1)
      theta <- rnorm(vo$n, vo$i.rho, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                  subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      if (M$Model == "BSS") {
        X3s <- try(predict.gam(gam3, newdata = data[, 
                                                    -dim(data)[2]], type = "lpmatrix"), silent = TRUE)
        if (any(class(X3s) == "try-error")) 
          stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?SemiParBIV for more information.")
      }
      l.sp3 <- length(gam3$sp)
      if (l.sp3 != 0) {
        ngc <- 2
        while (any(round(summary(gam3)$edf, 1) > 1)) {
          gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                        1, subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp3 != 0) 
        sp3 <- gam3$sp
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf
      if (M$Model == "BSS") 
        start.v <- c(gam1$coefficients, c.gam2, gam3$coefficients)
      if (M$Model != "BSS") 
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients)
    }
    if (M$l.flist == 4) {
      formula.eq3 <- formula[[3]]
      formula.eq4 <- formula[[4]]
      nad1 <- "sigma2"
      nad2 <- "theta"
      formula.eq3 <- as.formula(paste(nad1, "~", formula.eq3[2], 
                                      sep = ""))
      formula.eq4 <- as.formula(paste(nad2, "~", formula.eq4[2], 
                                      sep = ""))
      set.seed(1)
      sigma2 <- rnorm(vo$n, vo$log.sig2, 0.001)
      theta <- rnorm(vo$n, vo$i.rho, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                  knots = knots, drop.unused.levels = vo$drop.unused.levels)
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                  knots = knots, drop.unused.levels = vo$drop.unused.levels)
      l.sp3 <- length(gam3$sp)
      l.sp4 <- length(gam4$sp)
      if (l.sp3 != 0) {
        ngc <- 2
        while (any(round(summary(gam3)$edf, 1) > 1)) {
          gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                        1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp4 != 0) {
        ngc <- 2
        while (any(round(summary(gam4)$edf, 1) > 1)) {
          gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                        1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]
      if (l.sp3 != 0) 
        sp3 <- gam3$sp
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf
      if (l.sp4 != 0) 
        sp4 <- gam4$sp
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf
      start.v <- c(gam1$coefficients, gam2$coefficients, 
                   gam3$coefficients, gam4$coefficients)
    }
    if (M$l.flist == 5) {
      formula.eq3 <- formula[[3]]
      formula.eq4 <- formula[[4]]
      formula.eq5 <- formula[[5]]
      nad1 <- "sigma2"
      nad2 <- "nu"
      nad3 <- "theta"
      formula.eq3 <- as.formula(paste(nad1, "~", formula.eq3[2], 
                                      sep = ""))
      formula.eq4 <- as.formula(paste(nad2, "~", formula.eq4[2], 
                                      sep = ""))
      formula.eq5 <- as.formula(paste(nad3, "~", formula.eq5[2], 
                                      sep = ""))
      set.seed(1)
      sigma2 <- rnorm(vo$n, vo$log.sig2, 0.001)
      nu <- rnorm(vo$n, vo$log.nu, 0.001)
      theta <- rnorm(vo$n, vo$i.rho, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                  knots = knots, drop.unused.levels = vo$drop.unused.levels)
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                  knots = knots, drop.unused.levels = vo$drop.unused.levels)
      gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                  knots = knots, drop.unused.levels = vo$drop.unused.levels)
      l.sp3 <- length(gam3$sp)
      l.sp4 <- length(gam4$sp)
      l.sp5 <- length(gam5$sp)
      if (l.sp3 != 0) {
        ngc <- 2
        while (any(round(summary(gam3)$edf, 1) > 1)) {
          gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                        1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp4 != 0) {
        ngc <- 2
        while (any(round(summary(gam4)$edf, 1) > 1)) {
          gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                        1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp5 != 0) {
        ngc <- 2
        while (any(round(summary(gam5)$edf, 1) > 1)) {
          gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                        1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]
      X5 <- model.matrix(gam5)
      X5.d2 <- dim(X5)[2]
      if (l.sp3 != 0) 
        sp3 <- gam3$sp
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf
      if (l.sp4 != 0) 
        sp4 <- gam4$sp
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf
      if (l.sp5 != 0) 
        sp5 <- gam5$sp
      environment(gam5$formula) <- environment(gam2$formula)
      gp5 <- gam5$nsdf
      start.v <- c(gam1$coefficients, gam2$coefficients, 
                   gam3$coefficients, gam4$coefficients, gam5$coefficients)
    }
  }
  if (type == "copR") {
    BivD <- M$BivD
    if (M$surv == TRUE) 
      BivD <- "N"
    if (margins[1] %in% c(M$m2, M$m3) && margins[2] %in% 
        c(M$m2, M$m3) && BivD == "T") {
      if (margins[1] %in% c(M$m2) && margins[2] %in% c(M$m2)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        formula.eq6 <- formula[[6]]
        nad2.1 <- "sigma2.1"
        nad2.2 <- "sigma2.2"
        nad.3 <- "dof"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad.3, "~", 
                                        formula.eq5[2], sep = ""))
        formula.eq6 <- as.formula(paste(nad3, "~", formula.eq6[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        dof <- rnorm(vo$n, vo$dof.st, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        l.sp6 <- length(gam6$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp6 != 0) {
          ngc <- 2
          while (any(round(summary(gam6)$edf, 1) > 1)) {
            gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        if (l.sp6 != 0) 
          sp6 <- gam6$sp
        environment(gam6$formula) <- environment(gam2$formula)
        gp6 <- gam6$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                     gam6$coefficients)
      }
      if (margins[1] %in% c(M$m3) && margins[2] %in% c(M$m3)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        formula.eq6 <- formula[[6]]
        formula.eq7 <- formula[[7]]
        formula.eq8 <- formula[[8]]
        nad2.1 <- "sigma2.1"
        nad2.2 <- "sigma2.2"
        nad.1 <- "nu.1"
        nad.2 <- "nu.2"
        nad.3 <- "dof"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad.1, "~", 
                                        formula.eq5[2], sep = ""))
        formula.eq6 <- as.formula(paste(nad.2, "~", 
                                        formula.eq6[2], sep = ""))
        formula.eq7 <- as.formula(paste(nad.3, "~", 
                                        formula.eq7[2], sep = ""))
        formula.eq8 <- as.formula(paste(nad3, "~", formula.eq8[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        nu.1 <- rnorm(vo$n, vo$log.nu.1, 0.001)
        nu.2 <- rnorm(vo$n, vo$log.nu.2, 0.001)
        dof <- rnorm(vo$n, vo$dof.st, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam7 <- gam(formula.eq7, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam8 <- gam(formula.eq8, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        l.sp6 <- length(gam6$sp)
        l.sp7 <- length(gam7$sp)
        l.sp8 <- length(gam8$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp6 != 0) {
          ngc <- 2
          while (any(round(summary(gam6)$edf, 1) > 1)) {
            gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp7 != 0) {
          ngc <- 2
          while (any(round(summary(gam7)$edf, 1) > 1)) {
            gam7 <- gam(formula.eq7, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp8 != 0) {
          ngc <- 2
          while (any(round(summary(gam8)$edf, 1) > 1)) {
            gam8 <- gam(formula.eq8, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        X7 <- model.matrix(gam7)
        X7.d2 <- dim(X7)[2]
        X8 <- model.matrix(gam8)
        X8.d2 <- dim(X8)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        if (l.sp6 != 0) 
          sp6 <- gam6$sp
        environment(gam6$formula) <- environment(gam2$formula)
        gp6 <- gam6$nsdf
        if (l.sp7 != 0) 
          sp7 <- gam7$sp
        environment(gam7$formula) <- environment(gam2$formula)
        gp7 <- gam7$nsdf
        if (l.sp8 != 0) 
          sp8 <- gam8$sp
        environment(gam8$formula) <- environment(gam2$formula)
        gp8 <- gam8$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                     gam6$coefficients, gam7$coefficients, gam8$coefficients)
      }
      if (margins[1] %in% c(M$m2) && margins[2] %in% c(M$m3)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        formula.eq6 <- formula[[6]]
        formula.eq7 <- formula[[7]]
        nad2.1 <- "sigma2.1"
        nad2.2 <- "sigma2.2"
        nad.2 <- "nu.2"
        nad.3 <- "dof"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad.2, "~", 
                                        formula.eq5[2], sep = ""))
        formula.eq6 <- as.formula(paste(nad.3, "~", 
                                        formula.eq6[2], sep = ""))
        formula.eq7 <- as.formula(paste(nad3, "~", formula.eq7[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        nu.2 <- rnorm(vo$n, vo$log.nu.2, 0.001)
        dof <- rnorm(vo$n, vo$dof.st, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam7 <- gam(formula.eq7, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        l.sp6 <- length(gam6$sp)
        l.sp7 <- length(gam7$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp6 != 0) {
          ngc <- 2
          while (any(round(summary(gam6)$edf, 1) > 1)) {
            gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp7 != 0) {
          ngc <- 2
          while (any(round(summary(gam7)$edf, 1) > 1)) {
            gam7 <- gam(formula.eq7, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        X7 <- model.matrix(gam7)
        X7.d2 <- dim(X7)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        if (l.sp6 != 0) 
          sp6 <- gam6$sp
        environment(gam6$formula) <- environment(gam2$formula)
        gp6 <- gam6$nsdf
        if (l.sp7 != 0) 
          sp7 <- gam7$sp
        environment(gam7$formula) <- environment(gam2$formula)
        gp7 <- gam7$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                     gam6$coefficients, gam7$coefficients)
      }
      if (margins[1] %in% c(M$m3) && margins[2] %in% c(M$m2)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        formula.eq6 <- formula[[6]]
        formula.eq7 <- formula[[7]]
        nad2.1 <- "sigma2.1"
        nad2.2 <- "sigma2.2"
        nad.1 <- "nu.1"
        nad.3 <- "dof"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad.1, "~", 
                                        formula.eq5[2], sep = ""))
        formula.eq6 <- as.formula(paste(nad.3, "~", 
                                        formula.eq6[2], sep = ""))
        formula.eq7 <- as.formula(paste(nad3, "~", formula.eq7[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        nu.1 <- rnorm(vo$n, vo$log.nu.1, 0.001)
        dof <- rnorm(vo$n, vo$dof.st, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam7 <- gam(formula.eq7, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        l.sp6 <- length(gam6$sp)
        l.sp7 <- length(gam7$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp6 != 0) {
          ngc <- 2
          while (any(round(summary(gam6)$edf, 1) > 1)) {
            gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp7 != 0) {
          ngc <- 2
          while (any(round(summary(gam7)$edf, 1) > 1)) {
            gam7 <- gam(formula.eq7, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        X7 <- model.matrix(gam7)
        X7.d2 <- dim(X7)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        if (l.sp6 != 0) 
          sp6 <- gam6$sp
        environment(gam6$formula) <- environment(gam2$formula)
        gp6 <- gam6$nsdf
        if (l.sp7 != 0) 
          sp7 <- gam7$sp
        environment(gam7$formula) <- environment(gam2$formula)
        gp7 <- gam7$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                     gam6$coefficients, gam7$coefficients)
      }
    }
    else {
      if (margins[1] %in% c(M$m1d, M$bl) && margins[2] %in% 
          c(M$m1d, M$bl)) {
        formula.eq3 <- formula[[3]]
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad3, "~", formula.eq3[2], 
                                        sep = ""))
        set.seed(1)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients)
      }
      if (margins[1] %in% c(M$m1d) && margins[2] %in% 
          c(M$m2, M$m2d)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        nad2.2 <- "sigma2.2"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad3, "~", formula.eq4[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients)
      }
      if (margins[1] %in% c(M$m1d) && margins[2] %in% 
          c(M$m3)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        nad2.2 <- "sigma2.2"
        nad.2 <- "nu.2"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad3, "~", formula.eq5[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        nu.2 <- rnorm(vo$n, vo$log.nu.2, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients)
      }
      if (margins[1] %in% c(M$m2, M$m2d) && margins[2] %in% 
          c(M$m2, M$m2d)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        
        # new stuff to incorporate zero-inflation
        formula.eq6 <- formula[[6]]
        formula.eq7 <- formula[[7]]
        
        formula.eq6 <- as.formula(paste("y1 ~ ", formula.eq6[2], sep=""))
        formula.eq7 <- as.formula(paste("y2 ~ ", formula.eq7[2], sep=""))
        
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp6 <- length(gam6$sp)
        
        gam7 <- gam(formula.eq7, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp7 <- length(gam7$sp)
        
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        
        X7 <- model.matrix(gam7)
        X7.d2 <- dim(X7)[2]
        
        if (l.sp6 != 0) 
          sp6 <- gam6$sp
        environment(gam6$formula) <- environment(gam1$formula)
        gp6 <- gam6$nsdf
        
        if (l.sp7 != 0) 
          sp7 <- gam7$sp
        environment(gam7$formula) <- environment(gam2$formula)
        gp7 <- gam7$nsdf
        # end of new stuff
        
        
        nad2.1 <- "sigma2.1"
        nad2.2 <- "sigma2.2"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad3, "~", formula.eq5[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients,
                     # new stuff
                     gam6$coefficients, gam7$coefficients)
                     # end of new stuff
      }
      if (margins[1] %in% c(M$m2, M$m2d) && margins[2] %in% 
          c(M$bl)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        nad2.1 <- "sigma2.1"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad3, "~", formula.eq4[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients)
      }
      if (margins[1] %in% c(M$m3, M$m3d) && margins[2] %in% 
          c(M$bl)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        nad2.1 <- "sigma2.1"
        nad.1 <- "nu.1"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad.1, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad3, "~", formula.eq5[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        nu.1 <- rnorm(vo$n, vo$log.nu.1, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients)
      }
      if (margins[1] %in% c(M$m3, M$m3d) && margins[2] %in% 
          c(M$m3, M$m3d)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        formula.eq6 <- formula[[6]]
        formula.eq7 <- formula[[7]]
        nad2.1 <- "sigma2.1"
        nad2.2 <- "sigma2.2"
        nad.1 <- "nu.1"
        nad.2 <- "nu.2"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad.1, "~", 
                                        formula.eq5[2], sep = ""))
        formula.eq6 <- as.formula(paste(nad.2, "~", 
                                        formula.eq6[2], sep = ""))
        formula.eq7 <- as.formula(paste(nad3, "~", formula.eq7[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        nu.1 <- rnorm(vo$n, vo$log.nu.1, 0.001)
        nu.2 <- rnorm(vo$n, vo$log.nu.2, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam7 <- gam(formula.eq7, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        l.sp6 <- length(gam6$sp)
        l.sp7 <- length(gam7$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp6 != 0) {
          ngc <- 2
          while (any(round(summary(gam6)$edf, 1) > 1)) {
            gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp7 != 0) {
          ngc <- 2
          while (any(round(summary(gam7)$edf, 1) > 1)) {
            gam7 <- gam(formula.eq7, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        X7 <- model.matrix(gam7)
        X7.d2 <- dim(X7)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        if (l.sp6 != 0) 
          sp6 <- gam6$sp
        environment(gam6$formula) <- environment(gam2$formula)
        gp6 <- gam6$nsdf
        if (l.sp7 != 0) 
          sp7 <- gam7$sp
        environment(gam7$formula) <- environment(gam2$formula)
        gp7 <- gam7$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                     gam6$coefficients, gam7$coefficients)
      }
      if (margins[1] %in% c(M$m2, M$m2d) && margins[2] %in% 
          c(M$m3, M$m3d)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        formula.eq6 <- formula[[6]]
        nad2.1 <- "sigma2.1"
        nad2.2 <- "sigma2.2"
        nad.2 <- "nu.2"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad.2, "~", 
                                        formula.eq5[2], sep = ""))
        formula.eq6 <- as.formula(paste(nad3, "~", formula.eq6[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        nu.2 <- rnorm(vo$n, vo$log.nu.2, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        l.sp6 <- length(gam6$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp6 != 0) {
          ngc <- 2
          while (any(round(summary(gam6)$edf, 1) > 1)) {
            gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        if (l.sp6 != 0) 
          sp6 <- gam6$sp
        environment(gam6$formula) <- environment(gam2$formula)
        gp6 <- gam6$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                     gam6$coefficients)
      }
      if (margins[1] %in% c(M$m3, M$m3d) && margins[2] %in% 
          c(M$m2, M$m2d)) {
        formula.eq3 <- formula[[3]]
        formula.eq4 <- formula[[4]]
        formula.eq5 <- formula[[5]]
        formula.eq6 <- formula[[6]]
        nad2.1 <- "sigma2.1"
        nad2.2 <- "sigma2.2"
        nad.1 <- "nu.1"
        nad3 <- "theta"
        formula.eq3 <- as.formula(paste(nad2.1, "~", 
                                        formula.eq3[2], sep = ""))
        formula.eq4 <- as.formula(paste(nad2.2, "~", 
                                        formula.eq4[2], sep = ""))
        formula.eq5 <- as.formula(paste(nad.1, "~", 
                                        formula.eq5[2], sep = ""))
        formula.eq6 <- as.formula(paste(nad3, "~", formula.eq6[2], 
                                        sep = ""))
        set.seed(1)
        sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
        sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
        nu.1 <- rnorm(vo$n, vo$log.nu.1, 0.001)
        theta <- rnorm(vo$n, vo$i.rho, 0.001)
        rm(list = ".Random.seed", envir = globalenv())
        gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        gam6 <- gam(formula.eq6, data = data, gamma = ngc, 
                    knots = knots, drop.unused.levels = vo$drop.unused.levels)
        l.sp3 <- length(gam3$sp)
        l.sp4 <- length(gam4$sp)
        l.sp5 <- length(gam5$sp)
        l.sp6 <- length(gam6$sp)
        if (l.sp3 != 0) {
          ngc <- 2
          while (any(round(summary(gam3)$edf, 1) > 1)) {
            gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp4 != 0) {
          ngc <- 2
          while (any(round(summary(gam4)$edf, 1) > 1)) {
            gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp5 != 0) {
          ngc <- 2
          while (any(round(summary(gam5)$edf, 1) > 1)) {
            gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        if (l.sp6 != 0) {
          ngc <- 2
          while (any(round(summary(gam6)$edf, 1) > 1)) {
            gam6 <- gam(formula.eq6, data = data, gamma = ngc + 
                          1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
            ngc <- ngc + 1
            if (ngc > 5) 
              break
          }
        }
        X3 <- model.matrix(gam3)
        X3.d2 <- dim(X3)[2]
        X4 <- model.matrix(gam4)
        X4.d2 <- dim(X4)[2]
        X5 <- model.matrix(gam5)
        X5.d2 <- dim(X5)[2]
        X6 <- model.matrix(gam6)
        X6.d2 <- dim(X6)[2]
        if (l.sp3 != 0) 
          sp3 <- gam3$sp
        environment(gam3$formula) <- environment(gam2$formula)
        gp3 <- gam3$nsdf
        if (l.sp4 != 0) 
          sp4 <- gam4$sp
        environment(gam4$formula) <- environment(gam2$formula)
        gp4 <- gam4$nsdf
        if (l.sp5 != 0) 
          sp5 <- gam5$sp
        environment(gam5$formula) <- environment(gam2$formula)
        gp5 <- gam5$nsdf
        if (l.sp6 != 0) 
          sp6 <- gam6$sp
        environment(gam6$formula) <- environment(gam2$formula)
        gp6 <- gam6$nsdf
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     gam3$coefficients, gam4$coefficients, gam5$coefficients, 
                     gam6$coefficients)
      }
    }
  }
  if (type == "gaml") {
    sp2 <- NULL
    if (margins %in% c(M$m2, M$m2d)) {
      formula.eq2 <- formula[[2]]
      nad2.1 <- "sigma2"
      formula.eq2 <- as.formula(paste(nad2.1, "~", formula.eq2[2], 
                                      sep = ""))
      set.seed(1)
      sigma2 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam2 <- gam(formula.eq2, data = data, gamma = ngc, 
                  knots = knots, drop.unused.levels = vo$drop.unused.levels)
      if (M$sp.method != "perf") {
        gam2ff <- gam(formula.eq2, data = data, gamma = ngc, 
                      knots = knots, drop.unused.levels = vo$drop.unused.levels, 
                      fit = FALSE)
        Sl.sf2 <- Sl.setup(gam2ff)
        rm(gam2ff)
      }
      l.sp2 <- length(gam2$sp)
      if (l.sp2 != 0) {
        ngc <- 2
        while (any(round(summary(gam2)$edf, 1) > 1)) {
          gam2 <- gam(formula.eq2, data = data, gamma = ngc + 
                        1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X2 <- model.matrix(gam2)
      X2.d2 <- dim(X2)[2]
      if (l.sp2 != 0) 
        sp2 <- gam2$sp
      environment(gam2$formula) <- environment(gam1$formula)
      gp2 <- gam2$nsdf
      start.v <- c(gam1$coefficients, gam2$coefficients)
    }
    if (margins %in% c(M$m3, M$m3d)) {
      formula.eq2 <- formula[[2]]
      formula.eq3 <- formula[[3]]
      nad2.1 <- "sigma2"
      nad.1 <- "nu"
      formula.eq2 <- as.formula(paste(nad2.1, "~", formula.eq2[2], 
                                      sep = ""))
      formula.eq3 <- as.formula(paste(nad.1, "~", formula.eq3[2], 
                                      sep = ""))
      set.seed(1)
      sigma2 <- rnorm(vo$n, vo$log.sig2.1, 0.001)
      nu <- rnorm(vo$n, vo$log.nu.1, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam2 <- gam(formula.eq2, data = data, gamma = ngc, 
                  knots = knots, drop.unused.levels = vo$drop.unused.levels)
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                  knots = knots, drop.unused.levels = vo$drop.unused.levels)
      if (M$sp.method != "perf") {
        gam2ff <- gam(formula.eq2, data = data, gamma = ngc, 
                      knots = knots, drop.unused.levels = vo$drop.unused.levels, 
                      fit = FALSE)
        gam3ff <- gam(formula.eq3, data = data, gamma = ngc, 
                      knots = knots, drop.unused.levels = vo$drop.unused.levels, 
                      fit = FALSE)
        Sl.sf2 <- Sl.setup(gam2ff)
        Sl.sf3 <- Sl.setup(gam3ff)
        rm(gam2ff, gam3ff)
      }
      l.sp2 <- length(gam2$sp)
      l.sp3 <- length(gam3$sp)
      if (l.sp2 != 0) {
        ngc <- 2
        while (any(round(summary(gam2)$edf, 1) > 1)) {
          gam2 <- gam(formula.eq2, data = data, gamma = ngc + 
                        1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp3 != 0) {
        ngc <- 2
        while (any(round(summary(gam3)$edf, 1) > 1)) {
          gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                        1, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X2 <- model.matrix(gam2)
      X2.d2 <- dim(X2)[2]
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      if (l.sp2 != 0) 
        sp2 <- gam2$sp
      environment(gam2$formula) <- environment(gam1$formula)
      gp2 <- gam2$nsdf
      if (l.sp3 != 0) 
        sp3 <- gam3$sp
      environment(gam3$formula) <- environment(gam1$formula)
      gp3 <- gam3$nsdf
      start.v <- c(gam1$coefficients, gam2$coefficients, 
                   gam3$coefficients)
    }
  }
  if (type == "copSS") {
    if (M$l.flist == 3) {
      formula.eq3 <- formula[[3]]
      nad1 <- "theta"
      formula.eq3 <- as.formula(paste(nad1, "~", formula.eq3[2], 
                                      sep = ""))
      set.seed(1)
      theta <- rnorm(vo$n, vo$i.rho, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                  subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      X3s <- try(predict.gam(gam3, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X3s) == "try-error")) 
        stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.")
      l.sp3 <- length(gam3$sp)
      if (l.sp3 != 0) {
        ngc <- 2
        while (any(round(summary(gam3)$edf, 1) > 1)) {
          gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                        1, subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      if (l.sp3 != 0) 
        sp3 <- gam3$sp
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf
      start.v <- c(gam1$coefficients, c.gam2, gam3$coefficients)
    }
    if (M$l.flist == 4) {
      formula.eq3 <- formula[[3]]
      formula.eq4 <- formula[[4]]
      nad1 <- "sigma2"
      nad2 <- "theta"
      formula.eq3 <- as.formula(paste(nad1, "~", formula.eq3[2], 
                                      sep = ""))
      formula.eq4 <- as.formula(paste(nad2, "~", formula.eq4[2], 
                                      sep = ""))
      set.seed(1)
      sigma2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
      theta <- rnorm(vo$n, vo$i.rho, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                  subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                  subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      X3s <- try(predict.gam(gam3, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X3s) == "try-error")) 
        stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.")
      X4s <- try(predict.gam(gam4, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X4s) == "try-error")) 
        stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.")
      l.sp3 <- length(gam3$sp)
      l.sp4 <- length(gam4$sp)
      if (l.sp3 != 0) {
        ngc <- 2
        while (any(round(summary(gam3)$edf, 1) > 1)) {
          gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                        1, subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp4 != 0) {
        ngc <- 2
        while (any(round(summary(gam4)$edf, 1) > 1)) {
          gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                        1, subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]
      if (l.sp3 != 0) 
        sp3 <- gam3$sp
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf
      if (l.sp4 != 0) 
        sp4 <- gam4$sp
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf
      start.v <- c(gam1$coefficients, c.gam2, gam3$coefficients, 
                   gam4$coefficients)
    }
    if (M$l.flist == 5) {
      formula.eq3 <- formula[[3]]
      formula.eq4 <- formula[[4]]
      formula.eq5 <- formula[[5]]
      nad1 <- "sigma2"
      nad2 <- "nu"
      nad3 <- "theta"
      formula.eq3 <- as.formula(paste(nad1, "~", formula.eq3[2], 
                                      sep = ""))
      formula.eq4 <- as.formula(paste(nad2, "~", formula.eq4[2], 
                                      sep = ""))
      formula.eq5 <- as.formula(paste(nad3, "~", formula.eq5[2], 
                                      sep = ""))
      set.seed(1)
      sigma2 <- rnorm(vo$n, vo$log.sig2.2, 0.001)
      nu <- rnorm(vo$n, vo$log.nu.2, 0.001)
      theta <- rnorm(vo$n, vo$i.rho, 0.001)
      rm(list = ".Random.seed", envir = globalenv())
      gam3 <- gam(formula.eq3, data = data, gamma = ngc, 
                  subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      gam4 <- gam(formula.eq4, data = data, gamma = ngc, 
                  subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      gam5 <- gam(formula.eq5, data = data, gamma = ngc, 
                  subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
      X3s <- try(predict.gam(gam3, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X3s) == "try-error")) 
        stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.")
      X4s <- try(predict.gam(gam4, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X4s) == "try-error")) 
        stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.")
      X5s <- try(predict.gam(gam5, newdata = data[, -dim(data)[2]], 
                             type = "lpmatrix"), silent = TRUE)
      if (any(class(X5s) == "try-error")) 
        stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.")
      l.sp3 <- length(gam3$sp)
      l.sp4 <- length(gam4$sp)
      l.sp5 <- length(gam5$sp)
      if (l.sp3 != 0) {
        ngc <- 2
        while (any(round(summary(gam3)$edf, 1) > 1)) {
          gam3 <- gam(formula.eq3, data = data, gamma = ngc + 
                        1, subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp4 != 0) {
        ngc <- 2
        while (any(round(summary(gam4)$edf, 1) > 1)) {
          gam4 <- gam(formula.eq4, data = data, gamma = ngc + 
                        1, subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      if (l.sp5 != 0) {
        ngc <- 2
        while (any(round(summary(gam5)$edf, 1) > 1)) {
          gam5 <- gam(formula.eq5, data = data, gamma = ngc + 
                        1, subset = inde, knots = knots, drop.unused.levels = vo$drop.unused.levels)
          ngc <- ngc + 1
          if (ngc > 5) 
            break
        }
      }
      X3 <- model.matrix(gam3)
      X3.d2 <- dim(X3)[2]
      X4 <- model.matrix(gam4)
      X4.d2 <- dim(X4)[2]
      X5 <- model.matrix(gam5)
      X5.d2 <- dim(X5)[2]
      if (l.sp3 != 0) 
        sp3 <- gam3$sp
      environment(gam3$formula) <- environment(gam2$formula)
      gp3 <- gam3$nsdf
      if (l.sp4 != 0) 
        sp4 <- gam4$sp
      environment(gam4$formula) <- environment(gam2$formula)
      gp4 <- gam4$nsdf
      if (l.sp5 != 0) 
        sp5 <- gam5$sp
      environment(gam5$formula) <- environment(gam2$formula)
      gp5 <- gam5$nsdf
      start.v <- c(gam1$coefficients, c.gam2, gam3$coefficients, 
                   gam4$coefficients, gam5$coefficients)
    }
  }
  L <- list(start.v = start.v, X3 = X3, X4 = X4, X5 = X5, 
            X6 = X6, X7 = X7, X8 = X8, X9 = X9, X3.d2 = X3.d2, X4.d2 = X4.d2, 
            X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2, 
            X9.d2 = X9.d2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, 
            gp7 = gp7, gp8 = gp8, gp9 = gp9, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8, 
            gam9 = gam9, l.sp3 = l.sp3, l.sp4 = l.sp4, l.sp5 = l.sp5, 
            l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9, 
            sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, 
            sp8 = sp8, sp9 = sp9, X3s = X3s, X4s = X4s, X5s = X5s, 
            X6s = X6s, X7s = X7s, X8s = X8s, X9s = X9s, Sl.sf2 = Sl.sf2, 
            Sl.sf3 = Sl.sf3)
  if (type == "gaml") {
    L$X2 <- X2
    L$X2.d2 <- X2.d2
    L$gp2 <- gp2
    L$gam2 <- gam2
    L$l.sp2 <- l.sp2
    L$sp2 <- sp2
  }
  L
}













reflect_triangle <- function(m, from=c("lower", "upper")) {
  ix <- switch(match.arg(from), lower=upper.tri, upper=lower.tri)(m, diag=FALSE)
  m[ix] <- t(m)[ix]
  m
}











FAST.CoExpress.nb <- function(formula, data, copula) {
  
  if (copula == "Gaussian"){
    copula <- "N"
  } else if (copula == "Clayton"){
    copula <- "C0"
  } else if (copula == "Frank"){
    copula <- "F"
  } else if (copula == "Joe"){
    copula <- "J0"
  } else if (copula == "Gumbel"){
    copula <- "G0"
  } else {
    stop("copula needs to be one of \"Gaussian\", \"Clayton\", \"Frank\", \"Joe\", \"Gumbel\"")
  }
  

  weights <- NULL
  subset <- NULL
  copula2 <- "N"
  margins<-c("NBI","NBI")
  model <- "B"
  dof <- 3
  dof2 <- 3
  ordinal <- FALSE
  surv <- FALSE
  cens1 <- NULL
  cens2 <- NULL
  cens3 <- NULL
  dep.cens <- FALSE
  upperBt1 <- NULL
  upperBt2 <- NULL
  fp <- FALSE
  infl.fac <- 1
  gamlssfit <- FALSE
  rinit <- 1
  rmax <- 1000
  iterlimsp <- 50
  tolsp <- 1e-07
  gc.l <- FALSE
  extra.regI <- "t"
  k1.tvc <- 0
  k2.tvc <- 0 
  knots <- NULL
  penCor <- "unpen"
  parscale <- 1
  sp.penCor <- 3
  Chol <- FALSE
  gamma <- 1
  w.alasso <- NULL
  drop.unused.levels <- TRUE
  min.dn <- 1e-40
  min.pr <- 1e-16
  max.pr <- 0.999999

  BivD <- copula
  BivD2 <- copula2
  Model <- model
  if (missing(margins)) {
    stop("You must choose the margins' values.")
  }
  if (missing(Model)) {
    stop("You must choose a model type.")
  }
  if (margins[1] == "PH" && surv == TRUE) {
    margins[1] <- "cloglog"
  }
  if (margins[1] == "PO" && surv == TRUE) {
    margins[1] <- "logit"
  }
  if (margins[2] == "PH" && surv == TRUE) {
    margins[2] <- "cloglog"
  }
  if (margins[2] == "PO" && surv == TRUE) {
    margins[2] <- "logit"
  }
  bl <- c("probit", "logit", "cloglog")
  v.rB1 <- upperBt1
  v.rB2 <- upperBt2
  if (Model == "ROY") {
    L <- eval(substitute(SemiParROY(formula, data, weights,
                                    subset,
                                    BivD1 = BivD, BivD2, margins, dof1 = dof,
                                    dof2, gamlssfit, fp, infl.fac, rinit, rmax, iterlimsp,
                                    tolsp, gc.l, parscale, extra.regI, knots = knots,
                                    drop.unused.levels = drop.unused.levels, min.dn = min.dn,
                                    min.pr = min.pr, max.pr = max.pr
    ), list(weights = weights)))
  } else {
    if (surv == FALSE && ordinal == FALSE) {
      if ((margins[1] %in% bl && margins[2] %in% bl &&
           is.na(margins[3])) || (margins[1] %in% bl &&
                                  !(margins[2] %in% bl) && Model == "B" && is.na(margins[3]))) {
        L <- eval(substitute(SemiParBIV(formula, data,
                                        weights, subset, Model, BivD, margins, dof,
                                        gamlssfit, fp,
                                        hess = TRUE, infl.fac, rinit,
                                        rmax, iterlimsp, tolsp, gc.l, parscale, extra.regI,
                                        intf = TRUE, theta.fx = NULL, knots = knots,
                                        drop.unused.levels = drop.unused.levels, min.dn = min.dn,
                                        min.pr = min.pr, max.pr = max.pr
        ), list(weights = weights)))
      }
    }
    if (surv == FALSE && ordinal == TRUE) {
      if ((margins[1] %in% bl && margins[2] %in% bl &&
           is.na(margins[3])) || (margins[1] %in% bl &&
                                  !(margins[2] %in% bl) && is.na(margins[3]))) {
        L <- eval(substitute(CopulaCLM(formula, data,
                                       weights, subset, Model, BivD, margins, dof,
                                       gamlssfit, fp,
                                       hess = TRUE, infl.fac, rinit,
                                       rmax, iterlimsp, tolsp, gc.l, parscale, extra.regI,
                                       intf = TRUE, theta.fx = NULL, knots = knots,
                                       drop.unused.levels = drop.unused.levels, min.dn = min.dn,
                                       min.pr = min.pr, max.pr = max.pr
        ), list(weights = weights)))
      }
    }
    if (margins[1] %in% bl && !(margins[2] %in% bl) && surv ==
        FALSE && is.na(margins[3]) && Model == "BSS" &&
        ordinal == FALSE) {
      L <- eval(substitute(
        copulaSampleSel(formula, data,
                        weights, subset, BivD, margins, dof, fp, infl.fac,
                        rinit, rmax, iterlimsp, tolsp, gc.l, parscale,
                        extra.regI, knots,
                        drop.unused.levels = drop.unused.levels,
                        min.dn = min.dn, min.pr = min.pr, max.pr = max.pr
        ),
        list(weights = weights)
      ))
    }
    if (!is.na(margins[3])) {
      if (margins[1] %in% bl && margins[2] %in% bl &&
          margins[3] %in% bl && surv == FALSE && ordinal ==
          FALSE) {
        L <- eval(substitute(
          SemiParTRIV(formula, data,
                      weights, subset, Model, margins, penCor, sp.penCor,
                      approx = FALSE, Chol, infl.fac, gamma, w.alasso,
                      rinit, rmax, iterlimsp, tolsp, gc.l, parscale,
                      extra.regI, knots, drop.unused.levels = drop.unused.levels,
                      min.dn = min.dn, min.pr = min.pr, max.pr = max.pr
          ),
          list(weights = weights)
        ))
      }
    }
    if ((!(margins[1] %in% bl) || surv == TRUE) && ordinal ==
        FALSE) {
      robust <- FALSE
      t.c <- 3
      sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss1 <- gamlss2 <- gam1 <- gam2 <- y1m <- y2m <- indexTeq1 <- indexTeq2 <- NULL
      i.rho <- log.sig2.2 <- log.nu.2 <- log.nu.1 <- log.sig2.1 <- dof.st <- NULL
      end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- 0
      sp1 <- sp2 <- NULL
      sp3 <- gp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- gam9 <- NULL
      sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- sp9 <- NULL
      c11 <- c10 <- c01 <- c00 <- NA
      cens1Mix <- cens2Mix <- NULL
      Sl.sf <- NULL
      sp.method <- "perf"
      Xd1 <- Xd2 <- mono.sm.pos1 <- mono.sm.pos2 <- mono.sm.pos <- NULL
      surv.flex <- FALSE
      Deq1 <- pos.pbeq1 <- Deq2 <- pos.pbeq2 <- list()
      BivD2 <- c(
        "C0C90", "C0C270", "C180C90", "C180C270",
        "J0J90", "J0J270", "J180J90", "J180J270", "G0G90",
        "G0G270", "G180G90", "G180G270", "GAL0GAL90",
        "GAL0GAL270", "GAL180GAL90", "GAL180GAL270"
      )
      opc <- c(
        "N", "C0", "C90", "C180", "C270", "J0",
        "J90", "J180", "J270", "G0", "G90", "G180",
        "G270", "F", "AMH", "FGM", "T", "PL", "HO",
        "GAL0", "GAL90", "GAL180", "GAL270"
      )
      scc <- c(
        "C0", "C180", "GAL0", "GAL180", "J0", "J180",
        "G0", "G180", BivD2
      )
      sccn <- c(
        "C90", "C270", "GAL90", "GAL270", "J90",
        "J270", "G90", "G270"
      )
      m2 <- c(
        "N", "GU", "rGU", "LO", "LN", "WEI", "iG",
        "GA", "BE", "FISK", "GP", "GPII", "GPo"
      )
      m3 <- c("DAGUM", "SM", "TW")
      m1d <- c("PO", "ZTP", "DGP0")
      m2d <- c("NBI", "NBII", "PIG", "DGP", "DGPII")
      m3d <- c("DEL", "SICHEL")
      ct <- data.frame(c(opc), c(
        1:14, 55, 56, 57, 60,
        61, 62:65
      ))
      cta <- data.frame(c(opc), c(
        1, 3, 23, 13, 33, 6,
        26, 16, 36, 4, 24, 14, 34, 5, 55, 56, 2, 60,
        61, 62:65
      ))
      if (BivD %in% BivD2) {
        if (BivD %in% BivD2[1:4]) {
          BivDt <- "C0"
        }
        if (BivD %in% BivD2[5:12]) {
          BivDt <- "J0"
        }
        if (BivD %in% BivD2[13:16]) {
          BivDt <- "C0"
        }
        nC <- ct[which(ct[, 1] == BivDt), 2]
        nCa <- cta[which(cta[, 1] == BivDt), 2]
      }
      if (!(BivD %in% BivD2)) {
        nC <- ct[which(ct[, 1] == BivD), 2]
        nCa <- cta[which(cta[, 1] == BivD), 2]
      }
      if (!is.list(formula)) {
        stop("You must specify a list of equations.")
      }
      l.flist <- length(formula)
      form.check(formula, l.flist)
      cl <- match.call()
      mf <- match.call(expand.dots = FALSE)
      pred.varR <- pred.var(formula, l.flist)
      v1 <- pred.varR$v1
      v2 <- pred.varR$v2
      pred.n <- pred.varR$pred.n
      if (!is.null(v.rB1)) {
        pred.n <- c(pred.n, v.rB1)
      }
      if (!is.null(v.rB2)) {
        pred.n <- c(pred.n, v.rB2)
      }
      fake.formula <- paste(v1[1], "~", paste(pred.n,
                                              collapse = " + "
      ))
      environment(fake.formula) <- environment(formula[[1]])
      mf$formula <- fake.formula
      mf$upperBt1 <- mf$upperBt2 <- mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$dep.cens <- mf$ordinal <- mf$Model <- mf$model <- mf$knots <- mf$k1.tvc <- mf$k2.tvc <- mf$surv <- mf$BivD <- mf$copula <- mf$copula2 <- mf$margins <- mf$fp <- mf$dof <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL
      mf$drop.unused.levels <- drop.unused.levels
      if (surv == TRUE) {
        mf$na.action <- na.pass
      }
      mf[[1]] <- as.name("model.frame")
      data <- eval(mf, parent.frame())
      if (surv == TRUE) {
        if (!("(cens1)" %in% names(data)) && margins[1] %in%
            bl) {
          stop("You must provide both censoring indicators.")
        }
        if (!("(cens2)" %in% names(data)) && margins[2] %in%
            bl) {
          stop("You must provide both censoring indicators.")
        }
      }
      if (gc.l == TRUE) {
        gc()
      }
      if (!("(weights)" %in% names(data))) {
        weights <- rep(1, dim(data)[1])
        data$weights <- weights
        names(data)[length(names(data))] <- "(weights)"
      } else {
        weights <- data[, "(weights)"]
      }
      if (!("(cens1)" %in% names(data))) {
        cens1 <- rep(0, dim(data)[1])
        data$cens1 <- cens1
        names(data)[length(names(data))] <- "(cens1)"
      } else {
        cens1 <- data[, "(cens1)"]
      }
      if (!("(cens2)" %in% names(data))) {
        cens2 <- rep(0, dim(data)[1])
        data$cens2 <- cens2
        names(data)[length(names(data))] <- "(cens2)"
      } else {
        cens2 <- data[, "(cens2)"]
      }
      if (!("(cens3)" %in% names(data))) {
        cens3 <- rep(0, dim(data)[1])
        data$cens3 <- cens3
        names(data)[length(names(data))] <- "(cens3)"
      } else {
        cens3 <- data[, "(cens3)"]
      }
      if (surv == TRUE) {
        if (is.factor(cens1) && !is.factor(cens2)) {
          stop("Both censoring indicators have to be factor variables for mixed censoring case.")
        }
        if (!is.factor(cens1) && is.factor(cens2)) {
          stop("Both censoring indicators have to be factor variables for mixed censoring case.")
        }
      }
      if (surv == TRUE && is.factor(cens1) && !is.null(v.rB1)) {
        data[!(cens1 == "I"), v.rB1] <- data[!(cens1 ==
                                                 "I"), v1[1]]
      }
      if (surv == TRUE && is.factor(cens2) && !is.null(v.rB2)) {
        data[!(cens2 == "I"), v.rB2] <- data[!(cens2 ==
                                                 "I"), v2[1]]
      }
      if (surv == TRUE) {
        if (any(is.na(data[, v1[1]]) | is.na(data[
          ,
          v2[1]
        ]))) {
          stop("Time to event with NA's. Please check your time covariates.")
        }
        actual.NAs <- as.numeric(which(apply(apply(
          data,
          1, is.na
        ), 2, any)))
        data <- na.omit(data)
        if (length(actual.NAs) > 0) {
          cens1 <- cens1[-actual.NAs]
          cens2 <- cens2[-actual.NAs]
        }
      }
      n <- dim(data)[1]
      if (surv == TRUE && is.factor(cens1) && is.factor(cens2)) {
        cens1Mix <- cens1
        cens2Mix <- cens2
        cens1 <- cens2 <- rep(1, n)
      }
      M <- list(
        m1d = m1d, m2 = m2, m2d = m2d, m3 = m3,
        m3d = m3d, BivD = BivD, bl = bl, robust = robust,
        opc = opc, extra.regI = extra.regI, margins = margins,
        BivD2 = BivD2, dof = dof, surv = surv, c1 = cens1,
        c2 = cens2, c3 = cens3, dep.cens = dep.cens
      )
      M$K1 <- NULL
      M$type.cens1 <- M$type.cens2 <- "R"
      # pream.wm(formula, margins, M, l.flist)
      formula.eq1 <- formula[[1]]
      formula.eq2 <- formula[[2]]
      form.eq12R <- form.eq12(
        formula.eq1, data, v1, margins[1],
        m1d, m2d
      )
      formula.eq1 <- form.eq12R$formula.eq1
      formula.eq1r <- form.eq12R$formula.eq1r
      y1 <- form.eq12R$y1
      y1.test <- form.eq12R$y1.test
      y1m <- form.eq12R$y1m
      if (surv == FALSE) {
        gam1 <- eval(substitute(
          gam(formula.eq1,
              gamma = infl.fac,
              weights = weights, data = data, knots = knots,
              drop.unused.levels = drop.unused.levels
          ),
          list(weights = weights)
        ))
      }
      if (surv == TRUE && margins[1] %in% c(m2, m3) &&
          margins[2] %in% bl) {
        gam1 <- eval(substitute(
          gam(formula.eq1,
              gamma = infl.fac,
              weights = weights, data = data, knots = knots,
              drop.unused.levels = drop.unused.levels
          ),
          list(weights = weights)
        ))
      } else {
        if (surv == TRUE && !(margins[1] %in% bl)) {
          gam1 <- eval(substitute(
            gam(formula.eq1,
                gamma = infl.fac,
                weights = weights * cens1, data = data,
                knots = knots, drop.unused.levels = drop.unused.levels
            ),
            list(weights = weights, cens1 = cens1)
          ))
        }
      }
      if (surv == TRUE && margins[1] %in% bl) {
        surv.flex <- TRUE
        f.eq1 <- form.eq12R$f.eq1
        data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
        tempb <- eval(substitute(
          gam(f.eq1,
              family = cox.ph(),
              data = data, weights = cens1, drop.unused.levels = drop.unused.levels
          ),
          list(cens1 = cens1)
        ))
        data$Sh <- as.vector(mm(predict(tempb, type = "response"),
                                min.pr = min.pr, max.pr = max.pr
        ))
        cens11 <- ifelse(cens1 == 0, 1e-07, cens1)
        gam1 <- eval(substitute(
          scam(formula.eq1,
               gamma = infl.fac,
               weights = weights * cens11, data = data
          ),
          list(weights = weights, cens11 = cens11)
        ))
        lsgam1 <- length(gam1$smooth)
        if (lsgam1 == 0) {
          stop("You must use at least a monotonic smooth function of time in the first equation.")
        }
        clsm <- ggr <- NA
        for (i in 1:lsgam1) {
          clsm[i] <- class(gam1$smooth[[i]])[1]
        }
        if (sum(as.numeric(clsm[1] %in% c("mpi.smooth"))) ==
            0) {
          stop("You must have a monotonic smooth of time and it has to be the first to be included.")
        }
        l.sp1 <- length(gam1$sp)
        if (l.sp1 != 0) {
          sp1 <- gam1$sp
        }
        sp1[1] <- 1
        gam.call <- gam1$call
        gam.call$sp <- sp1
        gam1 <- eval(gam.call)
        j <- 1
        for (i in 1:lsgam1) {
          if (max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$term))) !=
              0 && clsm[i] == "mpi.smooth") {
            mono.sm.pos1 <- c(mono.sm.pos1, c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para))
          }
        }
        X1 <- predict(gam1, type = "lpmatrix")
        if (!is.null(indexTeq1) && k1.tvc != 0) {
          if (range(X1[, indexTeq1])[1] < 0) {
            stop("Check design matrix for smooth(s) of tvc term(s) in eq. 1.")
          }
        }
        Xd1 <- Xdpred(gam1, data, v1[1])
        gam1$y <- data[, v1[1]]
        st.v1 <- c(gam1$coefficients)
        if (!is.null(indexTeq1)) {
          st.v1[mono.sm.pos1] <- exp(st.v1[mono.sm.pos1])
          while (range(Xd1 %*% st.v1)[1] < 0) {
            st.v1[indexTeq1] <- 0.999 *
              st.v1[indexTeq1]
          }
          gam1$coefficients <- gam1$coefficients.t <- st.v1
          gam1$coefficients.t[mono.sm.pos1] <- exp(gam1$coefficients.t[mono.sm.pos1])
        }
      }
      gam1$formula <- formula.eq1r
      lsgam1 <- length(gam1$smooth)
      y1 <- y1.test
      if (margins[1] %in% c("LN")) {
        y1 <- log(y1)
      }
      attr(data, "terms") <- NULL
      if (!(surv == TRUE && margins[1] %in% bl)) {
        names(gam1$model)[1] <- as.character(formula.eq1r[2])
        X1 <- predict(gam1, type = "lpmatrix")
        l.sp1 <- length(gam1$sp)
        sp1 <- gam1$sp
      }
      gp1 <- gam1$nsdf
      X1.d2 <- dim(X1)[2]
      form.eq12R <- form.eq12(
        formula.eq2, data, v2, margins[2],
        m1d, m2d
      )
      formula.eq2 <- form.eq12R$formula.eq1
      formula.eq2r <- form.eq12R$formula.eq1r
      y2 <- form.eq12R$y1
      y2.test <- form.eq12R$y1.test
      y2m <- form.eq12R$y1m
      if (surv == FALSE) {
        gam2 <- eval(substitute(
          gam(formula.eq2,
              gamma = infl.fac,
              weights = weights, data = data, knots = knots,
              drop.unused.levels = drop.unused.levels
          ),
          list(weights = weights)
        ))
      }
      if (surv == TRUE && !(margins[2] %in% bl)) {
        gam2 <- eval(substitute(
          gam(formula.eq2,
              gamma = infl.fac,
              weights = weights * cens2, data = data, knots = knots,
              drop.unused.levels = drop.unused.levels
          ),
          list(weights = weights, cens2 = cens2)
        ))
      }
      if (surv == TRUE && margins[2] %in% bl) {
        surv.flex <- TRUE
        f.eq2 <- form.eq12R$f.eq1
        data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
        tempb <- eval(substitute(
          gam(f.eq2,
              family = cox.ph(),
              data = data, weights = cens2, drop.unused.levels = drop.unused.levels
          ),
          list(cens2 = cens2)
        ))
        data$Sh <- as.vector(mm(predict(tempb, type = "response"),
                                min.pr = min.pr, max.pr = max.pr
        ))
        cens22 <- ifelse(cens2 == 0, 1e-07, cens2)
        gam2 <- eval(substitute(
          scam(formula.eq2,
               gamma = infl.fac,
               weights = weights * cens22, data = data
          ),
          list(weights = weights, cens22 = cens22)
        ))
        lsgam2 <- length(gam2$smooth)
        if (lsgam2 == 0) {
          stop("You must use at least a monotonic smooth function of time in the second equation.")
        }
        clsm <- ggr <- NA
        for (i in 1:lsgam2) {
          clsm[i] <- class(gam2$smooth[[i]])[1]
        }
        if (sum(as.numeric(clsm[1] %in% c("mpi.smooth"))) ==
            0) {
          stop("You must have a monotonic smooth of time and it has to be the first to be included.")
        }
        l.sp2 <- length(gam2$sp)
        if (l.sp2 != 0) {
          sp2 <- gam2$sp
        }
        sp2[1] <- 1
        gam.call <- gam2$call
        gam.call$sp <- sp2
        gam2 <- eval(gam.call)
        j <- 1
        for (i in 1:lsgam2) {
          if (max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$term))) !=
              0 && clsm[i] == "mpi.smooth") {
            mono.sm.pos2 <- c(mono.sm.pos2, c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para))
          }
        }
        X2 <- predict(gam2, type = "lpmatrix")
        if (!is.null(indexTeq2) && k2.tvc != 0) {
          if (range(X2[, indexTeq2])[1] < 0) {
            stop("Check design matrix for smooth(s) of tvc term(s) in eq. 2.")
          }
        }
        Xd2 <- Xdpred(gam2, data, v2[1])
        gam2$y <- data[, v2[1]]
        st.v2 <- c(gam2$coefficients)
        if (!is.null(indexTeq2)) {
          st.v2[mono.sm.pos2] <- exp(st.v2[mono.sm.pos2])
          while (range(Xd2 %*% st.v2)[1] < 0) {
            st.v2[indexTeq2] <- 0.999 *
              st.v2[indexTeq2]
          }
          gam2$coefficients <- gam2$coefficients.t <- st.v2
          gam2$coefficients.t[mono.sm.pos2] <- exp(gam2$coefficients.t[mono.sm.pos2])
        }
      }
      gam2$formula <- formula.eq2r
      lsgam2 <- length(gam2$smooth)
      y2 <- y2.test
      if (margins[2] %in% c("LN")) {
        y2 <- log(y2)
      }
      attr(data, "terms") <- NULL
      if (!(surv == TRUE && margins[2] %in% bl)) {
        names(gam2$model)[1] <- as.character(formula.eq2r[2])
        X2 <- predict(gam2, type = "lpmatrix")
        l.sp2 <- length(gam2$sp)
        sp2 <- gam2$sp
      }
      gp2 <- gam2$nsdf
      X2.d2 <- dim(X2)[2]
      res1 <- residuals(gam1)
      res2 <- residuals(gam2)
      ass.s <- cor(res1, res2, method = "kendall")
      ass.s <- sign(ass.s) * ifelse(abs(ass.s) > 0.9,
                                    0.9, abs(ass.s)
      )
      i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)
      dof.st <- log(dof - 2)
      names(dof.st) <- "dof.star"
      if (!(margins[1] %in% c(m1d, bl))) {
        start.snR <- startsn(margins[1], y1)
        log.sig2.1 <- start.snR$log.sig2.1
        names(log.sig2.1) <- "sigma1.star"
        if (margins[1] %in% c(m3)) {
          log.nu.1 <- start.snR$log.nu.1
          names(log.nu.1) <- "nu.1.star"
        }
      }
      if (!(margins[2] %in% c(m1d, bl))) {
        start.snR <- startsn(margins[2], y2)
        log.sig2.2 <- start.snR$log.sig2.1
        names(log.sig2.2) <- "sigma2.star"
        if (margins[2] %in% c(m3)) {
          log.nu.2 <- start.snR$log.nu.1
          names(log.nu.2) <- "nu.2.star"
        }
      }
      vo <- list(
        gam1 = gam1, gam2 = gam2, i.rho = i.rho,
        log.sig2.2 = log.sig2.2, log.nu.2 = log.nu.2,
        log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1,
        dof.st = dof.st, n = n, drop.unused.levels = drop.unused.levels
      )
      start.v <- overall.sv(margins, M, vo)
      if (l.flist > 2) {
        overall.svGR <- overall.svG(formula, data,
                                    ngc = 2,
                                    margins, M, vo, gam1, gam2, knots = knots
        )
        start.v <- overall.svGR$start.v
        X3 <- overall.svGR$X3
        X4 <- overall.svGR$X4
        X5 <- overall.svGR$X5
        X6 <- overall.svGR$X6
        X7 <- overall.svGR$X7
        X8 <- overall.svGR$X8
        X3.d2 <- overall.svGR$X3.d2
        X4.d2 <- overall.svGR$X4.d2
        X5.d2 <- overall.svGR$X5.d2
        X6.d2 <- overall.svGR$X6.d2
        X7.d2 <- overall.svGR$X7.d2
        X8.d2 <- overall.svGR$X8.d2
        gp3 <- overall.svGR$gp3
        gp4 <- overall.svGR$gp4
        gp5 <- overall.svGR$gp5
        gp6 <- overall.svGR$gp6
        gp7 <- overall.svGR$gp7
        gp8 <- overall.svGR$gp8
        gam3 <- overall.svGR$gam3
        gam4 <- overall.svGR$gam4
        gam5 <- overall.svGR$gam5
        gam6 <- overall.svGR$gam6
        gam7 <- overall.svGR$gam7
        gam8 <- overall.svGR$gam8
        l.sp3 <- overall.svGR$l.sp3
        l.sp4 <- overall.svGR$l.sp4
        l.sp5 <- overall.svGR$l.sp5
        l.sp6 <- overall.svGR$l.sp6
        l.sp7 <- overall.svGR$l.sp7
        l.sp8 <- overall.svGR$l.sp8
        sp3 <- overall.svGR$sp3
        sp4 <- overall.svGR$sp4
        sp5 <- overall.svGR$sp5
        sp6 <- overall.svGR$sp6
        sp7 <- overall.svGR$sp7
        sp8 <- overall.svGR$sp8
      }
      GAM <- list(
        gam1 = gam1, gam2 = gam2, gam3 = gam3,
        gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7,
        gam8 = gam8, gam9 = gam9
      )
      aregams <- FALSE
      if ((l.sp1 != 0 || l.sp2 != 0 || l.sp3 != 0 || l.sp4 !=
           0 || l.sp5 != 0 || l.sp6 != 0 || l.sp7 != 0 ||
           l.sp8 != 0) && fp == FALSE) {
        L.GAM <- list(
          l.gam1 = length(gam1$coefficients),
          l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients),
          l.gam4 = length(gam4$coefficients), l.gam5 = length(gam5$coefficients),
          l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients),
          l.gam8 = length(gam8$coefficients), l.gam9 = 0
        )
        L.SP <- list(
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3,
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,
          l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9
        )
        sp <- c(
          sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8,
          sp9
        )
        qu.mag <- S.m(GAM, L.SP, L.GAM)
        aregams <- TRUE
      }
      if (missing(parscale)) {
        parscale <- 1
      }
      respvec <- respvec2 <- respvec3 <- list(
        y1 = y1,
        y2 = y2, y1.y2 = NULL, y1.cy2 = NULL, cy1.y2 = NULL,
        cy1.cy2 = NULL, cy1 = NULL, cy = NULL, univ = 0
      )
      my.env <- new.env()
      my.env$signind <- 1
      lsgam3 <- length(gam3$smooth)
      lsgam4 <- length(gam4$smooth)
      lsgam5 <- length(gam5$smooth)
      lsgam6 <- length(gam6$smooth)
      lsgam7 <- length(gam7$smooth)
      lsgam8 <- length(gam8$smooth)
      lsgam9 <- length(gam9$smooth)
      indUR <- indUL <- indUI <- indUU <- indRR <- indRL <- indRI <- indRU <- indLR <- indLL <- indLI <- indLU <- indIR <- indIL <- indII <- indIU <- rep(
        0,
        n
      )
      if (surv == TRUE && dep.cens == FALSE) {
        if ((surv == TRUE && margins[1] %in% bl && margins[2] %in%
             bl && !is.factor(cens1) && !is.factor(cens2)) ||
            (surv == TRUE && margins[1] %in% m2 && margins[2] %in%
             m2)) {
          c11 <- cens1 * cens2
          c10 <- cens1 * (1 - cens2)
          c01 <- (1 - cens1) * cens2
          c00 <- (1 - cens1) * (1 - cens2)
        }
        if (surv == TRUE && margins[1] %in% c(m2, m3) &&
            margins[2] %in% bl) {
          c11 <- cens2
          c10 <- 1 - cens2
          c01 <- NULL
          c00 <- NULL
        }
        if (!is.null(cens1Mix) && !is.null(cens2Mix)) {
          if (surv == TRUE && margins[1] %in% bl &&
              margins[2] %in% bl && is.factor(cens1Mix) &&
              is.factor(cens2Mix)) {
            gamlssfit <- TRUE
            indUR <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix ==
                                                                "R")
            indUL <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix ==
                                                                "L")
            indUI <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix ==
                                                                "I")
            indUU <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix ==
                                                                "U")
            indRR <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix ==
                                                                "R")
            indRL <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix ==
                                                                "L")
            indRI <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix ==
                                                                "I")
            indRU <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix ==
                                                                "U")
            indLR <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix ==
                                                                "R")
            indLL <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix ==
                                                                "L")
            indLI <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix ==
                                                                "I")
            indLU <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix ==
                                                                "U")
            indIR <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix ==
                                                                "R")
            indIL <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix ==
                                                                "L")
            indII <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix ==
                                                                "I")
            indIU <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix ==
                                                                "U")
          }
        }
      }
      if (surv == TRUE && dep.cens == TRUE) {
        c11 <- NULL
        c10 <- cens1
        c01 <- cens2
        c00 <- cens3
      }
      my.env$k1 <- k1.tvc
      my.env$k2 <- k2.tvc
      VC <- list(
        lsgam1 = lsgam1, indexTeq1 = indexTeq1,
        indexTeq2 = indexTeq2, lsgam2 = lsgam2, Deq1 = Deq1,
        pos.pbeq1 = pos.pbeq1, Deq2 = Deq2, pos.pbeq2 = pos.pbeq2,
        lsgam3 = lsgam3, robust = FALSE, sp.fixed = NULL,
        lsgam4 = lsgam4, Sl.sf = Sl.sf, sp.method = sp.method,
        lsgam5 = lsgam5, K1 = NULL, lsgam6 = lsgam6,
        lsgam7 = lsgam7, lsgam8 = lsgam8, lsgam9 = lsgam9,
        X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5,
        X6 = X6, X7 = X7, X8 = X8, X1.d2 = X1.d2, X2.d2 = X2.d2,
        X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2,
        X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2,
        gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4,
        gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8,
        l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3,
        l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,
        l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = 0, my.env = my.env,
        infl.fac = infl.fac, weights = weights, fp = fp,
        gamlssfit = gamlssfit, hess = NULL, Model = "CC",
        univ.gamls = FALSE, model = model, end = end,
        BivD = BivD, nCa = nCa, copula = copula, copula2 = copula2,
        nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
        parscale = parscale, margins = margins, Cont = "YES",
        ccss = "no", m2 = m2, m3 = m3, m1d = m1d, m2d = m2d,
        m3d = m3d, bl = bl, triv = FALSE, y1m = y1m,
        y2m = y2m, tc = t.c, i.rho = i.rho, dof = dof,
        dof.st = dof.st, BivD2 = BivD2, cta = cta, ct = ct,
        zerov = -10, c11 = c11, c10 = c10, c01 = c01,
        c00 = c00, indUR = indUR, indUL = indUL, indUI = indUI,
        indUU = indUU, indRR = indRR, indRL = indRL,
        indRI = indRI, indRU = indRU, indLR = indLR,
        indLL = indLL, indLI = indLI, indLU = indLU,
        indIR = indIR, indIL = indIL, indII = indII,
        indIU = indIU, surv = surv, Xd1 = Xd1, Xd2 = Xd2,
        mono.sm.pos1 = mono.sm.pos1, mono.sm.pos2 = mono.sm.pos2,
        surv.flex = surv.flex, mono.sm.pos = mono.sm.pos,
        gp2.inf = NULL, informative = "no", zero.tol = 0.01,
        min.dn = min.dn, min.pr = min.pr, max.pr = max.pr
      )
      if (gc.l == TRUE) {
        gc()
      }
      if (gamlssfit == TRUE) {
        type.cens1 <- type.cens2 <- "R"
        surv1 <- surv2 <- surv
        form.gamlR <- form.gaml(formula, l.flist, M)
        if (surv == TRUE && margins[1] %in% c(m2, m3) &&
            margins[2] %in% bl) {
          surv1 <- FALSE
        }
        if (surv == TRUE && margins[1] %in% bl && margins[2] %in%
            bl && is.factor(cens1Mix) && is.factor(cens2Mix)) {
          cens1 <- cens1Mix
          cens2 <- cens2Mix
          type.cens1 <- type.cens2 <- "mixed"
          M$type.cens1 <- type.cens1
          M$type.cens2 <- type.cens2
        }
        gamlss1 <- eval(substitute(
          gamlss(form.gamlR$formula.gamlss1,
                 data = data, weights = weights, subset = subset,
                 margin = margins[1], surv = surv1, cens = cens1,
                 type.cens = type.cens1, upperB = upperBt1,
                 infl.fac = infl.fac, rinit = rinit, rmax = rmax,
                 iterlimsp = iterlimsp, tolsp = tolsp, gc.l = gc.l,
                 parscale = 1, extra.regI = extra.regI, k.tvc = k1.tvc,
                 drop.unused.levels = drop.unused.levels
          ),
          list(weights = weights, cens1 = cens1)
        ))
        gamlss2 <- eval(substitute(
          gamlss(form.gamlR$formula.gamlss2,
                 data = data, weights = weights, subset = subset,
                 margin = margins[2], surv = surv2, cens = cens2,
                 type.cens = type.cens2, upperB = upperBt2,
                 infl.fac = infl.fac, rinit = rinit, rmax = rmax,
                 iterlimsp = iterlimsp, tolsp = tolsp, gc.l = gc.l,
                 parscale = 1, extra.regI = extra.regI, k.tvc = k2.tvc,
                 drop.unused.levels = drop.unused.levels
          ),
          list(weights = weights, cens2 = cens2)
        ))
        SP <- list(
          sp1 = sp1, sp2 = sp2, sp3 = sp3,
          sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7,
          sp8 = sp8
        )
        gamls.upsvR <- gamls.upsv(gamlss1, gamlss2,
                                  margins, M, l.flist,
                                  nstv = names(start.v),
                                  VC, GAM, SP
        )
        sp <- gamls.upsvR$sp
        start.v <- gamls.upsvR$start.v
        VC$X1 <- gamlss1$VC$X1
        VC$Xd1 <- gamlss1$VC$Xd1
        VC$X1.2 <- gamlss1$VC$X2
        VC$X2 <- gamlss2$VC$X1
        VC$Xd2 <- gamlss2$VC$Xd1
        VC$X2.2 <- gamlss2$VC$X2
        rangeSurv1 <- gamlss1$rangeSurv
        rangeSurv2 <- gamlss2$rangeSurv
      }


      func.opt <- bdiscrdiscr_ARB

      if (aregams) {
        qu.mag$kappad <- TRUE
      }

      result <- try({
        # Code that might generate an error
        kappa1start <- unname(pscl::zeroinfl(respvec$y1 ~ VC$X1-1 | VC$X6-1, dist="negbin")$coefficients$zero)
        # kappa1start <- aaaaaaaaaaaaaaaaaaa()
      }, silent = TRUE)
        
      # You can check if there was an error
      if (inherits(result, "try-error")) {
        kappa1start <- rep(0, NCOL(VC$X6))
        kappa1start[1] <- logit(0.5*mean(respvec$y1==0) + 0.001)
      } 

      result <- try({
        # Code that might generate an error
        kappa2start <- unname(pscl::zeroinfl(respvec$y2 ~ VC$X2-1 | VC$X7-1, dist="negbin")$coefficients$zero)
        # kappa2start <- aaaaaaaaaaaaaaaaaaa()
      }, silent = TRUE)
        
      # You can check if there was an error
      if (inherits(result, "try-error")) {
        kappa2start <- rep(0, NCOL(VC$X7))
        kappa2start[1] <- logit(0.5*mean(respvec$y2==0) + 0.001)
      } 

      start.v_use <- start.v
      
      start.v_use[(length(start.v)-(VC$X6.d2 + VC$X7.d2 - 1)):length(start.v)] <- c(kappa1start, kappa2start)

      SemiParFit <- SemiParBIV.fit(
        func.opt = func.opt,
        start.v = start.v_use, rinit = rinit, rmax = rmax,
        iterlim = 10000, iterlimsp = iterlimsp, tolsp = tolsp,
        respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag
      )
      SemiParFit.p <- copulaReg.fit.post(
        SemiParFit = SemiParFit,
        VC = VC, GAM
      )
      y1.m <- y1
      if (margins[1] == "LN") {
        y1.m <- exp(y1)
      }
      y2.m <- y2
      if (margins[2] == "LN") {
        y2.m <- exp(y2)
      }
      SemiParFit <- SemiParFit.p$SemiParFit
      if (gc.l == TRUE) {
        gc()
      }
      
      
      cov.c(SemiParFit)
      
      # get parameter ests
      ests <- SemiParFit$fit$argument

      # get varcov mat
      varcovmat <- SemiParFit.p$Vb
      
      # calculate std errors
      ses <- sqrt(diag(varcovmat))
      
      # calculate p-values
      pvalz <- ifelse(pnorm(ests/ses) > 0.5, 2*(1-pnorm(ests/ses)), 2*pnorm(ests/ses))
      
      output.mat <- cbind(coefficients=ests, std.errors=ses, p.values=pvalz)
      
    }
  }
  list(output.mat=output.mat, varcovmat=varcovmat, loglik=SemiParFit.p$logLik)
}







bdiscrdiscr_ARB <- function(params, respvec, VC, ps, AT = FALSE) {
  
  p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA
  eta1 <- VC$X1 %*% params[1:VC$X1.d2]
  eta2 <- VC$X2 %*% params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
  etad <- etas1 <- etas2 <- l.ln <- NULL
  if (is.null(VC$X3)) {
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 +
                                     1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 +
                                     2)]
    teta.st <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 3)]

    # new stuff for zinf
    zeta1.st <- params[(VC$X1.d2 + VC$X2.d2 + 4)]

    zeta2.st <- params[(VC$X1.d2 + VC$X2.d2 + 5)]
    # end of new stuff for zinf


  }
  if (!is.null(VC$X3)) {
    sigma21.st <- etas1 <- VC$X3 %*% params[(VC$X1.d2 +
                                               VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- VC$X4 %*% params[(VC$X1.d2 +
                                               VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 +
                                                                           VC$X3.d2 + VC$X4.d2)]
    teta.st <- etad <- VC$X5 %*% params[(VC$X1.d2 + VC$X2.d2 +
                                           VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 +
                                                                       VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
    # new stuff for zinf
    zeta1.st <- VC$X6 %*% params[(VC$X1.d2 + VC$X2.d2 +
                                  VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + 1):(VC$X1.d2 + VC$X2.d2 +
                                                                         VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2)]

    zeta2.st <- VC$X7 %*% params[(VC$X1.d2 + VC$X2.d2 +
                                     VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + 1):(VC$X1.d2 + VC$X2.d2 +
                                                                                       VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2)]
    # end of new stuff for zinf

  }

  # new stuff for zinf
  zinf.tr <- function(vrb.st, margin)
  {
    vrb.st <- ifelse(vrb.st > 3, 3, vrb.st)
    vrb.st <- ifelse(vrb.st < -18, -18, vrb.st)
    vrb <- sigmoid(vrb.st)
    list(vrb = vrb, vrb.st = vrb.st)
  }
  # end of new stuff for zinf


  sstr1 <- esp.tr(sigma21.st, VC$margins[1])
  sstr2 <- esp.tr(sigma22.st, VC$margins[2])
  sigma21.st <- sstr1$vrb.st
  sigma22.st <- sstr2$vrb.st
  sigma21 <- sstr1$vrb
  sigma22 <- sstr2$vrb

  # new stuff for zinf
  zinftr1 <- zinf.tr(zeta1.st, VC$margins[1])
  zinftr2 <- zinf.tr(zeta2.st, VC$margins[2])
  zeta1.st <- zinftr1$vrb.st
  zeta2.st <- zinftr2$vrb.st
  pr1 <- c(zinftr1$vrb)
  pr2<- c(zinftr2$vrb)
  # end of new stuff for zinf

  eta1 <- eta.tr(eta1, VC$margins[1])
  eta2 <- eta.tr(eta2, VC$margins[2])
  resT <- teta.tr(VC, teta.st)
  teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
  teta1 <- teta2 <- teta <- resT$teta
  Cop1 <- Cop2 <- VC$BivD
  nC1 <- nC2 <- VC$nC
  teta.ind1 <- as.logical(c(1, 0, round(runif(VC$n - 2))))
  teta.ind2 <- teta.ind1 == FALSE
  if (!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1) {
    teta.st1 <- teta.st[teta.ind1]
    teta.st2 <- teta.st[teta.ind2]
    teta1 <- teta[teta.ind1]
    teta2 <- teta[teta.ind2]
  }
  if (VC$BivD %in% VC$BivD2) {
    if (VC$BivD %in% VC$BivD2[c(1:4, 13:16)]) {
      teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov),
                          TRUE, FALSE
      )
    }
    if (VC$BivD %in% VC$BivD2[5:12]) {
      teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov) +
                            1, TRUE, FALSE)
    }
    teta.ind2 <- teta.ind1 == FALSE
    VC$my.env$signind <- ifelse(teta.ind1 == TRUE, 1, -1)
    teta1 <- teta[teta.ind1]
    teta2 <- -teta[teta.ind2]
    teta.st1 <- teta.st[teta.ind1]
    teta.st2 <- teta.st[teta.ind2]
    if (length(teta) == 1) {
      teta.ind2 <- teta.ind1 <- rep(TRUE, VC$n)
    }
    Cop1Cop2R <- Cop1Cop2(VC$BivD)
    Cop1 <- Cop1Cop2R$Cop1
    Cop2 <- Cop1Cop2R$Cop2
    nC1 <- VC$ct[which(VC$ct[, 1] == Cop1), 2]
    nC2 <- VC$ct[which(VC$ct[, 1] == Cop2), 2]
  }
  dHs1 <- distrHsDiscr(respvec$y1, eta1, sigma21, sigma21.st,
                       nu = 1, nu.st = 1, margin2 = VC$margins[1], naive = FALSE,
                       y2m = VC$y1m, min.dn = VC$min.dn, min.pr = VC$min.pr,
                       max.pr = VC$max.pr
  )
  dHs2 <- distrHsDiscr(respvec$y2, eta2, sigma22, sigma22.st,
                       nu = 1, nu.st = 1, margin2 = VC$margins[2], naive = FALSE,
                       y2m = VC$y2m, min.dn = VC$min.dn, min.pr = VC$min.pr,
                       max.pr = VC$max.pr
  )
  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2
  p1 <- dHs1$p2
  p2 <- dHs2$p2
  C11 <- C01 <- C10 <- C00 <- NA
  if (length(teta1) != 0) {
    C11[teta.ind1] <- mm(BiCDF(
      p1[teta.ind1], p2[teta.ind1],
      nC1, teta1, VC$dof
    ), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C01[teta.ind1] <- mm(BiCDF(
      mm(p1[teta.ind1] - pdf1[teta.ind1],
         min.pr = VC$min.pr, max.pr = VC$max.pr
      ), p2[teta.ind1],
      nC1, teta1, VC$dof
    ), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C10[teta.ind1] <- mm(BiCDF(
      p1[teta.ind1], mm(p2[teta.ind1] -
                          pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr),
      nC1, teta1, VC$dof
    ), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C00[teta.ind1] <- mm(BiCDF(
      mm(p1[teta.ind1] - pdf1[teta.ind1],
         min.pr = VC$min.pr, max.pr = VC$max.pr
      ), mm(p2[teta.ind1] -
              pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr),
      nC1, teta1, VC$dof
    ), min.pr = VC$min.pr, max.pr = VC$max.pr)
  }
  if (length(teta2) != 0) {
    C11[teta.ind2] <- mm(BiCDF(
      p1[teta.ind2], p2[teta.ind2],
      nC2, teta2, VC$dof
    ), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C01[teta.ind2] <- mm(BiCDF(
      mm(p1[teta.ind2] - pdf1[teta.ind2],
         min.pr = VC$min.pr, max.pr = VC$max.pr
      ), p2[teta.ind2],
      nC2, teta2, VC$dof
    ), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C10[teta.ind2] <- mm(BiCDF(
      p1[teta.ind2], mm(p2[teta.ind2] -
                          pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr),
      nC2, teta2, VC$dof
    ), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C00[teta.ind2] <- mm(BiCDF(
      mm(p1[teta.ind2] - pdf1[teta.ind2],
         min.pr = VC$min.pr, max.pr = VC$max.pr
      ), mm(p2[teta.ind2] -
              pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr),
      nC2, teta2, VC$dof
    ), min.pr = VC$min.pr, max.pr = VC$max.pr)
  }
  E <- mm(C11 - C01 - C10 + C00, min.pr = VC$min.pr, max.pr = VC$max.pr)
  l.par <- VC$weights * log(E)
  if (length(teta1) != 0) {
    dHC11F <- copgHs(p1[teta.ind1], p2[teta.ind1],
                     eta1 = NULL,
                     eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn,
                     min.pr = VC$min.pr, max.pr = VC$max.pr
    )
    dHC01F <- copgHs(
      mm(p1[teta.ind1] - pdf1[teta.ind1],
         min.pr = VC$min.pr, max.pr = VC$max.pr
      ), p2[teta.ind1],
      eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1,
      VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr,
      max.pr = VC$max.pr
    )
    dHC10F <- copgHs(p1[teta.ind1], mm(p2[teta.ind1] - pdf2[teta.ind1],
                                       min.pr = VC$min.pr, max.pr = VC$max.pr
    ),
    eta1 = NULL,
    eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn,
    min.pr = VC$min.pr, max.pr = VC$max.pr
    )
    dHC00F <- copgHs(
      mm(p1[teta.ind1] - pdf1[teta.ind1],
         min.pr = VC$min.pr, max.pr = VC$max.pr
      ), mm(p2[teta.ind1] -
              pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr),
      eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1,
      VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr,
      max.pr = VC$max.pr
    )
  }
  if (length(teta2) != 0) {
    dHC11S <- copgHs(p1[teta.ind2], p2[teta.ind2],
                     eta1 = NULL,
                     eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn,
                     min.pr = VC$min.pr, max.pr = VC$max.pr
    )
    dHC01S <- copgHs(
      mm(p1[teta.ind2] - pdf1[teta.ind2],
         min.pr = VC$min.pr, max.pr = VC$max.pr
      ), p2[teta.ind2],
      eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2,
      VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr,
      max.pr = VC$max.pr
    )
    dHC10S <- copgHs(p1[teta.ind2], mm(p2[teta.ind2] - pdf2[teta.ind2],
                                       min.pr = VC$min.pr, max.pr = VC$max.pr
    ),
    eta1 = NULL,
    eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn,
    min.pr = VC$min.pr, max.pr = VC$max.pr
    )
    dHC00S <- copgHs(
      mm(p1[teta.ind2] - pdf1[teta.ind2],
         min.pr = VC$min.pr, max.pr = VC$max.pr
      ), mm(p2[teta.ind2] -
              pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr),
      eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2,
      VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr,
      max.pr = VC$max.pr
    )
  }
  derC11.derp1 <- derC01.derp1 <- derC10.derp1 <- derC00.derp1 <- derC11.derp2 <- derC01.derp2 <- derC10.derp2 <- derC00.derp2 <- derC11.derthet <- derC01.derthet <- derC10.derthet <- derC00.derthet <- derteta.derteta.st <- NA
  if (length(teta1) != 0) {
    derC11.derp1[teta.ind1] <- dHC11F$c.copula.be1
    derC01.derp1[teta.ind1] <- dHC01F$c.copula.be1
    derC10.derp1[teta.ind1] <- dHC10F$c.copula.be1
    derC00.derp1[teta.ind1] <- dHC00F$c.copula.be1
    derC11.derp2[teta.ind1] <- dHC11F$c.copula.be2
    derC01.derp2[teta.ind1] <- dHC01F$c.copula.be2
    derC10.derp2[teta.ind1] <- dHC10F$c.copula.be2
    derC00.derp2[teta.ind1] <- dHC00F$c.copula.be2
    derC11.derthet[teta.ind1] <- dHC11F$c.copula.thet
    derC01.derthet[teta.ind1] <- dHC01F$c.copula.thet
    derC10.derthet[teta.ind1] <- dHC10F$c.copula.thet
    derC00.derthet[teta.ind1] <- dHC00F$c.copula.thet
    derteta.derteta.st[teta.ind1] <- dHC11F$derteta.derteta.st
  }
  if (length(teta2) != 0) {
    derC11.derp1[teta.ind2] <- dHC11S$c.copula.be1
    derC01.derp1[teta.ind2] <- dHC01S$c.copula.be1
    derC10.derp1[teta.ind2] <- dHC10S$c.copula.be1
    derC00.derp1[teta.ind2] <- dHC00S$c.copula.be1
    derC11.derp2[teta.ind2] <- dHC11S$c.copula.be2
    derC01.derp2[teta.ind2] <- dHC01S$c.copula.be2
    derC10.derp2[teta.ind2] <- dHC10S$c.copula.be2
    derC00.derp2[teta.ind2] <- dHC00S$c.copula.be2
    derC11.derthet[teta.ind2] <- dHC11S$c.copula.thet
    derC01.derthet[teta.ind2] <- dHC01S$c.copula.thet
    derC10.derthet[teta.ind2] <- dHC10S$c.copula.thet
    derC00.derthet[teta.ind2] <- dHC00S$c.copula.thet
    derteta.derteta.st[teta.ind2] <- dHC11S$derteta.derteta.st
  }
  derpdf1.dereta1 <- dHs1$derpdf2.dereta2
  derp1.dereta1 <- dHs1$derp2.dereta2
  derp1m1.dereta1 <- derp1.dereta1 - derpdf1.dereta1
  derpdf1.dersigma21.st <- dHs1$derpdf2.dersigma2.st
  derp1.dersigma21.st <- dHs1$derp2.dersigma.st
  derp1m1.dersigma21.st <- derp1.dersigma21.st - derpdf1.dersigma21.st
  derpdf2.dereta2 <- dHs2$derpdf2.dereta2
  derp2.dereta2 <- dHs2$derp2.dereta2
  derp2m1.dereta2 <- derp2.dereta2 - derpdf2.dereta2
  derpdf2.dersigma22.st <- dHs2$derpdf2.dersigma2.st
  derp2.dersigma22.st <- dHs2$derp2.dersigma.st
  derp2m1.dersigma22.st <- derp2.dersigma22.st - derpdf2.dersigma22.st
  fE1 <- (derC11.derp1 - derC10.derp1) * derp1.dereta1 - (derC01.derp1 -
                                                            derC00.derp1) * derp1m1.dereta1
  fE2 <- (derC11.derp2 - derC01.derp2) * derp2.dereta2 - (derC10.derp2 -
                                                            derC00.derp2) * derp2m1.dereta2
  fEt <- (derC11.derthet - derC01.derthet - derC10.derthet +
            derC00.derthet) * derteta.derteta.st
  fE1s <- (derC11.derp1 - derC10.derp1) * derp1.dersigma21.st -
    (derC01.derp1 - derC00.derp1) * derp1m1.dersigma21.st
  fE2s <- (derC11.derp2 - derC01.derp2) * derp2.dersigma22.st -
    (derC10.derp2 - derC00.derp2) * derp2m1.dersigma22.st
  dl.dbe1 <- VC$weights * fE1 / E
  dl.dbe2 <- VC$weights * fE2 / E
  dl.dsigma21.st <- VC$weights * fE1s / E
  dl.dsigma22.st <- VC$weights * fE2s / E
  dl.dteta.st <- VC$weights * fEt / E
  der2C11.derp1p1 <- der2C01.derp1p1 <- der2C10.derp1p1 <- der2C00.derp1p1 <- der2C11.derp2p2 <- der2C01.derp2p2 <- der2C10.derp2p2 <- der2C00.derp2p2 <- der2C11.derp1p2 <- der2C01.derp1p2 <- der2C10.derp1p2 <- der2C00.derp1p2 <- der2C11.derp1t <- der2C01.derp1t <- der2C10.derp1t <- der2C00.derp1t <- der2C11.derp2t <- der2C01.derp2t <- der2C10.derp2t <- der2C00.derp2t <- der2C11.derthet2 <- der2C01.derthet2 <- der2C10.derthet2 <- der2C00.derthet2 <- der2teta.derteta.stteta.st <- NA
  if (length(teta1) != 0) {
    der2C11.derp1p1[teta.ind1] <- dHC11F$c.copula2.be1
    der2C01.derp1p1[teta.ind1] <- dHC01F$c.copula2.be1
    der2C10.derp1p1[teta.ind1] <- dHC10F$c.copula2.be1
    der2C00.derp1p1[teta.ind1] <- dHC00F$c.copula2.be1
    der2C11.derp2p2[teta.ind1] <- dHC11F$c.copula2.be2
    der2C01.derp2p2[teta.ind1] <- dHC01F$c.copula2.be2
    der2C10.derp2p2[teta.ind1] <- dHC10F$c.copula2.be2
    der2C00.derp2p2[teta.ind1] <- dHC00F$c.copula2.be2
    der2C11.derp1p2[teta.ind1] <- dHC11F$c.copula2.be1be2
    der2C01.derp1p2[teta.ind1] <- dHC01F$c.copula2.be1be2
    der2C10.derp1p2[teta.ind1] <- dHC10F$c.copula2.be1be2
    der2C00.derp1p2[teta.ind1] <- dHC00F$c.copula2.be1be2
    der2C11.derp1t[teta.ind1] <- dHC11F$c.copula2.be1t
    der2C01.derp1t[teta.ind1] <- dHC01F$c.copula2.be1t
    der2C10.derp1t[teta.ind1] <- dHC10F$c.copula2.be1t
    der2C00.derp1t[teta.ind1] <- dHC00F$c.copula2.be1t
    der2C11.derp2t[teta.ind1] <- dHC11F$c.copula2.be2t
    der2C01.derp2t[teta.ind1] <- dHC01F$c.copula2.be2t
    der2C10.derp2t[teta.ind1] <- dHC10F$c.copula2.be2t
    der2C00.derp2t[teta.ind1] <- dHC00F$c.copula2.be2t
    der2C11.derthet2[teta.ind1] <- dHC11F$bit1.th2ATE
    der2C01.derthet2[teta.ind1] <- dHC01F$bit1.th2ATE
    der2C10.derthet2[teta.ind1] <- dHC10F$bit1.th2ATE
    der2C00.derthet2[teta.ind1] <- dHC00F$bit1.th2ATE
    der2teta.derteta.stteta.st[teta.ind1] <- dHC11F$der2teta.derteta.stteta.st
  }
  if (length(teta2) != 0) {
    der2C11.derp1p1[teta.ind2] <- dHC11S$c.copula2.be1
    der2C01.derp1p1[teta.ind2] <- dHC01S$c.copula2.be1
    der2C10.derp1p1[teta.ind2] <- dHC10S$c.copula2.be1
    der2C00.derp1p1[teta.ind2] <- dHC00S$c.copula2.be1
    der2C11.derp2p2[teta.ind2] <- dHC11S$c.copula2.be2
    der2C01.derp2p2[teta.ind2] <- dHC01S$c.copula2.be2
    der2C10.derp2p2[teta.ind2] <- dHC10S$c.copula2.be2
    der2C00.derp2p2[teta.ind2] <- dHC00S$c.copula2.be2
    der2C11.derp1p2[teta.ind2] <- dHC11S$c.copula2.be1be2
    der2C01.derp1p2[teta.ind2] <- dHC01S$c.copula2.be1be2
    der2C10.derp1p2[teta.ind2] <- dHC10S$c.copula2.be1be2
    der2C00.derp1p2[teta.ind2] <- dHC00S$c.copula2.be1be2
    der2C11.derp1t[teta.ind2] <- dHC11S$c.copula2.be1t
    der2C01.derp1t[teta.ind2] <- dHC01S$c.copula2.be1t
    der2C10.derp1t[teta.ind2] <- dHC10S$c.copula2.be1t
    der2C00.derp1t[teta.ind2] <- dHC00S$c.copula2.be1t
    der2C11.derp2t[teta.ind2] <- dHC11S$c.copula2.be2t
    der2C01.derp2t[teta.ind2] <- dHC01S$c.copula2.be2t
    der2C10.derp2t[teta.ind2] <- dHC10S$c.copula2.be2t
    der2C00.derp2t[teta.ind2] <- dHC00S$c.copula2.be2t
    der2C11.derthet2[teta.ind2] <- dHC11S$bit1.th2ATE
    der2C01.derthet2[teta.ind2] <- dHC01S$bit1.th2ATE
    der2C10.derthet2[teta.ind2] <- dHC10S$bit1.th2ATE
    der2C00.derthet2[teta.ind2] <- dHC00S$bit1.th2ATE
    der2teta.derteta.stteta.st[teta.ind2] <- dHC11S$der2teta.derteta.stteta.st
  }
  der2pdf1.dereta1 <- dHs1$der2pdf2.dereta2
  der2p1.dereta1eta1 <- dHs1$der2p2.dereta2eta2
  der2p1m1.dereta1eta1 <- der2p1.dereta1eta1 - der2pdf1.dereta1
  der2pdf2.dereta2 <- dHs2$der2pdf2.dereta2
  der2p2.dereta2eta2 <- dHs2$der2p2.dereta2eta2
  der2p2m1.dereta2eta2 <- der2p2.dereta2eta2 - der2pdf2.dereta2
  der2pdf1.dersigma21.st2 <- dHs1$der2pdf2.dersigma2.st2
  der2p1.dersigma21.st2 <- dHs1$der2p2.dersigma2.st2
  der2p1m1.dersigma21.st2 <- der2p1.dersigma21.st2 - der2pdf1.dersigma21.st2
  der2pdf2.dersigma22.st2 <- dHs2$der2pdf2.dersigma2.st2
  der2p2.dersigma22.st2 <- dHs2$der2p2.dersigma2.st2
  der2p2m1.dersigma22.st2 <- der2p2.dersigma22.st2 - der2pdf2.dersigma22.st2
  der2pdf1.dereta1dersigma21.st <- dHs1$der2pdf2.dereta2dersigma2.st
  der2p1.dereta1dersigma21.st <- dHs1$der2p2.dereta2dersigma2.st
  der2p1m1.dereta1dersigma21.st <- der2p1.dereta1dersigma21.st -
    der2pdf1.dereta1dersigma21.st
  der2pdf2.dereta2dersigma22.st <- dHs2$der2pdf2.dereta2dersigma2.st
  der2p2.dereta2dersigma22.st <- dHs2$der2p2.dereta2dersigma2.st
  der2p2m1.dereta2dersigma22.st <- der2p2.dereta2dersigma22.st -
    der2pdf2.dereta2dersigma22.st
  d2l.be1.be1 <- -VC$weights * ((E * (der2C11.derp1p1 * derp1.dereta1^2 +
                                        derC11.derp1 * der2p1.dereta1eta1 - (der2C01.derp1p1 *
                                                                               derp1m1.dereta1^2 + derC01.derp1 * der2p1m1.dereta1eta1) -
                                        (der2C10.derp1p1 * derp1.dereta1^2 + derC10.derp1 *
                                           der2p1.dereta1eta1) + (der2C00.derp1p1 * derp1m1.dereta1^2 +
                                                                    derC00.derp1 * der2p1m1.dereta1eta1)) - fE1^2) / E^2)

  d2l.be2.be2 <- -VC$weights * ((E * (der2C11.derp2p2 * derp2.dereta2^2 +
                                        derC11.derp2 * der2p2.dereta2eta2 - (der2C01.derp2p2 *
                                                                               derp2.dereta2^2 + derC01.derp2 * der2p2.dereta2eta2) -
                                        (der2C10.derp2p2 * derp2m1.dereta2^2 + derC10.derp2 *
                                           der2p2m1.dereta2eta2) + (der2C00.derp2p2 * derp2m1.dereta2^2 +
                                                                      derC00.derp2 * der2p2m1.dereta2eta2)) - fE2^2) / E^2)
  d2l.sigma21.sigma21 <- -VC$weights * ((E * (der2C11.derp1p1 *
                                                derp1.dersigma21.st^2 + derC11.derp1 * der2p1.dersigma21.st2 -
                                                (der2C01.derp1p1 * derp1m1.dersigma21.st^2 + derC01.derp1 *
                                                   der2p1m1.dersigma21.st2) - (der2C10.derp1p1 * derp1.dersigma21.st^2 +
                                                                                 derC10.derp1 * der2p1.dersigma21.st2) + (der2C00.derp1p1 *
                                                                                                                            derp1m1.dersigma21.st^2 + derC00.derp1 * der2p1m1.dersigma21.st2)) -
                                           fE1s^2) / E^2)
  d2l.sigma22.sigma22 <- -VC$weights * ((E * (der2C11.derp2p2 *
                                                derp2.dersigma22.st^2 + derC11.derp2 * der2p2.dersigma22.st2 -
                                                (der2C01.derp2p2 * derp2.dersigma22.st^2 + derC01.derp2 *
                                                   der2p2.dersigma22.st2) - (der2C10.derp2p2 * derp2m1.dersigma22.st^2 +
                                                                               derC10.derp2 * der2p2m1.dersigma22.st2) + (der2C00.derp2p2 *
                                                                                                                            derp2m1.dersigma22.st^2 + derC00.derp2 * der2p2m1.dersigma22.st2)) -
                                           fE2s^2) / E^2)
  d2l.rho.rho <- -VC$weights * ((E * (der2C11.derthet2 * derteta.derteta.st^2 +
                                        derC11.derthet * der2teta.derteta.stteta.st - (der2C01.derthet2 *
                                                                                         derteta.derteta.st^2 + derC01.derthet * der2teta.derteta.stteta.st) -
                                        (der2C10.derthet2 * derteta.derteta.st^2 + derC10.derthet *
                                           der2teta.derteta.stteta.st) + (der2C00.derthet2 *
                                                                            derteta.derteta.st^2 + derC00.derthet * der2teta.derteta.stteta.st)) -
                                   fEt^2) / E^2)
  d2l.be1.be2 <- -VC$weights * ((E * (der2C11.derp1p2 * derp1.dereta1 *
                                        derp2.dereta2 - der2C01.derp1p2 * derp1m1.dereta1 *
                                        derp2.dereta2 - der2C10.derp1p2 * derp1.dereta1 * derp2m1.dereta2 +
                                        der2C00.derp1p2 * derp1m1.dereta1 * derp2m1.dereta2) -
                                   fE1 * fE2) / E^2)
  d2l.be1.sigma22 <- -VC$weights * ((E * (der2C11.derp1p2 *
                                            derp1.dereta1 * derp2.dersigma22.st - der2C01.derp1p2 *
                                            derp1m1.dereta1 * derp2.dersigma22.st - der2C10.derp1p2 *
                                            derp1.dereta1 * derp2m1.dersigma22.st + der2C00.derp1p2 *
                                            derp1m1.dereta1 * derp2m1.dersigma22.st) - fE1 * fE2s) / E^2)
  d2l.be1.rho <- -VC$weights * ((E * (der2C11.derp1t * derp1.dereta1 *
                                        derteta.derteta.st - der2C01.derp1t * derp1m1.dereta1 *
                                        derteta.derteta.st - der2C10.derp1t * derp1.dereta1 *
                                        derteta.derteta.st + der2C00.derp1t * derp1m1.dereta1 *
                                        derteta.derteta.st) - fE1 * fEt) / E^2)
  d2l.be1.sigma21 <- -VC$weights * ((E * (der2C11.derp1p1 *
                                            derp1.dereta1 * derp1.dersigma21.st + derC11.derp1 *
                                            der2p1.dereta1dersigma21.st - (der2C01.derp1p1 * derp1m1.dereta1 *
                                                                             derp1m1.dersigma21.st + derC01.derp1 * der2p1m1.dereta1dersigma21.st) -
                                            (der2C10.derp1p1 * derp1.dereta1 * derp1.dersigma21.st +
                                               derC10.derp1 * der2p1.dereta1dersigma21.st) + (der2C00.derp1p1 *
                                                                                                derp1m1.dereta1 * derp1m1.dersigma21.st + derC00.derp1 *
                                                                                                der2p1m1.dereta1dersigma21.st)) - fE1 * fE1s) / E^2)
  d2l.be2.sigma22 <- -VC$weights * ((E * (der2C11.derp2p2 *
                                            derp2.dereta2 * derp2.dersigma22.st + derC11.derp2 *
                                            der2p2.dereta2dersigma22.st - (der2C01.derp2p2 * derp2.dereta2 *
                                                                             derp2.dersigma22.st + derC01.derp2 * der2p2.dereta2dersigma22.st) -
                                            (der2C10.derp2p2 * derp2m1.dereta2 * derp2m1.dersigma22.st +
                                               derC10.derp2 * der2p2m1.dereta2dersigma22.st) +
                                            (der2C00.derp2p2 * derp2m1.dereta2 * derp2m1.dersigma22.st +
                                               derC00.derp2 * der2p2m1.dereta2dersigma22.st)) -
                                       fE2 * fE2s) / E^2)
  d2l.be2.sigma21 <- -VC$weights * ((E * (der2C11.derp1p2 *
                                            derp1.dersigma21.st * derp2.dereta2 - der2C01.derp1p2 *
                                            derp1m1.dersigma21.st * derp2.dereta2 - der2C10.derp1p2 *
                                            derp1.dersigma21.st * derp2m1.dereta2 + der2C00.derp1p2 *
                                            derp1m1.dersigma21.st * derp2m1.dereta2) - fE1s * fE2) / E^2)
  d2l.be2.rho <- -VC$weights * ((E * (der2C11.derp2t * derp2.dereta2 *
                                        derteta.derteta.st - der2C01.derp2t * derp2.dereta2 *
                                        derteta.derteta.st - der2C10.derp2t * derp2m1.dereta2 *
                                        derteta.derteta.st + der2C00.derp2t * derp2m1.dereta2 *
                                        derteta.derteta.st) - fE2 * fEt) / E^2)
  d2l.rho.sigma21 <- -VC$weights * ((E * (der2C11.derp1t *
                                            derp1.dersigma21.st * derteta.derteta.st - der2C01.derp1t *
                                            derp1m1.dersigma21.st * derteta.derteta.st - der2C10.derp1t *
                                            derp1.dersigma21.st * derteta.derteta.st + der2C00.derp1t *
                                            derp1m1.dersigma21.st * derteta.derteta.st) - fE1s *
                                       fEt) / E^2)
  d2l.rho.sigma22 <- -VC$weights * ((E * (der2C11.derp2t *
                                            derp2.dersigma22.st * derteta.derteta.st - der2C01.derp2t *
                                            derp2.dersigma22.st * derteta.derteta.st - der2C10.derp2t *
                                            derp2m1.dersigma22.st * derteta.derteta.st + der2C00.derp2t *
                                            derp2m1.dersigma22.st * derteta.derteta.st) - fE2s *
                                       fEt) / E^2)
  d2l.sigma21.sigma22 <- -VC$weights * ((E * (der2C11.derp1p2 *
                                                derp1.dersigma21.st * derp2.dersigma22.st - der2C01.derp1p2 *
                                                derp1m1.dersigma21.st * derp2.dersigma22.st - der2C10.derp1p2 *
                                                derp1.dersigma21.st * derp2m1.dersigma22.st + der2C00.derp1p2 *
                                                derp1m1.dersigma21.st * derp2m1.dersigma22.st) - fE1s *
                                           fE2s) / E^2)


  S.h <- ps$S.h
  if (length(S.h) != 1) {
    S.h1 <- 0.5 * crossprod(params, S.h) %*% params
    S.h2 <- S.h %*% params
  } else {
    S.h <- S.h1 <- S.h2 <- 0
  }


  ########## Incorporating kappas ##########
  Xd2_len <- VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2
  beta1_idcs <- 1:VC$X1.d2
  beta2_idcs <- (VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)
  alpha1_idcs <- (VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)
  alpha2_idcs <- (VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)
  tau_idcs <- (VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(Xd2_len)

  # new zinf stuff
  kappa1_idcs <- (Xd2_len+1):(Xd2_len + VC$X6.d2)
  kappa2_idcs <- (Xd2_len + VC$X6.d2 + 1):(Xd2_len + VC$X6.d2 + VC$X7.d2)
  # end of new zinf stuff



  y1zero <- ifelse(respvec$y1 == 0, 1, 0)
  y2zero <- ifelse(respvec$y2 == 0, 1, 0)
  bothzero <- ifelse((respvec$y1 == 0) & (respvec$y2 == 0), 1, 0)



  ########## likelihood (length n) vector for copula portion ##########
  dcop <- E

  ##########  likelihood (length n) vector for nb1 marginal portions ##########
  dnb1 <- c(dHs1$pdf2)
  dnb2 <- c(dHs2$pdf2)

  ##########  overall likelihood (length n) vector with kappas incorporated ##########
  doverall <- (1 - pr1) * (1 - pr2) * dcop + (1 - pr1) * pr2 * dnb1 * y2zero + pr1 * (1 - pr2) * dnb2 * y1zero + pr1 * pr2 * bothzero


  # kappa1 first derivatives
  doverall.dpr1 <- -(1 - pr2) * dcop - pr2 * dnb1 * y2zero + (1 - pr2) * dnb2 * y1zero + pr2 * bothzero
  dpr1.dzeta1 <- pr1 * (1 - pr1)
  doverall.dzeta1 <- doverall.dpr1 * dpr1.dzeta1
  
  dl.dkappa1 <- VC$X6 * c(1/doverall) * doverall.dzeta1
  
  
  # kappa2 first derivatives
  doverall.dpr2 <- -(1 - pr1) * dcop + (1 - pr1) * dnb1 * y2zero - pr1 * dnb2 * y1zero + pr1 * bothzero
  dpr2.dzeta2 <- pr2 * (1 - pr2)
  doverall.dzeta2 <- doverall.dpr2 * dpr2.dzeta2
  
  dl.dkappa2 <- VC$X7 * c(1/doverall) * doverall.dzeta2
  
  
  
  # beta1 first derivatives
  dcop.deta1 <- c(fE1)
  dnb1.deta1 <- c(derpdf1.dereta1)
  
  doverall.deta1 <- (1 - pr1) * (1 - pr2) * dcop.deta1 + (1 - pr1) * pr2 * dnb1.deta1 * y2zero
  dl.dbeta1 <- VC$X1 * c(1/doverall) * doverall.deta1
  
  
  # beta2 first derivatives
  dcop.deta2 <- c(fE2)
  dnb2.deta2 <- c(derpdf2.dereta2)
  
  doverall.deta2 <- (1 - pr1) * (1 - pr2) * dcop.deta2 + pr1 * (1 - pr2) * dnb2.deta2 * y1zero
  dl.dbeta2 <- VC$X2 * c(1/doverall) * doverall.deta2
  
  
  # alpha1 first derivatives
  dcop.dsigma21 <- c(fE1s)
  dnb1.dsigma21 <- c(derpdf1.dersigma21.st)
  
  doverall.dsigma21 <- (1 - pr1) * (1 - pr2) * dcop.dsigma21 + (1 - pr1) * pr2 * dnb1.dsigma21 * y2zero
  dl.dalpha1 <- VC$X3 * c(1/doverall) * doverall.dsigma21
  
  
  # alpha2 first derivatives
  dcop.dsigma22 <- c(fE2s)
  dnb2.dsigma22 <- c(derpdf2.dersigma22.st)
  
  doverall.dsigma22 <- (1 - pr1) * (1 - pr2) * dcop.dsigma22 + pr1 * (1 - pr2) * dnb2.dsigma22 * y1zero
  dl.dalpha2 <- VC$X4 * c(1/doverall) * doverall.dsigma22
  
  
  # tau first derivatives
  dcop.dteta <- c(fEt)
  
  doverall.dteta <- (1 - pr1) * (1 - pr2) * dcop.dteta
  dl.dtau <- VC$X5 * c(1/doverall) * doverall.dteta
  
  

  
  G_temp <- -colSums(cbind(cbind(dl.dbeta1, dl.dbeta2, dl.dalpha1, dl.dalpha2, dl.dtau, dl.dkappa1, dl.dkappa2)))
  

  # trying hessian for kappa1s
  hessmat <- matrix(0, nrow=kappa2_idcs[length(kappa2_idcs)], ncol=kappa2_idcs[length(kappa2_idcs)])
  
  hessmat[kappa1_idcs, kappa1_idcs] <- crossprod(VC$X6 * c(1/doverall) * doverall.dzeta1, VC$X6 * ((1-2*pr1)-c(1/doverall)*doverall.dzeta1))
  hessmat[kappa2_idcs, kappa2_idcs] <- crossprod(VC$X7 * c(1/doverall) * doverall.dzeta2, VC$X7 * ((1-2*pr2)-c(1/doverall)*doverall.dzeta2))
  hessmat[kappa1_idcs, kappa2_idcs] <- hessmat[kappa2_idcs, kappa1_idcs] <- crossprod(VC$X6 * c(1/doverall) * dpr1.dzeta1 * dpr2.dzeta2, VC$X7*((dcop - dnb1*y2zero - dnb2*y1zero + bothzero) - c(1/doverall)*doverall.dpr1*doverall.dpr2))

  hessmat[beta2_idcs, kappa1_idcs] <- crossprod(VC$X2 * c(1/doverall) * dpr1.dzeta1, VC$X6 * ((-(1-pr2)*dcop.deta2 + (1-pr2)*dnb2.deta2*y1zero)-c(1/doverall)*doverall.dpr1*doverall.deta2))
  hessmat[beta1_idcs, kappa1_idcs] <- crossprod(VC$X1 * c(1/doverall) * dpr1.dzeta1, VC$X6 * ((-(1-pr2)*dcop.deta1 -pr2*dnb1.deta1*y2zero)-c(1/doverall)*doverall.dpr1*doverall.deta1))
  hessmat[beta2_idcs, kappa2_idcs] <- crossprod(VC$X2 * c(1/doverall) * dpr2.dzeta2, VC$X7 * ((-(1-pr1)*dcop.deta2 - pr1*dnb2.deta2*y1zero)-c(1/doverall)*doverall.dpr2*doverall.deta2))
  hessmat[beta1_idcs, kappa2_idcs] <- crossprod(VC$X1 * c(1/doverall) * dpr2.dzeta2, VC$X7 * ((-(1-pr1)*dcop.deta1 + (1-pr1)*dnb1.deta1*y2zero)-c(1/doverall)*doverall.dpr2*doverall.deta1))
  
  hessmat[alpha2_idcs, kappa1_idcs] <- crossprod(VC$X4 * c(1/doverall) * dpr1.dzeta1, VC$X6 * ((-(1-pr2)*dcop.dsigma22 + (1-pr2)*dnb2.dsigma22*y1zero)-c(1/doverall)*doverall.dpr1*doverall.dsigma22))
  hessmat[alpha1_idcs, kappa1_idcs] <- crossprod(VC$X3 * c(1/doverall) * dpr1.dzeta1, VC$X6 * ((-(1-pr2)*dcop.dsigma21 -pr2*dnb1.dsigma21*y2zero)-c(1/doverall)*doverall.dpr1*doverall.dsigma21))
  hessmat[alpha2_idcs, kappa2_idcs] <- crossprod(VC$X4 * c(1/doverall) * dpr2.dzeta2, VC$X7 * ((-(1-pr1)*dcop.dsigma22 - pr1*dnb2.dsigma22*y1zero)-c(1/doverall)*doverall.dpr2*doverall.dsigma22))
  hessmat[alpha1_idcs, kappa2_idcs] <- crossprod(VC$X3 * c(1/doverall) * dpr2.dzeta2, VC$X7 * ((-(1-pr1)*dcop.dsigma21 + (1-pr1)*dnb1.dsigma21*y2zero)-c(1/doverall)*doverall.dpr2*doverall.dsigma21))
  
  hessmat[tau_idcs, kappa1_idcs] <- crossprod(VC$X5 * c(1/doverall) * dpr1.dzeta1, VC$X6 * ((-(1-pr2)*dcop.dteta) - c(1/doverall)*doverall.dpr1*doverall.dteta))
  hessmat[tau_idcs, kappa2_idcs] <- crossprod(VC$X5 * c(1/doverall) * dpr2.dzeta2, VC$X7 * ((-(1-pr1)*dcop.dteta) - c(1/doverall)*doverall.dpr2*doverall.dteta))
  
  
  
  hessmat_cop <- matrix(0, nrow=nrow(hessmat), ncol=ncol(hessmat))
  

  hessmat_cop[beta1_idcs, beta1_idcs] <- crossprod(VC$X1*(1-pr1)*(1-pr2), VC$X1 * c(fE1^2 / E - d2l.be1.be1 * E) * c(1 / doverall))
  hessmat_cop[beta2_idcs, beta2_idcs] <- crossprod(VC$X2*(1-pr1)*(1-pr2), VC$X2 * c(fE2^2 / E - d2l.be2.be2 * E) * c(1 / doverall))
  hessmat_cop[beta1_idcs, beta2_idcs] <- crossprod(VC$X1*(1-pr1)*(1-pr2), VC$X2 * c(fE1 * fE2 / E - d2l.be1.be2 * E) * c(1 / doverall))

  hessmat_cop[alpha1_idcs, alpha1_idcs] <- crossprod(VC$X3*(1-pr1)*(1-pr2), VC$X3 * c(fE1s^2 / E - d2l.sigma21.sigma21 * E) * c(1 / doverall))
  hessmat_cop[alpha2_idcs, alpha2_idcs] <- crossprod(VC$X4*(1-pr1)*(1-pr2), VC$X4 * c(fE2s^2 / E - d2l.sigma22.sigma22 * E) * c(1 / doverall))
  hessmat_cop[alpha1_idcs, alpha2_idcs] <- crossprod(VC$X3*(1-pr1)*(1-pr2), VC$X4 * c(fE1s * fE2s / E - d2l.sigma21.sigma22 * E) * c(1 / doverall))

  hessmat_cop[tau_idcs, tau_idcs] <- crossprod(VC$X5*(1-pr1)*(1-pr2), VC$X5 * c(fEt^2 / E - d2l.rho.rho * E) * c(1 / doverall))

  hessmat_cop[beta1_idcs, alpha1_idcs] <- crossprod(VC$X1*(1-pr1)*(1-pr2), VC$X3 * c(fE1 * fE1s / E - d2l.be1.sigma21 * E) * c(1 / doverall))
  hessmat_cop[beta1_idcs, alpha2_idcs] <- crossprod(VC$X1*(1-pr1)*(1-pr2), VC$X4 * c(fE1 * fE2s / E - d2l.be1.sigma22 * E) * c(1 / doverall))
  hessmat_cop[beta2_idcs, alpha1_idcs] <- crossprod(VC$X2*(1-pr1)*(1-pr2), VC$X3 * c(fE2 * fE1s / E - d2l.be2.sigma21 * E) * c(1 / doverall))
  hessmat_cop[beta2_idcs, alpha2_idcs] <- crossprod(VC$X2*(1-pr1)*(1-pr2), VC$X4 * c(fE2 * fE2s / E - d2l.be2.sigma22 * E) * c(1 / doverall))

  hessmat_cop[beta1_idcs, tau_idcs] <- crossprod(VC$X1*(1-pr1)*(1-pr2), VC$X5 * c(fE1 * fEt / E - d2l.be1.rho * E) * c(1 / doverall))
  hessmat_cop[beta2_idcs, tau_idcs] <- crossprod(VC$X2*(1-pr1)*(1-pr2), VC$X5 * c(fE2 * fEt / E - d2l.be2.rho * E) * c(1 / doverall))

  hessmat_cop[alpha1_idcs, tau_idcs] <- crossprod(VC$X3*(1-pr1)*(1-pr2), VC$X5 * c(fE1s * fEt / E - d2l.rho.sigma21 * E) * c(1 / doverall))
  hessmat_cop[alpha2_idcs, tau_idcs] <- crossprod(VC$X4*(1-pr1)*(1-pr2), VC$X5 * c(fE2s * fEt / E - d2l.rho.sigma22 * E) * c(1 / doverall))

  

  
  hessmat_nb1 <- matrix(0, nrow=nrow(hessmat), ncol=ncol(hessmat))
  hessmat_nb2 <- matrix(0, nrow=nrow(hessmat), ncol=ncol(hessmat))
  

  hessmat_nb1[beta1_idcs, beta1_idcs] <- crossprod(VC$X1*(1 - pr1)*pr2*y2zero, VC$X1 * c(der2pdf1.dereta1) * c(1 / doverall) * y2zero)
  hessmat_nb1[alpha1_idcs, alpha1_idcs] <- crossprod(VC$X3*(1 - pr1)*pr2*y2zero, VC$X3 * c(der2pdf1.dersigma21.st2) * c(1 / doverall) * y2zero)
  hessmat_nb1[beta1_idcs, alpha1_idcs] <- crossprod(VC$X1*(1 - pr1)*pr2*y2zero, VC$X3 * c(der2pdf1.dereta1dersigma21.st) * c(1 / doverall) * y2zero)

  hessmat_nb2[beta2_idcs, beta2_idcs] <- crossprod(VC$X2*pr1*(1 - pr2)*y1zero, VC$X2 * c(der2pdf2.dereta2) * c(1 / doverall) * y1zero)
  hessmat_nb2[alpha2_idcs, alpha2_idcs] <- crossprod(VC$X4*pr1*(1 - pr2)*y1zero, VC$X4 * c(der2pdf2.dersigma22.st2) * c(1 / doverall) * y1zero)
  hessmat_nb2[beta2_idcs, alpha2_idcs] <- crossprod(VC$X2*pr1*(1 - pr2)*y1zero, VC$X4 * c(der2pdf2.dereta2dersigma22.st) * c(1 / doverall) * y1zero)

  
  lik_grad <- cbind(cbind(dl.dbeta1, dl.dbeta2, dl.dalpha1, dl.dalpha2, dl.dtau))
  lik_grad <- cbind(lik_grad, matrix(0, nrow=nrow(lik_grad), ncol=length(kappa1_idcs)+length(kappa2_idcs)))
  
  piece1 <- crossprod(lik_grad, (-1)*lik_grad)# * c(-1 / doverall^2))

  piece2 <- hessmat + hessmat_cop + hessmat_nb1 + hessmat_nb2


  H_temp <- -1*reflect_triangle(piece1 + piece2, from="upper")
  
  
  l.par_temp <- c(log(doverall))
  res_temp <- -sum(l.par_temp)
  S.res_temp <- res_temp
  
  res_temp <- S.res_temp + S.h1
  if (length(S.h2) > 1){
    G_temp <- G_temp + S.h2
  }
  if (length(S.h) > 1){
    H_temp <- H_temp + S.h
  }
  
  

  list(
    value = res_temp, gradient = G_temp, hessian = H_temp, S.h = S.h,
    S.h1 = S.h1, S.h2 = S.h2, l = S.res_temp, l.ln = l.ln, l.par = l.par_temp,
    ps = ps, eta1 = eta1, eta2 = eta2, etad = etad, etas1 = etas1,
    etas2 = etas2, BivD = VC$BivD
  )
}


tmpfun <- get("gjrm", envir = asNamespace("GJRM"))

environment(FAST.CoExpress.nb) <- environment(tmpfun)




tmpfun <- get("overall.svG", envir = asNamespace("GJRM"))

environment(overall.svG) <- environment(tmpfun)

assignInNamespace("overall.svG", overall.svG, ns = "GJRM")

