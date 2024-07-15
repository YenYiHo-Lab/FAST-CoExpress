source("C:/Users/andyb/Downloads/FAST_CoExpress_functions.R")
library(dplyr)


n <- 2500

## true coefficient values. Note we will make marginal 1 be not zero inflated ##
# HSPA1B CFD
beta1_use <- c(1.59, -0.027, 0.097, -0.011)
beta2_use <- c(2.79, -0.017, 0.25, 0.077)
alpha1_use <- c(0.053, 0.311)
alpha2_use <- c(0.88, -0.067)
tau_use <- c(0.068, -0.0018, 0.076, -0.071)
kappa1_use <- c(-20, 0)
kappa2_use <- c(-2.29, .99)

paramvec <- c(beta1_use, beta2_use, alpha1_use, alpha2_use, tau_use, kappa1_use, kappa2_use)
names(paramvec) <- c("beta01", "beta11", "beta21", "beta31", "beta02", "beta12", "beta22", "beta32", "alpha01", "alpha11", "alpha02", "alpha12", "tau0", "tau1", "tau2", "tau3", "kappa01", "kappa11", "kappa02", "kappa12")


## specify how the different parameters depend on covariates ##
eq1 <- y1 ~ x1 + x2 + x3
eq2 <- y2 ~ x1 + x2 + x3
eq3 <- ~ x2
eq4 <- ~ x2
eq5 <- ~ x1 + x2 + x3
eqlist <- list(eq1, eq2, eq3, eq4, eq5)



## Run B iterations of robustness study ##
B <- 5 # 10000

mse <- mbe <- cilen <- covrg <- rep(0, length(paramvec))


for (b in 1:B){
  
  print(b)
  
  paramvec <- c(beta1_use, beta2_use, alpha1_use, alpha2_use, tau_use, kappa1_use, kappa2_use)
  names(paramvec) <- c("beta01", "beta11", "beta21", "beta31", "beta02", "beta12", "beta22", "beta32", "alpha01", "alpha11", "alpha02", "alpha12", "tau0", "tau1", "tau2", "tau3", "kappa01", "kappa11", "kappa02", "kappa12")
  

  ## generate covariates ##
  x2 <- rbinom(n, 1, 0.59)
  x1 <- rnbinom(n, mu=11-3.5*x2, size=4.8 - 0.3*x2)
  x3 <- x1*x2
  
  ## calculate parameter values ##
  mu1.true <- (cbind(1,x1,x2,x3) %*% beta1_use) %>% exp %>% c
  mu2.true <- (cbind(1,x1,x2,x3) %*% beta2_use) %>% exp %>% c
  sig1.true <- (cbind(1,x2) %*% alpha1_use) %>% exp %>% c
  sig2.true <- (cbind(1,x2) %*% alpha2_use) %>% exp %>% c
  rho.true <- (cbind(1,x1,x2,x3) %*% tau_use) %>% tanh %>% c
  p1.true <- (cbind(1,x2) %*% kappa1_use) %>% sigmoid %>% c
  p2.true <- (cbind(1,x2) %*% kappa2_use) %>% sigmoid %>% c

  ## simulate data from nbnb Gaussian copula model ##
  y1y2 <- copnb.sim(p1.true, p2.true, mu1.true, mu2.true, sig1.true, sig2.true, rho.true, "Gaussian")
  ourdat <- cbind(y1y2, x1, x2, x3) %>% as.data.frame
  
  ## run the fitting function ##
  out_obj <- FAST.CoExpress.nb(formula = eqlist,
                               data = ourdat,
                               copula = "Gaussian")
  
  ## find out if true values were contained in their CI ##
  ests <- out_obj$output.mat[,"coefficients"]
  ses <- out_obj$output.mat[,"std.errors"]
  varcovmat <- out_obj$varcovmat
  
  names(ests) <- names(ses) <- dimnames(varcovmat)[[1]] <- dimnames(varcovmat)[[2]] <- names(paramvec)
  
  # kappa11, alpha11 corresponds to x2 = 0
  ests["kappa11"] <- sum(ests[c("kappa01", "kappa11")])
  ests["alpha11"] <- sum(ests[c("alpha01", "alpha11")])
  paramvec["kappa11"] <- sum(paramvec[c("kappa01", "kappa11")])
  paramvec["alpha11"] <- sum(paramvec[c("alpha01", "alpha11")])
  
  
  ses["kappa11"] <- sqrt(t(c(1,1))%*%varcovmat[c("kappa01","kappa11"),c("kappa01","kappa11")]%*%c(1,1))
  ses["alpha11"] <- sqrt(t(c(1,1))%*%varcovmat[c("alpha01","alpha11"),c("alpha01","alpha11")]%*%c(1,1))
  
  
  # kappa12, alpha12 corresponds to x1 = 0
  ests["kappa12"] <- sum(ests[c("kappa02", "kappa12")])
  ests["alpha12"] <- sum(ests[c("alpha02", "alpha12")])
  paramvec["kappa12"] <- sum(paramvec[c("kappa02", "kappa12")])
  paramvec["alpha12"] <- sum(paramvec[c("alpha02", "alpha12")])
  
  ses["kappa12"] <- sqrt(t(c(1,1))%*%varcovmat[c("kappa02","kappa12"),c("kappa02","kappa12")]%*%c(1,1))
  ses["alpha12"] <- sqrt(t(c(1,1))%*%varcovmat[c("alpha02","alpha12"),c("alpha02","alpha12")]%*%c(1,1))
  
  
  low_est_high <- cbind(ests, ests, ests) + 1.96*cbind(-ses, 0, ses)
  
  low_est_high[c("kappa01","kappa11","kappa02","kappa12"),] <- sigmoid(low_est_high[c("kappa01","kappa11","kappa02","kappa12"),])
  paramvec[c("kappa01","kappa11","kappa02","kappa12")] <- sigmoid(paramvec[c("kappa01","kappa11","kappa02","kappa12")])
  
  low_est_high[c("alpha01","alpha11","alpha02","alpha12"),] <- exp(low_est_high[c("alpha01","alpha11","alpha02","alpha12"),])
  paramvec[c("alpha01","alpha11","alpha02","alpha12")] <- exp(paramvec[c("alpha01","alpha11","alpha02","alpha12")])
  
  
  cilen <- cilen + (low_est_high[,3] - low_est_high[,1])/B
  
  mse <- mse + (low_est_high[,2] - paramvec)^2/B 
  
  mbe <- mbe + abs(low_est_high[,2] - paramvec)/B 
  
  covrg <- covrg + 1*((low_est_high[,1] < paramvec) & (low_est_high[,3] > paramvec))/B
  
}

cbind(mse, mbe, cilen, covrg)[-which(names(paramvec) %in% c("kappa01", "kappa11")),]


