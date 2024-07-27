source("/path/to/FAST_CoExpress_functions.R")
library(dplyr)

n <- 2500

## true coefficient values ##
beta1_use <- c(3.16)
beta2_use <- c(2.44)
alpha1_use <- c(-0.23)
alpha2_use <- c(-0.50)
tau_use <- c(-1.62, 0.21)
kappa1_use <- c(-0.93)
kappa2_use <- c(-1.22)

paramvec <- c(beta1_use, beta2_use, alpha1_use, alpha2_use, tau_use, kappa1_use, kappa2_use)
names(paramvec) <- c("beta01", "beta02", "alpha01", "alpha02", "tau0", "tau1", "kappa01", "kappa02")


## specify how the different parameters depend on covariates ##
eq1 <- y1 ~ 1
eq2 <- y2 ~ 1
eq3 <- ~ 1
eq4 <- ~ 1
eq5 <- ~ x1
eq6 <- ~ 1
eq7 <- ~ 1
eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)



## Run B iterations of coverage study ##
B <- 5 # 10000

mse <- mbe <- cilen <- covrg <- rep(0, length(paramvec))


for (b in 1:B){
  
  print(b)

  ## generate covariates ##
  x1 <- runif(n, 0, 15)
  
  ## calculate parameter values ##
  mu1.true <- (rep(1,n) * beta1_use) %>% exp %>% c
  mu2.true <- (rep(1,n) * beta2_use) %>% exp %>% c
  sig1.true <- (rep(1,n) * alpha1_use) %>% exp %>% c
  sig2.true <- (rep(1,n) * alpha2_use) %>% exp %>% c
  rho.true <- (cbind(1,x1) %*% tau_use) %>% tanh %>% c
  p1.true <- (rep(1,n) * kappa1_use) %>% sigmoid %>% c
  p2.true <- (rep(1,n) * kappa2_use) %>% sigmoid %>% c

  ## simulate data from nbnb Gaussian copula model ##
  y1y2 <- copnb.sim(p1.true, p2.true, mu1.true, mu2.true, sig1.true, sig2.true, rho.true, "Gaussian")
  ourdat <- cbind(y1y2, x1) %>% as.data.frame
  
  ## run the fitting function ##
  out_obj <- FAST.CoExpress.nb(formula = eqlist,
                               data = ourdat,
                               copula = "Gaussian")
  
  ## find out if true values were contained in their CI ##
  ests <- out_obj$output.mat[,"coefficients"]
  ses <- out_obj$output.mat[,"std.errors"]
  varcovmat <- out_obj$varcovmat
  
  names(ests) <- names(ses) <- dimnames(varcovmat)[[1]] <- dimnames(varcovmat)[[2]] <- names(paramvec)
  
  low_est_high <- cbind(ests, ests, ests) + 1.96*cbind(-ses, 0, ses)
  
  covvec <- ((low_est_high[,1] < paramvec) & (low_est_high[,3] > paramvec))
  
  retmat <- cbind(low_est_high[,1], paramvec, low_est_high[,c(2,3)], covvec)
  colnames(retmat) <- c("lower", "true", "est", "upper", "covered")
  
  retmat[c("alpha01", "alpha02"), c("lower","true","est","upper")] <- exp(retmat[c("alpha01", "alpha02"), c("lower","true","est","upper")])
  retmat[c("kappa01", "kappa02"), c("lower","true","est","upper")] <- sigmoid(retmat[c("kappa01", "kappa02"), c("lower","true","est","upper")])
  
  mse <- mse + (retmat[,"est"] - retmat[,"true"])^2/B
  mbe <- mbe + abs(retmat[,"est"] - retmat[,"true"])/B
  cilen <- cilen + (retmat[,"upper"] - retmat[,"lower"])/B
  covrg <- covrg + retmat[, "covered"]/B 
  
  
}

cbind(mse, mbe, cilen, covrg)


