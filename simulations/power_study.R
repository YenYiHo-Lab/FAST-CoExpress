source("C:/Users/andyb/Downloads/FAST_CoExpress_functions.R")
library(MASS)
library(truncnorm)
library(rCoCoA)
library(purrr)
library(geepack)
library(LiquidAssociation)
library(dplyr)
library(tidyr)


n <- 2500

## true coefficient values ##
beta1_use <- c(3.16)
beta2_use <- c(2.44)
alpha1_use <- c(-0.23)
alpha2_use <- c(-0.50)
tau_use <- c(-1.62, 0)
kappa1_use <- c(-0.93)
kappa2_use <- c(-1.22)

paramvec <- c(beta1_use, beta2_use, alpha1_use, alpha2_use, tau_use, kappa1_use, kappa2_use)
names(paramvec) <- c("beta01", "beta02", "alpha01", "alpha02", "tau0", "tau1", "kappa01", "kappa02")
which_tau1 <- which(names(paramvec) == "tau1")


## specify how the different parameters depend on covariates ##
eq1 <- y1 ~ 1
eq2 <- y2 ~ 1
eq3 <- ~ 1
eq4 <- ~ 1
eq5 <- ~ x1
eqlist <- list(eq1, eq2, eq3, eq4, eq5)



## Run B iterations of power study ##
B <- 1 # 1000

tau1vec <- seq(0, 0.175, length=3)

powmat <- matrix(0, ncol=4, nrow=length(tau1vec))
colnames(powmat) <- c("FAST.CoExpress", "CoCoA", "CNM", "LA")


for (j in 1:length(tau1vec)){

  tau_use[2] <- tau1vec[j]
  
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
    
    ## fit FAST-CoExpress ##
    out_obj <- FAST.CoExpress.nb(formula = eqlist,
                                 data = ourdat,
                                 copula = "Gaussian")
    
    outmat <- out_obj$output.mat
    pval <- outmat[which_tau1, "p.values"]
    powmat[j, "FAST.CoExpress"] <- powmat[j, "FAST.CoExpress"]  + (pval < 0.05)/B

    
    ## We must normalize y1, y2 for use in LA, CNM.full, and CoCoA ##
    marg1fit <- glm.nb(y1 ~ 1, data=ourdat)

    mu1_est <- c(exp(rep(1,n)*marg1fit$coefficients))
    size1_est <- rep(marg1fit$theta, n)

    marg2fit <- glm.nb(y2 ~ 1, data=ourdat)

    mu2_est <- c(exp(rep(1,n)*marg2fit$coefficients))
    size2_est <- rep(marg2fit$theta, n)

    ourdat$y1.norm <- apply(cbind(ourdat$y1, mu1_est, size1_est), 1,
                       function(a) ifelse(a[1]==0, rtruncnorm(n=1, a=-Inf, b=qnorm(p=pnbinom(q=a[1], mu=a[2], size=a[3]))),
                                          runif(n=1, min=qnorm(p=pnbinom(q=a[1]-1, mu=a[2], size=a[3])), max=qnorm(p=pnbinom(q=a[1], mu=a[2], size=a[3])))))

    ourdat$y2.norm <- apply(cbind(ourdat$y2, mu2_est, size2_est), 1,
                       function(a) ifelse(a[1]==0, rtruncnorm(n=1, a=-Inf, b=qnorm(p=pnbinom(q=a[1], mu=a[2], size=a[3]))),
                                          runif(n=1, min=qnorm(p=pnbinom(q=a[1]-1, mu=a[2], size=a[3])), max=qnorm(p=pnbinom(q=a[1], mu=a[2], size=a[3])))))


    ## fit CoCoA using REML ##
    cocoa_df = data.frame(X = ourdat$y1.norm, Y = ourdat$y2.norm, Z=ourdat$x1)
    cocoa_res <- get_params_reml(cocoa_df)$out2
    
    powmat[j, "CoCoA"] <- powmat[j, "CoCoA"] + cocoa_res[nrow(cocoa_res),"wald"]/B
    



    # ## fit CNM using GEE ##
    ourdat$rownum <- 1:nrow(ourdat)

    ourdat_long <- ourdat %>%
      pivot_longer(cols = c("y1.norm", "y2.norm"),
                   names_to = "Margin",
                   values_to = "Value") %>%
      mutate(Margin = ifelse(Margin == "y1.norm", "Margin 1", "Margin 2"))

      zcor_matrix <- model.matrix(~ x1, data = ourdat)

      cnm_res <- summary(geese(Value ~ Margin,
                         id = rownum,
                         data = ourdat_long,
                         family = gaussian,
                         corstr = "exchangeable",
                         sformula = ~ Margin,
                         zcor = zcor_matrix,
                         cor.link = "fisherz"))$correlation
      
      powmat[j, "CNM"] <- powmat[j, "CNM"] + 1*(cnm_res[nrow(cnm_res), "p"]<0.05)/B
      


    ## fit LA ##
    LA_res <- getsLA(as.matrix(cocoa_df))
    
    powmat[j, "LA"] <- powmat[j, "LA"] + 1*(LA_res[2]<0.05)/B

    
  }
  
}

powmat

