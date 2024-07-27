<h1 align="center" style="font-weight: bold;">FAST-CoExpress</h1>

<p align="center">
  <align="center"><ins>F</ins>ast, <ins>A</ins>daptive, <ins>S</ins>calable <ins>T</ins>ool for differential <ins>Co</ins>-<ins>Express</ins>ion analyses
</p>

## The Model

Let $Y_1, Y_2$ be random variables representing the expression levels of two genes in a cell, and let $\mathbf{x} = \[x_1, \dots, x_s\]^{\top}$ be a vector containing covariates from that cell. The FAST-CoExpress framework models the joint distribution of $Y_1, Y_2$ conditional on $\mathbf{x}$ using a zero-inflated bivariate copula model where all parameters can be made dependent on subsets of covariates. 

The joint distribution of $Y_1, Y_2$ is given by 

$$
\begin{aligned}
&f_{m_1, m_2}(y_{1}, y_{2} ; \  p_{1}, p_{2}, \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) = \\
&\begin{cases} 
    (1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) & \text{if } y_{1}\neq 0 \text{, } y_{2}\neq 0 \\
    (1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) + (1-p_{1})p_{2}f_{m_1}(y_{1}; \ \mu_1, \sigma_1)  & \text{if } y_{1}\neq 0 \text{, } y_{2}= 0 \\
    (1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) + p_{1}(1-p_{2})f_{m_2}(y_{2}; \ \mu_2, \sigma_2) & \text{if } y_{1}= 0 \text{, } y_{2}\neq 0 \\
    (1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) + p_{1}p_{2} & \text{if } y_{1}= 0 \text{, } y_{2}=  0   
  \end{cases}
\end{aligned}
$$

where $f_1, f_2$ are the marginal distributions of $Y_1, Y_2$, $g$ is the distribution function from a bivariate copula with association parameter $\rho$, and $p_1$, $p_2$ are parameters representing the probability of an observation from their respective marginal being zeroed out by a drop-out event. Parameters can be made dependent on covariates in the following way 

$$
\begin{align}
 \eta_1(\mu_{1}) &= \mathbf{x_{s_1}}^{\top} \mathbf{\beta_1} \\
 \eta_2(\mu_{2}) &= \mathbf{x_{s_2}}^{\top} \mathbf{\beta_2} \\
 \eta_3(\sigma_{1}) &= \mathbf{x_{s_3}}^{\top} \mathbf{\alpha_1} \\
 \eta_4(\sigma_{2}) &= \mathbf{x_{s_4}}^{\top} \mathbf{\alpha_2} \\
 \eta_5(\rho) &= \mathbf{x_{s_5}}^{\top} \mathbf{\tau} \\
 \text{logit}(p_{1}) &= \mathbf{x_{s_6}}^{\top} \mathbf{\kappa_1} \\
 \text{logit}(p_{2}) &=  \mathbf{x_{s_7}}^{\top} \mathbf{\kappa_2} \\
 \end{align}
$$

where $\eta_k$ is a link function ensuring the $k^{\text{th}}$ parameter does not stray outside its domain and $S_k$ is a subset of $\{1,\dots, s\}$ reflecting the elements of $\boldsymbol{x}$ which the $k^{\text{th}}$ parameter is to be made dependent on. 


## Fitting Function

The fitting function for the FAST-CoExpress model with negative binomial marginals is called `FAST.CoExpress.nb`. For the negative binomial case, we have  $\eta_k = \log$ for $k=1,\dots,4$. 

The function works as follows:

```{r}
source("/path/to/FAST_CoExpress_functions.R")

FAST.CoExpress.nb(formula, copula, data)
```
Parameters:
* `formula`: A list of seven formula objects specifying the covariate-dependence of the different parameters.
  * `[[1]]`: Covariate-dependence of $\mu_1$
  * `[[2]]`: Covariate-dependence of $\mu_2$
  * `[[3]]`: Covariate-dependence of $\sigma_1$
  * `[[4]]`: Covariate-dependence of $\sigma_2$
  * `[[5]]`: Covariate-dependence of $\rho$
  * `[[6]]`: Covariate-dependence of $p_1$
  * `[[7]]`: Covariate-dependence of $p_2$
* `copula`: A string specifying one of the following copulas: "Gaussian", "Frank", "Gumbel", "Joe", "Clayton".
* `data`: A data.frame whose column names correspond to the variables referenced in the formula list.



The output is a list containing coefficient estimates, standard errors, and p-values, along with the value of the log-likelihood of the fitted model and the variance-covariance matrix of the coefficients.

## Simulation Studies

We conducted simulation studies on the coverage, power, and robustness of the FAST-CoExpress method with negative binomial marginals and a Gaussian copula. These studies can be found in the `simulations` folder of this repository.

## Real Data Analysis

We performed differential co-expression analysis on breast cancer data from gene expression omnibus GSE176078. Code for the analysis can be found in the simulations folder of this repository, and the data used can be found at the following link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5354529.

## Example Usage

```{r}
source("/path/to/FAST_CoExpress_functions.R")
library(dplyr)

n <- 5000

## true coefficient values ##
beta1_use <- c(2, 0.1, -0.2, -0.05)
beta2_use <- c(1, -0.05, 0.3, -0.1)
alpha1_use <- c(-0.3, 0.6)
alpha2_use <- c(-1.1, 0.7)
tau_use <- c(-0.4, 0.2, 1, -0.25)
kappa1_use <- c(-2, 1)
kappa2_use <- c(-1, -2)

## generate covariates ##
x1 <- runif(n, 0, 5)
x2 <- rbinom(n=n, size=1, prob=0.59)
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

## specify how the different parameters depend on covariates ##
eq1 <- y1 ~ x1 + x2 + x3
eq2 <- y2 ~ x1 + x2 + x3
eq3 <- ~ x2
eq4 <- ~ x2
eq5 <- ~ x1 + x2 + x3
eq6 <- ~ x2
eq7 <- ~ x2
eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)

## run the fitting function ##
out_obj <- FAST.CoExpress.nb(formula = eqlist,
                             data = ourdat,
                             copula = "Gaussian")

round(out_obj$output.mat, 4)
```



