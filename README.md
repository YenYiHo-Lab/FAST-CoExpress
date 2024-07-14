<h1 align="center" style="font-weight: bold;">FAST-CoExpress</h1>

<p align="center">
  <align="center"><ins>F</ins>ast, <ins>A</ins>daptive, <ins>S</ins>calable <ins>T</ins>ool for differential <ins>Co</ins>-<ins>Express</ins>ion analyses
</p>

## The Model

Let $Y_1, Y_2$ be random variables representing the expression levels of two genes in a cell, and let $\mathbf{x} = \[x_1, \dots, x_s\]^{\top}$ be a vector containing covariates from that cell. The FAST-CoExpress framework models the joint distribution of $Y_1, Y_2$ conditional on $\mathbf{x}$ using a zero-inflated bivariate copula model where all parameters can be made dependent on subsets of covariates. 

The joint distribution of $Y_1, Y_2$ is given by 

$$
\begin{aligned}
&f_{m_1, m_2}(y_{1}, y_{2} ; \  p_{1}, p_{2}, \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) =
\begin{cases} 
    (1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) & \text{if } y_{1}\neq 0 \text{, } y_{2}\neq 0 \\
    (1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) + (1-p_{1})p_{2}f_{m_1}(y_{1}; \ \mu_1, \sigma_1)  & \text{if } y_{1}\neq 0 \text{, } y_{2}= 0 \\
    (1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) + p_{1}(1-p_{2})f_{m_2}(y_{2}; \ \mu_2, \sigma_2) & \text{if } y_{1}= 0 \text{, } y_{2}\neq 0 \\
    (1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) + p_{1}p_{2} & \text{if } y_{1}= 0 \text{, } y_{2}=  0   
  \end{cases}
\end{aligned}
$$

where $f_1, f_2$ are the marginal distributions of $Y_1, Y_2$ and $g$ is a bivariate copula with association parameter $\rho$. Parameters can be made dependent on covariates in the following way 

$$
\begin{align}
 \text{logit}(p_{1}) &= \mathbf{x_{s_1}}^{\top} \mathbf{\kappa_1} \\
 \text{logit}(p_{2}) &=  \mathbf{x_{s_2}}^{\top} \mathbf{\kappa_2} \\
 \log(\mu_{1}) &= \mathbf{x_{s_3}}^{\top} \mathbf{\beta_1} \\
 \log(\mu_{2}) &= \mathbf{x_{s_4}}^{\top} \mathbf{\beta_2} \\
 \log(\sigma_{1}) &= \mathbf{x_{s_5}}^{\top} \mathbf{\alpha_1} \\
 \log(\sigma_{2}) &= \mathbf{x_{s_6}}^{\top} \mathbf{\alpha_2} \\
 \text{atanh}(\rho) &= \mathbf{x_{s_7}}^{\top} \mathbf{\tau_1} \\
 \end{align}
$$

where $S_k$ is a subset of $\{1,\dots, s\}$ reflecting the elements of $\boldsymbol{x}$ which the $k^{\text{th}}$ parameter is to be made dependent on. 


## Usage




```{r}

source("/path/to/FAST_CoExpress_functions.R")

```
Parameters:
* `marginals`: The two marginals. Options are NB, ZINB, GA, ZIGA, Beta, ZIBEta
* `x`: The vector (or matrix) containing the covariate values to be regressed for mean and rho parameters.
* `eta1.true`: The coefficients of the 1st marginal's zero-inflation parameter. 
* `eta2.true`: The coefficients of the 2nd marginal's zero-inflation parameter. 
* `beta1.true`: The coefficients of the 1st marginal's mean parameter. 
* `beta2.true`: The coefficients of the 2nd marginal's mean parameter. 
* `alpha1.true`: The coefficient of the 1st marginal's second parameter. 
* `alpha2.true`: The coefficient of the 2nd marginal's second parameter.
* `tau.true`: The coefficients of the correlation parameter. 
* `w`: A vector (or matrix) containing the covariate values to be regressed for zero-inflation parameters.

This will simulate a 2-column matrix of `NROW(x)` rows of observations from the scdeco.cop model.





