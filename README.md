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

where $f_1, f_2$ are the marginal distributions of $Y_1, Y_2$ and $g$ is the distribution function from a bivariate copula with association parameter $\rho$. Parameters can be made dependent on covariates in the following way 

$$
\begin{align}
 \eta_1(\mu_{1}) &= \mathbf{x_{s_3}}^{\top} \mathbf{\beta_1} \\
 \eta_2(\mu_{2}) &= \mathbf{x_{s_4}}^{\top} \mathbf{\beta_2} \\
 \eta_3(\sigma_{1}) &= \mathbf{x_{s_5}}^{\top} \mathbf{\alpha_1} \\
 \eta_4(\sigma_{2}) &= \mathbf{x_{s_6}}^{\top} \mathbf{\alpha_2} \\
 \eta_5(\rho) &= \mathbf{x_{s_7}}^{\top} \mathbf{\tau_1} \\
 \text{logit}(p_{1}) &= \mathbf{x_{s_1}}^{\top} \mathbf{\kappa_1} \\
 \text{logit}(p_{2}) &=  \mathbf{x_{s_2}}^{\top} \mathbf{\kappa_2} \\
 \end{align}
$$

where $\eta_k$ is a link function ensuring the $k^{\text{th}}$ parameter does not stray outside its domain and $S_k$ is a subset of $\{1,\dots, s\}$ reflecting the elements of $\boldsymbol{x}$ which the $k^{\text{th}}$ parameter is to be made dependent on. 


## Usage

The fitting function for the FAST-CoExpress model with negative binomial marginals is called `FAST.CoExpress.nb` and is contained in the `FAST_CoExpress_NB.R` folder of this repository. For the negative binomial case, we have  $\eta_k = \log$ for $k=1,\dots,4$. 

The function works as follows:

```{r}
FAST.CoExpress.nb(formula,
                  copula,
                  data)
```
Parameters:
* `formula`: A list of five formulas specifying the covariate-dependence of the different parameters.
* `copula`: A string specifying one of the following copulas: "Gaussian", "Frank", "Gumbel", "Joe", "Clayton".
* `data`: A data.frame whose column names correspond to the variables referenced in the formula list.

where the 5 elements of the `formula` list are:
* `eq1`: Formula object specifying the covariate-dependence of $\mu_1$.
* `eq2`: Formula object specifying the covariate-dependence of $\mu_2$.
* `eq3`: Formula object specifying the covariate-dependence of $\sigma_1$ and $p_1$.
* `eq4`: Formula object specifying the covariate-dependence of $\sigma_2$ and $p_2$.
* `eq5`: Formula object specifying the covariate-dependence of $\rho$.

The output is a list containing coefficient estimates, standard errors, and p-values, along with the value of the log-likelihood of the fitted model and the variance-covariance matrix of the coefficients.

## Simulation Studies

We conducted simulation studies on the coverage, power, and robustness of the FAST-CoExpress method with negative binomial marginals and a Gaussian copula. These studies can be found in the `simulations` folder of this repository.







