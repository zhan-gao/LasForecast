# LasForecast

This package develops a framework for economic forecasting in a data-rich environment with a particular emphasis on linear predictive regression. The goal is to automate the processes of parameter tuning, rolling window forecasting and backtesting, and visualization in a unified framework.

The package covers the following methods.

|                                       | Functions                                  | Reference                                                    |
| ------------------------------------- | ------------------------------------------ | ------------------------------------------------------------ |
| Lasso                                 | `train_lasso()`, `glmnet::glmnet()`        | Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. *Journal of the Royal Statistical Society: Series B (Methodological)*, *58*(1), 267-288. |
| Adaptive Lasso (Alasso)               | `train_lasso(ada = TRUE)`, `adalasso()`    | Zou, H. (2006). The adaptive lasso and its oracle properties. *Journal of the American statistical association*, *101*(476), 1418-1429<br /> Medeiros, M. C., & Mendes, E. F. (2016). ℓ1-regularization of high-dimensional time-series models with non-Gaussian and heteroskedastic errors. *Journal of Econometrics*, *191*(1), 255-271. |
| Twin Adaptive Lasso (TAlasso)         | `train_replasso(ada = TRUE)`, `replasso()` | Lee, J. H., Shi, Z., & Gao, Z. (2022). On LASSO for predictive regression. *Journal of Econometrics*, *229*(2), 322-349. |
| Best Subset Selection                 | `bss.bic()`, `bss()`                       | Bertsimas, D., King, A., & Mazumder, R. (2016). Best subset selection via a modern optimization lens. *The annals of statistics*, *44*(2), 813-852. |
| Complete Subset Regression (CSR)      | `csr.bic()`, `csr()`                       | Elliott, G., Gargano, A., & Timmermann, A. (2013). Complete subset regressions. *Journal of Econometrics*, *177*(2), 357-373.<br />Elliott, G., Gargano, A., & Timmermann, A. (2015). Complete subset regressions with large-dimensional sets of predictors. *Journal of Economic Dynamics and Control*, *54*, 86-110. |
| L2-relaxation Forecast Combination | `l2relax()`, `train_l2_relax()`            | [Shi, Z., Su, L., & Xie, T. (2024)](https://direct.mit.edu/rest/article-abstract/doi/10.1162/rest_a_01261/113783/2-Relaxation-With-Applications-to-Forecast?redirectedFrom=fulltext). : With Applications to Forecast Combination and Portfolio Analysis, *The Review of Economics and Statistics* |

The package exploits the advantage of well-established packages like `glmnet` and model training framework `caret`. 

We can run backtesting and compare the forecasting performance among different methods based on rolling windows and forecasting horizons using the `roll_predict` function. For example, we can replicate the empirical results in [Lee, Shi and Gao (2022)](https://www.sciencedirect.com/science/article/pii/S030440762100049X) as in [this script](https://github.com/zhan-gao/Alasso_Predictive_Regression/blob/master/Welch_Goyal/master_rolling.R).

### Installation
```r
install.packages("devtools")
devtools::install_github("zhan-gao/LasForecast")
```

### Usage

