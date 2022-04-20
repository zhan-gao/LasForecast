# LasForecast

This package intends to develop a framework for economic forecasting in a data-rich envrionment with a particular emphasis on Lasso-type linear models. The scope of methods is expected to expand to cover forecast combination methods machine learning technologies. The goal is to automate the processes of parameter tuning, rolling window forecasting and performance visualization in a unified framework.

The package exploits the advantage of well-established packages like `glmnet` and model training framework `caret`.

### Installation
```{r}
install.packages("devtools")
devtools::install_github("chenyang45/BoostedHP")
library("bHP")
```