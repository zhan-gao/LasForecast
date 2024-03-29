% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bss.R
\name{bss}
\alias{bss}
\title{Implement best subset selection via Mixed Integer Optimization (MIO)
for given k}
\usage{
bss(
  y,
  X,
  k,
  intercept = TRUE,
  b0 = NULL,
  tau = 2,
  tol = 1e-04,
  MaxIter = 10000,
  polish = TRUE,
  time.limit = 1200
)
}
\description{
Implement best subset selection via Mixed Integer Optimization (MIO)
for given k
}
