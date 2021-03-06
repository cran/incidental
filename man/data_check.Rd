% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_processing.R
\name{data_check}
\alias{data_check}
\title{Input data check}
\usage{
data_check(reported, delay_dist)
}
\arguments{
\item{reported}{An integer vector of reported cases.}

\item{delay_dist}{A positive vector that sums to one, which describes the delay distribution.}
}
\description{
Check input data for: \itemize{
  \item minimum length of reported
  \item integer for reported
  \item positivity for delay_dist and reported
  \item sums to 1 for delay_dist}
Throw an error if any conditions are violated.
}
