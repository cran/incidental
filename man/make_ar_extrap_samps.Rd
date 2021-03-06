% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_processing.R
\name{make_ar_extrap_samps}
\alias{make_ar_extrap_samps}
\title{Make AR samples for extrapolation past end point}
\usage{
make_ar_extrap_samps(
  reported,
  num_ar_steps = 10,
  num_ar_samps = 50,
  seed = 1,
  linear_tail = 14,
  extrapolation_prior_precision = 2
)
}
\arguments{
\item{reported}{An integer vector of reported cases.}

\item{num_ar_steps}{An integer number of AR steps after last observation.}

\item{num_ar_samps}{An integer number of AR samples.}

\item{seed}{Seed for RNG.}

\item{linear_tail}{An integer number of days used to fit linear model on tail to be used as a mean
for AR extrapolation.}

\item{extrapolation_prior_precision}{A positive scalar for extrapolation slope shrinkage prior precision.}
}
\value{
A matrix of size (num_ar_samps x n + num_ar_steps)
}
\description{
Make auto-regressive (AR) samples for extrapolation past end point to help with right-censoring problems.
}
