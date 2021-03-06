% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spline_model.R
\name{compute_log_incidence}
\alias{compute_log_incidence}
\title{Compute log likelihood of incidence model}
\usage{
compute_log_incidence(beta, Q, Tobs)
}
\arguments{
\item{beta}{parameter vector of num_params}

\item{Q}{spline basis matrix, of size Tmod x num_params}

\item{Tobs}{maximum observed time point}
}
\value{
I Tobs-length vector that models log incidence curve
}
\description{
This function computes log likelihood of incidence model given parameters and observations.
}
