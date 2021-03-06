% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{incidence_to_df}
\alias{incidence_to_df}
\title{Export incidence model to data frame}
\usage{
incidence_to_df(x, times = NULL, low_quantile = 0.05, high_quantile = 0.95)
}
\arguments{
\item{x}{An "incidence_spline_model" output from \code{\link{fit_incidence}}.}

\item{times}{An optional vector of time indices.}

\item{low_quantile}{A scalar that specifies the low quantile value for the output CI.}

\item{high_quantile}{A scalar that specifies the high quantile value for the output CI.}
}
\value{
A data frame with the following entries: \itemize{
        \item Time -- a time index; if `ts` is `NULL` it is the observation number;
        \item Reported -- the value of `reported`;
        \item Ihat -- MAP incidence curve estimate;
        \item Chat -- expected cases given MAP incidence curve estimate;
        \item LowCI  -- lower pointwise credible interval bands around the incidence curve; and
        \item HighCI  -- higher pointwise credible interval bands around the incidence curve.}
}
\description{
Export the output of \code{\link{fit_incidence}} to a data frame with an optional
addition of a time index.
}
\examples{
\donttest{
indiana_model <- fit_incidence(
                  reported = spanish_flu$Indiana, 
                  delay_dist = spanish_flu_delay_dist$proportion)
 indiana_df <- incidence_to_df(indiana_model, times = spanish_flu$Date)
}
}
