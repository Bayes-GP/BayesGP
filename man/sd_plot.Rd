% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_post_fit.R
\name{sd_plot}
\alias{sd_plot}
\title{Plot the posterior density of a variance parameter in the fitted model}
\usage{
sd_plot(object, component = NULL, h = NULL, theta_logprior = NULL)
}
\arguments{
\item{object}{The fitted object from the function `model_fit`.}

\item{component}{The component of the variance parameter that you want to show. By default this is `NULL`, indicating the family.sd is of interest.}

\item{h}{For PSD, the unit of predictive step to consider, by default is set to `NULL`, indicating the result is using the same `h` as in the model fitting.}

\item{theta_logprior}{The log prior function used on the selected variance parameter. By default is `NULL`, and the current Exponential prior will be used.}
}
\description{
Plot the posterior density of a variance parameter in the fitted model
}
