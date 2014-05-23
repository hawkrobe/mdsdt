\name{gaussian_fit}
\alias{gaussian_fit}
\title{Uses Newton-Raphson gradient descent to fit a Gaussian model}
\usage{
  gaussian_fit(freq)
}
\arguments{
  \item{freq}{4x4 confusion matrix containing counts.
  (Assumes row/col order: aa, ab, ba, bb).}
}
\value{
  List containing the following elements:
  \item{parameters}{The fitted means and correlations for
  all 4 response distributions} \item{prob}{The predicted
  response probabilities according to the fit}
  \item{info_mat}{The information matrix (inverse of Fisher
  Information matrix) for all parameters} \item{nll}{The
  negative log likelihood} \item{aic}{Akaike information
  criterion, for model comparison} \item{bic}{Bayesian
  information criterion, for model comparison}
  \item{icomp}{Bozdogan's information complexity (ICOMP),
  for model comparison}
}
\description{
  Uses Newton-Raphson gradient descent to fit a Gaussian
  model
}
\examples{
gaussian_fit(observerB)
}

