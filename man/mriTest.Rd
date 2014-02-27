\name{mriTest}
\alias{mriTest}
\title{Tests marginal response invariance at both levels on each dimension}
\usage{
  mriTest(x)
}
\arguments{
  \item{x}{four-by-four confusion matrix}
}
\value{
  data frame containing z-scores and p-values for all four
  tests
}
\description{
  Tests marginal response invariance at both levels on each
  dimension
}
\examples{
mri_test(observerA)
}

