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
\details{
  If the p value for either level of the x dimension is
  significant, we are justified in rejecting the null
  hypothesis of perceptual separability on the x dimension.
  Similarly for the y dimension.
}
\examples{
data(thomasA)
mriTest(thomasA)
}

