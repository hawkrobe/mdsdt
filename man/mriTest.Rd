\name{mriTest}
\alias{mriTest}
\title{Tests marginal response invariance at both levels on each dimension}
\source{
  Ashby, F. G., & Townsend, J. T. (1986). Varieties of
  perceptual independence. Psychological review, 93(2),
  154.

  Thomas, R. D. (2001).Perceptual interactions of facial
  dimensions in speeded classification and identification.
  Perception \& Psychophysics, 63(4), 625--650.

  Silbert, N. H., & Thomas, R. D. (2013). Decisional
  separability, model identification, and statistical
  inference in the general recognition theory framework.
  Psychonomic bulletin & review, 20(1), 1-20.
}
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

  The estimator is derived in a footnote of Thomas (2001).
}
\examples{
data(thomasA)
mriTest(thomasA)
}

