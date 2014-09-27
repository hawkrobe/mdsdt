\name{riTest}
\alias{riTest}
\title{Tests report independence for each stimulus response distribution}
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
  riTest(x)
}
\arguments{
  \item{x}{four-by-four confusion matrix}
}
\value{
  data frame containing z-scores and p-values for all four
  tests
}
\description{
  Tests report independence for each stimulus response
  distribution
}
\details{
  If p value is sufficiently low, we're justified in
  rejecting the null hypothesis of sampling within that
  factor. p values come from a chi-squared test on the
  confusion matrix, as explaned in a footnote of Thomas
  (2001).
}
\examples{
data(thomasA)
riTest(thomasA)
}

