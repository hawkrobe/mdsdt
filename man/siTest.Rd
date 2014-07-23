\name{siTest}
\alias{siTest}
\title{Tests sampling independence for each stimulus response distribution}
\usage{
  siTest(x)
}
\arguments{
  \item{x}{four-by-four confusion matrix}
}
\value{
  data frame containing z-scores and p-values for all four
  tests
}
\description{
  Tests sampling independence for each stimulus response
  distribution
}
\details{
  If p value is sufficiently low, we're justified in
  rejecting the null hypothesis of sampling independence.
}
\examples{
data(thomasA)
siTest(thomasA)
}

