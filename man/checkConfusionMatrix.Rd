\name{checkConfusionMatrix}
\alias{checkConfusionMatrix}
\title{Checks whether the input is a valid confusion matrix}
\usage{
  checkConfusionMatrix(x)

  checkConfusionMatrix(x)
}
\arguments{
  \item{x:}{four-by-four confusion matrix of counts of
  probabilities}

  \item{x:}{four-by-four confusion matrix of counts of
  probabilities}
}
\value{
  Boolean indicating whether input is valid

  Boolean indicating whether input is valid
}
\description{
  Checks whether the input is a valid confusion matrix

  Checks whether the input is a valid confusion matrix
}
\examples{
checkConfusionMatrix(matrix(c(3,1,1,4,3,3),3))
checkConfusionMatrix(matrix(c(.5,.2,.2,.2,
                              .2,.5,.1,.2,
                              .1,.2,.5,.1,
                              .2,.2,.3,.5), 4))
checkConfusionMatrix(matrix(c(3,1,1,4,3,3),3))
checkConfusionMatrix(matrix(c(.5,.2,.2,.2,
                              .2,.5,.1,.2,
                              .1,.2,.5,.1,
                              .2,.2,.3,.5), 4))
}

