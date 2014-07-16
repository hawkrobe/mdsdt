\name{fit.grt}
\alias{fit.grt}
\title{Fit full Gaussian GRT model}
\usage{
  fit.grt(freq, PS_x = FALSE, PS_y = FALSE, PI = "none")
}
\arguments{
  \item{freq}{Can be entered in two ways: 1) a 4x4
  confusion matrix containing counts, with each row
  corresponding to a stimulus and each column corresponding
  to a response. row/col order must be aa, ab, ba, bb. 2) A
  three-way 'xtabs' table with the stimuli as the third
  index and the NxN possible responses as the first two
  indices.}

  \item{PS_x}{if TRUE, will fit model with assumption of
  perceptual separability on the x dimension}

  \item{PS_y}{if TRUE, will fit model with assumption of
  perceptual separability on the y dimension}

  \item{PI}{'none' by default, imposing no restrictions and
  fitting different correlations for all distributions. If
  'same_rho', will constrain all distributions to have same
  correlation parameter. If 'all', will constain all
  distribution to have 0 correlation.}
}
\value{
  An S3 \code{grt} object
}
\description{
  Use Newton-Raphson gradient descent to fit the mean and
  covariance of a bivariate Gaussian distribution for each
  stimulus class, subject to given constraints. Standard
  case uses confusion matrix from a 2x2 full-report
  identification experiment, but will also work in designs
  with N levels of confidence associated with each
  dimension (e.g. in Wickens, 1992).
}
\examples{
# Fit unconstrained model
data(thomas01);
grt_obj <- grt.fit(observerB);

# Use standard S3 generics to examine
print(grt_obj);
summary(grt_obj);
plot(grt_obj);

# Fit model with assumption of perceptual separability on both dimensions
grt_obj_PS <- grt.fit(observerB, PS_x = TRUE, PS_y = TRUE);
summary(grt_obj_PS);
plot(grt_obj_PS);

# Compare models
GOF(grt_obj, teststat = 'AIC');
GOF(grt_obj_PS, teststat = 'AIC');
}

