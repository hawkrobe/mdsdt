% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/grt_base.R
\name{anova.grt}
\alias{anova.grt}
\title{Compare nested GRT models}
\usage{
\method{anova}{grt}(object, ...)
}
\arguments{
\item{object}{A fitted GRT model returned by fit.grt}

\item{...}{A larger GRT model, with model1 nested inside}
}
\description{
Conducts a likelihood-ratio G-test on nested GRT models. Currently only accepts pairs of nested models, not arbitrary sequences.
}

