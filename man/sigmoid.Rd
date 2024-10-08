% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{sigmoid}
\alias{sigmoid}
\title{Sigmoid function}
\usage{
sigmoid(x, r = 1, s = 0)
}
\arguments{
\item{x}{numeric scalar/vector/matrix}

\item{r}{inverse scale parameter of the sigmoid function}

\item{s}{midpoint parameter of the sigmoid function}
}
\value{
A numeric scalar/vector/matrix of the same dimensions containing
  the transformed values.
}
\description{
Sigmoid function
}
\examples{

p <- sigmoid(0.5)

}
