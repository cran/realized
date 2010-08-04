\name{rMarginal}
\alias{rMarginal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maginal Contribution to Realized Estimate }
\description{
Plots the marginal contribution to the realized estimate. 
}
\usage{
rMarginal(x, y = NULL, period, align.period = 1, plotit = FALSE, cts = TRUE,makeReturns = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{  RealizedObject or TimeSeries for S+}
  \item{y}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{align.period}{ Align the returns to this period first}
  \item{plotit}{ T for plot }
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
}
\details{
Plots the marginal contribution to the realized estimate.  This is a good tool to determine what obersations are adding (possibly subtracting for covariance) to the estimate.  For version 0.7 this is only implemented for the naive estimators, in 1.0 it will be implemented generically.
}
\value{
Marginal contribution vector if plotit = F
}
\references{
S. W. Payseur.  A One Day Comparison of Realized Variance and Covariance Estimators. \emph{Working Paper: University 
of Washington}, 2007
}
\author{ Scott Payseur <spayseur@u.washington.edu>}
\seealso{ \code{\link{rAccumulation}}, \code{\link{rRealizedVariance}} }

\examples{
data(msft.real.cts)
data(dates.example)

par(mfrow=c(2,1))
plot(rCumSum(msft.real.cts[[1]], period=10, align.period=60), xlab="", ylab="Cumulative Ruturns", main="MSFT", sub=dates.example[[1]], type="p")
barplot(rMarginal(msft.real.cts[[1]], period=10, align.period=60)$y, main="Marginal Contribution Plot") 



}
\keyword{methods}
