\name{rAccumulation}
\alias{rAccumulation}

\title{ Realized Accumulation Plot }
\description{
Plots the realized estimate as it accumulates over a time interval.  
}
\usage{
rAccumulation(x, period = 1, y = NULL, align.period = 1, plotit = FALSE, cts = TRUE, makeReturns = FALSE)
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
Plots the realized estimate as it accumulates over a time interval.  This is a good tool to determine what obersations are adding (possibly subtracting for covariance) to the estimate.  For version 0.7 this is only implemented for the naive estimators, in 1.0 it will be implemented generically.
}
\value{
Realized accumulation vector if plotit = F
}
\references{ 
S. W. Payseur.  A One Day Comparison of Realized Variance and Covariance Estimators. \emph{Working Paper: University 
of Washington}, 2007
}
\author{ Scott Payseur <spayseur@u.washington.edu>}

\seealso{ \code{\link{rMarginal}}, \code{\link{rRealizedVariance}} }

\examples{
data(msft.real.cts)
data(dates.example)

cumm <- list() 
cumm[[1]] <- rCumSum(msft.real.cts[[1]], period=1, align.period=60) 
cumm[[2]] <- rCumSum(msft.real.cts[[1]], period=10, align.period=60) 
cumm[[3]] <- rCumSum(msft.real.cts[[1]], period=20, align.period=60) 
cumm[[4]] <- rCumSum(msft.real.cts[[1]], period=30, align.period=60) 
accum <- list() 
accum[[1]] <- rAccumulation(msft.real.cts[[1]], period=10, align.period=60) 
accum[[2]] <- rAccumulation(msft.real.cts[[1]], period=20, align.period=60) 
accum[[3]] <- rAccumulation(msft.real.cts[[1]], period=30, align.period=60)

par(mfrow=c(2,1)) 
plot(cumm[[1]], xlab="", ylab="Cumulative Ruturns", main="MSFT", sub=dates.example[[1]], type="p", col=16, lwd=2) 
lines(cumm[[2]], col=2, lwd=2) 
lines(cumm[[3]], col=3, lwd=2) 
lines(cumm[[4]], col=4, lwd=2) 
plot(accum[[1]], xlab="", ylab="Realized Accumulation", type="l",main="MSFT", sub=dates.example[[1]], col=2, lwd=2) 
lines(accum[[2]], col=3, lwd=2) 
lines(accum[[3]], col=4, lwd=2) 

}



\keyword{methods}
