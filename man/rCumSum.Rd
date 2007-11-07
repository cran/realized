\name{rCumSum}
\alias{rCumSum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot cummulative returns }
\description{
Plots cummulative returns at a certain alignment given a return series.
}
\usage{
rCumSum(x, period = 1, align.period = 1, plotit = FALSE, type = "l", cts = TRUE, makeReturns = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{align.period}{ Align the returns to this period first}
  \item{plotit}{ T for plot }
    \item{type}{ Line or points}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}

  
}
\value{
Cummulative return vector if plotit = F
}

\author{ Scott Payseur <spayseur@u.washington.edu>}

\examples{

data(msft.real.cts)
data(dates.example)

cumm <- list() 
cumm[[1]] <- rCumSum(msft.real.cts[[1]], period=1, align.period=60) 
cumm[[2]] <- rCumSum(msft.real.cts[[1]], period=10, align.period=60) 
cumm[[3]] <- rCumSum(msft.real.cts[[1]], period=20, align.period=60) 
cumm[[4]] <- rCumSum(msft.real.cts[[1]], period=30, align.period=60) 
plot(cumm[[1]], xlab="", ylab="Cumulative Ruturns", main="MSFT", sub=dates.example[[1]], type="p", col=16, lwd=2) 
lines(cumm[[2]], col=2, lwd=2) 
lines(cumm[[3]], col=3, lwd=2) 
lines(cumm[[4]], col=4, lwd=2)
}
\keyword{methods}