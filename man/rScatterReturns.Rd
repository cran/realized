\name{rScatterReturns}
\alias{rScatterReturns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Scatterplot of aligned returns }
\description{
Creates a scatterplot of cross returns.
}
\usage{
rScatterReturns(x,y, period, align.period=1,numbers=FALSE,xlim= NULL, ylim=NULL, plotit=TRUE, pch=NULL, cts=TRUE, makeReturns=FALSE, scale.size=0, col.change=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{  RealizedObject or TimeSeries for S+}
  \item{y}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
  \item{plotit}{ T for plot }
  \item{numbers}{ T for count }
  \item{pch}{type of point}
  \item{ylim}{ylimit}
  \item{xlim}{xlimit}
  \item{scale.size}{.}
  \item{col.change}{.}
    \item{...}{ ...}
}
\details{
Scatterplot of returns.
}

\references{
S. W. Payseur.  A One Day Comparison of Realized Variance and Covariance Estimators. \emph{Working Paper: University 
of Washington}, 2007
}
\author{ Scott Payseur <spayseur@u.washington.edu>}


\examples{


data(msft.real.cts)
data(ge.real.cts)
par(mfrow=c(2,1))
rScatterReturns(msft.real.cts[[1]],y=ge.real.cts[[1]], period=1, align.period=20,ylab="GE",xlab="MSFT",numbers=FALSE) 
rScatterReturns(msft.real.cts[[1]],y=ge.real.cts[[1]], period=1, align.period=20,ylab="GE",xlab="MSFT",numbers=TRUE) 


}
\keyword{methods}