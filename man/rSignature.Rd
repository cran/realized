\name{rSignature}
\alias{rSignature}

\title{ Signature Plots}
\description{
Creates realized variance, covariance, and correlation plots for one or multiple days for each type of estimator.
}
\usage{
rSignature(range, x, y = NULL, type = "naive", cor = FALSE, args = list(), xscale = 1, iteration.funct = "", iterations = NULL, plotit = FALSE, cts = TRUE, makeReturns = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 
  \item{range}{  x axis of signature plot (inputs for realized estimator)}
 \item{x}{  RealizedObject or TimeSeries for S+}
  \item{y}{  RealizedObject or TimeSeries for S+}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
  \item{type}{ Type of realized estimator to use, a rv. or rc. is appended to this value and that function is called }
  \item{cor}{ T for correlation }
  \item{args}{ List of extra parameters to pass into rv.* or rc.* }
  \item{xscale}{ value to multiply range by for x axis (see below) }
  \item{iteration.funct}{ Name of function to be used for average signature plot }
  \item{iterations}{ Numeric vector to pass into iteration.funct }
  \item{plotit}{ T for plot }
}
\details{
Creates a list that contains an x and y value that is used for plotting a signature plot.  
}
\value{
List containing:
  \item{x }{  x axis of signature plot (range * xscale)}
  \item{y }{  Realized variance or covariance}
  \item{type }{  type of estimator used for calculation}
  \item{cor }{  is this a correlation signture}
  \item{cov }{  it this realized covariance}
}
\references{
T.G. Andersen, T. Bollerslev, F.X. Diebold, and P. Labys. Great realizations. \emph{Risk}, 13:105 108, 2000.

Y. Fang. Volatility modeling and estimation of high-frequency data with gaussian noise. \emph{unpublished doctoral
thesis, MIT, Sloan School of Management}, 1996.

J. E. Griffen and R. C. A. Oomen. Covariance measurement in the presence of non-synchronous trading and
market microstructure noise. \emph{Working Paper}, June 27th, 2006.

S. W. Payseur.  A One Day Comparison of Realized Variance and Covariance Estimators. \emph{Working Paper: University 
of Washington}, 2007
}

\author{ Scott Payseur <spayseur@u.washington.edu>}

\seealso{ \code{\link{rRealizedVariance}} }
\examples{
data(msft.real.cts)
data(ge.real.cts)

# see users manual for more examples

test.cov <- rSignature(1:1200,x=msft.real.cts[[1]], y=ge.real.cts[[1]], xscale=1/60) 
test.rect <- rSignature(1:600,msft.real.cts[[1]], ge.real.cts[[1]],type="kernel",args=list(type="rectangular"), xscale=1/30)
test.mth <- rSignature(1:600,msft.real.cts[[1]], ge.real.cts[[1]],type="kernel",args=list(type="mth"), xscale=1/30)
plot(test.cov, ylab="Realized Covariance", xlab="Minutes", main="GE | MSFT") 
lines(test.rect, col=3, lwd=1) 
lines(test.mth, col=4, lwd=2) 
axis(3, c(0,(1:5)*4), c("Lags:",as.character((1:5)*120)))
legend(13,.00015,c("Rectangular", "Mod TH"), lwd=c(1,2), col=c(3,4)) 

}


\keyword{methods}