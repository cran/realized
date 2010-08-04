\name{rRealizedVariance}
\alias{rRealizedVariance}
\alias{rRealisedVariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate realized variance, covariance, or correlation.}
\description{
Calculate realized variance, covariance, correlation, covariance matrix, or correlation matrix.
}
\usage{
rRealizedVariance(x, y = NULL, type = "naive", period = 1, lags = 1, cor = FALSE, rvargs = list(), cts = TRUE, makeReturns = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  RealizedObject or TimeSeries for S+}
  \item{y}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{type}{ Type of realized estimator to use, a rv. or rc. is appended to this value and that function is called }
  \item{lags}{ Number of lags or subgrids for kernel and subsample based estimators }
  \item{cor}{ T for correlation }
  \item{rvargs}{ List of extra parameters to pass into rv.* or rc.* }
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
}
\details{
Calculate realized variance, covariance, correlation, covariance matrix, or correlation matrix.
}
\value{
A single numeric value or a matrix if x is multicolumn matrix.
}
\author{ Scott Payseur <spayseur@u.washington.edu>}


\seealso{ \code{\link{rc.avg}}, \code{\link{rc.kernel}}, \code{\link{rc.naive}}, \code{\link{rc.timescale}}, \code{\link{rv.avg}}, \code{\link{rv.kernel}}, \code{\link{rv.naive}}, \code{\link{rv.timescale}}  }
\examples{
data(msft.real.cts)
data(ge.real.cts)

# Traditional Estimate at highest frequency
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=1)

# Traditional Estimate at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=1, rvargs=list(align.period=60)) 

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=10, rvargs=list(align.period=60)) 

# Bartlett Kernel Estimate with minute aligned data at 20 lags 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=20, rvargs=list(align.period=60, type="Bartlett"))

# Cubic Kernel Estimate with second aligned data at 400 lags 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=400, rvargs=list(type="Cubic")) 

# Lead-Lag with one lag at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=1, rvargs=list(align.period=60)) 
 
# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="avg", period=600) 


# Traditional Estimate at highest frequency 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="naive", period=1)

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="naive", period=10, rvargs=list(align.period=60))

# Lead-Lag with one lag at one minute frequency> 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="kernel", lags=1, rvargs=list(align.period=60))

# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="avg", period=600)
 
 
 
 
# Traditional Estimate at highest frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=1, cor=TRUE) 

# Traditional Estimate at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=1, rvargs=list(align.period=60), cor=TRUE)

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=10, rvargs=list(align.period=60), cor=TRUE) 

# Bartlett Kernel Estimate with minute aligned data at 20 lags 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=20, rvargs=list(align.period=60, type="Bartlett"), cor=TRUE)

# Cubic Kernel Estimate with second aligned data at 400 lags 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=400, rvargs=list(type="Cubic"), cor=TRUE) 

# Lead-Lag with one lag at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=1, rvargs=list(align.period=60), cor=TRUE) 

# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="avg", period=600, cor=TRUE) 


# Correlation Matrices
# Traditional Estimate at highest frequency 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="naive", period=1, cor=TRUE)

# Traditional Estimate at 10 minute frequency
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="naive", period=10, rvargs=list(align.period=60), cor=TRUE)

# Lead-Lag with one lag at one minute frequency 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="kernel", lags=1, rvargs=list(align.period=60), cor=TRUE) 

# Subsample Average Estimate with second aligned data at 600 subgrids >
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="avg", period=600, cor=TRUE) 


}

\keyword{methods}