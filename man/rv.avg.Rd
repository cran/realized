\name{rv.avg}
\alias{rv.avg}

\title{Realized Variance: Average Subsample}
\description{
  Realized Variance using average subsample.
}
\usage{
rv.avg(x,period, align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}

\arguments{
   \item{x}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
    \item{...}{ ...}
}
\value{
 Realized variance using average subsample.
}
\references{
L. Zhang, P.A Mykland, and Y. Ait-Sahalia. A tale of two time scales: Determining integrated volatility
with noisy high-frequency data. \emph{Journal of the American Statistical Association}, 2005.

}
\author{ Scott Payseur <spayseur@u.washington.edu> }
\seealso{\code{\link{rv.avg}}, \code{\link{rRealizedVariance}} }
\examples{
 
data(msft.real.cts)


#
# Average subsampled realized variance for CTS aligned at one second returns at 
# 600 subgrids (10 minutes).
#
rv.avg(x = msft.real.cts[[1]], period = 600)

#
# Average subsampled realized variance for CTS aligned at one minute returns at 
# 5 subgrids (5 minutes).
#
rv.avg(x = msft.real.cts[[1]], period = 5, align.period=60)


}
\keyword{methods}