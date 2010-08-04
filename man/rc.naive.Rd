\name{rc.naive}
\alias{rc.naive}
\title{ Realized Covariance }
\description{
  Traditional realized covariance estimator.
}
\usage{
rc.naive(x, y, period, align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{  RealizedObject or TimeSeries for S+}
  \item{y}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
    \item{...}{ ...}
  }

\value{
Sum of cross products of high frequency returns.
}
\references{
T.G. Andersen, T. Bollerslev, F.X. Diebold, and P. Labys. The distribution of exchange rate volatility.
\emph{Journal of the American Statistical Association}, 96:42-55, 2001.

}
\author{Scott Payseur <spayseur@u.washington.edu> }

\seealso{\code{\link{rv.naive}}, \code{\link{rSignature}}, \code{\link{rRealizedVariance}}}
\examples{


data(msft.real.cts)
data(ge.real.cts)

# Traditional RC aligned to 60 seconds
#
rc.naive(x = msft.real.cts[[1]], y = ge.real.cts[[1]], period = 60)
}
\keyword{methods}
