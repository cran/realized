\name{rv.naive}
\alias{rv.naive}
\title{ Realized Variance }
\description{
  Traditional realized variance estimator.
}
\usage{
rv.naive(x, period, align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
    \item{...}{ ...}
  }

\value{
Sum of products of high frequency returns.
}
\references{

O. E. Barndorff-Nielson and N. Sheppard. Econometric analysis of realised volatility and its use in estimating
stochastic volatility models. \emph{Journal of the Royal Statistical Society}, Series B 64:253-280, 2002

T.G. Andersen, T. Bollerslev, F.X. Diebold, and P. Labys. The distribution of exchange rate volatility.
\emph{Journal of the American Statistical Association}, 96:42-55, 2001.

}
\author{Scott Payseur <spayseur@u.washington.edu> }

\seealso{\code{\link{rc.naive}}, \code{\link{rSignature}}, \code{\link{rRealizedVariance}}}
\examples{


data(msft.real.cts)
data(ge.real.cts)

# Traditional RV aligned to 60 seconds
#
rv.naive(x = msft.real.cts[[1]], period = 60)

}
\keyword{methods}