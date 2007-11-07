\name{rc.timescale}
\alias{rc.timescale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Realized Covariance: Two Timescales}
\description{
Realized Covariance using a generalization of the popular two timescale variance method.}
\usage{
rc.timescale(x, y, period, align.period = 1, adj.type = "classic", cts = TRUE, makeReturns = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  RealizedObject or TimeSeries for S+}
  \item{y}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
  \item{adj.type}{"classic", "adj" or "aa"}
    \item{...}{ ...}
}
\details{
Realized Covariance using two timescale method.  

}
\value{
Realized covariance using two timescale method}

\references{
L. Zhang, P.A Mykland, and Y. Ait-Sahalia. A tale of two time scales: Determining integrated volatility
with noisy high-frequency data. \emph{Journal of the American Statistical Association}, 2005.

Michiel de Pooter, Martin Martens, and Dick van Dijk. Predicting the daily covariance matrix for sp100
stocks using intraday data - but which frequency to use? \emph{Working Paper}, October 2005.
}
\author{ Scott Payseur <spayseur@u.washington.edu> }
\seealso{\code{\link{rv.timescale}}, \code{\link{rRealizedVariance}} }
\examples{
data(msft.real.cts)
data(ge.real.cts)



rc.timescale(x = msft.real.cts[[1]], y = ge.real.cts[[1]], period = 60, adj.type="aa")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
