\name{rv.timescale}
\alias{rv.timescale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Realized Variance: Two Timescales}
\description{
Realized variance using the two timescale variance method.
}
\usage{
rv.timescale(x, period, align.period = 1, adj.type = "classic", cts = TRUE, makeReturns = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{  RealizedObject or TimeSeries for S+}
  \item{period}{ Sampling period }
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
  \item{adj.type}{"classic", "adj" or "aa"}
    \item{...}{ ...}
}
\details{
 Realized variance using two timescale method.  


}
\value{
  Realized variance using two timescale method
}
\references{
L. Zhang, P.A Mykland, and Y. Ait-Sahalia. A tale of two time scales: Determining integrated volatility with noisy high-frequency data. \emph{Journal of the American Statistical Association}, 2005.


}
\author{ Scott Payseur <spayseur@u.washington.edu> }
\seealso{\code{\link{rc.timescale}}, \code{\link{rRealizedVariance}} }
\examples{

data(msft.real.cts)

rv.timescale(x = msft.real.cts[[1]], period = 60, adj.type="aa")

}
\keyword{methods}