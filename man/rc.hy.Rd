\name{rc.hy}
\alias{rc.hy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hayashi-Yoshida }
\description{
Hayashi-Yoshida covariance estimator
}
\usage{
rc.hy(x, y, period = 1, align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
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

\references{ T. Hayashi and N. Yoshida. On covariance estimation of non-synchronously observed diffusion processes.
\emph{Bernoulli}, 11:359-379, 2005. }
\author{ Scott Payseur }
\seealso{ rRealizedVariance }
\examples{
data(msft.real.cts)
data(ge.real.cts)

# Hayashi-Yoshida realized covariance for CTS aligned at one second returns at 
# 600 subgrids (10 minutes).
#
rc.hy(x = msft.real.cts[[1]], y = ge.real.cts[[1]], period = 600)

}
\keyword{methods}
