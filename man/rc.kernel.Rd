\name{rc.kernel}
\alias{rc.kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Realized Covariance: Kernel  }
\description{
 Realized covariance calculation using a kernel estimator. 
}
\usage{
rc.kernel(x, y, q, align.period = 1, adj = TRUE, type = 0, cts = TRUE, makeReturns = FALSE,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  RealizedObject or TimeSeries for S+}
  \item{y}{  RealizedObject or TimeSeries for S+}
  \item{q}{Number of lags}
  \item{align.period}{ Align the returns to this period first}
  \item{cts}{ Create calendar time sampling if a non realizedObject is passed}
  \item{makeReturns}{ Prices are passed make them into log returns}
  \item{adj}{T to use dof adjustment}
  \item{type}{0-11 or a character string}
    \item{...}{ ...}
}

\details{
The different types of kernels can be found using rKernel.available().
}
\value{
Kernel estimate of realized covariance.
}
\references{
Ole E. Barndorff-Nielsen, Peter Reinhard Hansen, Asger Lunde, and Neil Shephard. Regular and modified kernel-based estimators of integrated variance: The case with independent noise. \emph{Working Paper}, 2004.

Michiel de Pooter, Martin Martens, and Dick van Dijk. Predicting the daily covariance matrix for sp100
stocks using intraday data - but which frequency to use? \emph{Working Paper}, October 2005.

J. E. Griffen and R. C. A. Oomen. Covariance measurement in the presence of non-synchronous trading and
market microstructure noise. \emph{Working Paper}, June 27th, 2006.
}
\author{Scott Payseur <spayseur@u.washington.edu>}

\seealso{\code{\link{rRealizedVariance}},\code{\link{rv.kernel}}, \code{\link{rKernel.available}},\code{\link{rKernel}}}
\examples{

data(msft.real.cts)
data(ge.real.cts)

# kernel realized covariance for CTS aligned at one minute returns
#
rc.kernel(x = msft.real.cts[[1]], y = ge.real.cts[[1]], q=1, type="bartlett", align.period=60)
}
\keyword{methods}
