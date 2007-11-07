\name{rv.kernel}
\alias{rv.kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Realized Variance: Kernel }
\description{
Realized variance calculation using a kernel estimator. 
}
\usage{
rv.kernel(x, q, align.period = 1, adj = TRUE, type = 0, cts = TRUE, makeReturns = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{  RealizedObject or TimeSeries for S+}
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
Kernel estimate of realized variance.
}
\references{
Ole E. Barndorff-Nielsen, Peter Reinhard Hansen, Asger Lunde, and Neil Shephard. Regular and modified kernel-based estimators of integrated variance: The case with independent noise. \emph{Working Paper}, 2004.

B. Zhou. High-frequency data and volatility in foreign-exchange rates. \emph{Journal of Buiness & Economic
Statistics}, 14:45-52, 1996.

P. Hansen and A. Lunde. Realized variance and market microstructure noise. \emph{Journal of Business and
Economic Statistics}, 24:127-218, 2006.
}
\author{Scott Payseur <spayseur@u.washington.edu>}

\seealso{\code{\link{rRealizedVariance}},\code{\link{rc.kernel}}, \code{\link{rKernel.available}},\code{\link{rKernel}}}
\examples{

data(msft.real.cts)
data(ge.real.cts)

# kernel realized variance for CTS aligned at one minute returns
#
rv.kernel(x = msft.real.cts[[1]], q=1, type="bartlett", align.period=60)

}
\keyword{methods}
