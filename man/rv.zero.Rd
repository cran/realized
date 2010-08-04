\name{rv.zero}
\alias{rv.zero}
\title{ Calculates the percentage of zero returns at a specified sampling period }
\description{
Calculates the percentage of zero returns at a specified sampling period.
}
\usage{
rv.zero(x, period, align.period = 1, cts = TRUE, makeReturns = FALSE, ...)
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
Percentage of zero returns.
}
\references{
S. W. Payseur.  A One Day Comparison of Realized Variance and Covariance Estimators. \emph{Working Paper: University 
of Washington}, 2007
}
\author{ Scott Payseur <spayseur@u.washington.edu>}

\seealso{ \code{\link{rSignature}}, \code{\link{rc.zero}} }


\examples{

data(msft.real.cts)

#
# Plot the percentage of zero returns for sampling periods of 1 second
# to 600 seconds.
#
plot(rSignature(1:600, x=msft.real.cts[[1]], type="zero", xscale=1/60), main="Percentage of co-zeros", ylab="Percent Zero", xlab= "minutes")

}
\keyword{methods}