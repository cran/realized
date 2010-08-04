\name{rKernel}
\alias{rKernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Kernel Function }
\description{
  The kernel function that is used for the realized variance kernel estimators.
}
\usage{
rKernel(x, type = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Value between [0,1] to calculate kernel weight }
  \item{type}{ Type of kernel 0-11, or values from rKernel.available() }
}
\value{
Kernel value
}
\references{
Ole E. Barndorff-Nielsen, Peter Reinhard Hansen, Asger Lunde, and Neil Shephard. Regular and modified kernel-based estimators of integrated variance: The case with independent noise. \emph{Working Paper}, 2004.

 }
\author{ Scott Payseur <spayseur@u.washington.edu> }
\seealso{\code{\link{rKernel.available}},\code{\link{rv.kernel}},\code{\link{rc.kernel}}}
\examples{
#
# Plot all available kernels
#
par(mfrow=c(3,4))
x <- (0:100)*.01
for(i in 1:length(rKernel.available()))
    plot(x=x,y=sapply(x, FUN="rKernel", type=rKernel.available()[i]), xlab="", ylab="", main=rKernel.available()[i],ylim=c(0,1))


}
\keyword{methods}