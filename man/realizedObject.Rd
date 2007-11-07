\name{realizedObject}
\alias{realizedObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ realizedObject }
\description{
Creates an object of type realizedObject to be passed to most of the functions from the realized library
}
\usage{
realizedObject(x, cts = TRUE, makeReturns = FALSE, millisstart = 34200000, millisend = 57600000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Either an S+ timeSeries or a list with $data that is prices or returns, and $ milliseconds }
  \item{cts}{ Calendar time sampling or tick time sampling }
  \item{makeReturns}{ If you pass prices you should have this set to T }
  \item{millisstart}{ When does your integrated variance start 34200000 is 09:30 }
  \item{millisend}{ When does your integrated variance end 57600000 is 16:00 }
}
\details{
See user's manual.
}
\value{
Object of type realizedObject
}
\author{ Scott Payseur <spayseur@u.washington.edu> }
\examples{
# see users manual

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
