\name{sjedist}
\title{Spatial distributions}
\description{
  top level library for analysing retinal mosaics.
}
\usage{
}
\arguments{
}
\value{
}
\details{
  
  Follow the Rodieck prescription for finding the effective radius and
  the packing factor.  Note this method is very close to that of the
  Ripley method for computing Khat.

  \code{drp.makestf} takes a set of data points and will output a text
  file (with Macintosh end of lines) that can then be input into the
  MacDRP program from Bob Rodieck.  This serves as a useful comparison.
}
\references{
  Rodieck RW (1991). The density recovery profile: A method for the
  analysis of points in the plane applicable to retinal studies. 
  Visual Neuroscience 6: 95-111.
}
\seealso{
  khat() and other functions by Ripley and Diggle.  See also
  `rodieckach' for the data and code to generate the DRP from Rodieck
  (1992).
}


\author{Stephen Eglen}
\keyword{spatial}
