\name{ssToroid}
\alias{ssToroid}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Toroidal boundary correction for \dQuote{sampSurf} Objects
  
}

\description{

  This routine takes a \dQuote{sampSurf} object and makes a larger
  \dQuote{\code{\linkS4class{Raster}}} object with copies of the original all around it
  using toroidal correction, more commonly called periodization or circular correction for
  the boundary.

  Toroidal correction is specifically used for correction of \dQuote{edge effect} in the
  application of a \acronym{MODWT} wavelet filter as discussed in detail in Lark and
  Webster (2004); see details below.

}

\usage{
ssToroid(ss, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ss}{ An object of class \dQuote{\code{\linkS4class{sampSurf}}}. }

  \item{\dots}{ Gobbled. }

} %arguments

\details{

  This method is meant to be used stand-alone as an illustration for what happens in
  periodic correction. An example is given below, and a more extensive discussion can be
  found in the package vignette.

  It should be kept in mind that the \code{ssMODWT} constructor uses the \pkg{waveslim}
  code for the decomposition. This code \emph{always} applies a periodic correction to the
  image. Therefore, there is no need to use this routine for this type of correction when
  creating objects of class \dQuote{ssMODWT}.

  I suppose this could be used in \pkg{sampSurf} proper as a boundary correction technique
  for inclusion zones. Perhaps some version of it will make its way into that package in
  the future to compliment those boundary correction methods already available.

} %details

\value{

  A \dQuote{Raster} image padded in all directions using the toriodal correction.

}

\references{

\describe{
 
  \item{}{R. M. Lark and R. Webster. 2004. Analysing soil variation in two dimensions with
    the discrete wavelet transform. \emph{European Journal of Soil Science},
    \bold{55}:777--797.}

} %describe
}

\author{
Jeffrey H. Gove
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ssReflect}}
}

\examples{
#
# creates a sampSurf object with horizontal point sampling
# then does a toroidal correction...
#
tr = Tract(c(x = 64, y = 64), cellSize = 1) #square tract ~0.5ha
btr = bufferedTract(10, tr)
ag3m = angleGauge(3)        #metric BAF
sshps = sampSurf(10, btr, iZone = 'horizontalPointIZ', angleGauge = ag3m,
        topDiam = c(0,0), startSeed = 123)
torhps = ssToroid(sshps)    #raster object
\dontrun{
plot(torhps, col = palMODWT(100, range=cellStats(torhps, range)))
plot(perimeter(sshps), add = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ function }% use one of  RShowDoc("KEYWORDS")

