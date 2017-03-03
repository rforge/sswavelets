\name{ssReflect}
\alias{ssReflect}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Reflection boundary correction for \dQuote{sampSurf} Objects

}

\description{

  This routine takes a \dQuote{sampSurf} object and makes a larger
  \dQuote{\code{\linkS4class{Raster}}} object with copies of the original all around it
  using reflection.

  Reflection is specifically used for correction of \dQuote{edge effect} in the
  application of a \acronym{MODWT} wavelet filter as discussed in detail in Lark and
  Webster (2004); see details below.

}

\usage{
ssReflect(ss, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{ss}{ An object of class \dQuote{\code{\linkS4class{sampSurf}}}. }

  \item{\dots}{ Gobbled. }

}

\details{

  This routine may be used stand-alone, but its main use it to be called from within
  \code{\link{ssMODWT}} when constructing objects of class
  \dQuote{\code{\linkS4class{ssMODWT}}}. In this case, the filter is applied to the
  extended raster scene encompassing the reflected replicates and then clipped back to the
  original extents of the \dQuote{sampSurf} image for further processing.

  It should be kept in mind that the \code{ssMODWT} constructor uses the \pkg{waveslim}
  code for the decomposition. This code \emph{always} applies a periodic correction to the
  image. Therefore, if one chooses the reflection as well, the periodic correction is
  applied to the extended relected image set in the decomposition.

} %details

\value{

  A \dQuote{Raster} image padded in all directions using the reflection correction.

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
\code{\link{ssToroid}}
}
\examples{
#
# creates a sampSurf object with horizontal point sampling
# then applies reflection...
#
tr = Tract(c(x = 64, y = 64), cellSize = 1) #square tract ~0.5ha
btr = bufferedTract(10, tr)
ag3m = angleGauge(3)        #metric BAF
sshps = sampSurf(10, btr, iZone = 'horizontalPointIZ', angleGauge = ag3m,
        topDiam = c(0,0), startSeed = 123)
refhps = ssReflect(sshps)    #raster object
\dontrun{
plot(refhps, col = palMODWT(100, range=cellStats(refhps, range)))
plot(perimeter(sshps), add = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

