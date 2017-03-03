\name{ssEnergy}
\alias{ssEnergy}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Calcluate the \dQuote{Average Energy} from a \dQuote{sampSurf} Image
}

\description{

  This will calculate the \dQuote{average energy} surface for the \dQuote{sampSurf} object
  slot in an \dQuote{\code{\linkS4class{ssMODWT}}} or
  \dQuote{\code{\linkS4class{ssCovMODWT}}} object. In the case case of covariance, this is
  simply the average covariance energy.

}

\usage{
ssEnergy(ssMODWT, showPlot = TRUE, showIZs = TRUE, boxCol = "gray", col = NA, title = NA, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ssMODWT}{ Either a valid \dQuote{ssMODWT} or \dQuote{ssCovMODWT}
	object. }

  \item{showPlot}{ \code{TRUE}: display the surface; \code{FALSE}: no display. }

  \item{showIZs}{ \code{TRUE}: display the inclusion zones; \code{FALSE}: don't display. }

  \item{boxCol}{ Outline color for perimeter box. }

  \item{col}{ A vector palette of colors for the surface; see, e.g.,
      \code{\link{palMODWT}}. }

  \item{title}{ If \code{showPlot}, a title, or \code{NA} for nothing. }

  \item{\dots}{ Gobbled. }

} %arguments

\details{

  The \dQuote{average energy} in a scene is essentially the mean sums of squares that is
  probably more familiar to most users. The former terminology is used by physicists.

  The package vignette in the references below gives more information on this subject.

}

\value{
A \dQuote{\code{\linkS4class{Raster}}} object corresponding to the average energy.
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Jeffrey H. Gove
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{covMODWT}}
}
\examples{
#
# 1. creates a sampSurf object with horizontal point sampling
# 2. creates a J_0 = 3-level MODWT decomposition object
# 3. calculates the average energy of the surface...
#
tr = Tract(c(x = 64, y = 64), cellSize = 1) #square tract ~0.5ha
btr = bufferedTract(10, tr)
ag3m = angleGauge(3)        #metric BAF
sshps = sampSurf(10, btr, iZone = 'horizontalPointIZ', angleGauge = ag3m,
        topDiam = c(0,0), startSeed = 123)
modwt.hps = ssMODWT(sshps, J = 3)
ae = ssEnergy(modwt.hps, showPlot = FALSE) #showPlot=TRUE to display
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ function }% use one of  RShowDoc("KEYWORDS")

