\name{plotLevel2D}
\alias{plotLevel2D}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

Plot a Full Set of Images for a \dQuote{sampSurf} Wavelet Decomposition

}

\description{

  This routine will plot the full set of isotropic and anisotropic images for a given
  \code{level} of decomposition either as the raw wavelet coefficients or various wavelet
  variances. Both \dQuote{ssMODWT} and covariance objects, \dQuote{ssCovMODWT}, are
  supported.

  The \code{level = j} decompositions are plotted in a 2x2 matrix with\dots
 
  \tabular{lc}{ 
     1st row: \tab  Horizontal | Diagonal \cr
     2nd row: \tab  Isotropic   | Vertical \cr
  }

  At \code{level = J0}, there are five subfigures as the smooth, \sQuote{LLJ}, is also
  included centered in the 3rd row.


}

\usage{
plotLevel2D(ssMODWT,
            level = 1,
            decompType = c("modwt", "mra"),
            type = c("raw", "var"),
            isoSmooth = TRUE,
            showPlot = TRUE,
            showIZs = TRUE,
            col = NA,
            boxCol = "gray",
            title = NA,
            showLegend = TRUE,
            runQuiet = FALSE,
            ...
           )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ssMODWT}{ An object of class \dQuote{\code{\linkS4class{ssMODWT}}} for raw or
    variance display, or an object of class \dQuote{\code{\linkS4class{ssCovMODWT}}} for
    covariance display. }

  \item{level}{ The decomposition level desired, \eqn{j = 1,\ldots,J_0}{j =
        1,\ldots,J0}. }
  

  \item{decompType}{ Either \dQuote{MODWT} or \dQuote{MRA} decomposition may be displayed
      for \code{type = 'raw'}. Only the former is available for variance or covariance
      displays. }

  \item{type}{ \sQuote{raw} corresponds to the raw wavelet coefficients. \sQuote{var} is
      for a plot of the variance. }

  \item{isoSmooth}{ \code{TRUE}: include \sQuote{LLJ} in the result for isotropic at level
      \eqn{J}{J} (i.e., \eqn{J_0}{J0}); \code{FALSE}: do not include the smooth
      component. }

  \item{showPlot}{ \code{TRUE}: display the plot; \code{FALSE}: no display. }

  \item{showIZs}{ \code{TRUE}: display the inclusion zones for the \dQuote{Stem} objects;
      \code{FALSE}: no display. }

  \item{col}{ A vector palette of colors for the surface; if this is \code{NA}, the
      default, then \code{\link{palMODWT}} is used. }

  \item{boxCol}{ Outline color for perimeter box. }

  \item{title}{ Add a title to the plot if desired. The default (\code{NA}) is a title
      (and subtitles) constructed from the arguments passed. Otherwise, one can pass their
      own title or \code{title = \dQuote{}} for no title. }


  \item{showLegend}{ \code{TRUE}: Show the \sQuote{legend} identifiers around the outside
      of the maps; \code{FALSE}: suppress the legend. }

  \item{runQuiet}{ \code{TRUE}: No feedback when run; \code{FALSE}: Some summary results. }

  \item{\dots}{ Gobbled. }

} %arguments

\details{

  This routine uses \code{\link{plotMODWT2D}} to generate each of the individual
  subfigures, which can be altered using the arguments above that are in common to
  \code{plotMODWT2D}.

  Please see the details in \code{plotMODWT2D} and the vignette for more information and
  examples. 

}

\value{

The final \dQuote{\code{\linkS4class{RasterLayer}}} object invisibly. 

}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Jeffrey H. Gove
}

%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotMODWT2D}}
}

\examples{
#
# 1. creates a sampSurf object with horizontal point sampling
# 2. creates a J_0 = 3-level MODWT decomposition object
# 3. a basic plot...
#
tr = Tract(c(x = 64, y = 64), cellSize = 1) #square tract ~0.5ha
btr = bufferedTract(10, tr)
ag3m = angleGauge(3)        #metric BAF
sshps = sampSurf(10, btr, iZone = 'horizontalPointIZ', angleGauge = ag3m,
        topDiam = c(0,0), startSeed = 123)
modwt.hps = ssMODWT(sshps, J = 3)
\dontrun{
r2 = plotLevel2D(modwt.hps, type = 'var')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ function }% use one of  RShowDoc("KEYWORDS")

