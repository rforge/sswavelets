\name{plotMODWT2D}
\alias{plotMODWT2D}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  
Plot Two-Dimensional Wavelet Decomposition on \dQuote{sampSurf} Objects 

}

\description{

  This routine will plot any of the individual \dQuote{ssMODWT} decompositions below
  (\code{waveType}) at any \code{level} either as the raw wavelet coefficients or various
  wavelet variances (see Details for some restrictions). Covariance objects,
  \dQuote{ssCovMODWT}, are also supported.

}

\usage{
plotMODWT2D(ssMODWT, 
            level = 1, 
            waveType = c("ISO", "LH", "HL", "HH", "LL"), 
            decompType = c("modwt", "mra"), 
            type = c("raw", "var"), 
            isoSmooth = TRUE, 
            showPlot = TRUE, 
            showIZs = TRUE, 
            addLattice = FALSE, 
            col = NA, 
            boxCol = "gray", 
            title = NA, 
            ...
           )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ssMODWT}{ An object of class \dQuote{\code{\linkS4class{ssMODWT}}} for raw or
    variance display, or an object of class \dQuote{\code{\linkS4class{ssCovMODWT}}} for
    covariance display. }

  \item{level}{ The decomposition level desired, \eqn{j = 1,\ldots,J_0}{j = 1,\ldots,J0}. 
	Note that for \code{waveType = 'ISO'} and \code{type = 'var'}, it is
	legal to specify \code{j = 0}; please see the details below for an explanation. }

  \item{waveType}{ \sQuote{ISO} is for isotropic; the next three are the anisotropic
      decompositions: \sQuote{LH} horizontal; \sQuote{HL} vertical; \sQuote{HH}
      diagonal. Finally, the \sQuote{LL} is the smooth at the highest scale (level
      \eqn{J_0}{J0}). Please see the vignette for more details.}

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

  \item{addLattice}{ \code{TRUE}: create a \pkg{lattice} version of the graph (see
      details); \code{FALSE}: no lattice version. }

  \item{col}{ A vector palette of colors for the surface; if this is \code{NA}, the
      default, then \code{\link{palMODWT}} is used. }

  \item{boxCol}{ Outline color for perimeter box. }

  \item{title}{ Add a title to the plot if desired. The default (\code{NA}) is a title
      (and subtitles) constructed from the arguments passed. Otherwise, one can pass their
      own title or \code{title = \dQuote{}} for no title. }

  \item{\dots}{ Gobbled. }

} %arguments

\details{

  This routine has a fair bit of functionality and it is the basis for the
  \code{\link{plotLevel2D}} function as well. The vignette has a number of examples that
  can be reviewed for more information.

  Note that not all combinations of argument options are available. For example, variance
  is calculated only for \code{decompType = \dQuote{MODWT}} as noted above. Also,
  \code{level = 0} is only meaningful at this point for \code{waveType = \dQuote{ISO}} and
  \code{type = \dQuote{var}}, which presupposes \code{decompType = \dQuote{MODWT}}. Also,
  only \code{type = \dQuote{var}} is allowed when the \code{ssMODWT} argument passed is an
  object of class \dQuote{ssCovMODWT}. This latter restriction is because there the
  covariance is composed of two \dQuote{ssMODWT} objects; therefore, \code{type = \sQuote{raw}} %
  would be ambiguous. One can always plot the raw wavelets for the individual
  \dQuote{ssMODWT} objects that were used to create the covariance object.

  With \code{addLattice = TRUE}, the actual graphical object is created using
  \link[rasterVis]{levelplot} in the \pkg{rasterVis} package. Therefore, the latter must
  be available to use this option. The \dQuote{lattice} graph is returned and may be
  plotted as usual.


}

\value{

  A list with\dots

  \item{r: }{ The final \dQuote{\code{\linkS4class{RasterLayer}}} object. }

  \item{plt: }{ The \dQuote{lattice} plot object if desired; \code{NULL} otherwise. }

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
\code{\link{plotLevel2D}}
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
r1 = plotMODWT2D(modwt.hps, waveType = 'ISO', type = 'var', level = 1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ function }% use one of  RShowDoc("KEYWORDS")

