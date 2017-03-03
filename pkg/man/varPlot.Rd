\name{varPlot}
\alias{varPlot}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Summary Wavelet Variance Plot

}

\description{

  This routine creates a plot of the marginal isotropic variances against distance/scale
  conditioned on either \emph{(i)} each of the different sampling methods, one panel for
  each method+size (i.e., baf, plot size) combination: \code{groupSM = FALSE}; or
  \emph{(ii)} just the \code{sampMethod} level with each of the \sQuote{sizes} grouped:
  \code{groupSM = TRUE}. Note that using \emph{(i)} could conceivably get rather busy
  panel-wise.

}

\usage{
varPlot(..., 
        sampMeth = c("Meth.A", "Meth.B"), 
        groupSM = TRUE, 
        units = c("metric", "English"), 
        isoSmooth = TRUE, 
        showPlot = TRUE, 
        fileName = "", 
        ylab = "Isotropic Variance", 
        xlab = "Distance", 
        scales = "same", 
        type = "b", 
        pch = 19, 
        as.table = TRUE, 
        theme = c("plain", "custom", "ggplot", "economist"), 
        runQuiet = FALSE
       )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{\dots}{ Any number of \emph{long-form} data frames as generated from
    \code{\link{hfsMODWT}}. }

  \item{sampMeth}{ A character vector of sampling method identifiers that will give a label
     to each of the data frame objects in \sQuote{\ldots} that is useful for
     grouping. For example, if the first data frame is from fixed-area plot sampling
     and the second from critial height sampling, one might use
     \code{sampMeth = c("FAP","CHS")} for the identifiers of these methods. }
  
  \item{groupSM}{ \code{TRUE}: condition on \code{sampMeth} and group on the \code{id}
     column in the data frames; \code{FALSE}: condition on \code{id} and do not group. }

  \item{units}{ The appropriate units that match the underlying sampling surface. Please
      be mindful that no checking is possible for the correct units as there is no connection
      between the original \dQuote{sampSurf} objects and the data frames used for input. }

  \item{isoSmooth}{ \code{TRUE}: include \sQuote{LLJ} (smooth) in the result for the
      isotropic variance at level \eqn{J}{J}; \code{FALSE}: do not include the smooth
      component. }

  \item{showPlot}{ \code{TRUE}: display the plot; \code{FALSE}: no display. }

  \item{fileName}{ The file name for a \sQuote{hard copy}; \code{\dQuote{}} for no hard copy file
      of the plot. }

  \item{ylab}{ See \code{\link{par}}. }

  \item{xlab}{ See \code{\link{par}}. }

  \item{scales}{ See \code{\link{xyplot}}. }

  \item{type}{ See \code{\link{par}}. }

  \item{pch}{ See \code{\link{par}}. }

  \item{as.table}{ See \code{\link{xyplot}}. }

  \item{theme}{ One of the lattice themes listed in the argument default; see the styles
      in \pkg{latticeExtra}. }

  \item{runQuiet}{ \code{TRUE}: No feedback when run; \code{FALSE}: Some summary results. }

} %arguments

\details{

  The \code{sampMeth} argument is helpful for distinguishing between sampling methods if
  more than one different method is used in the data frame(s) passed. Simply use a vector
  of blank strings if nothing is desired.

  A hard copy to a file may be \code{print}ed if \code{hardcopyLattice} is available (it
  resides in another of the author's packages and is available on request). Otherwise,
  simply \code{print}/\code{plot} the \code{plt} component of the return list to a
  \code{\link{trellis.device}} as usual for \dQuote{lattice} objects.

  Note that if \code{groupSM = TRUE} then we will have as many panels as sampling methods;  if it is
  \code{FALSE}, then the number of panels is as above in \emph{(i)}.

  If we want to look at a plot of just the isotropic wavelet variances without the
  smooth/scale component (\code{LL.var} in the \dQuote{ssMODWT} object) at level
  \eqn{J_0}{J0}, then one can use \code{isoSmooth = FALSE}. 

  Details of its use are given in the vignette along with examples, which are somewhat
  lengthy to set up, therefore there are no examples presented below.

} %details

\value{

  A list with\dots

  \item{df: }{ The concatenated data frames from \dQuote{\dots}. }

  \item{plt: }{ The \dQuote{lattice} plot object. }

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
\code{\link{hfsMODWT}}, \code{\link{hfsPlot}}
}

\examples{
#
# please see the examples in the vignette.
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ plot }% use one of  RShowDoc("KEYWORDS")

