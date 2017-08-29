\name{hfsPlot}
\alias{hfsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

H. F. Smith Plot 

}

\description{

  This routine will generate a lattice plot for the marginal isotropic wavelet variance in
  the form of an \sQuote{H. F. Smith} plot; that is, variance versus the average inclusion
  zone area increasing on the \eqn{x}{x}-axis for the different sampling methods.

  The routine assumes the existence of one or more data frames from
  \code{\link{hfsMODWT}}, which takes an \dQuote{\code{\linkS4class{ssMODWT}}} object set
  and exports the appropriate data frame---in this case the \code{hfsMODWT} output should
  be in the \emph{long} form.

}

\usage{
hfsPlot(...,
        sampMeth = c("Meth.A", "Meth.B"),
        conditionOn = c("iso.j", "tau.j", "j"),
        units = c("metric", "English"),
        showPlot = TRUE,
        fileName = "",
        ylab = "Isotropic Variance",
        xlab = "Average Inclusion Zone Area",
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

  \item{conditionOn}{ One of the factor or character columns in the data frame; for columns
     \code{tau.j} and \code{j}, they will be converted appropriately. }

  \item{units}{ The appropriate units that match the underlying sampling surface. Please
      be mindful that no checking is possible for the correct units as there is no connection
      between the original \dQuote{sampSurf} objects and the data frames used for input. }

  \item{showPlot}{ \code{TRUE}: display the plot; \code{FALSE}: no display. }

  \item{fileName}{ The file name for a \sQuote{hard copy}; \dQuote{\code{}} for no hard copy file
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

  Details of its use are given in the vignette along with examples, which are somewhat
  lengthy to set up, therefore there are no examples presented below.

}

\value{

  A list with\dots

  \item{df: }{ The concatenated data frames from \dQuote{\dots}. }

  \item{plt: }{ The \dQuote{lattice} plot object. }

}

\references{

\describe{
 
  \item{}{H. F. Smith. 1938. An empirical law describing heterogeneity in the yields of
    agricultural crops. \emph{Journal of Agricultural Science}, \bold{28}:1--23 }

} %describe

} %references

\author{
Jeffrey H. Gove
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hfsMODWT}}, \code{\link{varPlot}}
}

\examples{
#
# please see the examples in the vignette.
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ function }% use one of  RShowDoc("KEYWORDS")

