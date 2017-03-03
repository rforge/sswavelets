\name{hfsMODWT}
\alias{hfsMODWT}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

Prepare \dQuote{ssMODWT} Objects for Further Summary Analysis

}

\description{

  This function takes a number of \dQuote{ssMODWT} result objects and sets up a data frame
  with the isotropic variances for each level in the decomposition.  It also includes
  total variance, average energy and surface mean from the decompositions, among others.

  This is helpful in further analysis of the H. F. Smith type (using the companion
  \code{hfsPlot}), either by model estimation, plotting, or both.

}

\usage{
hfsMODWT(..., ids = NA, long = TRUE, runQuiet = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{\dots}{ Any number of \dQuote{\code{\linkS4class{ssMODWT}}} objects. }

  \item{ids}{ A character vector of simplified ids for each of the \dQuote{\dots}
      objects. The names of the objects will be used by default (\code{ids = NA}). The
      identifying names (\code{ids}) should correspond to some measure of the size of the
      inclusion zone; e.g., for circular plots, these might be \code{ids = c('cp.8',
        'cp.10')} for radii of 8 and 10m (ft), respectively. }

  \item{long}{ \code{TRUE}: The results are a data frame in \dQuote{long} format (more
      useful for, e.g., plotting with \pkg{lattice}); \code{FALSE}: the results are
      returned in \dQuote{wide} format. }

  \item{runQuiet}{ \code{TRUE}: No feedback when run; \code{FALSE}: Some summary results. }
}

\details{

  Normally one run of this routine might correspond to a given sampling method, but with
  different design parameters like plot size or basal area factor. The reason for this is
  that then multiple methods (data frames from this function) are passed to one of the
  plotting methods (see below) and are assigned secondary identifiers to distinguish
  between sets of sampling methods (i.e., circular plots, horizontal point sampling,
  etc.). This seemed to be the simplest way to structure it at the time.

  Note that one must not mix units of measure (\dQuote{English} and \dQuote{metric'}) in
  the \dQuote{ssMODWT} objects passed. No formal check is made for this, so be good.

  The average inclusion zone area is returned in the data frame as the unit of
  \sQuote{effort} corresponding to the various variances. For fixed-area plots, this is
  constant, for inclusion zones with area that is variable based on some size attribute of
  the standing trees or downed logs, this is the average over all inclusion areas in the
  population.

  Details of its use are given in the vignette along with examples, which are somewhat
  lengthy to set up, therefore there are no examples presented below.

}

\value{

  A data frame in either wide or long format. The columns for both
  options are discussed in the vignette.

}

\references{

\describe{
 
  \item{}{H. F. Smith. 1938. An empirical law describing heterogeneity in the yields of
    agricultural crops. \emph{Journal of Agricultural Science}, \bold{28}:1--23 }

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
\code{\link{hfsPlot}}, \code{\link{varPlot}}
}

\examples{
#
# please see the examples in the vignette.
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ function }% use one of  RShowDoc("KEYWORDS")

