\name{befp34}
\alias{befp34}

\docType{data}

\title{

Plot 34 Bartlett Density Study
  
}
\description{

  This data set comes from 1989 remeasurements on the \sQuote{Density Study} at the
  Bartlett Experimental Forest, Bartlett, NH, USA. The Bartlett Experimental Forest is
  part of the White Mountain National Forest.

}

\usage{data("befp34")}

\format{
  The format is:
 chr "befp34"
}

\details{

  The \sQuote{Density Study} is one of many silvicultural treatments in the Bartlett
  Experiment Forest, Bartlett, NH, USA. The study was designed with the idea of applying
  different treatments to existing even-aged stands with the intent of eventual conversion
  to uneven-aged. More details may be found in Leak and Solomon (1975) or Leak and Gove
  (2008). The original treatment for plot 34 was cutting to 40 square feet of basal area
  per acre with 30 percent of the residual basal area in sawtimber trees over 10.5 inches
  dbh.

  Note that this study was conducted in \sQuote{English} (U.S. Customary) units, therefore
  all locations and measurements are in this system. It is simple enough to convert to
  metric if desired. These data have left it in \sQuote{English} because it gives better
  detail in the \dQuote{sampSurf} and subsequent \dQuote{ssMODWT} image results.

  The columns or variables in the data frame are those required to create an object of
  class \code{\linkS4class{standingTrees}}, plus a few extras\ldots

\describe{

  \item{plotNo}{The plot number.}

  \item{treeNo}{The tree number within this plot.}

  \item{x,y}{The tree location coordinates in feet. The northwest corner of the physical
      plot is at \eqn{(0,0)}. }

  \item{spp.num}{The species number corresponding to the codes below.}

  \item{species}{The species code. Legal species codes are\ldots

    \tabular{ll}{
      \sQuote{HE} \tab Eastern hemlock (\emph{Tsuga canadensis}) \cr
      \sQuote{BE} \tab American beech (\emph{Fagus grandifolia}) \cr
      \sQuote{PB} \tab Paper birch (\emph{Betula papyrifera}) \cr
      \sQuote{RM} \tab Red maple (\emph{Acer rubrum}) \cr
      \sQuote{RS} \tab Red spruce (\emph{Picea rubens}) \cr
      \sQuote{SM} \tab Sugar maple (\emph{Acer saccharum}) \cr
      \sQuote{WA} \tab White ash (\emph{Fraxinus americana}) \cr
      \sQuote{YB} \tab Yellow birch (\emph{Betula alleghaniensis}) \cr
      \sQuote{OT} \tab Other species \cr
    } %tabular species
  } %item

  \item{dbh}{Diameter at breast height in inches.}

  \item{height}{A modeled height in feet; the \sQuote{all species} equation in Fast and
      Ducey (2011) was used.}

  \item{solidType}{A randomly generated solid type in the range \eqn{[2,4]}{[2,4]}; see
      the details for the slot of the same name in \code{\linkS4class{standingTree}}.}

  \item{topDiam}{The top diameter in inches. These should all be zero, tapering to a
      \sQuote{tip} for all species.}

} %describe

} %details

\source{

The full data set will eventually be available online. It is currently available from the
package author.

}

\references{

Please refer to the package vignette for examples of the use of this data set.

\describe{

  \item{}{The \href{https://www.nrs.fs.fed.us/ef/locations/nh/bartlett/}{Bartlett}
      Experimental Forest.}

  \item{}{The White Mountain National Forest, \href{https://www.fs.usda.gov/r9}{Eastern
        Region, R9}.}

  \item{}{Leak, W. B. and D. S. Solomon. 1975. Influence of Residual Stand Density on
      Regeneration of Northern Hardwoods. USDA Forest Service Research Paper
      NE-310. Northeastern Forest Experiment Station, Upper Darby, PA.
      \url{https://www.treesearch.fs.fed.us/pubs/15417}. }

  \item{}{Leak, W. B. and J. H. Gove. 2008. Growth of Northern Hardwoods in New England: a
      25-Year Update. \emph{Northern Journal of Applied Forestry} \bold{25}(2):103--105.}

  \item{}{Fast, A. J. and M. J. Ducey. 2011. Height-Diameter Equations for Select New
      Hampshire Tree Species. \emph{Northern Journal of Applied Forestry} \bold{28}(3):157--160.}

} %describe
} %references

\examples{
#
#LazyData, so just use it; e.g.,...
#
dim(befp34)
head(befp34)
}
\keyword{datasets}
