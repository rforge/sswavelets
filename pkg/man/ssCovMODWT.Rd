\name{ssCovMODWT}
\alias{ssCovMODWT}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Generate Objects of Class \dQuote{\code{\linkS4class{ssCovMODWT}}}

}

\description{

  This is the generic function for class \dQuote{ssCovMODWT}. Please see the
  associated method in \code{\link{ssCovMODWT-methods}} for more details.

} %description


\usage{
ssCovMODWT(ssMODWT.a, ssMODWT.b, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ssMODWT.a}{ First signature object. Because there is only one extant method, this
    argument name is derived from the class \dQuote{\code{\linkS4class{ssMODWT}},} to
    which it applies. }

  \item{ssMODWT.b}{ Second signature object. Again, this argument name is derived from the
    class \dQuote{\code{\linkS4class{ssMODWT}},} to which it applies. }

  \item{\dots}{  See the method(s) for details. }

}

\details{

  Only one method is available for this generic function as noted above. Please see the
  link above for other arguments that are possible in the constructor.

} %details


\value{
  A valid object of class \dQuote{\code{\linkS4class{ssCovMODWT}}.}
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
Jeffrey H. Gove
}

\note{

  The example below creates an object of class \dQuote{ssCovMODWT} in such a way that the
  two sampling methods used are comparable. For example, the same tract and population of
  standing trees is used for both sampSurf runs. Both the horizontal point and critical
  height methods use the same angle gauge, though this is not necessary, but nice for
  illustration. In addition, to make sure the inclusion zones for critical height sampling
  match those for horizontal point sampling, notice that the reference height is requested
  to be at \acronym{DBH} rather than ground level (the default), though again this is not
  strictly necessary, depending on the intent of the simulations. Finally, the
  \code{ssMODWT} decompositions share the same total decomposition level \eqn{J}{J}. The
  points above that should be routinely followed in all comparisons are to use \emph{(i)}
  the same tract, \emph{(ii)} the same population of trees and \emph{(iii)} the same
  decomposition level. Note that it is perfectly acceptable to compare two different basal
  area factors in \acronym{HPS}, for example, or two different plot sizes for fixed area
  plot sampling. The results of such comparisons might surprise you.

} %note

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\dQuote{\code{\linkS4class{ssMODWT}}}, \dQuote{\code{\linkS4class{sampSurf}}}.
}

\examples{
#
# creates a sampSurf object with horizontal point sampling BAF 5,
# then creates a J_0 = 4-level MODWT decomposition object...
#
tr = Tract(c(x = 64, y = 64), cellSize = 1) #square tract ~0.5ha
btr = bufferedTract(10, tr)
ag5m = angleGauge(5)
strees = standingTrees(10, btr, dbhs=c(15,25), topDiam=c(0,0), startSeed = 123)
strees.hps = standingTreeIZs(strees, 'horizontalPointIZ', angleGauge = ag5m)
ss.hps = sampSurf(strees.hps, btr)
modwt.hps = ssMODWT(ss.hps, J = 4)
#
# creates a sampSurf object with critical height sampling BAF 5,
# then creates a J_0 = 4-level MODWT decomposition object...
#
strees.chs = standingTreeIZs(strees, 'criticalHeightIZ', angleGauge = ag5m,
                             referenceHeight='dbh')
ss.chs = sampSurf(strees.chs, btr)
modwt.chs = ssMODWT(ss.chs, J = 4)
#
# now both HPS and CHS are to the same reference height with the 
# same population on the same tract, so we can decompose 
# the covariance...
#
modwt.cov = ssCovMODWT(modwt.hps, modwt.chs)
#
# take a look at the level j=3 surface...
#
\dontrun{
plotMODWT2D(modwt.cov, level=3, type='var')
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ generic }% use one of  RShowDoc("KEYWORDS")

