\name{ssMODWT}
\alias{ssMODWT}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Generate Objects of Class \dQuote{\code{\linkS4class{ssMODWT}}}

}

\description{

  This is the generic function for class \dQuote{ssMODWT}. Please see the
  associated method in \code{\link{ssMODWT-methods}} for more details.

} %description

\usage{
ssMODWT(ss, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{ss}{ Signature object. Because there is only one extant method, this argument name
    is derived from the class \dQuote{\code{\linkS4class{sampSurf}},} to which it applies. }

  \item{\dots}{  See the method(s) for details. }

 } %arguments
 
\details{

  Only one method is available for this generic function as noted above. Please see the
  link above for other arguments that are possible in the constructor.

} %details

\value{
  A valid object of class \dQuote{\code{\linkS4class{ssMODWT}}.}
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
\dQuote{\code{\linkS4class{ssMODWT}}}, \dQuote{\code{\linkS4class{sampSurf}}}.
}

\examples{
#
# creates a sampSurf object with horizontal point sampling
# then creates a J_0 = 5-level MODWT decomposition object...
#
tr = Tract(c(x = 64, y = 64), cellSize = 1) #square tract ~0.5ha
btr = bufferedTract(10, tr)
ag3m = angleGauge(3)        #metric BAF
sshps = sampSurf(10, btr, iZone = 'horizontalPointIZ', angleGauge = ag3m,
        topDiam = c(0,0), startSeed = 123)
modwt.hps = ssMODWT(sshps, J = 5)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ generic }% use one of  RShowDoc("KEYWORDS")

