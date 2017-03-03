\name{palMODWT}
\alias{palMODWT}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Simple Palette Function for \dQuote{ssWavelets} and \dQuote{sampSurf} 
}

\description{
This simply gives a pleasing palette with reds for negative (if any)
and blues up to brownish for positive; zero should always be white-ish,
which was chosen for the background.
}

\usage{
palMODWT(n, bias = 5, range = NA, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{n}{ The length of the palette vector returned. }

  \item{bias}{ Please see \code{\link{colorRampPalette}}. }

  \item{range}{ The range of the variable in question to be displayed. }

  \item{\dots}{ Passed to \code{colorRampPalette}. }

}

%\details{
%}

\value{
  A vector of length \code{n} of colors in the palette.
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
\code{\link{colorRampPalette}}
}
\examples{
#
# display a simple palette...
#
n = 100
pie(rep(1,n), col=palMODWT(n, range=c(0,100)))
pie(rep(1,n),
    labels = as.character(round(seq(-10, 10, len = n), 1)),
    col = palMODWT(n, range = c(-10, 10)))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ function }% use one of  RShowDoc("KEYWORDS")

