\name{covMODWT}
\alias{covMODWT}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Function to Calculate (Co-) Variances for MODWT Decompositions

}

\description{

  This routine calculates various variances \emph{for} objects of class
  \dQuote{\code{\linkS4class{ssMODWT}}}, or covariances \emph{for} objects of class
  \dQuote{\code{\linkS4class{ssCovMODWT}}}. It is rather complex, so please see the
  \emph{Details} section and the vignette referenced below for more information. Please
  also see the caveat in \emph{Details} concerning the general use of this function.

}

\usage{
covMODWT(ss.modwt.x, ss.modwt.y, ...)
}


\arguments{

  \item{ss.modwt.x}{ An list object as returned from \code{waveslim::modwt.2d}. }

  \item{ss.modwt.y}{ A second list object from \code{waveslim::modwt.2d} for covariance
      calculations. If variance calculations are desired, then leave this argument out as
      the routine tests for this argument as \code{missing} and makes a second copy of it
      from the first argument for variance calculation. }

  \item{\dots}{ Gobbled. }
} %arguments

\details{

   This routine acts much like the \code{stats::cov} routine in the sense that if
   the second argument, \code{ss.modwt.y} is \code{missing}, then it calculates the
   variance instead of the covariance. The return list object from this routine is quite
   long and complicated and is explained in adequate detail in the vignette.

   Please note that this routine is an integral part of the constructor methods for both
   \dQuote{ssMODWT} and \dQuote{ssCovMODWT} objects, where the results are automatically
   incorporated into the correct slot in these objects. Because of this, there is very
   little need to call this routine on your own, in fact, it is discouraged.

}

\value{

This is a very succinct summary, please see the vignette and/or code for details. A list
with components\ldots\cr

\item{isCovariance: }{\code{TRUE}: covariances; \code{FALSE}: variances}

\item{summary: }{a list of marginal total (co-) variances (vectors)}

\item{total: }{a list of surface decomposition total (co-) variances (scalars) as well as
    the smooth (surface) mean}

\item{image: }{a list of matrix/image (co-) variances}

}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Jeffrey H. Gove
}

\note{

Again, there should be little need to use this outside of the constructor functions. If
one is extending the classes in this package and requires the calculation of (co-)
variances, then please see the extant constructors for detailed examples. For this reason,
no examples are presented below.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\dQuote{\code{\linkS4class{ssMODWT}}}, \dQuote{\code{\linkS4class{ssCovMODWT}}}.
}

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ function }% use one of  RShowDoc("KEYWORDS")

