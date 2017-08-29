\name{ssCovMODWT-methods}
\docType{methods}
\alias{ssCovMODWT-methods}
%\alias{ssCovMODWT,ssMODWT,ssMODWT-method}

\title{

  Methods for \dQuote{ssCovMODWT} object construction in Package \pkg{ssWavelets}

}

\description{

  This is the single method for construction of \dQuote{\code{\linkS4class{ssCovMODWT}}}
  objects via the generic \code{\link{ssCovMODWT}}. Simply pass two valid \dQuote{ssMODWT}
  objects to the method with other appropriate arguments listed below as desired. If the
  first and second objects are \code{identical} one should simply get the variance
  results back (as a check).

}

\section{Methods}{
  
  \describe{
	

  \item{\strong{\code{signature(ssMODWT.a = "ssMODWT", ssMODWT.b = "ssMODWT")}}}{}


  \item{}{\bold{usage\ldots}
    \preformatted{ssCovMODWT(ssMODWT.a,
           ssMODWT.b, 
           description = 'sampSurf Covariance MODWT wavelet decomposition object',
           runQuiet = FALSE,
           ...) }
	   
    \itemize{
			
      \item{\code{ssMODWT.a}:}{ The first object of class
        \dQuote{\code{\linkS4class{ssMODWT}}} for the covariance decomposition. }
			
      \item{\code{ssMODWT.b}:}{ The second object of class
        \dQuote{\code{\linkS4class{ssMODWT}}} for the covariance decomposition. }
			
      \item{\code{description}:}{ A character description as desired. }
			
      \item{\code{runQuiet}:}{ \code{TRUE}: no feedback; \code{FALSE}: some results
          printed. }

	  \item{\code{\dots}:}{ Arguments currently gobbled. }
		
    } %itemize
  } %item

  } %describe
} %methods

\keyword{methods}

