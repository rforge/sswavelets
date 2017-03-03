\name{ssMODWT-methods}
\docType{methods}
\alias{ssMODWT-methods}
\alias{ssMODWT,sampSurf-method}

\title{ Methods for \dQuote{ssMODWT} object construction in Package \pkg{ssWavelets} }

\description{

  This is the single method for construction of \dQuote{\code{\linkS4class{ssMODWT}}}
  objects via the generic \code{\link{ssMODWT}}. Simply pass a valid \dQuote{sampSurf}
  object to the method with other appropriate arguments listed below as desired.

}

\section{Methods}{

\describe{

  \item{\strong{\code{\code{signature(ss = "sampSurf")}}}}{}


  \item{}{\bold{usage\ldots}
    \preformatted{ssMODWT(ss,
        J = NA,  
        wfName = c('haar'),  
        reflect = FALSE,
        shift = TRUE,        
        trimRight = TRUE,    
        description = 'sampSurf MODWT wavelet decomposition object',
        runQuiet = FALSE,
        ...) }
	   
    \itemize{
			
      \item{\code{ss}:}{ An object of class \dQuote{\code{\linkS4class{sampSurf}}} on
        which the MODWT decomposition is to be made.}
			
      \item{\code{J}:}{ Desired maximum level for the MODWT decomposition, often denoted as
          \eqn{J_0}{J0}.}
			
      \item{\code{wfName}:}{The wavelet filter name. Ostensibly, this should match one of
          the filters allowed in \code{waveslim::wave.filter}. See the class definition
          for those currently supported.}
			
      \item{\code{reflect}:}{ \code{TRUE}: reflect the sampSurf results (note: this is in
          addition to periodic correction, which is always done); \code{FALSE}: don't
          reflect.}
			
      \item{\code{shift}:}{ \code{TRUE}: apply \code{waveslim::shift.2d} to
          rectify/realign the \acronym{MODWT} images with the original at all levels;
          \code{FALSE}: leave the alignment as is.}
			
      \item{\code{trimRight}:}{ When \code{shift=TRUE} occasionally the image has come out
          larger at the first decomposition level that the original image, this will trim
          the offending extra on the right side and bottom; \code{FALSE}: leave as is. }
			
      \item{\code{description}:}{ A character description as desired. }
			
      \item{\code{runQuiet}:}{ \code{TRUE}: no feedback; \code{FALSE}: some results
          printed. }

	  \item{\code{\dots}:}{ Arguments currently gobbled. }
		
    } %itemize
  } %item


  } %describe
} %methods

\keyword{methods}

