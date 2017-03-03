\name{ssMODWT-class}
\Rdversion{1.1}
\docType{class}
\alias{ssMODWT-class}
\alias{ssCovMODWT,ssMODWT,ssMODWT-method}

\title{Class \code{"ssMODWT"}}

\description{
  
  A subclass of virtual class \dQuote{\code{\linkS4class{ssWavelet}}} that can be used to
  create objects of that class. These objects contain the results from a two-dimensional
  maximal overlap discrete wavelet transform (\acronym{MODWT}) that has been applied to
  the \code{"sampSurf"} object stored in the object's \code{ss} slot.

  Examples and information concerning this class are found in the package vignette in the
  references below. This vignette should be consulted for more details on creating and
  using objects from this class.

}

\section{Objects from the Class}{

  Objects are created using the the constructor method \code{\link{ssMODWT}}.

}


\section{Slots}{

  In addition to the slots provided by the virtual superclass
  \dQuote{\code{\linkS4class{ssWavelet}}}, the following slots are
  represented\ldots
  
  \describe{

  \item{\code{ss.modwt}:}{Object of class \code{"list"}: The raw wavelet coefficient
      results from running a 2-D \acronym{MODWT} wavelet transformation using
      \code{\link{modwt.2d}} in package \pkg{waveslim}. }
	
    \item{\code{vars.modwt}:}{Object of class \code{"list"}: The variance results. This is
      a rather long list containing marginal summary, total, and image variances. Details
      are provided in the vignette.}
	
    \item{\code{ss.mra}:}{Object of class \code{"list"}: The results of applying a 2-D
        multi-resolution analysis (\acronym{MRA}) to the \code{ss.modwt} decomposition. }
	
    \item{\code{levels}:}{Object of class \code{"list"}: A list with information on the
        image dimensions and levels of the decomposition. Details are provided in the
        vignette.}
	
    \item{\code{description}:}{Object of class \code{"character"}: See
        \dQuote{\code{\linkS4class{ssWavelet}}}. }
	
    \item{\code{wfName}:}{Object of class \code{"character"}: See
        \dQuote{\code{\linkS4class{ssWavelet}}}. }
	
    \item{\code{ss}:}{Object of class \code{"sampSurf"}: See
        \dQuote{\code{\linkS4class{ssWavelet}}}. }
  }
  
} %slots

\section{Extends}{
Class \code{"\linkS4class{ssWavelet}"}, directly.
}

\section{Methods}{

  Methods other than \code{\link{ssMODWT}} for object creation that are
  applicable to this class include\ldots

  \describe{

    \item{ssCovMODWT}{\code{signature(ssMODWT.a = "ssMODWT", ssMODWT.b = "ssMODWT")}: 
        Covariance analysis of two \dQuote{ssMODWT} wavelet decompositions. }

    \item{summary}{\code{signature(object = "ssMODWT")}: Summary of the
      object. }

	 }
} %Methods

\references{
%%  ~~put references to the literature/web site here~~
}

\author{
Jeffrey H. Gove
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

\dQuote{\code{\linkS4class{ssWavelet}}}, \dQuote{\code{\linkS4class{ssCovMODWT}}},
and package references: \code{\link[ssWavelets:ssWavelets-package]{ssWavelets}} and
 \code{\link[sampSurf:sampSurf-package]{sampSurf}}

}

\examples{
showClass("ssMODWT")
}

\keyword{classes}
