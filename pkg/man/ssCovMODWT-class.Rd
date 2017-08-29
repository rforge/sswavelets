\name{ssCovMODWT-class}
\Rdversion{1.1}
\docType{class}
\alias{ssCovMODWT-class}
%\alias{show,ssCovMODWT-method}
\alias{summary,ssCovMODWT-method}

\title{Class \code{"ssCovMODWT"}}

\description{

This class encompasses the essentials of two \code{\linkS4class{ssMODWT}} classes that share the same
decomposition attributes as found in the \code{levels} slot. The two surfaces produce a
set of covariances that are also contained in the object.

The \dQuote{ssMODWT} objects used to create the class must be from the same underlying
tract and population of trees or logs in the sense of a \pkg{sampSurf} simulation. The
difference in the two lies in applying different sampling methods or, say, plot sizes for
circular plot sampling as an example. The wavelet covariance decomposition is then based
on the two wavelet decompositions of the \dQuote{sampSurf} objects.

Examples and information concerning this class are found in the package vignette in the
references below. This vignette should be consulted for more details on creating and using
objects from this class.

} %description

\section{Objects from the Class}{

  Objects are created using the the constructor method \code{\link{ssCovMODWT}}.


}

\section{Slots}{ 

  Two \dQuote{ssMODWT} objects are used to construct this object. They are 
  variously referred to as the primary or first object, which appears in the
  \dQuote{default} fields corresponding to those in \dQuote{ssMODWT}
  below, and the second(ary) object. However, the
  primary often appears with an appended \sQuote{\code{a}} in the code to distinguish it from the
  second object, which always appears with an appended \sQuote{\code{b}} as seen below. I
  have tried to use \code{[a,b]} to distinguish the objects in the code and have stayed
  away from using \code{[x,y]} so as not confuse the latter's use in referring to the
  dimensions of the images; the exception being in \code{\link{covMODWT}} and some of the
  text in the vignette.

 
  \describe{
	
    \item{\code{ss.b}:}{Object of class \code{"sampSurf"}: The second sampling surface. }
%	
    \item{\code{ss.modwt.b}:}{Object of class \code{"list"}: The raw wavelet coefficient
        results from running a 2-D \acronym{MODWT} wavelet transformation using
        \code{\link{modwt.2d}} in package \pkg{waveslim} on the \code{ss.b} surface. }
%	
    \item{\code{ss.mra.b}:}{Object of class \code{"list"}: The results of applying a 2-D
        multi-resolution analysis (\acronym{MRA}) to the \code{ss.modwt.b}
        decomposition. }
%	
    \item{\code{covStats}:}{Object of class \code{"list"}: A succinct summary of the
        overall covariances and correlations. These are used when \code{show}ing or
        performing a \code{\link{summary}} on the object. }
%	
    \item{\code{ss.modwt}:}{Object of class \code{"list"}: See
        \dQuote{\code{\linkS4class{ssMODWT}}}; this is for the first surface in the
        \code{ss} slot.  }
%	
    \item{\code{vars.modwt}:}{Object of class \code{"list"}: The covariance results. This
        is a rather long list containing marginal summary, total, and image
        covariances. Details are provided in the vignette. Please note that this slot does
        not contain variances as in \dQuote{\code{\linkS4class{ssMODWT}}}, though the two
        have the same structure.  }
%	
    \item{\code{ss.mra}:}{Object of class \code{"list"}: See
        \dQuote{\code{\linkS4class{ssMODWT}}}; this is for the primary object.. }
%	
    \item{\code{levels}:}{Object of class \code{"list"}: See
        \dQuote{\code{\linkS4class{ssMODWT}}}. The information here should match exactly
        as checked in the constructor for both decomposition objects and is thus the
        shared result. }
%	
    \item{\code{description}:}{Object of class \code{"character"}: See
        \dQuote{\code{\linkS4class{ssWavelet}}}. }
%	
    \item{\code{wfName}:}{Object of class \code{"character"}: See
        \dQuote{\code{\linkS4class{ssWavelet}}}. }
%	
    \item{\code{ss}:}{Object of class \code{"sampSurf"}: See
        \dQuote{\code{\linkS4class{ssWavelet}}}. }
	
    } %describe
  } %slots

\section{Extends}{

Class \code{"\linkS4class{ssMODWT}"}, directly. \cr
Class \code{"\linkS4class{ssWavelet}"}, by class "ssMODWT", distance 2.

}

\section{Methods}{

No methods defined with class "ssCovMODWT" in the signature. However the
\sQuote{summary} method is available as a subclass of \dQuote{ssMODWT}.

}

%\references{
%%  ~~put references to the literature/web site here~~
%}

\author{
Jeffrey H. Gove
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

\dQuote{\code{\linkS4class{ssWavelet}}}, \dQuote{\code{\linkS4class{ssMODWT}}},
and package references: \code{\link[ssWavelets:ssWavelets-package]{ssWavelets}} and
 \code{\link[sampSurf:sampSurf-package]{sampSurf}}

}
\examples{
showClass("ssCovMODWT")
}
\keyword{classes}
