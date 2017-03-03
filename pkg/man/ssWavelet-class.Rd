\name{ssWavelet-class}
\Rdversion{1.1}
\docType{class}
\alias{ssWavelet-class}
%\alias{show,ssWavelet-method}
%\alias{summary,ssWavelet-method}

\title{Class \code{"ssWavelet"}}

\description{
This is a virtual class that is the base class for other wavelet classes
based on \pkg{sampSurf}. It is very simply by design, and can easily be
extended to subclasses.
}

\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}

\section{Slots}{
  \describe{
    \item{\code{description}:}{Object of class \code{"character"}: Some
	  form of identifier. }

    \item{\code{wfName}:}{Object of class \code{"character"}: The
        wavelet filter name. This should be a legal name from
        \code{waveslim::wave.filter}; though for now it is restricted to
        \dQuote{haar} wavelets. }

    \item{\code{ss}:}{Object of class \code{"sampSurf"}: the decomposition is based on
        this object. }

  }
} %slots

\section{Methods}{
Virtual class: No methods defined with class "ssWavelet" in the signature.
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
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}

\examples{
showClass("ssWavelet")
}

\keyword{classes}
