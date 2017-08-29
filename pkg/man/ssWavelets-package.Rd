\name{ssWavelets-package}
\alias{ssWavelets-package}
\alias{ssWavelets}
\docType{package}

\title{
\packageTitle{ssWavelets}
}

\description{
  %\packageDescription{ssWavelets}

  This package adds several classes, generics and associated methods as well as a few
  various functions to help with wavelet decomposition of sampling surfaces generated using
  \pkg{sampSurf}. As such, it can be thought of as an extension to \pkg{sampSurf} for
  wavelet analysis.

  At present, the maximal overlap discrete wavelet transform (\acronym{MODWT}) as
  implemented in the \R package \pkg{waveslim} is currently supported. Because objects of
  class \dQuote{\code{\linkS4class{sampSurf}}} primarily contain results that are simply
  raster images, the two-dimensional routines from \pkg{waveslim} are used exclusively
  (though the 2-D pyramid algorithm employs the 1-D routines behind the scenes). Also,
  because one is interested in working with the full image rather than compression,
  methods such as the discrete wavelet transform (\acronym{DWT}) seem to make little sense
  on these data. However, other wavelet filters may prove useful in the future, and the
  package has been designed for easy extension to other methods.

  Much of the package is geared towards decomposition of the sample variance; that is, the
  total variance of the sampSurf raster object. The idea is to decompose the sample
  (sampSurf) variance by scale (distance) to see how sampling methods compare at different
  scales. This distance-based decomposition gives an idea of the variation at scales that
  is unavailable in the overall sampSurf variance.

  The \dQuote{raw} wavelet decomposition is available in the objects, but its use seems
  limited in the simulated sampling context at this point. Plotting methods are available
  for both, but analysis functions are geared towards variance decomposition.

  A rudimentary S4 class structure has been established for extending the current classes
  should one desire to do so.

} %description


\details{

%The DESCRIPTION file:
%\packageDESCRIPTION{ssWavelets}
%\packageIndices{ssWavelets}
  
\tabular{ll}{
Package: \tab ssWavelets \cr
Type: \tab Package \cr
Version: \tab 0.1-2 \cr
Date: \tab 29-Aug-2017 \cr
License: \tab GPL (>= 3) \cr
LazyData: \tab TRUE \cr
}


A list of the main resources are described in the following sections. Please refer to the
links below for more information. The package vignette provides further details and
examples on how to use this package.

} %details


%--------------------------------------------------------------------------
%
%  S4 classes section...
%
\section{Classes For Use In \pkg{ssWavelets}\ldots}{ The following S4
  classes are available.


%  
  \subsection{The \dQuote{ssWavelet} class\ldots}{

  Objects of the non-virtual classes below can be created using
  constructor functions of the same name, see the \emph{Object
  Constructors} section for details.

  \tabular{ll}{
	\code{\linkS4class{ssWavelet}} \tab Virtual base class for the following\cr
	\code{\linkS4class{ssMODWT}} \tab A class for \acronym{MODWT} decomposition\cr
	\code{\linkS4class{ssCovMODWT}} \tab A class for covariance decomposition\cr
  }
  
  } %ssWavelet

} %S4 Classes





%--------------------------------------------------------------------------
%  constructors section...
%
\section{Object Constructors}{
  
  For each of the non-virtual classes defined in the table above, we must be able to
  create objects that can be used in \R. This is done using class-specific
  \dQuote{constructor} methods that take the drudgery away from creating what can often be
  somewhat complicated (with all the graphical components) new object
  instances. Eventually, there may be more than one constructor for a given class of
  object, and these are differentiated by the method signature; see the links provided
  below for more details.

%
  \subsection{\dQuote{ssWavelet} class constructors\ldots}{

  \tabular{ll}{
	\code{\link{ssMODWT}} \tab Constructor for individual
	\code{\linkS4class{ssMODWT}} objects\cr
%
	\code{\link{ssCovMODWT}} \tab Constructor for individual
	   \code{\linkS4class{ssCovMODWT}} objects\cr   
  }

  } %Constructors subsection
} %Constructors section






%--------------------------------------------------------------------------
%  summary & plotting section...
%
\section{Summary and Plotting Methods}{

  The objects created above have graphical content made possible by
  using classes from the \code{\link{sp}} or \code{\link{raster}}
  packages. Therefore, functions have been developed to allow for
  graphical display of these objects to show different aspects of the
  wavelet decomposition (i.e., raw coefficients, variance, covariance).

  In addition, \code{\link{summary}} and \code{\link{show}} methods are
  also available (sometimes producing the same result) for printing
  information within the object\ldots

  \tabular{ll}{
	
	\code{\link{plotMODWT2D}} \tab  Plots the \acronym{MODWT} object at
	a given level and anisotropic decomposition\cr
%	
	\code{\link{plotLevel2D}} \tab  Plots the set of \acronym{MODWT}
	object decompositions at chosen level\cr
%	
%	\code{\link{plot3D}} \tab  To plot \dQuote{sampSurf} objects using
%	package \pkg{rgl}; can also be used for \dQuote{ssWavelet} object images\cr
%	
	\code{\link{show}} \tab  To succinctly print one of the package objects\cr

	\code{\link{summary}} \tab  To summarize one of the package objects\cr
  } %tabular

  Note that \code{\link{plot3D}} is used to plot \dQuote{raster} and \dQuote{sampSurf} objects using
  package \pkg{rgl}. However, it can also be used for \dQuote{ssWavelet} object images
  by first \sQuote{casting} the desired matrix representation to an object of class 
  \dQuote{\code{\linkS4class{RasterLayer}}}. Casting is necessary as there is no method
  written specifically for objects of class \dQuote{ssWavelets}, since the number of
  possible images in these objects is large. Examples of this procedure are given in the
  vignette. 
  
} %summary and plotting







%--------------------------------------------------------------------------
%  miscellaneous functions section...
%
\section{Other Functions}{
  A few other functions that will not normally be required listed
  below. These include some plotting methods in addition to those listed
  above\ldots

  \tabular{ll}{
	\code{\link{covMODWT}} \tab Much like the \R \pkg{stats} routine \link{cov}, but for wavelets\cr
%	
	\code{\link{hfsMODWT}} \tab For comparing different sampling methods
	by scale\cr
%	
	\code{\link{hfsPlot}} \tab Plots the results of
	\code{\link{hfsMODWT}}: An H. F. Smith plot\cr
%	
	\code{\link{varPlot}} \tab A variance plot using results from \code{\link{hfsMODWT}}\cr
%	
	\code{\link{palMODWT}} \tab A nice palette function for display\cr
%	
	\code{\link{ssEnergy}} \tab Average \dQuote{energy} for a
	\dQuote{sampSurf} object\cr
%	
	\code{\link{ssToroid}} \tab Returns the toroidal boundary correction object\cr
%	
	\code{\link{ssReflect}} \tab Returns the reflection (periodization)
	boundary correction object\cr
  } %tabular
  
} %other functions





%--------------------------------------------------------------------------
%  data sets section
%
\section{Data Sets}{
  Data sets included in the package\ldots

  \tabular{ll}{
	\code{\link{befp34}} \tab Mapped tree data from the Bartlett
	Experimental Forest. \cr

    } %tabular
  
} %data sets






\author{
  %\packageAuthor{ssWavelets}
Author & Maintainer: Jeffrey H. Gove, <jgove@fs.fed.us>

%Maintainer: \packageMaintainer{ssWavelets} %works in html, not pdf!
%Maintainer: Jeffrey H. Gove, <jgove@fs.fed.us>
}

\references{
\describe{
%
\item{}{Gove, J. H. 2017. A User's Guide to the ssWavelets Package.
  \url{http://sswavelets.r-forge.r-project.org/}.}
%
\item{}{Gove, J. H. 2017. Some Refinements on the Comparison of Areal
Sampling Methods via Simulation. \emph{Forests}, Submitted.}
} %describe  

} %references

%see hsearch_db_keywords() for a list of keywords...
\keyword{ package }

\seealso{
 \code{\link[sampSurf:sampSurf-package]{sampSurf}}
}

%\examples{
%Please see the links above and the vignette for examples.
%}
