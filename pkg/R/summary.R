#---------------------------------------------------------------------------
#
#   The following are the summary methods for classes...
#
#   1. 'ssWavelets' virtual class
#   2. 'ssMODWT' subclass
#   3. 'ssCovMODWT' subclass
#
#   Returns...
#     summary information invisibly
#
#Author...									Date: 5-Oct-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#



#================================================================================
#  1. method for class ssWavelet & subclasses...
#
setMethod('summary',
          signature(object = 'ssWavelet'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of common items...
#------------------------------------------------------------------------------
    cat('\nObject of class:', class(object))
    .StemEnv$underLine(60)
    if(!is.na(object@description))
      cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')

    cat('\nWavelet filter:', object@wfName)


    units = object@ss@tract@units
    cat('\nUnits =', units)

    cat('\n')
    
    return(invisible())
}   #summary for 'ssWavelet'
) #setMethod




#================================================================================
#  2. method for subclass "ssMODWT"...
#
setMethod('summary',
          signature(object = 'ssMODWT'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'ssWavelet' method for 'ssMODWT'...
#------------------------------------------------------------------------------
    callNextMethod()

    cat('\nWavelet decomposition levels =', object@levels$J)
    cat('\nRows x columns =', object@levels$N, 'x', object@levels$M)
    
#
#   perhaps a little on total variance...
#
    if(!object@vars.modwt$isCovariance) {
      cat('\nsampSurf summary statistics...')
      cat('\n  Surface mean =', object@ss@surfStats$mean)
      cat('\n  Surface variance =', object@ss@surfStats$var)
    
      cat('\nWavelet (MODWT) summary statistics...')
      cat('\n  Wavelet mean =', object@vars.modwt$total$LL.mean)
      cat('\n  Wavelet sample variance =', object@vars.modwt$total$modwt.var)
    }
    cat('\n')
    
    return(invisible())
}   #summary for 'ssMODWT'
) #setMethod
    





#================================================================================
#  3. method for subclass "ssCovMODWT"...
#
setMethod('summary',
          signature(object = 'ssCovMODWT'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'ssMODWT' method for 'ssCovMODWT'...
#------------------------------------------------------------------------------
    callNextMethod()
    
    cat('\nsampSurf summary statistics...')
    cat('\n  Covariance =', object@covStats$ssCov)
    cat('\n  Correlation =', object@covStats$ssCor)
    cat('\n  --Surface "A" mean =', object@ss@surfStats$mean)
    cat('\n  --Surface "A" variance =', object@ss@surfStats$var)
    cat('\n  --Second "B" Surface mean =', object@ss.b@surfStats$mean)
    cat('\n  --Surface "B" variance =', object@ss@surfStats$var)
    
    cat('\nWavelet (MODWT) summary statistics...')
    cat('\n  Wavelet sample covariance =', object@vars.modwt$total$modwt.var)
    cat('\n  Wavelet sample correlation =', object@covStats$modwtCor)
    cat('\n  --Wavelet mean "A" =', object@vars.modwt$total$LL.mean)
    cat('\n  --Wavelet mean "B" =', object@vars.modwt$total$LL.mean.y)
    cat('\n')

    
    return(invisible())
}   #summary for 'ssCovMODWT'
) #setMethod
    
    
    
    

