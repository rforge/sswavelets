#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   'ssMODWT' class...
#
#Author...									Date: 17-Jan-2017
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jgove@fs.fed.us or jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   generic definition...
#
setGeneric('ssCovMODWT',  
           function(ssMODWT.a, ssMODWT.b, ...) standardGeneric('ssCovMODWT'),
             signature = c('ssMODWT.a', 'ssMODWT.b')
          )



          
#================================================================================
#
setMethod('ssCovMODWT',
          signature(ssMODWT.a = 'ssMODWT', ssMODWT.b = 'ssMODWT'),
function(ssMODWT.a,
         ssMODWT.b,          #must be identical to varMODWT -- a check
         description = 'sampSurf Covariance MODWT wavelet decomposition object',
         runQuiet = FALSE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   Creates a sample covariance object from two ssMODWT objects.
#
#   Note that we could add correlations here. Not so much on the images,
#   per se, but the summary and totals would be valuable. The totals
#   for correlation are computed below, these could also be done for summary.
#
#   Prototype.
#
#**>Note: Adapted to S4 25-Jan-2017, JHG.
#
#Author...									Date: 7-Sept-2016
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   check for missing and arguments' classes...
#
 #   if( missing(ssMODWT.a) || missing(ssMODWT.b) )
 #     stop('One or both ssMODWT objects is missing!')
#
  #  if(!is(ssMODWT.a, 'ssMODWT') || !is(ssMODWT.b, 'ssMODWT'))
  #    stop("Arguments ssMODWT.a & ssMODWT.b must be of class 'ssMODWT'!")
      
#
#   check for congruent dimensions, level of decompostion, etc....
#
    if(!all.equal(ssMODWT.a@levels, ssMODWT.b@levels))
      stop('The two "ssMODWT" objects must have the exact same @levels components')
      
    #all okay--dimensions...
    NM = ssMODWT.a@ss@surfStats$nc
      
    #raster objects, just multiply like a matrix & sum...
    ssCov = cellStats(ssMODWT.a@ss@tract * ssMODWT.b@ss@tract, sum)/NM - 
            ssMODWT.a@ss@surfStats$mean*ssMODWT.b@ss@surfStats$mean
    ssCor = ssCov/(ssMODWT.a@ss@surfStats$stDev * ssMODWT.b@ss@surfStats$stDev)

    covs = covMODWT(ssMODWT.a@ss.modwt, ssMODWT.b@ss.modwt)
    
    modwtCov = covs$total$modwt.var
    modwtCor = covs$total$modwt.var/(sqrt(ssMODWT.a@vars.modwt$total$modwt.var) *
                                     sqrt(ssMODWT.b@vars.modwt$total$modwt.var))
 
    if(!runQuiet) {   
      cat('\nsampSurf covariance =', ssCov)
      cat('\nsampSurf correlation =', ssCor)
      cat('\nWavelet sample covariance =', modwtCov)
      cat('\nWavelet sample correlation =', modwtCor, '(total)' )
    }
    cat('\n')  
    
    summaryStats = list(ssCov = ssCov,
                        ssCor = ssCor,
                        modwtCov = modwtCov, #also in vars.modwt$total$modwt.var
                        modwtCor = modwtCor
                       )
    
    #mimic an ssMODWT object for now...
    ssCov = new('ssCovMODWT',
                description = description,
                ss = ssMODWT.a@ss,
                ss.b = ssMODWT.b@ss,
                ss.modwt = ssMODWT.a@ss.modwt,
                vars.modwt = covs,
                ss.modwt.b = ssMODWT.b@ss.modwt,
                ss.mra = ssMODWT.a@ss.mra,
                ss.mra.b = ssMODWT.b@ss.mra,
                wfName = ssMODWT.a@wfName,
                levels = ssMODWT.a@levels,
                covStats = summaryStats
               )
    
    return(invisible(ssCov))

}   #ssCovMODWT constructor
)   #setMethod

                

