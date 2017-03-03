ssEnergy = function(ssMODWT,
                    showPlot = TRUE,
                    showIZs = TRUE,
                    boxCol = 'gray', 
                    col = NA,  
                    title = NA,  
                    ...
                   )
{
#---------------------------------------------------------------------------
#
#   This calculates and returns the "average energy" sampSurf surface from 
#   the sampSurf MODWT object. This can be compared against the wavelet  
#   decomposition as, e.g.,...
#
#     M = nrow(ssMODWT@ss@tract)
#     N = ncol(ssMODWT@ss@tract)
#     plot3D(raster(ssMODWT@vars.modwt$image$isoSurf.var, 0,M, 0,N), 
#            col=palMODWT(100))
#
#**>Note: I added the ability to handle covariance objects. 31-Oct-2016
#
#**>Note: Updated to S4 ssWavelets package, 20-Jan-2017, JHG.
#
#   Arguments...
#     ssMODWT = an object from ssMODWT() or descendent class ssCovMODWT()
#     showPlot = TRUE: display the plot; FALSE: don't
#     showIZs = TRUE: display the inclusion zones on the horizontal map; FALSE: none
#     boxCol = outline color for perimeter "boxes" at each level
#     col = see below, must be a vector palette of colors
#     title = use NA for default titles; '' for no title
#     ... = gobbled
#
#   Returns...
#     the raster average energy image
#
#Author...									Date: 25-Oct-2016
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   check for compatibility...
#
    if(!is(ssMODWT@ss, 'sampSurf'))     #ssMODWT object passed
      stop('ssMODWT@ss must be of class "sampSurf"')
      
    isCovariance = ssMODWT@vars.modwt$isCovariance
    if(isCovariance) {
      if(!is(ssMODWT@ss.b, 'sampSurf'))     #two ssMODWT objects passed
        stop('ssMODWT@ss.b must be of class "sampSurf"')
      ss.b = ssMODWT@ss.b
    }
      
    ss = ssMODWT@ss

   
    #sampSurf inclusion zone class ==> type of sampling method used...
    izClass = class(ss@izContainer@iZones[[1]])[1]
      
    M = nrow(ss@tract)
    N = ncol(ss@tract)
    NM = N*M
   
#
#   average energy -- all there is to it!...    
#
    if(!isCovariance) {
      ae = ss@tract^2/NM                       #a raster object
      xtraTitle = ''
    }
    else {
      ae = ss@tract * ss.b@tract / NM          #a raster object
      xtraTitle = '\n(covariance)'
    }
    
#
#   plot it if desired...
#
    if(showPlot) {
      if(any(is.na(col))) {
        #dwt.colors = colorRampPalette(c('whitesmoke', 'lightskyblue2', 'saddlebrown'), bias=5)
        #col = dwt.colors(100)
        #we could put n.pal=100 and bias=5 as arguments above...  <********* consider this <**************
        col = palMODWT(100, bias=5, range=cellStats(ae, range))
      }

      plot(ae, col=col, asp=1)
      if(showIZs) { #add the inclusion zones if desired...
        plot(ssMODWT@ss@izContainer, add=TRUE, izColor=NA)
        if(isCovariance)
          plot(ss.b@izContainer, add=TRUE, izColor=NA, lty='dashed') #for now, dash to distinguish
      }
      
      #same dimension as the original, so take the easy way out for bbox...
      plot(perimeter(ssMODWT@ss), add=TRUE, border = boxCol)

      if(all(is.na(title))) {
        subTitle = paste('(',izClass,': ',deparse(substitute(ssMODWT)),')',sep='')
        theTitle = paste('sampSurf Average Energy', xtraTitle)
        title(theTitle, sub = subTitle)
      }
      else
        title(title)

    } #showPlot
    
    return(invisible(ae))
}   #ssEnergy
        
      
      

