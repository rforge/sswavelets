plotMODWT2D = function(ssMODWT,
                       level = 1,
                       waveType = c('ISO', 'LH', 'HL', 'HH', 'LL'),
                       decompType = c('modwt', 'mra'),
                       type = c('raw', 'var'), #default is raw coefficients, variance is alternative
                       isoSmooth = TRUE,       #include LLJ in iso?
			           showPlot = TRUE,
 			           showIZs = TRUE,
 			           addLattice = FALSE,
                       col = NA,          #color scheme for surfaces
                       boxCol = 'gray',   #outline color for "boxes" at each level
                       title = NA,        #use NA for default titles; '' for no title
#                       runQuiet = FALSE,
                       ...
                      )
{
#---------------------------------------------------------------------------
#
#   This routine will plot any of the individual decompositions below at
#   any level either as the raw wavelet coefficients or the variance. It is
#   for the objects generated from ssMODWT().
#
#   Note that in order to make life much much simpler, the results from ssMODWT
#   are converted to raster objects to take advantage of its algebra and plotting.
#
#   Note that 'var'+'mra' is not allowed as a valid combination.
#
#   Note that level=0, waveType="ISO", type="var" will plot the isometric 
#   variance surface matrix in isoSurf.var. JHG 19-July-2016
#
#   Note that type="var" calculates E[X^2] for any matrix (including isoSurf.var)
#   and does not subtract off the mean, or negative variances would result because
#   the individual components at each pixel are so small relative to the mean.
#
#   Raw MRA gives exact alginment with the original ('raw'+'mra'). This is most
#   conveniently demostrated at level one, where the wavelet filter picks out the 
#   boundaries of the inclusion zones. With raw+mra, they will align very
#   well with the actual plotted inclusion zones, with raw+modwt, you will
#   notice the misalignment. 
# 
#**>Note: we have the following for waveType...
#
#     LH = scaling-wavelet matrices V   (horizontal filter)
#     HH = wavelet-wavelet matrices W   (diagonal filter)
#     HL = wavelet-scaling matrices U   (vertical filter)
#     LL = scaling-scaling matrix Z     (residual average at J)
#     ISO = isotropic: LH? + HL? + LH?  (+LLJ if level==J) [computed here]
#
#   Arguments...
#     ssMODWT = an object from ssMODWT()
#     level = the decomposition level 1:J; zero means plot the aggregate variance surface
#     waveType = see above
#     decompType = 'modwt' or 'mra' (multi-resolution via modwt)
#     type = 'raw' = the wavelet coefficients; 'var' = wavelet variance
#     isoSmooth = TRUE: include LLJ in the result for isotropic at level J; FALSE: do
#                 not include the smooth component
#     showPlot = TRUE: display the plot; FALSE: don't
#     showIZs = TRUE: display the inclusion zones on the horizontal map; FALSE: none
#     addLattice = TRUE: adds a levelplot() to the return list; FALSE: NULL
#     col = see below, must be a function
#     boxCol = outline color for perimeter "boxes" at each level
#     title = use NA for default titles; '' for no title
#     runQuiet = currently not used
#     ... = gobbled
#
#   Returns...
#     a list invisibly with...
#     -- the final raster object from the display
#     -- a levelplot lattice object version if addLattice=TRUE
#
#**>Note that it would be better to use a logarithmic color scale than
#   to take the logarithm of the variance values. See if there is this
#   option in raster or levelplot.
#
#   Update: This is version 2, and reflects the changes made to ssMODWT when
#           varMODWT was written to pre-compute the variances; 25-Aug-2016, JHG
##
#Author...									Date: 8-Apr-2016
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   check the first argument's class...
#
    if(!is(ssMODWT, 'ssMODWT'))
      stop("Argument ssMODWT must be of class 'ssMODWT'!")

    #sampSurf inclusion zone class ==> type of sampling method used...
    izClass = class(ssMODWT@ss@izContainer@iZones[[1]])[1]

    type = match.arg(type)
    waveType = match.arg(waveType)
    decompType = match.arg(decompType)

    if(decompType == 'mra' && type == 'var')
      stop('Variance plot only valid with modwt, not mra')
      
    isCovariance = ssMODWT@vars.modwt$isCovariance
    if(type == 'raw' && isCovariance)
      stop("type='raw' is not available for covariance objects.") 


    j = ssMODWT@levels$j
    J = ssMODWT@levels$J
    tau.j = ssMODWT@levels$tau.j
    M = ssMODWT@levels$M
    N = ssMODWT@levels$N
    Lj = ssMODWT@levels$Lj
    
    NM = N*M              #ncells in the raster
    
    x.min = xmin(ssMODWT@ss@tract) 
    x.max = xmax(ssMODWT@ss@tract)
    y.min = ymin(ssMODWT@ss@tract)
    y.max = ymax(ssMODWT@ss@tract)

#
#   allow level==0 for ISO surface variance...
#
    if(level == 0 && (waveType != "ISO" || type != "var"))
      stop("Total surface variance (at level 0) is only available for 'ISO' & 'var'!")

    if(level > J || level < 0)
      stop('level = ',level,' must be in [0, J=',J,']')

#
#   MODWT decomposition on the sampling surface or MODWT+MRA decompostion...
#
    if(decompType == 'modwt')    
      ss.modwt = ssMODWT@ss.modwt
    else
      ss.modwt = ssMODWT@ss.mra
      
#
#   pull out the image variances slot to make things simpler...
#
    im.vars = ssMODWT@vars.modwt$image                       


#
#   put this into raster() format for ease of further manupulation (plotting, etc.);
#   we must specify the extents of the images, otherwise using simply, e.g., 
#   r.hl = raster(ss.modwt[[HL]]) will scale the image to [0,1]...
#
    if(level > 0 && type == 'raw') {
      LH = paste('LH', level, sep='')
      r.lh = raster(ss.modwt[[LH]], x.min, x.max, y.min, y.max)
      HL = paste('HL', level, sep='')
      r.hl = raster(ss.modwt[[HL]], x.min, x.max, y.min, y.max)
      HH = paste('HH', level, sep='')
      r.hh = raster(ss.modwt[[HH]], x.min, x.max, y.min, y.max)
    
      #include the smoothing from the last level if need be...
      if(level == J) {
        LL = paste('LL', level, sep='')
        r.ll = raster(ss.modwt[[LL]], x.min, x.max, y.min, y.max)
      }
    }
    else  #works with covMODWT() results...
      r.isoSurf = raster(im.vars$isoSurf.var, x.min, x.max, y.min, y.max)

    iso.j = paste('iso', level, sep='')
    
#
#   set up the raster for the plot...
#
    r = switch(waveType,
               ISO = {
                      if(type == 'raw') {                    #raw
                        r = r.lh + r.hl + r.hh
                        if(level == J && isoSmooth)          #include the smooth 
                          r = r + r.ll
                        else                                 #if() returns NULL as default else
                          r
                      }
                      else {                                 #variance
                        if(level == 0) {                     #total surface variance
                          r = r.isoSurf
                        }
                        else {                               
                          r = raster(im.vars$iso.var[[iso.j]], x.min, x.max, y.min, y.max)
                        }
                        #subtract out the smooth if desired from levels that contain it...
                        if(level %in% c(0,J) && !isoSmooth)
                          r = r - raster(im.vars$LL.var, x.min, x.max, y.min, y.max)
                        else                                 #if() returns NULL as default else
                          r
                      } #variance
                    },
               LH = {
                      if(type == 'raw')
                        r = r.lh
                      else              #variance
                        r = raster(im.vars$LH.var[[level]], x.min, x.max, y.min, y.max)
                    },
               HL = {
                      if(type == 'raw')
                        r = r.hl
                      else              #variance
                        r = raster(im.vars$HL.var[[level]], x.min, x.max, y.min, y.max)
                    },
               HH = {
                      if(type == 'raw')
                        r = r.hh
                      else              #variance
                        r = raster(im.vars$HH.var[[level]], x.min, x.max, y.min, y.max)
                    },
               LL = {
                     if(level < J)
                       stop('LL is only available at level J = ',J)
                     if(type == 'raw')
                        r = r.ll
                      else              #variance
                        r = raster(im.vars$LL.var, x.min, x.max, y.min, y.max)
                    }
             )

  

#----------------------------------------
#   now plot the object if desired...
#

    if(any(is.na(col)))
      #we could put n.pal=100 and bias=5 as arguments above...  <********* consider this <**************
      col = palMODWT(100, bias=5, range=cellStats(r, range))
    
    if(showPlot) {
      #plot(r, col=.StemEnv$blue.colors(1000), asp=1) #raster needs aspect here to make it square
      #
      #plot it...
      suppressWarnings({      #to suppress the '"add" is not a graphical parameter' warnings
        plot(r, col=col, asp=1) #raster needs aspect here to make it square

        if(showIZs) { #add the inclusion zones if desired...
          plot(ssMODWT@ss@izContainer, add=TRUE, izColor=NA)
          if(isCovariance)
            plot(ssMODWT@ss.b@izContainer, add=TRUE, izColor=NA, lty='dashed') #for now, dash to distinguish
        }
          

        #same dimension as the original, so take the easy way out for bbox...
        plot(perimeter(ssMODWT@ss), add=TRUE, border = boxCol)
        if(all(is.na(title))) {
          subTitle = paste('(',izClass,': ',deparse(substitute(ssMODWT)),')',sep='')
          if(level %in% c(0,J) && waveType == 'ISO')                 #show only for levels 0 or J
            subTitle = paste(subTitle, ' (isoSmooth = ',isoSmooth,')', sep='')
          if(level > 0)
            lower = paste('(level j=', level,  ', J=',J, ', tau.j=', tau.j[level], ')', sep='')
          else
            lower = paste('(level j=', level,  ', J=',J, ')', sep='')
          if(type == 'raw' )
            title.type = ' Raw coefficients\n'
          else if(!isCovariance)
            title.type = ' Wavelet variance\n'
           else
            title.type = ' Wavelet covariance\n'
            
          title = paste(toupper(decompType), ' ', waveType, title.type, lower, sep='')  
          title(title, sub = subTitle)
        }
        else
          title(title)
      }) #suppressWarnings
    } #showPlot

#
#   we can always return the level plot version; note that just the main title is used here...
#
    if(addLattice && requireNamespace('rasterVis'))
      plt = rasterVis::levelplot(r, col.regions = col, main = title)
    else
      plt = NULL
    
       
    return(invisible(list(r = r,
                          plt = plt
                         )
                    )
          )

}   #plotMODWT2D
                         
