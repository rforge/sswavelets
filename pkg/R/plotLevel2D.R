plotLevel2D = function(ssMODWT,
                       level = 1,
                       decompType = c('modwt', 'mra'),
                       type = c('raw', 'var'), #default is raw coefficients, variance is alternative
                       isoSmooth = TRUE,       #include LLJ in iso?
			           showPlot = TRUE,
			           showIZs = TRUE,
                       col = NA,          #color scheme for surfaces
                       boxCol = 'gray',   #outline color for "boxes" at each level
                       title = NA,        #use NA for default titles; '' for no title
                       showLegend = TRUE,
                       runQuiet = FALSE,
                       ...
                      )
{
#---------------------------------------------------------------------------
#
#   This routine will plot the MODWT or MODWT-MRA analysis as a matrix...
#     1st row:   LH|HH
#     2nd row:  ISO|HL
#   in one 2x2 figure. It does this using raster to piece the individual 
#   decompositions together. At level=J, there are five subfigures as LLJ is
#   also included centered in the 3rd row.
#
#   Note that 'var'+'mra' is not allowed as a valid combination--see
#   plotMODWR2D(), which is the workhorse here.
#
#**>Note: this uses plotMODWT2D for display of each of the decompositions...
#
#     LH = scaling-wavelet matrices V   (horizontal filter)
#     HH = wavelet-wavelet matrices W  (diagonal filter)
#     HL = wavelet-scaling matrices U  (vertical filter)
#     LL = residual average at J
#     ISO = isotropic: LH? + HL? + LH? (+LLJ if level==J)
#
#   Arguments...
#     ssMODWT = an object from ssMODWT()
#     level = the decomposition level 1:J
#     decompType = 'modwt' or 'mra' (multi-resolution via modwt)
#     type = 'raw' = the wavelet coefficients; 'var' = wavelet variance
#     isoSmooth = TRUE: include LLJ in the result for isometric at level J; FALSE: do
#                 not include the smooth component
#     showPlot = TRUE: display the plot; FALSE: don't
#     showIZs = TRUE: display the inclusion zones on the horizontal map; FALSE: none
#     col = see below, must be a function
#     boxCol = outline color for perimeter "boxes" at each level
#     title = use NA for default titles; '' for no title
#     showLegend = TRUE: display key around the maps; FALSE: no key
#     runQuiet = currently not used
#     ... = gobbled
#
#   Returns...
#     A raster object of the combined figure invisibly
#
#   Note that the legend is displayed on the screen as there is room to the outside
#   between the plot and the axes box. However, when exporting to a pdf file, this 
#   extra room disappears, so the resulting plot will have no labels.
#
#**>Aside: This is the first routine I have used withCallingHandlers() to suppress
#          warning messages. Evidently, suppressWarnings() does not work (I have used
#          it elsewhere with the same problem) when one is calling other routines
#          within the code block--we still get the message about the number of 
#          warnings generated in the calls, and this is irritating since it is always
#          the `"add" is not a graphical parameter' warning and is harmless. 
#          Note that returning invisible() does not work, but NULL does.
#**>Aside Revisited: Nope, withCallingHandlers does not work. And take a look
#                    at the code for suppressWarnings, it sets options(warn=-1)
#                    and uses withCallingHandlers() !
#
#Author...									Date: 9-May-2016
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

    decompType = match.arg(decompType)
    type = match.arg(type)

    #sampSurf inclusion zone class ==> type of sampling method used...
    izClass = class(ssMODWT@ss@izContainer@iZones[[1]])[1]

#
#   get the correct spatial extents of the sampling surface...
#
	ex = extent(ssMODWT@ss@tract)
    xmin = ex@xmin
    xmax = ex@xmax
    ymin = ex@ymin
    ymax = ex@ymax
    
    j = ssMODWT@levels$j
    J = ssMODWT@levels$J
    if(level < j[1] || level > J) #could allow level 0 if we make cumulative LH, HL, HH too sometime
      stop('Illegal level = ', level, ' for this routine!')
    tau.j = ssMODWT@levels$tau.j
  
#
#   the individual decompositions...
#
    r.lh = plotMODWT2D(ssMODWT, level=level, waveType='LH', decompType=decompType, type=type, 
                       isoSmooth=isoSmooth, showPlot = FALSE)
    r.hl = plotMODWT2D(ssMODWT, level=level, waveType='HL', decompType=decompType, type=type, 
                       isoSmooth=isoSmooth, showPlot = FALSE)
    r.hh = plotMODWT2D(ssMODWT, level=level, waveType='HH', decompType=decompType, type=type, 
                       isoSmooth=isoSmooth, showPlot = FALSE)
    r.iso = plotMODWT2D(ssMODWT, level=level, waveType='ISO', decompType=decompType, type=type, 
                        isoSmooth=isoSmooth, showPlot = FALSE)
    
#
#   merge them all into one figure compliments of raster...
#
    r = merge(r.lh$r, shift(r.hh$r, x = xmax - xmin))               #top row: horizontal | diagonal
    r = merge(r, shift(r.hl$r, x = xmax - xmin, y = -ymax + ymin))  #bottom right: vertical
    r = merge(r, shift(r.iso$r, y = -ymax + ymin))                  #bottom left: isotropic

    if(level == J) {
      r.ll = plotMODWT2D(ssMODWT, level, type=type, wave='LL', decompType=decompType, showPlot = FALSE)
      #need as.integer below to keep the origin correct in compareRaster (called from merge)...
      r = merge(r, shift(r.ll$r, x = xmin + as.integer((xmax - xmin)/2), y = -2*ymax + 2*ymin))
    }    

 
#
#   now plot the object if desired...
#
    if(any(is.na(col)))
      #we could put n.pal=100 and bias=5 as arguments above...  <********* consider this <**************
      col = palMODWT(100, bias=5, range=cellStats(r, range))
    
    if(showPlot) {
      withCallingHandlers({       #to suppress the '"add" is not a graphical parameter' warnings
      #suppressWarnings({          #to suppress the '"add" is not a graphical parameter' warnings
      plot(r, col=col)
      
      if(showIZs)
        plot(ssMODWT@ss@izContainer, add=TRUE, izColor=NA)

      #same dimension as the original, so take the easy way out for bbox...
      r.bb = bbox(r)
      rownames(r.bb) = c('x','y')             #or it will fail in bboxCheck() below
      plot(bboxToPoly(r.bb), add = TRUE)      #full bbox
      lines(c(xmin, 2*xmax - xmin), c(ymin, ymin))   #horizontal divider
      lines(c(xmax, xmax), c(-ymax+2*ymin, ymax))    #vertical divider
      #and if it's the last level, we have a little more work to do...
      if(level == J) {
        lines(c(xmin, 2*xmax - xmin), c(-ymax + 2*ymin, -ymax + 2*ymin))     #horizontal divider
        l.xmin = (xmax - xmin)/2 + xmin
        lines(c(l.xmin, l.xmin), c(-2*ymax + 3*ymin, -ymax + 2*ymin))        #left vertical divider
        r.xmin = l.xmin + xmax - xmin
        lines(c(r.xmin, r.xmin), c(-2*ymax + 3*ymin, -ymax + 2*ymin))        #right vertical divider
      }

      if(all(is.na(title))) {
        lower = paste('(level j=', level,  ', J=',J, ', tau.j=', tau.j[level], ')', sep='')
        subTitle = paste('(',izClass,': ',deparse(substitute(ssMODWT)),')',sep='')
        if(level %in% c(0,J))                 #show only for levels 0 or J
          subTitle = paste(subTitle, ' (isoSmooth = ',isoSmooth,')', sep='')
        if(type == 'raw' )
          title(paste(toupper(decompType), ' Raw Coefficients\n', lower, sep=''), sub=subTitle)
        else if(!ssMODWT@vars.modwt$isCovariance)
          title(paste(toupper(decompType), ' Wavelet Variance\n', lower, sep=''), sub=subTitle)
        else
          title(paste(toupper(decompType), ' Wavelet Covariance\n', lower, sep=''), sub=subTitle)
      }
      
      if(showLegend) {
        x.left = xmin
        x.right = 2*xmax - xmin 
        y.top = ymax/2
        y.bot = ymin - ymax/2
        if(level < J) {
          text(x.left, y.top, 'LH', pos = 2)          
          text(x.right, y.top, 'HH', pos = 4)
          text(x.left, y.bot, 'iso', pos = 2)
          text(x.right, y.bot, 'HL', pos = 4)
        }
        else {
          x.smooth = xmin + (xmax - xmin)/2
          y.smooth = ymin - 3*ymax/2
          text(x.left, y.top, 'horizontal', pos = 2)
          text(x.right, y.top, 'diagonal', pos = 4)
          text(x.left, y.bot, 'isotropic', pos = 2)
          text(x.right, y.bot, 'vertical', pos = 4)
          text(x.smooth, y.smooth, 'smooth', pos = 2)
        }
      }
      #}) #suppressWarnings
      }, warning = function(cond) return(NULL)) #withCallingHandlers
    } #showPlot

#
    if(!runQuiet) {
      cat('\nTop left: horizontal')
      cat('\nTop right: diagonal')
      cat('\nBottom right: vertical')
      cat('\nBotton left: isotropic')
      if(level == J)
        cat('\nBottom middle: final smooth')
      cat('\n')
    }
       
    return(invisible(r)
          )

}   #plotLevel2D
      
