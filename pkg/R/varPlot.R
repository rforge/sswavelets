varPlot = function(...,
                   sampMeth = c('Meth.A', 'Meth.B'),
                   groupSM = TRUE,
                   units = c('metric', 'English'),
                   isoSmooth = TRUE,       #include LLJ.var in iso?
                   showPlot = TRUE,
                   fileName = '',
                   #lattice...
                   ylab = 'Isotropic Variance',
                   xlab = 'Distance', 
                   scales = 'same',
                   type = 'b',
                   pch = 19,
                   as.table = TRUE,
                   theme = c('plain','custom','ggplot','economist'),
                   #other...
                   runQuiet = FALSE
                  )
{
#---------------------------------------------------------------------------
#
#   This routine creates a plot of the marginal isotropic variances against
#   distance/scale conditioned on...
#       (i) each of the different sampling methods, one panel for each 
#           method+size (i.e., baf, plot size) combination: groupSM=FALSE
#       (ii) conditioned on just the sampMethod level with each of the 
#           "sizes" grouped: groupSM=TRUE
#   So using (i) could conceivably get rather busy panel-wise.
#
#**>Note: If groupSM=TRUE then we will have as many panels as sampling methods;
#         if it is FALSE, then the number of panels is as above in (i).
#
#**>Note: If we want to look at a plot of just the isotropic wavelet variances
#         without the smooth/scale component (LL.var) at level J_0, then one
#         can use isoSmooth=FALSE. This will be useful in calculating the 
#         fractal dimension if I ever get that far. (19-Dec-2016, JHG)
#
#   Arguments...
#     ... = any number of data frames from hfsMODWT
#     sampMeth = a vector of sampling method identifiers in whatever form
#                that will five a label useful for conditioning to each
#                of the data frame objects in ...
#     groupSM = TRUE: condition on sampMeth and group on id; FALSE: condition on
#               id and do not group
#     units = the appropriate units that match the underlying sampling surface
#     isoSmooth = TRUE: include LLJ in the result for isotropic at level J; FALSE: do
#                 not include the smooth component
#     showPlot = TRUE: display the plot; FALSE: just return &/or file a copy
#     fileName = '': no output file; character==>hardcopy
#     ylab, xlab, scales, type, pch, as.table are all either graphics parameters
#       or lattice options--see some examples below
#     theme = see the themes in latticeExtra
#     runQuiet = TRUE: sssssh; FALSE: verbose!
#
#   Returns...
#     a list invisibly with...
#       -- the concatenated data frames
#       -- the lattice object
#
#   Lattice options...
#
#   For example, we can make the scales logarithmic using either...
#       scales=list(x = list(log = 2), y=list(log=2)) 
#   for both, or
#       scales=list( y=list(log=2))
#   for just y
#
#Author...									Date: 6-Oct-2016
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   a few checks and set units...
#
    units = match.arg(units)
    
    if(!runQuiet) {
      cat('\nunits =', units)
      cat('\nisoSmooth =', isoSmooth)
    }
    if(units == 'metric')
      dist = 'm'
    else
      dist = 'ft'
      
      

#
#   how many hfsMOWDT lists were passed? ...
#
    args = list(...)
    n.args = length(args)
    arg.names = sapply(substitute(list(...)), deparse)[-1]

#
#   now assign the correct sampling method to each data frame and concatenate...
#
    if(length(sampMeth) < n.args)
      stop('number of ... args > sampMeth length!')
      
    df = NULL
    for(i in 1:n.args) {
      x = args[[i]]
      x$sm = sampMeth[i]
      df = rbind(df, x)
    }
    n = nrow(df)
    rownames(df) = 1:n
    
#
#   if we want to look at just the isotropic wavelet variance, we must subtract
#   total$LL.var off from the isoJ variance...
#
    if(!isoSmooth) {
      J = max(df$j)
      isoJ = paste('iso',J, sep='')
      idx = df[, 'iso.j'] == isoJ 
      LL.var = df[idx, 'avg.E'] - df[idx, 'waveVar']
      isoJ.var = df[idx, 'iso.var'] - LL.var
      df[idx, 'iso.var'] = isoJ.var
    }
    
      
#
#   set the theme as desired...
#
    #require(latticeExtra, quietly=TRUE)
    theme = match.arg(theme)
    lattice.options = lattice.options()
    switch(theme,
           plain = {par.settings = trellis.par.get()},
           custom = {par.settings = custom.theme.2()},
           ggplot = {par.settings = ggplot2like()
                     lattice.options = ggplot2like.opts()
                    },
           economist = {par.settings = theEconomist.theme()
                        lattice.options = theEconomist.opts()
                       }
          )

#
#   finally the plot...
#
    if(groupSM) {
      theFormula = formula(iso.var ~ tau.j | sm)
      groups = with(df, id)
      key = list(x = .05, y = .8, corner = c(0, 0))
    }
    else {
      theFormula = formula(iso.var ~ tau.j | id)
      groups = NULL
      key = NULL
    }
    
    plt = xyplot(theFormula, data=df, scales=scales, 
                 groups = groups, auto.key = key,
                 type=type,
                 as.table = as.table,
                 ylab = ylab,
                 xlab = xlab,
                 pch = pch, 
                 par.settings = par.settings,
                 lattice.options = lattice.options
                 #...
                )

    if(showPlot)
      print(plt) 
      
    if(nchar(fileName) > 0) {
      if(exists(hardcopyLattice))
        hardcopyLattice(plt, fileName, ...) 
      else 
        message('Function hardcopyLattice unavailable, please print from the returned object.')
    }
                 

    if(!runQuiet)
      cat('\n')           
    
    return(invisible(list(df=df,
                          plt=plt
                         )
                    )
          )
}   #varPlot
    
    

