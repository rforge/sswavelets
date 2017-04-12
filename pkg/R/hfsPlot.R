hfsPlot = function(...,
                   sampMeth = c('Meth.A','Meth.B'),
                   conditionOn = c('iso.j', 'tau.j', 'j'),
                   units = c('metric', 'English'),
                   showPlot = TRUE,
                   fileName = '',
                   ylab = 'Isotropic Variance',
                   xlab = 'Average Inclusion Zone Area', 
                   scales = 'same',
                   type = 'b',
                   pch = 19,
                   as.table = TRUE,
                   theme = c('plain','custom','ggplot','economist'),
                   runQuiet = FALSE
                  )
{
#---------------------------------------------------------------------------
#
#   This routine will generate a lattice plot for the marginal isotropic
#   wavelet variance in the form of a H. F. Smith plot; that is, with
#   the average inclusion zone area increasing on the x axis for the different
#   sampling methods.
#
#   The routine assumes that you have already got one or more data frames
#   from hfsMODWT, which takes an ssMODWT object (list) and exports the 
#   appropriate data frame--in this case using the long form.
#
#**>Note: we should probably add a key to the plot for the groups variable;
#         i.e., sampling method.
#
#   Arguments...
#     ... = any number of data frames from hfsMODWT
#     sampMeth = a vector of sampling method identifiers in whatever form
#                that will give a label useful for conditioning to each
#                of the data frame objects in ...
#     conditionOn = one of the factor or character columns in the data
#                   frame; for tau.j and j, they will be converted appropriately
#     units = the appropriate units that match the underlying sampling surface
#     showPlot = TRUE: display the plot; FALSE: just return &/or file a copy
#     fileName = '': no output file; character==>hardcopy
#     ylab, xlab, scales, type, pch, as.table are all either graphics parameters
#       or lattice options
#     theme = see the themes in latticeExtra
#     runQuiet = TRUE: sssssh; FALSE: verbose!
#
#   Returns...
#     a list invisibly with...
#       -- the concatenated data frames
#       -- the lattice object
#
#Author...									Date: 5-Oct-2016
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
    conditionOn = match.arg(conditionOn)
    
    if(!runQuiet) {
      cat('\nunits =', units)
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
    
    if(conditionOn == 'tau.j') {
      #df$tau.j = paste(df$tau.j, dist, sep='')
      tau.j = df$tau.j
      tau.labs = unique(paste(sort(tau.j), dist, sep=''))
      tau.j = paste(tau.j, dist, sep='')
      df$tau.j = factor(tau.j, unique(tau.j), labels=tau.labs, ordered=TRUE)
    }
    else if (conditionOn == 'j')
      df$j = paste('j', df$j, sep='.')
      #df$j = factor(df$j, unique(df$j), ordered=TRUE) #works fine too
      
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
    theFormula = formula(paste('iso.var ~ izArea|', conditionOn))
    plt = xyplot(theFormula, data = df, scales = scales, 
                 groups = sm, auto.key = list(x = .85, y = .85, corner = c(0, 0)),  
                 type = type,
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
      
#
#   hard copy only in PDF if no hardcopyLattice available...
#
    if(nchar(fileName) > 0) {
      if(exists(hardcopyLattice))
        hardcopyLattice(plt, fileName) 
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
}   #hfsPlot
    
    

