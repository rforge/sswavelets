hfsMODWT = function(...,
                    ids = NA,
                    long = TRUE,
                    runQuiet = TRUE
                   )
{
#---------------------------------------------------------------------------
#
#   This routine takes a number of ssMODWT result objects and sets up a
#   data frame with the isotropic variances for each level in the decomposition.
#   It also includes total variance and total energy and surface mean, all
#   from the decomposition.
#
#   This is helpful in further analysis of the H. F. Smith type, either by
#   model estimation, plotting, or both.
#
#   Selecting the inclusion zone area only works with non-pps methods, so we
#   need to take an average over all trees/logs perhaps for pps methods?
#   Or, perhaps the sum of the inclusion areas over all stems, since this
#   reduces correctly in the non-pps methods.
#
#   Update: Note that the average inclusion zone area is used by Wensel
#           and John. 1969. A statistical procedure for combining... FS:307-317
#           so that is what I adopted below (2-May-2016).
#
#   Update: This is version 2, and reflects the changes made to ssMODWT when
#           varMODWT was written to pre-compute the variances; 25-Aug-2016, JHG
#
#   Update: I made some changes to the format, including adding the ids and long
#           arguments to facilitate use with lattice plotting; 29-Sept-2016, JHG
#           See also hfsPlot(); 5-Oct-2016, JHG
#
#**>Note: Updated to S4 ssWavelets package, 23-Jan-2017, JHG.
#
#   Arguments...
#     ... = any number of ssMODWT objects for different plot sizes;
#     ids = simplified ids for each of the ... objects, otherwise, their
#           names will be used; these will correspond to some measure of the
#           size of the inclusion zone (e.g., for circular plots, these might
#           be c('cp.8', 'cp.10') for radii of 8 & 10m (ft)
#     long = TRUE: return a df in long format; FALSE: wide format
#     runQuiet = TRUE: no feedback; FALSE: a short report on the results
#
#   Returns...
#     a data frame invisibly
#
#   Examples...
#
#     1. Long form conditional plots from lattice...
#
#        These plots condition on the "iso.j" column, which is the names for the
#        different "iso" levels. Rather than this, one could use the actual levels
#        j, or the scales, tau.j, would even be better. Note that as of 5-Oct-2016
#        the function hfsPlot() can be used to make plots of the following types
#        with other options.
#
#     1.1 The first uses 'circularPlotIZ' with different (metric) radii...
#
#          sscp.8 = makeSS('ss.3tree', plotRadius=8)
#          sscp.10 = makeSS('ss.3tree', plotRadius=10)
#          sscp.16 = makeSS('ss.3tree', plotRadius=16)
#          mocp.8 = ssMODWT(sscp.8, J=4)
#          mocp.10 = ssMODWT(sscp.10, J=4)
#          mocp.16 = ssMODWT(sscp.16, J=4)
#          jl = hfsMODWT(mocp.10, mocp.16, mocp.8, ids=c('10','16','8'), long=TRUE)
#          xyplot(iso.var ~ izArea|iso.j, jl, scales='free')
#   
#        where the first 3 arguments in hfsMODWT correspond to ssMODWT objects
#        with different plot radii. The second argument gives simple ids. 
#        Remember that these are all in metric units, so a 16m plot radius is
#        essentially a 1/5th-acre plot.
#
#     1.2 This is a continuation of 1.1 above, it adds hps as a comparison...
#
#          sshps.3 = makeSS('ss.3tree', 'horizontalPointIZ', angleGauge=ag.3m)
#          sshps.4 = makeSS('ss.3tree', 'horizontalPointIZ', angleGauge=ag.4m)
#          sshps.5 = makeSS('ss.3tree', 'horizontalPointIZ', angleGauge=ag.5m)
#          mohps.3 = ssMODWT(sshps.3, J=4)
#          mohps.4 = ssMODWT(sshps.4, J=4)
#          mohps.5 = ssMODWT(sshps.5, J=4)
#          hps = hfsMODWT(mohps.3, mohps.4, mohps.5, ids=NA, long=TRUE)
#          hps$sm = 'hps'
#          jl$sm = 'cp'
#          j = rbind(jl,hps)
#          xyplot(iso.var ~ izArea|iso.j, j, scales='free', groups=sm)
#
#       This use of hfsMODWT as ids=NA, so the object names are used as ids.
#
#       Remember, the above are all volume sampling surfaces, this same exercise
#       can be done with any attibute that SS can handle. The angleGauge objects
#       are metric BAFs.
#
#
#Author...									Date: 27-Apr-2016
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   need sampSurf for area methods...
#
    #require(sampSurf, quietly = TRUE)
    
#
#   how many ssMODWT objects were passed? ...
#
    args = list(...)
    n.args = length(args)
    arg.names = sapply(substitute(list(...)), deparse)[-1]
    
#
#   make sure all "..." are ssMODWT objects...
#
    if(any(!sapply(args, is, 'ssMODWT')))
      stop('All objects passed in "..." must be of class "ssMODWT" or subclass!')

#
#   assumes all objects were decomposed to the same level... (might want a check here? <*******)
#
    j = args[[1]]@levels$j
    J = args[[1]]@levels$J
    tau.j = args[[1]]@levels$tau.j
    
    if(any(sapply(args, function(x) x@levels$J) != J))
      stop('Number of levels must be the same in each object: J = ',J)

#
#   create a wide format data frame for the results; assumes all methods have same J...
#
    ncol = J+7
    df = data.frame(matrix(NA, nrow = n.args, ncol = ncol))
    names(df) = c('id', names(args[[1]]@vars.modwt$summary$iso.var),
                  'mean', 'scaleVar', 'waveVar', 'waveSampVar', 'avg.E', 'izArea')

#
#   add the sampling method size ids, they default to the argument names themselves...
#
    if(any(is.na(ids))) {                    #default to the object names that were passed--promised
      ids = arg.names
      rownames(df) = paste(1:n.args, sep='') #no need to have row names the same
    }
    else
      rownames(df) = c(arg.names)
    df$id = ids
    rowNames = rownames(df)

#
#   just loop through each object, collecting the iso.var's, etc....
#
    for(i in seq_len(n.args)) {
      md = args[[i]]
      if(md@levels$J != J)
        stop('Number of levels must be the same in each object: J = ',J)
      #an = arg.names[i] 
      an = rowNames[i]
      df[an, 2:(J+1)] = md@vars.modwt$summary$iso.var                    #assumes all J's are equal!
      df[i, 'avg.E'] = sum(md@vars.modwt$summary$iso.var)                #E[x^2] == total energy
      df[i, 'waveVar'] = df[i, 'avg.E'] - md@vars.modwt$total$LL.var     #emp. wavelet variance, no smooth
      df[i,'mean'] = md@vars.modwt$total$LL.mean                         #same as the surfStats mean
      df[i, 'scaleVar'] = with(md@vars.modwt$total, LL.var - LL.mean^2)  #empirical scaling variance
      #df[i, 'izArea'] = area(md@ss@izContainer@iZones[[1]])   #IZ area
      df[i, 'izArea'] = mean(sapply(md@ss@izContainer@iZones, area))     #average IZ area
    }

#
#   this will match the sampSurf object's surfStats variance as in ssMODWT...
#
    df[, 'waveSampVar'] = df[, 'avg.E'] - df[, 'mean']^2

#
#   sort by ascending inclusion zone area...
#
    df = df[order(df$izArea), ]

#
#   reshape for long format if desired...    
#
    if(long) {
      isoCols = colnames(df)[grep('iso', colnames(df))]
      df = reshape(df, varying = isoCols,             #columns to move into a long column
                       v.names = 'iso.var',           #name for the varying column
                       timevar = 'iso.j',             #name of the times column in long
                       times = isoCols,               #the times ids in timevar
                       new.row.names = paste(1:(J*n.args),sep=''),  #best to supply these
                       direction = 'long'
                  )
      df$tau.j = as.vector(sapply(tau.j, rep, n.args))
      df$j = as.vector(sapply(j, rep, n.args))
    } 
    
    if(!runQuiet){ 
      cat('\nNumber of input objects =', n.args)
      cat('\nOutput format =', ifelse(long, 'long', 'wide'))
      cat('\n')
    }                 

    return(invisible(df))
}   #hfsMODWT

