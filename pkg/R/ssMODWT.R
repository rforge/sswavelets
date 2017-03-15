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
setGeneric('ssMODWT',  
           function(ss, ...) standardGeneric('ssMODWT'),
             signature = c('ss')
            )


          
#================================================================================
#
setMethod('ssMODWT',
          signature(ss = 'sampSurf'),
function(ss,
         J = NA,               #NA or NULL gives full decomposition
         wfName = c('haar'),   #just allow haar for now
         reflect = FALSE,
         shift = TRUE,         #shift back w/r to boundary correction
         trimRight = TRUE,     #trim off N+1,M+1
         description = 'sampSurf MODWT wavelet decomposition object',
         runQuiet = FALSE,
         ...
        )
{
#
#---------------------------------------------------------------------------
#
#   All decompositions are done with package waveslim.
#
#   This version of ssMODWT is the S4 version for package ssWavelts.
#
#   This routine does a maximal overlap discrete wavelet transform, MODWT
#   decomposition, on a sampSurf object. The main interest is the actual
#   MODWT results and several different versions of the variance, which are
#   conserved over the different scales in the MODWT analysis. The variances
#   and wavelet smooth (LL) are compared with the results from sampSurf
#   and printed as part of the output if !runQuiet. The routine modwt.2d 
#   is used for the decomposition. Note that we can apply shift.2d to the
#   MODWT matrices to re-align them to the original image, otherwise they
#   are shifted to the right consecutively more for each level j and hence
#   become more out of alignment with the original. Shifting realigns so
#   that they look like the results of the MRA analysis.
#
#   In addition, a MODWT multiresolution analysis (MRA) is also done here.
#   The routine mra.2d is used for this decomposition. MRA is, by definition,
#   aligned properly with the original, so it does not require shifting 
#   like the MODWT does.
#
#--------------------------
#   Note that modwt.2d returns a list with levels 1,...,J including the
#   following decomposed matrices...
#
#     LH = upper left == wavelet-scaling matrices V   (horizontal filter)
#     HH = upper right == wavelet-wavelet matrices W  (diagonal filter)
#     HL = lower right == scaling-wavelet matrices U  (vertical filter)
#
#   At level J, the final "smooth" (scaling-scaling) result is returned in
#   LLJ, where J is the last level (e.g., LL5 if J=5).
#--------------------------
#
#   Results can pe presented graphically using plotMODWT2D() and plotLevel2D().
#
#   The variances all work out to the actual surface variance given by sampSurf, 
#   for details see...
#
#   D. Mondal and D. B. Percival (2012), Wavelet Variance Analysis for Random 
#   Fields on a Regular Lattice, IEEE Transactions on Image Processing, 21(2), 
#   537-549.
#
#   Arguments...
#     ss = a valid sampSurf object; note that the internal raster sampling surface
#          must be square and dyadic in extent
#     J = the highest level for decomposition, which can be less than log(N,2),
#         where N is the number of rows and columns. If J=NA, then the default
#         is the full J=log(N,2) levels of decomposition
#     wfName = a wavelet filter name recognized by waveslim--only 'haar' for now
#     reflect = TRUE: use reflection+periodization via ssReflect; FALSE: periodization
#     shift = TRUE: apply shift.2d to rectify the MODWT images with the original;
#             FALSE: leave the alignment as is
#     trimRight = TRUE: if(shift) then this will trim off the N+1,M+1 row and column;
#                 FALSE: if(shift) this will trim the 1,1 row and column off
#     runQuiet = TRUE: no feedback; FALSE: a short report on the results
#     ... = gobbled
#
#   Note that trimRight is only used on matrices whose dimension is larger than the
#        original. This sometimes seems to happen on the level one set of matrices
#        returned from modwt.2d for some reason.
#
#   Returns...
#     an ssMODWT object invisibly with...
#     -- the input sampSurf object
#     -- the MODWT results
#     -- a list of the variances, see varMODWT
#     -- the results of the multiresolution analyis (MRA) using MODWT
#     -- the wavelet filter name--'haar' for now
#     -- a list of other information, see code
#
#   ++ Adapted to S4 for the ssWavelets package 17-Jan-2017, JHG.
#   ++ Removed trim in-line code and replaced it with the Trim function so that
#      it would also work on the full Reflected images rather than just the
#     image that is the size of the ss@tract. 28-Feb-2017, JHG.
#
#**>Note that the date below is for the original code that was essentially
#        the development version of this S4 version in the ssWavelets package.
#
#Author...									Date: 21-Apr-2016
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   the usual few checks--not very comprehensive...
#
    wfName = match.arg(wfName)
   
    #see the vignettes about this... 
    if(xres(ss@tract != 1) && shift)
      warning('\nResolutions other than 1 do not work with waveslim::shift.2d!')

    #inclusion zone class ==> type of sampling method used...
    izClass = class(ss@izContainer@iZones[[1]])[1]
    
#
#   dimensions, scales, etc...
#
    N = ncol(ss@tract)               #number of cells in x
    M = nrow(ss@tract)               #and y
    if(M != N)
      stop('For now, ss tract must be square')
     
    #maximum level is half the smallest tract dimension in N or M;
    #this is my constraint and not necessary--put an override in?...
    min.NM = min(N,M)
    if(is.na(J) || is.null(J))
      J = log(min.NM, 2) - 1
    if(2^J > min.NM/2) {
      J = log(min.NM, 2) - 1
      cat('\nJ too large, reset to J =', J,'\n')
    }      
      
#***
#***>need to check Mj and Nj: do they assume index begin at zero? <***
#***
    L = 2                            #number of coefficients for Haar
    j = 1:J                          #level index for scale
    tau.j = 2^(j-1)                  #the scale
    Lj = (2*tau.j - 1)*(L-1) + 1     #the filter width at level j, or scale 2*tau.j
    Nj = N - Lj + 1                  #number of non-boundary jth-level coefficients in x
    Mj = M - Lj + 1                  #number of non-boundary jth-level coefficients in y
    
#
#   trim function: compare the dimension of each of the return matrices from
#   shift.2d to those of the base image (extent) as the former will pad an
#   extra row and column for some reason on j=1 only...
#
    Trim = function(ss.modwt, ext) {
      #assume min extents != 0...
      N = ext@xmax - ext@xmin     #cols
      M = ext@ymax - ext@xmin     #rows
      #loop through each modwt.2d object matrix to check for trimming...
      for(i in seq_len(length(ss.modwt))) {
        m.d = ss.modwt[[i]] 
        if(!all(dim(m.d) == c(M,N)))               #trim it back to correct dimension
          if(trimRight)
            ss.modwt[[i]] = ss.modwt[[i]][1:N, 1:M]          #trim right: N+1,M+1
          else
            ss.modwt[[i]] = ss.modwt[[i]][2:(N+1), 2:(M+1)]  #trim left: 1,1
      }
      return(ss.modwt)
    } #Trim

#
#   convert the raster results to a matrix for use in modwt.2d...
#
    m = as.matrix(ss@tract)
	ex.tr = extent(ss@tract)
    #r = raster(m, ex@xmin, ex@xmax, ex@ymin, ex@ymax)
    NM = length(m)                            #ncells or simply N*M

#
#   perform the MODWT and MRA analysis...
#
    if(reflect) {                                      #Reflection + Periodic
      cat('\nReflection + Periodic')
      r.m = ssReflect(ss)                              #returns raster
      exr = extent(r.m)                                #get these expanded extents for cropping below
      m = as.matrix(r.m)                               #must cast to matrix...
      #MODWT analysis...
      ss.modwt = modwt.2d(m, wfName, J=J)              #MODWT -- only periodic boundary is supported
      att.ss.modwt = attributes(ss.modwt)
      if(shift) {                                      #shift the entire enlarged scene
        ss.modwt = shift.2d(ss.modwt)                  #shifting MUST precede cropping <****
        ss.modwt = Trim(ss.modwt, exr)                 #trim to the full reflected image
      }
      #need to crop each of the modwt matrices back to original extents; these must be two steps--don't combine!...
      tmp = sapply(ss.modwt, function(x) crop(raster(x, exr@xmin, exr@xmax, exr@ymin, exr@ymax), ex.tr))
      ss.modwt = lapply(tmp, as.matrix)                #convert back to a list of matrices
      attributes(ss.modwt) = att.ss.modwt              #reset attributes for a recognizable object
      #multi-resolution analysis...
      ss.mra = mra.2d(m, wfName, J=J, method='modwt')  #MRA for MODWT -- only periodic boundary is supported
      att.ss.mra = attributes(ss.mra)
      tmp = sapply(ss.mra, function(x) crop(raster(x, exr@xmin, exr@xmax, exr@ymin, exr@ymax), ex.tr)) #crop--see above
      ss.mra = lapply(tmp, as.matrix)
      attributes(ss.mra) = att.ss.mra                  #reset attributes for a recognizable object
    }
    else {                                             #Periodic only
      ss.modwt = modwt.2d(m, wfName, J=J)              #MODWT -- only periodic boundary is supported
      if(shift) {
        ss.modwt = shift.2d(ss.modwt)
        ss.modwt = Trim(ss.modwt, ex.tr)               #trim the images back to tract extents
      }
      ss.mra = mra.2d(m, wfName, J=J, method='modwt')  #MRA for MODWT -- only periodic boundary is supported
    }

#
#   a few corrections for modwt.2d and mra.2d under waveslim v1.7.5...
#
    names(ss.mra)[3*J + 1] = paste('LL', J, sep='')     #correct for mra--wrongly named 'LL1'
    #ss.modwt = fixMODWT(ss.modwt, FALSE)
    #ss.mra = fixMODWT(ss.mra, FALSE)


#
#   get the variances...
#
    ##vars = varMODWT(ss.modwt)
    vars = covMODWT(ss.modwt)

    if(!runQuiet) {
      cat('\nInclusion zone class:', izClass)
      modwt.var = vars$total$modwt.var
      isoSurf.var = vars$image$isoSurf.var
      LL.mean = vars$total$LL.mean
      cat('\nSurface variance =', ss@surfStats$var)
      cat('\nWavelet variance =', modwt.var, '(MODWT)')
      cat('\nSurf/Wave var ratio =', ss@surfStats$var/modwt.var)
      cat('\nCheck: surface var matrix =', sum(isoSurf.var)- ss@surfStats$mean^2, 
          '= E[X^2] - E[X]^2')   
      cat('\nSurface mean =', ss@surfStats$mean)
      cat('\nWavelet mean =', LL.mean, '(MODWT)')    
      cat('\nWavelet mean =', mean(ss.mra[[grep('LL', names(ss.mra))]]), '(MRA)') 
      cat('\n\n')
    }

 
    ssMODWT = new('ssMODWT',
                  description = description,
                  ss = ss,
                  ss.modwt = ss.modwt,
                  vars.modwt = vars,
                  ss.mra = ss.mra,
                  wfName = wfName,
                  levels = list(J = J,
                                N = N,
                                M = M,
                                L = L,
                                j = j,
                                tau.j = tau.j,
                                Lj = Lj,
                                Nj = Nj,
                                Mj = Mj
                               )
                 )
    return(invisible(ssMODWT))

}   #ssMODWT constructor
)   #setMethod

