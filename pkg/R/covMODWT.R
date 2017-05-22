covMODWT = function(ss.modwt.x, #if these two are the same, the results
                    ss.modwt.y, #must be identical to varMODWT -- a check
                    ...
                   )
{
#---------------------------------------------------------------------------
#
#   This routine will calculate the decomposition of the sample (co-)variance
#   that is applicable to sampSurf images. It includes...
#
#       -- individual matrix (image) components/contributions
#       -- marginal variances for each level and filter type (see below)
#       -- totals, scalar variances and surface mean
#
#   This should work on a modwt and mra+modwt list returned from waveslim.
#   Note that it is up to the user to determine whether the latter variance
#   calculation is applicable (i.e., Haar), as no information about the 
#   type of wavelet filter is required here.
#
#   The following applies to waveslim v 1.7.5...
#
#   modwt.2d and mra.2d both return a list with levels j=1,...,J including
#   the following decomposed matrices...
#
#     LHj = upper left == wavelet-scaling matrices V   (horizontal filter)
#     HHj = upper right == wavelet-wavelet matrices W  (diagonal filter)
#     HLj = lower right == scaling-wavelet matrices U  (vertical filter)
#
#   At level J, the final "smooth" (scaling-scaling) result is returned in
#   LLJ, where J is the last level (e.g., LL5 if J=5).
#
#   Note: All isotropic variances at level J include the LLJ smooth component
#         variances, whether totals, summary or image.
#
#   References...
#       M&P = D. Mondal and D. B. Percival. 2012. Wavelet variance analysis
#             for random fields on a regular lattice. IEEE Transactions on
#             Image Processing, 21(2):537-549.
#       GPS = M. Geilhufe, D. B. Percival, and H. L. Stern. 2013.
#             Two-dimensional wavelet variance estimation with
#             application to sea ice sar images. Computers and Geosciences, 
#             54:351-360.
#       Univariate but still applicable...
#       P&M = D. B. Percival & D. Mondal. 2012. A wavelet variance primer. Handbook
#             of statistics, Vol. 30, chap. 22.
#
#   Arguments...
#     ss.modwt.x = a list as returned from modwt.2d() or mra.wd() with the caveat
#                  given above about variances -- for the first sampling method
#     ss.modwt.y = same for the second sampling method; if it is missing, it will
#                  be set to ss.modwt.x and the variance will be calculated
#     ... = gobbled
#
#   Returns...
#     a list with the components described above in...
#       -- summary: marginal total variances
#       -- image: matrix/image variances
#       -- total: surface decomposition totals
#
#Author...									Date: 6-Sept-2016
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   do we have covariance or variance?...
#
    if(missing(ss.modwt.y)) {    #if missing, then simple variance
      ss.modwt.y = ss.modwt.x    #set y==x
      isCovariance = FALSE
    }
    else {
      isCovariance = TRUE
    }
    
      
#
#   some checks on the two objects passed...
#
      
    Js = function(x)                  #calculate J on the object
          return( (length(x) - 1)/3 )
          
    J.x = Js(ss.modwt.x)
    J.y = Js(ss.modwt.y)
    if(J.x != J.y)
      stop('Both objects must have the same number of levels J') 
      
    N.x = ncol(ss.modwt.x[[1]])                
    M.x = nrow(ss.modwt.x[[1]])                
    N.y = ncol(ss.modwt.y[[1]])                
    M.y = nrow(ss.modwt.y[[1]])
    if(N.x!=N.y || M.x!=M.y)
      stop('Dimensions of the MODWT objects are not comparible!')                

    x.names = names(ss.modwt.x)
    y.names = names(ss.modwt.y)
    if(!all(x.names == y.names))
      stop('Names must be the same for the two list objects!')
    

#
#   a few common/shared metrics about the surfaces and decomposition levels...
#
    J = J.x 
    j = 1:J                               #level index for scale

    N = N.x                               #number of cells in u (x) across columns
    M = M.x                               #and v (y) along rows
    NM = N*M                              #total: ncells() from sampSurf image
    
    w.names = x.names                     #names of the wavelet matrices j=1:J
    
#
#**>Note: The following gives a biased estimate of the wavelet variance because it includes
#         all N*M cells; i.e., those with periodic boundary correction. The unbiased
#         estimator uses N_j*M_j, which do not include the boundary cells; i.e., it excludes
#         these cells from the calculation. See M&P and P&M for more details. However,
#         I believe it is okay for the sample variance breakdowns, though not certain.
#         In addition, if we use a large enough external buffer that little or nothing
#         overlaps with periodic correction, then I don't believe it matters; therefore, 
#         what we have below should be fine.
#
#   everything is based on this result; it is essentially E[X^2_{u,v}] for each
#   decomposition component and includes: all: LHj, HHj, HLj, LLJ, j=1:J,
#   and is a list of E[X^2_{u,v}] matrices of the above names at each level,
#   where (u,v) is an individual spatial location/cell u=1:N, v=1:M
#
#   See M&P, specifically the "contributions to the sample variance" given on p.542
#   and 543, where the last equation in the section is the sample variance
#   decomposition when j=j' ==> the scales are the same in each direction: rows
#   & columns -- (12) is more general and allows differing (cross-variance) scales
#
#   the energy matrix/image would not include the divisor below; note that it would
#   be nice not to divide by NM here, but to do it later with the summary variances,
#   but there seems no better place and it would not make that much difference in
#   execution time...
#
    #w.var = lapply(ss.modwt, function(x) (x*x)/NM)      #over all: LHj, HHj, HLj, LLJ, j=1:J
    #over all: LHj, HHj, HLj, LLJ, j=1:J...
    w.cov = mapply(function(xx,yy) (xx*yy)/NM, ss.modwt.x, ss.modwt.y, SIMPLIFY=FALSE)
    

#
#   LLJ is the final residual/smooth matrix--taking the mean of LLJ gives us the sampSurf mean...
#
    LL.x = ss.modwt.x[[grep('LL', w.names)]]                #LL == scaling-scaling average matrix at J
    LL.y = ss.modwt.y[[grep('LL', w.names)]]                #LL == scaling-scaling average matrix at J
    LL.mean.x = mean(LL.x)                                  #a scalar
    LL.mean.y = mean(LL.y)                                  #a scalar
    #LL.mean = LL.mean.x * LL.mean.y                     #for sample covariance    
    LL.wcov = w.cov[[grep('LL', w.names)]]              #an NxM matrix
    LL.var = sum(LL.wcov)                               #scalar E[LL^2]

    
#
#   calculate the (anisotropic) variances for different filters at all levels
#   summing over the (u,v) -- the results are a vector -- these are the 
#   marginal component wavelet variances...
#
    LH.wcov = w.cov[grep('LH', w.names)]                #LH matrices
    LH.var = sapply(LH.wcov, function(x) sum(x))        #i.e., E[LH^2]
    
    HL.wcov = w.cov[grep('HL', w.names)]                #HL matrices
    HL.var = sapply(HL.wcov, function(x) sum(x))        #i.e., E[HL^2]
    
    HH.wcov = w.cov[grep('HH', w.names)]                #HH matrices
    HH.var = sapply(HH.wcov, function(x) sum(x))        #i.e., E[HH^2]

#
#   this computes the marginal isotropic "sample variance" vector from the above
#   and includes the smooth component at the Jth level...
#
    iso.var = LH.var + HL.var + HH.var
    iso.var[[J]] = iso.var[[J]] + LL.var          #add in the trend variance component for last level
    names(iso.var) = paste('iso', j, sep='')      #otherwise names are LH_j
    
#
#   a list of simple vectors...
#
    summary = list(LH.var = LH.var, HL.var = HL.var, HH.var = HH.var, 
                   iso.var = iso.var)
 
    
#
#   The surface variance should be equal to the sum of the "sample variances" 
#   (including the smooth, LLJ) minus the squared mean:  (coefficients^2) - mean^2
#   i.e., E[X^2] - E[X]^2, where the X are the different decomposition coefficients;
#   this is (12) [diagonal: j==j'] in Mondal & Percival (2012); this will be equal 
#   (very close) to the sampSurf variance in ss@surfStats$var...
#
#   The same idea for a covariance: E[XY] - E[X]E[Y]...
#
    modwt.var = sum(iso.var) - LL.mean.x * LL.mean.y      #== LL.mean^2 for variance
    
#
#   a list of scalar totals variances, including the surface mean from the decomposition...
#
    total = list(LL.mean = LL.mean.x, 
                 LL.var = LL.var, 
                 modwt.var = modwt.var
                )
    if(isCovariance)
      total$LL.mean.y = LL.mean.y               
 
#
#   the isotropic spatial surface variance matrices for each level, j=1:J;
#   note that the smooth variance is included here as with both iso.var and 
#   isoSurf.var for level J..
#
    Iso.var = vector('list', J)   
    for(i in j) {
      Iso.var[[i]] = LH.wcov[[i]] + HL.wcov[[i]] + HH.wcov[[i]]
    }
    Iso.var[[J]] = Iso.var[[J]] + LL.wcov       #include LLJ at level J
    names(Iso.var) = paste('iso', j, sep='') 

   
#
#   lastly, the cumulative isotropic spatial surface variance matrix over all levels, 
#   the sum of this matrix will also equal the total surface variance in modwt.var
#   as above when we subtract ss@surfStats$mean^2 off the sum since it includes
#   the smooth component; however, subtracting the LL.mean^2 may drive many of
#   the (u,v) variance components in the matrix negative because their
#   individual contribution to the variance is very small comparatively; so
#   for now, we leave this off; i.e., the sample variance image + LL.mean^2...
#
    isoSurf.var = matrix(0, nrow = M, ncol=N)   #==> modwt.var = sum(isoSurf.var) - LL.mean^2
    for( i in j)
      isoSurf.var = isoSurf.var + Iso.var[[i]]  #LL.wcov is already in Iso.var[[J]]

#
#   a list of list of variance matrices at each level...
#
    image = list(LH.var = LH.wcov,              #list of matrices j=1:J
                 HL.var = HL.wcov,              #list of matrices j=1:J 
                 HH.var = HH.wcov,              #list of matrices j=1:J 
                 LL.var = LL.wcov,              #a matrix: j=J 
                 iso.var = Iso.var,             #list of matrices j=1:J 
                 isoSurf.var = isoSurf.var      #single cumulative matrix level J
                )
    
    return(invisible(list(isCovariance = isCovariance,
                          summary = summary,
                          total = total, 
                          image = image 
                         )
                    )
          )
}   #covMODWT
    
    
    

