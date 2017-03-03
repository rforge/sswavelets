ssToroid = function(ss,
                     ...
                   )
{
#---------------------------------------------------------------------------
#
#   This routine takes a sampSurf object and makes a larger raster object with
#   copies of the original all around it using toroidal correction, more
#   commonly called periodization or circular correction for the boundary.
#
#   This routine is specifically used for correction of "edge effect" in
#   the application of a MODWT wavelet filter as discussed in detail in
#   Lark and Webster (2004) Eur. J. Soil Sci. (L&W).
#
#   Also see ssReflect() for reflection boundary correction.
#
#   Here we have done complete periodization in all directions to pad the
#   larger image. L&W used a cropped padding on the left because of the
#   actual wavelet used. This can easily be done if desired using the
#   raster crop() function on the returned image.
#
#   Arguments...
#     ss = a valid 'sampSurf' object
#     ... = gobbled
#
#   Returns...
#     a raster imaged padded in all directions using periodization as
#     described above.
#
#Author...									Date: 14-July-2016
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   first make sure raster is loaded and this is a sampSurf object...
#
    #require(raster)
    if(!is(ss, 'sampSurf'))
      stop("Argument 'ss' must be of class 'sampSurf'.")

#
#   get the correct spatial extents of the sampling surface...
#	
	ex = extent(ss@tract)
    xmin = ex@xmin
    xmax = ex@xmax
    ymin = ex@ymin
    ymax = ex@ymax
    
#
#   and the offsets for reflection...
#
    xOff = xmax - xmin
    yOff = ymax - ymin

#
#   make it a matrix, then a rasterLayer with correct extents...
#
    m = as.matrix(ss@tract)
	r = raster(m, xmin, xmax, ymin, ymax)

#
#   slide horizontally -- columns...
#
    r.h = merge(shift(r, x = -xOff), merge(r, shift(r, x = xOff))) #shift to left and right

#
#   slide vertically -- rows...
#    
    r.v = merge(shift(r, y = -yOff), merge(r.h, shift(r, y = yOff))) #shift to top and bottom

#
#   now add the diagonals/corners (bottom-left, top-right) -- shift in both directions...
#
    r.d = merge(shift(r, x = -xOff, y = -yOff), merge(r.v, shift(r, x = xOff, y = yOff)))
    
#
#   now add the diagonals/corners (top-left, bottom-right)...
#
    r2 = merge(shift(r, x = -xOff, y = yOff), merge(r.d, shift(r, x = xOff, y = -yOff)))  
    
    return(r2)
}   #ssToroid

