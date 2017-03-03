ssReflect = function(ss,
                     ...
                    )
{
#---------------------------------------------------------------------------
#
#   This routine really shows how the functions in raster shine when you
#   need to manipulate an image. Here we are reflecting the original image
#   to each side, top and bottom, and on the diagonals. To do this we
#   must flip the image various ways, and tile (merge) it all together.
#
#   This routine is specifically used for correction of "edge effect" in
#   the application of a MODWT wavelet filter as discussed in detail in
#   Lark and Webster (2004) Eur. J. Soil Sci. (L&W).
#
#   Also see ssToroid() for circular boundary correction.
#
#   Here we have done complete reflections in all directions to pad the
#   larger image. L&W used a cropped padding on the left because of the
#   actual wavelet used. This can easily be done if desired using the
#   raster crop() function on the returned image.
#
#   Arguments...
#     ss = a valid 'sampSurf' object
#     ... = gobbled
#
#   Returns...
#     a raster imaged padded in all directions using reflection as
#     described above.
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
#   make it a matrix, then a rasterLayer with correct extents...
#    
    m = as.matrix(ss@tract)
	r = raster(m, xmin, xmax, ymin, ymax)

#
#   reflect horizontally -- columns...
#
    xOff = xmax - xmin
    yOff = ymax - ymin
	rf.h = flip(r, 'x')
    r.h = merge(shift(rf.h, x = -xOff), merge(r, shift(rf.h, x = xOff))) #reflect to left and right

#
#   reflect vertically -- rows...
#    
    rf.v = flip(r, 'y')
    r.v = merge(shift(rf.v, y = -yOff), merge(r.h, shift(rf.v, y = yOff))) #reflect to top and bottom

#
#   now add the diagonals/corners (bottom-left, top-right) -- reflect in both directions...
#
    rf.hv = flip(rf.h, 'y')  
    rd = merge(shift(rf.hv, x = -xOff, y = -yOff), merge(r.v, shift(rf.hv, x = xOff, y = yOff)))
    
#
#   now add the diagonals/corners (top-left, bottom-right)...
#
    r2 = merge(shift(rf.hv, x = -xOff, y = yOff), merge(rd, shift(rf.hv, x = xOff, y = -yOff)))  
    
    return(r2)
}   #ssReflect

