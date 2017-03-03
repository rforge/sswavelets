palMODWT = function(n,
                    bias = 5,
                    range = NA,
                    ...
                   )
{
#---------------------------------------------------------------------------
#
#   A very simple but effective and nice looking palatte function. 
#   This is biased to the higher end of the scale so that the whitesmoke 
#   only occupies a small part of the range (see Note 2 below).
#
#   Arguments...
#     n = the number of colors (≥ 1) to be in the palette
#     bias = see colorRampPalette
#     range = the range of the variable in question
#     ... = passed to colorRampPalette
#
#**>Note 1: the palette(s) chosen here yield a white-ish around zero, 
#           and then heavily weight the other colors to more area,  
#           essentially squashing it near zero so not much other than zero 
#           will be shown in that color.
#
#**>Note 2: Actually, the "zero" color begins at what ever the minimum
#           is for all positive, or the maximum, for all negative. This
#           could be changed to make them both relative to zero if desired
#           in the future. JHG, 19-Jan-2017
#
#   Returns...
#     a vector of colors in the palette
#
#   An example to illustrate how it looks might be...
#
#       n=100; pie(rep(1,n), col=palMODWT(n, range=c(10,100)))
#
#   simply change the range to different values to test the palette results.
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
#   a little check...
#
    if(length(range) != 2L || any(is.na(range)))
      stop('Argument "range" must be length 2 with no NAs!')
    if(range[2] < range[1])
      stop('Non-standard range: range[2] < range[1]!')

   
    #is the range zero; both could be constant or zero, adjust if so...
    zeroRange = if(identical(range[1], range[2])) TRUE else FALSE
    if(zeroRange)
      if(identical(range[1], 0))              #both zero
        range = c(-1, 1)
      else if(range[1] > 0)                   #both positive
        range = c(0.9*range[1], 1.1*range[2])
      else                                    #both negative
        range = c(1.1*range[1], 0.9*range[2])
      
#
#   positive palette...
#
    pal =  colorRampPalette(c('whitesmoke', 'lightskyblue2', 'saddlebrown'), 
                            bias = bias, ...)
    #see if the range is all nonnegative...                           
    allPos = if(all(range >= 0)) TRUE else FALSE
    if(allPos) 
      return(pal(n))                            


#
#   check if all negative...
#
    #neg.colors = c('whitesmoke', 'snow', 'peachpuff', 'orangered2')
    neg.colors = c('whitesmoke', 'peachpuff', 'orangered2') 
    allNeg = if(all(range <= 0)) TRUE else FALSE   #entire range negative?...
    if(allNeg) {
      neg.pal = colorRampPalette(neg.colors, bias = bias, ...)
      return(rev(neg.pal(100)))
    }

    #both positive and negative?...
    posNeg = if(!allPos && !allNeg) TRUE else FALSE
    
        
    #equal range above and below, we've already caught zeroRange...
    zeroDiff = identical(diff(abs(range)), 0.0)  #equal pos&neg

#
#   for covariance + raw we can generate two, one for positive and one for negative, 
#   but then rev() the negative one so that the bias is on the large negative end
#
    if(posNeg) {
      if(zeroDiff)
        prop = 0.5 #1.0
      else 
        prop = abs(range[1]/sum(abs(range)))           #let's not divide by zero here, see above
      npal.neg = n * prop
      neg.bias = bias*prop
      neg.pal = colorRampPalette(neg.colors, bias = neg.bias, ...)
        tpal = c(rev(neg.pal(npal.neg)), pal(n-npal.neg))      
      return(tpal)
    }
    else {
      return(pal(n))
    }

}   #palMOWDT


