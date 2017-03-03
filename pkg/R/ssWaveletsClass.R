#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the "ssWavelts" virtual 
#   and subclass(es)
#
#   Note that we could extend sampSurf here, but I have elected to rather
#   provide a ss as a slot.
#
#   Classes and subclasses in this file include...
#
#     1. ssWavelet  --  VIRTUAL
#     2. ssMODWT    --  subclass of ssWavelet
#     3. ssCovMODWT --  ssMODWT covariance subclass   ****this is pairMODWT****
#
#   Note that it would be a simple matter to include other wavelet filters
#   under the above setup. 
#
#Author...									Date: 13-Jan-2017
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#




#=================================================================================================
#
#  1. define the ssWavelet class...
#
setClass('ssWavelet',
#
#  slots for the class and its subclasses...
#
    representation(description = 'character',
                   wfName = 'character',                       #wavelet filter name
                   ss = 'sampSurf'                             #always works on sampSurf objects
                   ),
                   
    contains = 'VIRTUAL',                                      #note!
    
    #some defaults for validity checking...
    prototype = list(description = 'sampSurf wavelet object', 
                     wfName = '', 
                     ss = new('sampSurf')                      #dummy, generates an invalid object
                    ),
                    
    validity = function(object) {

#              no checks for now...
                   
                 return(TRUE)
               } #validity check
) #class ssWavelet 





#=================================================================================================
#
#  2. the ssMODWT class is a direct descendant/subclass of 'ssWavelet'...
#
#
setClass('ssMODWT',
#
#  slots for the class and its subclasses...
#
    representation(ss.modwt = 'list',                          #modwt wavelet decomposition
                   vars.modwt = 'list',                        #modwt variances -- lots!
                   ss.mra = 'list',                            #modwt mra
                   levels = 'list'                             #levels, dimensions, etc.
                  ),
                   
    contains = 'ssWavelet',                                    #descendent of ssWavelet
    
    #some defaults for validity checking...
    prototype = list(ss.modwt = list(),
                     vars.modwt = list(),
                     ss.mra = list(),
                     levels = list()
                    ),
                    
    validity = function(object) {
    
                 Js = function(x)                  #calculate J on the object
                      return( (length(x) - 1)/3 )
                 #ss.modwt...
                 J1 = Js(object@ss.modwt)
                 if(J1 != object@levels$J)
                   return('Levels in ss.modwt slot does not match levels$J slot!')
                 #ss.modwt...
                 J2 = Js(object@ss.mra)
                 if(J2 != object@levels$J)
                   return('Levels in ss.mra slot does not match levels$J slot!')

#                check correct set of names in the modwt & mra objects...

                 dnames = c(as.vector(t(sapply(c('LH','HL','HH'), paste, 1:J1, sep=''))),
                            paste('LL', J1, sep=''))
                 if(!all(dnames == names(object@ss.modwt)))
                   return('Incorrect names for ss.modwt decomposition matrices')                            
                 if(!all(dnames == names(object@ss.mra)))
                   return('Incorrect names for ss.mra decomposition matrices') 
                   
#                check for haar filter, could have a list of others sometime...
                 if(!(object@wfName == 'haar'))
                   return('Only the haar wavelet filter is supported right now')
                   
                 return(TRUE)
               } #validity check
) #class ssMODWT




#*****>>>> add correlation slots below sometime????  <<<<<***********************

#=================================================================================================
#
#  3. the ssCovMODWT class is a direct descendant/subclass of 'ssMODWT'...
#
#
setClass('ssCovMODWT',
#
#  slots for the class and its subclasses...
#
    representation(ss.b = 'sampSurf',                         #the second sampSurf object
                   ss.modwt.b = 'list',
                   ss.mra.b = 'list',
                   covStats = 'list'
                  ),
                   
    contains = 'ssMODWT',                                     #descendent of ssMODWT
    
    #some defaults for validity checking...
    prototype = list(ss = new('sampSurf'),                    #dummy, generates an invalid object
                     ss.modwt.b = list(),
                     ss.mra.b = list(),
                     covStats = list()
                    ),

    #all validity checks should be in the superclasses...                    
    validity = function(object) {

                 return(TRUE)
               } #validity check
) #class ssCovMODWT



