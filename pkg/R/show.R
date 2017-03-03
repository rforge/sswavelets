#---------------------------------------------------------------------------
#
#   Methods for generic show() for class...
#     (1) ssWavelet and subclasses (e.g., ssMODWT, ssCovMODWT)
#
#   For now, these all just call the corresponding summary method.
#
#Author...									Date: 25-Jan-2017
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#

#================================================================================
#  method for base class ssWavelet & subclasses...
#
setMethod('show',
          signature(object = 'ssWavelet'),
function(object)
{
    return(summary(object))
}   #show for 'ssWavelet'
) #setMethod


