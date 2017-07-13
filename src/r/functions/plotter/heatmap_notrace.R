# R function
#
# Author: baderd
###############################################################################


heatmap_notrace= function( x, 
    denscol='green', col=bluered(50), key.par=list(las=1), keysize=1
    ,...
){
    require(gplots)
    heatmap.2(
        x=x, trace='none', 
        key.ylab='', key.title='', keysize=keysize, key.par=key.par,
        col=col, denscol=denscol
        ,...
    )
}

