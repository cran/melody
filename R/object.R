


require(methods)
setClass("sg", representation("list"))

sg <- function(x, ...){

  if(!is.matrix(x))
    stop('x must be a matrix')
  
  obj <- new('sg', list(x=x, ...))

  obj$width <- ncol(obj$x)  
  obj$height <- nrow(obj$x)


  return(obj)
}
 


   
