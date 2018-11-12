ortho <-
function(x0, x) {for(i in 1:ncol(x0)) x <- x - crossprod(x0[,i],x)/crossprod(x0[,i],x0[,i])*x0[,i]; x}
