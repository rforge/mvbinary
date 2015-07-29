setMethod(
  f="plot",
  signature = c("MvBinaryResult"),
  definition = function(x){
    op <- par(no.readonly = TRUE)
    colo <- rainbow( max(x@bestparam[,4]))
    if (length(colo)>5)  colo[c(2,5)] <- colo[c(5,2)]
    palette(colo)
    plot(x@bestparam[,1:2], pch=x@bestparam[,4]+15,
         col=x@bestparam[,4], xlab=expression(alpha), ylab=expression(epsilon),
         xlim=c(min(x@bestparam[,1])-0.05, max(x@bestparam[,1])))
    text(x@bestparam[,1], x@bestparam[,2], abbreviate(rownames(x@bestparam),3), pos = 2, col=x@bestparam[,4])
    
    par(op)
  }
)

plotblock <- function(x){
  colo <- rainbow( max(x@bestparam[,4]))
  if (length(colo)>5)  colo[c(2,5)] <- colo[c(5,2)]
  palette(colo)
  
  alpha <- eps <- rep(NA, max(x@bestparam[,4]))
  for (k in 1:max(x@bestparam[,4])){
    alpha[k] <- mean(x@bestparam[which(x@bestparam[,4]==k),1])
    eps[k] <- mean(x@bestparam[which(x@bestparam[,4]==k),2])
  }
  plot(alpha, eps, pch=(1:length(alpha))+15,
       xlab=expression(alpha), ylab=expression(epsilon),
       xlim=c(min(alpha)-0.05, max(alpha)), col=1:length(alpha))
  text(alpha, eps, paste("Block",1:length(alpha),sep=""), pos = 2, col=1:length(alpha))
}

plotdep <- function(x){
  plot(x@dependencies[which(x@dependencies[,4]!=0),3], x@dependencies[which(x@dependencies[,4]!=0),4], xlab="Empiric Cramer's V", ylab="Estimated Cramer's V")
  lines(c(0,1),c(0,1), lty=2)
}

setMethod(
  f="plot",
  signature = c("MvBinaryResult", "character"),
  definition = function(x, y){
    if (y == "variables"){
      plot(x)
    }else if (y=="blocks"){
      plotblock(x)
    }else if (y=="dependencies"){
      plotdep(x)
    }else{
      stop("Argument y must be equal to \"variables\" or \"blocks\" or \"dependencies\" ")
    } 
  }
)
