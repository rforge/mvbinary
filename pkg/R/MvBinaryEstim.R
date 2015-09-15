##' MvBinary a package for Multivariate Binary data
##'
##' MvBinary is a tool for fitting distribution on correlated multivariate binary data.
##'
##' \tabular{ll}{
##'   Package: \tab MvBinary\cr 
##'   Type: \tab Package\cr 
##'   Version: \tab 1.0.0\cr
##'   Date: \tab 2015-07-29\cr 
##'   License: \tab GPL-2\cr 
##'   LazyLoad: \tab yes\cr
##' }
##'
##'
##' @name MvBinary-package
##' @aliases MvBinary
##' @rdname MvBinary-package
##' @docType package
##' @keywords package
##' @import mgcv
##' @import parallel
##' @import graphics
##' @export MvBinaryEstimCAH
##' @exportClass MvBinaryResult
##'
##' @author
##' Author: Marbac M., and Sedki S.
##'
##' @references Todo.
##'
##' @examples
##' rm(list=ls())
##' require(MvBinary)
##' data(MvBinaryExample)
##' test <- MvBinaryEstimCAH(MvBinaryExample)
##'
NULL

MvBinaryEstimCAH <- function(x, nbcores=1, tol=0.01, nbiter=10){
  if (is.null(colnames(x))) colnames(x) <- paste("x",1:ncol(x), sep="")
  # Computation of the Cramer's V 
  VcramerEmpiric <- matrix(0, ncol(x), ncol(x))
  rownames(VcramerEmpiric) <- colnames(VcramerEmpiric) <- colnames(x)
  alpha <- colMeans(x)
  for (j in 1:ncol(x)){
    obs <- rbind((t(x[,j])%*%(x)), (t(x[,j])%*%(1-x)), (t(1-x[,j])%*%(x)), (t(1-x[,j])%*%(1-x)))/nrow(x)
    th <- rbind(alpha[j]*alpha, alpha[j] * (1-alpha), (1-alpha[j])*alpha, (1-alpha[j])*(1-alpha))
    VcramerEmpiric[j,] <- sqrt(colSums((obs - th)**2 / th))
  }
  tree <- hclust(as.dist(1-VcramerEmpiric), method="ward")
  models <- list(); for (k in 1:ncol(x)) models[[k]] <- cutree(tree, k)
  # Inference for the competiting models
  nb.cpus <- min(detectCores(all.tests = FALSE, logical = FALSE), nbcores)
  if ((nbcores>1)&&(Sys.info()["sysname"] != "Windows")){
      reference <- mclapply(X = models, FUN = XEMmodel, dataset=x, alpha=alpha, tol=tol, nbiter=nbiter, mc.cores = nb.cpus, mc.preschedule = TRUE, mc.cleanup = TRUE)
  }else{
    reference <- list(); for (loc in 1:length(models)) reference[[loc]] <- XEMmodel(x, alpha, tol, nbiter, models[[loc]])
  }
  # Design outputs
  allBIC <- rep(NA, length(reference))
  allModels <- matrix(NA, length(reference), ncol(x))
  for (loc in 1:length(reference)){
    allBIC[loc] <- reference[[loc]]$bic
    allModels[loc,] <- models[[loc]]
  }
  Best <- reference[[which.max(allBIC)]]
  names(reference[[which.max(allBIC)]]$epsilon) <-  names(reference[[which.max(allBIC)]]$delta) <- colnames(x)
  return( new("MvBinaryResult", 
              alpha=alpha, 
              epsilon=reference[[which.max(allBIC)]]$epsilon, 
              delta=reference[[which.max(allBIC)]]$delta, 
              blocks=reference[[which.max(allBIC)]]$model,
              nbparam=reference[[which.max(allBIC)]]$nbparam,
              loglike=reference[[which.max(allBIC)]]$loglike,
              bic=reference[[which.max(allBIC)]]$bic)
  )
}
