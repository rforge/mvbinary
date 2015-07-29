MvBinaryEstim <- function(x, models, nbcores=1, tol=0.01, nbinit=10){
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
  
  # If a list of models if not specified, a model collection is given by a CAH based on the Cramer's V matrix of similarities
  if (missing(models)){
    tree <- hclust(as.dist(1-VcramerEmpiric), method="ward")
    models <- list(); for (k in 1:ncol(x)) models[[k]] <- cutree(tree, k)
  }

  # Inference for the competiting models
  nb.cpus <- min(detectCores(all.tests = FALSE, logical = FALSE), nbcores)
  if ((nbcores>1)&&(Sys.info()["sysname"] != "Windows"))
      reference <- mclapply(X = models, FUN = AnalyseMultiBlock, x=x, tol=tol, nbinit=nbinit, mc.cores = nb.cpus, mc.preschedule = TRUE, mc.cleanup = TRUE)
  else
    reference <- list(); for (loc in 1:length(models)) reference[[loc]] <- AnalyseMultiBlock(x, models[[loc]], tol, nbinit)

  
  # Design outputs
  allBIC <- rep(NA, length(reference))
  allModels <- matrix(NA, length(reference), ncol(x))
  error <- FALSE
  for (loc in 1:length(reference)){
    allBIC[loc] <- reference[[loc]]$BIC
    allModels[loc,] <- models[[loc]]
    if (reference[[loc]]$error==TRUE) error <- TRUE
  }
  Bestparam <- reference[[which.max(allBIC)]]$param
  
  CramerModel <- matrix(0, ncol(x), ncol(x))
  for (j in 1:ncol(x)){
    inf <- Bestparam[,1]<Bestparam[j,1]
    alpha1 <- Bestparam[,1]*inf + Bestparam[j,1]*(1-inf)
    alpha2 <- Bestparam[,1]*(1-inf) + Bestparam[j,1]*(inf)
    if (any(Bestparam[,4]==Bestparam[j,4])){
      keep <- which(Bestparam[,4]==Bestparam[j,4])
      CramerModel[keep,j] <- (Bestparam[j,2]*Bestparam[,2] * sqrt(alpha1*(1-alpha2)/(alpha2*(1-alpha1))))[keep]
    }   
  }
  colnames(CramerModel) <- rownames(CramerModel) <- rownames(Bestparam)
  # Permutation of the columns and rows of VcramerEmpiric
  or <- rep(NA, ncol(x)); for (it in 1:ncol(x)) or[it] <- which(colnames(x)==colnames(CramerModel)[it])
  VcramerEmpiric <- VcramerEmpiric[,or]
  VcramerEmpiric <- VcramerEmpiric[or,]
  
  # Backup of the dependencies (observed and modelled)
  CramerModel[upper.tri(CramerModel, diag = TRUE)] <- VcramerEmpiric[upper.tri(VcramerEmpiric, diag = TRUE)] <- 0  
  mat1 <- as.character(matrix(rownames(Bestparam), ncol(x), ncol(x)))
  mat2 <- as.character(matrix(rownames(Bestparam), ncol(x), ncol(x), byrow = TRUE))
  dependencies <- data.frame(variable1=mat1, variable2=mat2, empiric=as.numeric(VcramerEmpiric), model=as.numeric(CramerModel), stringsAsFactors = FALSE)  
  dependencies <- dependencies[-which(dependencies$empiric==0),]
  dependencies <- dependencies[order(dependencies$model, decreasing = TRUE),]
  return( new("MvBinaryResult", 
              allModels=allModels, 
              allBIC=allBIC, 
              dependencies=dependencies, 
              bestparam=Bestparam, 
              bestinfo=list(loglike=reference[[which.max(allBIC)]]$loglike, 
                            nbparam=reference[[which.max(allBIC)]]$nbparam,
                            BIC=reference[[which.max(allBIC)]]$BIC, 
                            groupe=reference[[which.max(allBIC)]]$groupe),
              error=error)
  )
}