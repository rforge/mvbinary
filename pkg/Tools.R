Generator <- function(n, d, param=SampleParam(d)){
  x <- matrix(NA, n, d)
  u <- runif(n)
  for (j in 1:d){
    pj <- (1-param$epsilon[j]) * param$alpha[j] + param$epsilon[j] * (u < param$alpha[j])
    x[,j] <- (runif(n) < pj) * 1
  }
  return(x)
}

EM <- function(x, weight, tol=0.001, nbinit=5){
  # Initialization
  n <- sum(weight)
  if (ncol(x)==1){
    param <- matrix(c(sum(x*weight)/n, 0), 1, 2)
    colnames(param) <- c("alpha", "epsilon")
    rownames(param) <- colnames(x)
    loglike <- sum(weight*log(param[1,1]*x[,1] + (1-param[1,1])*(1-x[,1])))
    nbparam <- 1
    error <- 0
    output <- list(param=param, loglike=loglike, error=error, nbparam=nbparam)
  }else{
    output <- list(loglike=-Inf)
    for (it in 1:nbinit){
      param <- cbind(colSums(x*(weight))/n, runif(ncol(x)))
      colnames(param) <- c("alpha", "epsilon")
      rownames(param) <- colnames(x)
      # Tools for fij computation
      matInf <- matrix(0, length(param[,1]) , length(param[,1]) +1)
      matInf[upper.tri(matInf)] <- 1
      matSup <- 1-matInf
      # fij computation
      baseSup <- log(sweep(x, 2,(1-param[,2])*param[,1] + param[,2], "*") + sweep(1-x, 2, (1-param[,2])*(1-param[,1]), "*"))
      baseInf <- log(sweep(x, 2,(1-param[,2])*param[,1], "*") + sweep(1-x, 2, (1-param[,2])*(1-param[,1])  + param[,2], "*"))
      fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
      # Tools for loglikelihood
      repere <- (c(param[,1], 1)-c(0, param[,1]))
      loglikeprec <- -Inf
      # Loglikelihood computation
      probaInd <- fij %*% repere
      loglikeactu <- sum(log(probaInd)*weight)
      
      while ((loglikeactu - loglikeprec) > tol){
        # Estep
        tik <- sweep(sweep(fij,2,repere,"*"), 1, probaInd, "/")  
        # Mstep
        a <- colSums(weight*x*t(apply(tik[,-ncol(tik)],1, cumsum)))
        d <- colSums(weight*(1-x)*(1-t(apply(tik[,-ncol(tik)],1, cumsum))))
        e <- n - a - d
        param[,2] <- -(sqrt((4*param[,1]^2-4*param[,1]+1)*e^2+((4*param[,1]^2-2*param[,1])*d+4*a*param[,1]^2-6*a*param[,1]+2*a)*e+param[,1]^2*d^2+(2*a*param[,1]-2*a*param[,1]^2)*d+a^2*param[,1]^2-2*a^2*param[,1]+a^2)+
                         (-2*param[,1]^2+2*param[,1]-1)*e+(param[,1]-2*param[,1]^2)*d-2*a*param[,1]^2+3*a*param[,1]-a)/ ((2*param[,1]^2-2*param[,1])*e+(2*param[,1]^2-2*param[,1])*d+2*a*param[,1]^2-2*a*param[,1])
        if (any(param[,2]<0)) param[,2][which(param[,2]<0)] <- 0
        if (any(param[,2]>1)) param[,2][which(param[,2]>1)] <- 1
        if (any(is.na(param[,2]))) param[,2][which(is.na(param[,2]))] <- 0
        # fij computation
        baseSup <- log(sweep(x, 2,(1-param[,2])*param[,1] + param[,2], "*") + sweep(1-x, 2, (1-param[,2])*(1-param[,1]), "*"))
        baseInf <- log(sweep(x, 2,(1-param[,2])*param[,1], "*") + sweep(1-x, 2, (1-param[,2])*(1-param[,1])  + param[,2], "*"))
        fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
        # Loglikelihood computation
        probaInd <- fij %*% repere
        loglikeprec <- loglikeactu
        loglikeactu <- sum(weight*log(probaInd))
      }
      error <- (loglikeactu<loglikeprec)
      nbparam <- 2 * ncol(x)
      loglike=loglikeactu
      if (output$loglike<loglike){
        output <- list(param=param, loglike=loglike, error=error, nbparam=nbparam)
      }
    }
  }
  return(output)
}

FindSensDependencies <- function(x){
  delta <- rep(NA, ncol(x))
  a <- mean(x[,1]) 
  for (j in 1:ncol(x)){
    obs <- table(x[,1],x[,j])/nrow(x)
    b <- mean(x[,j])
    th <- matrix(c((1-a)*(1-b), a*(1-b), b*(1-a), b*a), 2, 2)
    delta[j]  <- unique(sign(diag(obs - th)))
  }  
  names(delta) <- colnames(x)
  return(delta)
}

# Attention ici x et weight sont des listes données par les groupes
AnalyseMultiBlock <- function(x, gr, tol=0.001, nbinit=10){
  loglike <- 0
  nbparam <- 0
  error <- 0
  for (k in 1:max(gr)){
    dataset <- matrix(x[,which(gr==k)], nrow(x), sum(gr==k))
    delta <- FindSensDependencies(dataset)
    if (any(delta==-1)) dataset[,which(delta==-1)] <- 1 - dataset[,which(delta==-1)] 
    ord <- order(colMeans(dataset))
    tmp <- uniquecombs(dataset)
    weight <- as.numeric(table(attr(tmp,"index")))
    dataset <- matrix(tmp, length(weight), ncol(dataset))
    colnames(dataset) <- colnames(x)[which(gr==k)]
    if (ncol(dataset)>1){
      dataset <- dataset[,ord]
      delta <- delta[ord]      
    }
    tmp <- EM(dataset, weight, tol, nbinit)
    nbparam <- tmp$nbparam + nbparam
    loglike <- tmp$loglike + loglike
    error <- tmp$error + error
    if (k==1){
      param <- cbind(tmp$param, delta, rep(k, length(delta)))
      colnames(param) <- c("alpha", "epsilon", "delta", "groupe")
    } 
    else
      param <- rbind(param, cbind(tmp$param, delta, rep(k, length(delta))))
  }
  
  output <- list(loglike=loglike, nbparam=nbparam, BIC=loglike-nbparam*0.5*log(nrow(x)), param=param, groupe=gr, error=(error>0))
  return(output)  
}





ProbaPostModel <- function(x, param){
  name <- colnames(x)[which(is.na(x[1,]))]
  groupe <- param[which(rownames(param)==name),4]
  param <- param[which(param[,4]==groupe),]
  buildx0 <- buildx1 <- matrix(NA, nrow(x), nrow(param))
  for (j in 1:ncol(buildx0))  buildx0[,j] <- buildx1[,j] <- x[, which(colnames(x)==rownames(param)[j])]
  if (any(param[,3]==-1)) buildx0[,which(param[,3]==-1)] <- buildx1[,which(param[,3]==-1)] <- 1-buildx1[,which(param[,3]==-1)]
  colnames(buildx0) <- colnames(buildx1) <- rownames(param)
  buildx0[, which(colnames(buildx0)==name)] <- 0
  buildx1[, which(colnames(buildx1)==name)] <- 1
  output <- ComputeProba(buildx1, param) / ComputeProba(buildx0, param)
  if (param[which(rownames(param)==name),3]==-1) output <- 1/output
  return(output)
}

ComputeProba <- function(x, param){
  matInf <- matrix(0, length(param[,1]) , length(param[,1]) +1)
  matInf[upper.tri(matInf)] <- 1
  matSup <- 1-matInf
  
  baseSup <- log(sweep(x, 2,(1-param[,2])*param[,1] + param[,2], "*") + sweep(1-x, 2, (1-param[,2])*(1-param[,1]), "*"))
  baseInf <- log(sweep(x, 2,(1-param[,2])*param[,1], "*") + sweep(1-x, 2, (1-param[,2])*(1-param[,1])  + param[,2], "*"))
  fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
  # Tools for loglikelihood
  repere <- (c(param[,1], 1)-c(0, param[,1]))
  loglikeprec <- -Inf
  # Loglikelihood computation
  probaInd <- fij %*% repere
  return( probaInd )
}
