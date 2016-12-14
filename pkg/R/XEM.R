feps <- function(eps, al, u1, u2, u3, u4)     u1*log((1-eps)*al + eps)               + u2*log((1-eps)*al)           + u3*log((1-eps)*(1-al))               + u4*log(1-al + eps*al)
dfeps <- function(eps, al, u1, u2, u3, u4)    u1*(1-al)/((1-eps)*al + eps)           - u2*al/((1-eps)*al)           + u3*(al-1)/((1-eps)*(1-al))           + u4*al/(1-al + eps*al)
ddfeps <- function(eps, al, u1, u2, u3, u4)  -u1*((1-al)**2)/(((1-eps)*al+ eps)**2)  - u2*(al**2)/(((1-eps)*al)**2) - u3*((al-1)**2)/(((1-eps)*(1-al))**2) - u4*(al**2)/((1-al + eps*al)**2)

NRepsfree <- function(eps,al, a, b, c, d){
  eps1 <- eps0 <- eps
  # car delta 1: mettre ordre a, b, c, d
  prec <- -Inf
  actu1 <- feps(eps1, al, a, b, c, d)
  while ((actu1-prec)>0.0001) {eps1 <- max(0.01, min(0.99, eps1 - dfeps(eps1, al, a, b, c, d)/ddfeps(eps1, al, a, b, c, d))); prec <- actu1; actu1 <-feps(eps1, al, a, b, c, d);}
  # car delta 1: mettre ordre b, a, d, c
  prec <- -Inf
  actu0 <- feps(eps0, al, b, a, d, c)
  while (max(actu0-prec)>0.0001) {eps0 <- max(0.01, min(0.99, eps0 - dfeps(eps0, al, b, a, d, c)/ddfeps(eps0, al, b, a, d, c))); prec <- actu0; actu0 <-feps(eps0, al, b, a, d, c);}
  list(delta=(actu1>actu0), eps=eps0*(actu0>=actu1)+eps1*(actu1>actu0))
}


NRepsequal <- function(eps1, al, a, b, c, d){
  # car delta 1: mettre ordre a, b, c, d
  prec <- -Inf
  actu1 <- sum(feps(eps1, al, a, b, c, d))
  while ((actu1-prec)>0.0001) {eps1 <- max(0.01, min(0.99, eps1 - sum(dfeps(eps1, al, a, b, c, d))/sum(ddfeps(eps1, al, a, b, c, d)))); prec <- actu1; actu1 <- sum(feps(eps1, al, a, b, c, d));}
  # car delta 1: mettre ordre b, a, d, c
  rep(eps1, length(al))
}


OneXemBlockfree <- function(x, weight, alpha, tol){
  db <- length(alpha)
  n <- sum(weight)
  # Initialization
  delta <- sample(0:1, db, replace = TRUE)
  epsilon <- runif(db)
  # Computation of the posterior probabilities: beta corresponds to the bounds of the integral where the likelihood is constant
  beta <- alpha * delta + (1-alpha) * (1-delta)
  ord <- order(beta, decreasing = FALSE)
  beta <- beta[ord]
  # fb sont les valeurs des constantes pour chaque individus
  matSup <- matrix(0, db , db+1)
  matSup[upper.tri(matSup)] <- 1
  matInf <- 1-matSup
  # Proba per individuals
  partinf <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * delta[ord] 
  partsup <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * (1-delta[ord])   
  baseSup <- log(sweep(x[,ord], 2, partsup, "*") + sweep(1-x[,ord], 2, 1-partsup, "*"))
  baseInf <- log(sweep(x[,ord], 2, partinf, "*") + sweep(1-x[,ord], 2, 1 - partinf, "*"))
  fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
  repere <- (c(beta, 1)-c(0, beta))
  probaInd <- fij %*% repere
  # Loglikelihood computation
  loglikeprec <- -Inf
  loglikeactu <- sum(log(probaInd)*weight)
  compteur <- 0
  if (is.na(loglikeactu)) loglikeactu <- loglikeprec <- 10**(-8)
  while ((loglikeactu - loglikeprec)>tol){
    compteur <- compteur + 1
    # E step
    tik <- sweep(sweep(fij,2,repere,"*"), 1, probaInd, "/")  
    # M step
    tmp2 <- lapply(1:db, function(j){
      if (j==1) tij <- tik[,1] else tij <- rowSums(tik[,1:j])
      NRepsfree(epsilon[ord[j]],alpha[ord[j]], sum(x[,ord[j]]*weight*tij), sum(x[,ord[j]]*weight*(1-tij)),  sum((1-x[,ord[j]])*weight*tij), sum((1-x[,ord[j]])*weight*(1-tij)))
    })
    epsilon[ord] <- unlist(lapply(1:db, function(b) tmp2[[b]]$eps))
    delta[ord] <- unlist(lapply(1:db, function(b) tmp2[[b]]$delta))
    # Computation of the posterior probabilities:
    beta <- alpha * delta + (1-alpha) * (1-delta)
    ord <- order(beta, decreasing = FALSE)
    beta <- beta[ord]
    # Proba per individuals
    partinf <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * delta[ord] 
    partsup <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * (1-delta[ord])   
    baseSup <- log(sweep(x[,ord], 2, partsup, "*") + sweep(1-x[,ord], 2, 1-partsup, "*"))
    baseInf <- log(sweep(x[,ord], 2, partinf, "*") + sweep(1-x[,ord], 2, 1 - partinf, "*"))
    fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
    repere <- (c(beta, 1)-c(0, beta))
    probaInd <- fij %*% repere
    # Loglikelihood computation
    loglikeprec <- loglikeactu
    loglikeactu <- sum(log(probaInd)*weight)
  }
  #print(loglikeactu-loglikeprec)
  if (delta[1]==0) delta <-  1- delta
  nbparam <- 2 * ncol(x)
  return(list(epsilon=epsilon, delta=delta, loglike=loglikeactu, nbparam=nbparam, bic=loglikeactu - nbparam * 0.5 * log(sum(weight))))
}




OneXemBlockEqual <- function(x, weight, alpha, tol){
  db <- length(alpha)
  n <- sum(weight)
  # Initialization
  delta <- rep(1,db)
  epsilon <- rep(runif(1), db)
  # Computation of the posterior probabilities: beta corresponds to the bounds of the integral where the likelihood is constant
  beta <- alpha * delta + (1-alpha) * (1-delta)
  ord <- order(beta, decreasing = FALSE)
  beta <- beta[ord]
  # fb sont les valeurs des constantes pour chaque individus
  matSup <- matrix(0, db , db+1)
  matSup[upper.tri(matSup)] <- 1
  matInf <- 1-matSup
  # Proba per individuals
  partinf <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * delta[ord] 
  partsup <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * (1-delta[ord])   
  baseSup <- log(sweep(x[,ord], 2, partsup, "*") + sweep(1-x[,ord], 2, 1-partsup, "*"))
  baseInf <- log(sweep(x[,ord], 2, partinf, "*") + sweep(1-x[,ord], 2, 1 - partinf, "*"))
  fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
  repere <- (c(beta, 1)-c(0, beta))
  probaInd <- fij %*% repere
  # Loglikelihood computation
  loglikeprec <- -Inf
  loglikeactu <- sum(log(probaInd)*weight)
  if (is.na(loglikeactu)) loglikeactu <- loglikeprec <- 10**(-8)
  while ((loglikeactu - loglikeprec)>tol){
    # E step
    tij <- t(apply(sweep(sweep(fij,2,repere,"*"), 1, probaInd, "/"), 1, cumsum))
    # M step
    epsilon <- rep(NRepsequal(epsilon[1], alpha[ord], colSums(x[,ord]*(tij[,1:db]*weight)), colSums(x[,ord]*((1-tij[,1:db])*weight)), colSums((1-x[,ord])*(tij[,1:db]*weight)), colSums((1-x[,ord])*((1-tij[,1:db])*weight))), db)
    # Computation of the posterior probabilities:
    beta <- alpha * delta + (1-alpha) * (1-delta)
    ord <- order(beta, decreasing = FALSE)
    beta <- beta[ord]
    # Proba per individuals
    partinf <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * delta[ord] 
    partsup <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * (1-delta[ord])   
    baseSup <- log(sweep(x[,ord], 2, partsup, "*") + sweep(1-x[,ord], 2, 1-partsup, "*"))
    baseInf <- log(sweep(x[,ord], 2, partinf, "*") + sweep(1-x[,ord], 2, 1 - partinf, "*"))
    fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
    repere <- (c(beta, 1)-c(0, beta))
    probaInd <- fij %*% repere
    # Loglikelihood computation
    loglikeprec <- loglikeactu
    loglikeactu <- sum(log(probaInd)*weight)
  }
  #print(loglikeactu-loglikeprec)
  if (delta[1]==0) delta <-  1- delta
  nbparam <- 1 + ncol(x)
  return(list(epsilon=epsilon, delta=delta, loglike=loglikeactu, nbparam=nbparam, bic=loglikeactu - nbparam * 0.5 * log(sum(weight))))
}


XEMblock <- function(x, weight, alpha, tol, nbiter){
  allblock <- c(lapply(as.list(1:nbiter), function(it) OneXemBlockfree(x, weight, alpha, tol)), lapply(as.list(1:nbiter), function(it) OneXemBlockEqual(x, weight, alpha, tol)))
  bic <- unlist(lapply(1:length(allblock), function(it)  allblock[[it]]$bic))
  return(allblock[[which.max(bic)]])
}


XEMmodel <- function(dataset, alpha, tol, nbinit.EM, model){
  resEM <- lapply(unique(model), function(b){
    if (sum(model==b)>1){
      tmp <- uniquecombs(dataset[,which(model==b)])
      tmp <- XEMblock(as.matrix(tmp), as.numeric(table(attr(tmp,"index"))), alpha[which(model==b)], tol, nbinit.EM)
    }else{
      tmp <- list(delta=1, epsilon=0, nbparam=1, loglike=log(alpha[which(model==b)])*sum(dataset[,which(model==b)]) + log(1-alpha[which(model==b)])*sum(1-dataset[,which(model==b)]))
    }
    tmp
  })
  epsilon <- delta <- rep(NA, ncol(dataset))
  loglike <- sum(unlist(lapply(unique(model), function(b) resEM[[b]]$loglike)))
  epsilon[order(model)] <- unlist(lapply(unique(model), function(b) resEM[[b]]$epsilon))
  delta[order(model)] <- unlist(lapply(unique(model), function(b) resEM[[b]]$delta))
  nbparam <- sum(unlist(lapply(unique(model), function(b) resEM[[b]]$nbparam)))
  return(list(delta=delta, epsilon=epsilon, alpha=alpha, loglike=loglike, nbparam=nbparam, bic=loglike - nbparam*0.5*log(nrow(dataset)), model=model))
}

XEMblock2 <- function(x, alpha, tol, nbinit.EM, model){
  resEM <- list()
  if (length(model)>1){
    tmp <- uniquecombs(x[,model])
    resEM <- XEMblock(as.matrix(tmp), as.numeric(table(attr(tmp,"index"))), alpha[model], tol, nbinit.EM)
  }else{
    loglike <- log(alpha[model])*sum(x[,model]) + log(1-alpha[model])*sum(1-x[,model])
    resEM <- list(delta=1, epsilon=0, nbparam=1, loglike=loglike, bic=loglike-0.5*log(nrow(x)))
  }
  resEM$model <- model
  return(resEM)
}