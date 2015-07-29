setMethod(
  f="boxplot",
  signature = c("MvBinaryResult"),
  definition = function(x){
    op <- par(no.readonly = TRUE)
    tmp <- data.frame(Cramer=x@dependencies[,3], Model=(x@dependencies[,4]!=0))
    tmp[which(tmp[,2]==FALSE), 2] <- "no"
    tmp[which(tmp[,2]==TRUE), 2] <- "yes"
    tmp[,2] <- as.factor(tmp[,2])
    boxplot(Cramer~Model, data=tmp, xlab="Dependency modelled", ylab="Empiric Cramer's V")
    par(op)
  }
)