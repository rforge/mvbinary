setClass(
  Class = "MvBinaryResult", 
  representation = representation(
    alpha="numeric",
    epsilon="numeric",
    delta="numeric",
    blocks="numeric",
    nbparam="numeric",
    loglike="numeric",
    bic="numeric"
  ), 
  prototype = prototype(
    alpha=numeric(),
    epsilon=numeric(),
    delta=numeric(),
    blocks=numeric(),
    nbparam=numeric(),
    loglike=numeric(),
    bic=numeric()
  )
)
