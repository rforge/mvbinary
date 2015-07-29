setClass(
  Class = "MvBinaryResult", 
  representation = representation(
    allModels="matrix",
    allBIC="numeric",
    dependencies="data.frame",
    bestparam="matrix",
    bestinfo="list",
    error="logical"
  ), 
  prototype = prototype(
    allModels=matrix(),
    allBIC=numeric(),
    dependencies=data.frame(),
    bestparam=matrix(),
    bestinfo=list(),
    error=logical()
  )
)
