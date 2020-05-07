library(rlang)
methods::setClass("simulation", 
         slots = c(
           "data" = "data.frame",
           "param" = "list",
           "type" = "vector"
         ))

# Getters and setters   
methods::setGeneric("data", function(x) standardGeneric("data"))
methods::setGeneric("data<-", function(x, value) standardGeneric("data<-"))
methods::setMethod("data", "simulation", function(x) x@data)
methods::setMethod("data<-", "simulation", function(x, value){
  x@data <- value
  methods::validObject(x)
  x
})


methods::setGeneric("param", function(x) standardGeneric("param"))
methods::setGeneric("param<-", function(x, value) standardGeneric("param<-"))
methods::setMethod("param", "simulation", function(x) x@param)
methods::setMethod("param<-", "simulation", function(x, value){
  x@param <- value
  methods::validObject(x)
  x
})

methods::setGeneric("type", function(x) standardGeneric("type"))
methods::setGeneric("type<-", function(x, value) standardGeneric("type<-"))
methods::setMethod("type", "simulation", function(x) x@type)
methods::setMethod("type<-", "simulation", function(x, value){
  x@param <- value
  methods::validObject(x)
})

# Validators 
methods::setValidity("simulation", function(object){
  if (is.list(simulation@data) == F | is.data.frame(simulation@data[[1]]) == F){
    "@data has to be a list of type data frame"
  } else if (is.list(simulation@param) == F){
    "@param has to be of type list"
  } else if (is.vector(simulation@type) == F){
    "@type has to be a character vector"
  } else if (is.numeric(simulation@type) == T){
    "@type has to be a character vector"
  } else {
    TRUE
  }
})

# show method
methods::setMethod("show", "simulation", function(object){
  cat(is(object)[[1]], "\n",
      "  Data:  ", length(object@data), " data sets of [", nrow(object@data[[1]]), " samples, ", ncol(object@data[[1]]), " variables]", "\n",
      "  Parameters:  ", length(object@param), " parameters", "\n",
      "  Type:  ", object@type, "\n",
      sep = ""
  )
})

# adding elements method
methods::setGeneric("adding_elements", function(object, value, slot) standardGeneric("update_elements"))
methods::setMethod("adding_elements", "simulation", function(object, slot, value){
  if (slot == "param"){
    if (is.list(value) != T | is.null(names(value))){
      stop("Updating @param needs a named list")
    }
    object@param <- rlist::list.append(object@param, value)
  } else if (slot == "type") {
    if (is.numeric(value)){
      stop("Updating @type needs a string")
    }
    object@type <- rlist::list.apend(object@type, value)
  } else if (slot == "data") {
    if (is.data.frame(value) != T){
      stop("Updating @data needs a data frame object")
    }
    object@data <- rlist::list.append(object@type, value)
  } else {
    stop("Unrecognized slot")
  }
})


# constructor for inital simulation object 
simulation <- function(data, param, type){
  obj <- new("simulation", data = data, param = param, type = type)
  return(obj)
}










