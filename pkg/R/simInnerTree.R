# simulation of the inner tree
simulate <- function(model, eventlog) {
  comps <- model$get.compartments()
  compnames <- model$get.names(comps)
  
  extant_l <- model$get.extant_lineages()
  not_extant_l <- list()
  not_yet_sampled_l <- list()
  
  extant_c <- unique(sapply(extant_l,function(x){x$get.location()}))
  not_extant_c <- list()
  not_yet_sampled_c <- list()
}