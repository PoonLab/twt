require(twt)
setwd('~/git/treeswithintrees')
source('pkg/R/loadInputs.R')
settings <- yaml.load_file('tests/fixtures/test.yaml')
input <- MODEL$new(settings)


timeline <- 0
comps <- input$get.compartments()

for (comp in comps){
  extant <- comp$get.extant_lineages() # need to incorporate extant lineages in each compartment
  if (length(extant)>1){
    # get popn.growth.dynamics
    popn.growth <- comp$get.type()$get.popn.growth.dynamics()
    for (piece in popn.growth){
      # draw waiting time for one piece of population growth dynamics
      wait.time <- function(k, alpha, beta){
        u <- runif(1, 0, 1)
        wait <- (1-(1-u)^(beta / choose (extant,2))) * alpha / beta
      }
      timeline = timeline+wait.time
      if (timeline > piece$end){
        timeline <- piece$end
        # and draw waiting time for next piece
      }
      else {
        timeline = timeline+wait.time
      }
    }
  }
}

