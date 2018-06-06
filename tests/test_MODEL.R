## test script
require(twt)
setwd('~/git/treeswithintrees')
settings <- yaml.load_file('tests/fixtures/example2.yaml')
test <- MODEL$new(settings)
e <- EventLogger$new()
e <- sim.outer.tree(test, e)
transm.tree <- .to.transmission.tree(e)



tips.n.heights <- init.fixed.samplings(test)
init.branching.events(test, e)    # applies only to example1.yaml for now, since they provide a "host tree" w/ transmission events