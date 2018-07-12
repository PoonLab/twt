## test script
require(twt)
setwd('~/git/treeswithintrees')
settings <- yaml.load_file('tests/fixtures/example2.yaml')
test <- MODEL$new(settings)
e <- EventLogger$new()

o <- sim.outer.tree(test, e)
i <- sim.inner.tree(test, o)





transm.tree <- .to.transmission.tree(o)

tips.n.heights <- init.fixed.samplings(test)
init.branching.events(test, e)    # applies only to example1.yaml for now, since they provide a "host tree" w/ transmission events