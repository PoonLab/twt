## test script
require(twt)
set.seed(34)
setwd('~/git/treeswithintrees')
settings <- yaml.load_file('tests/fixtures/example4.yaml')
test <- MODEL$new(settings)
fixed.samplings <- init.fixed.samplings(test)

e <- EventLogger$new()
sim.outer.tree(test, e)
sim.inner.tree(test, e)

# transm.tree <- .outer.tree.to.phylo(o)
phy <- .inner.tree.to.phylo(e, fixed.samplings)
plot(phy)

# tips.n.heights <- init.fixed.samplings(test)
# init.branching.events(test, e)    # applies only to example1.yaml for now, since they provide a "host tree" w/ transmission events


