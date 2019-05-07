## test script for 'example4.yaml'
require(twt)
#set.seed(34)
setwd('~/git/treeswithintrees')
# setwd('C:/Users/tng92/git/treeswithintrees')

# load settings into a MODEL object, create an Eventlogger object
settings <- yaml.load_file('tests/fixtures/example4.yaml')
m <- MODEL$new(settings)
e <- EventLogger$new()

# simulation of the outer and inner tree (transmission, migrations, coalescent events)
sim.outer.tree(m, e)
sim.migrations(m, e)
sim.inner.tree(m, e)

# print, plot, and output tree in ape::phylo format
# default format excludes node labels and transmission and migration events
print(e)
plot(e)
phy <- .eventlogger.to.phylo(e)
write.tree(phy)




## test script for 'example1.yaml'
settings <- yaml.load_file('tests/fixtures/example1.yaml')
m <- MODEL$new(settings)
e <- EventLogger$new()

init.branching.events(m, e)   # user provides a "host tree" w/ transmission events
sim.inner.tree(m, e)

print(e)
plot(e)
phy <- .eventlogger.to.phylo(e)
write.tree(phy)

