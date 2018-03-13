## test script
require(twt)
setwd('~/git/treeswithintrees')
settings <- yaml.load_file('tests/fixtures/example2.yaml')
test <- MODEL$new(settings)
e <- EventLogger$new()
tips.n.heights <- init.fixed.samplings(test)
init.branching.events(test, e)    # applies only to example1.yaml for now, since they provide a "host tree" w/ transmission events
e <- generate.transmission.events(test, e)
transm.tree <- .to.transmission.tree(e)

e <- EventLogger$new()
sapply(1:nrow(file), function(x) {
  e$add.event(name=file[x,'event.type'], 
              time=file[x,'time'], 
              obj1=file[x,'lineage1'], 
              obj2=file[x,'compartment1'], 
              obj3=file[x,'compartment2'])
})
e.mixed <- EventLogger$new()
sapply(1:nrow(mixed.file), function(x) {
  e.mixed$add.event(name=mixed.file[x,'event.type'], 
                    time=mixed.file[x,'time'], 
                    obj1=mixed.file[x,'lineage1'], 
                    obj2=mixed.file[x,'compartment1'], 
                    obj3=mixed.file[x,'compartment2'], 
                    cumulative=F)
})


x <- EventLogger$new()
x$add.event('transmission',  2.597914e-04, NA, 'US_host_8', 'US_host_2')
x$add.event('transmission',  8.423922e-04, NA, 'US_host_10', 'US_host_8')
