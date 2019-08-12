require(twt)
#setwd('~/git/treeswithintrees')
#source('pkg/R/classes.R')
settings <- yaml.load_file('test.yaml')
test <- MODEL$new(settings)

# create an EventLogger object for testing
e <- EventLogger$new()
e$add.event("transmission", 4, "NA", "I_95", "I_63")
e$add.event("transmission", 6, "NA", "I_73", "I_95")
e$add.event("transmission", 3, "NA", "I_20", "I_73")
e$add.event("transmission", 2, "NA", "I_94", "I_20")
e$add.event("transmission", 1, "NA", "I_97", "I_20")

test.get.all.events <- function(){
  result.noncumul <- e$get.all.events(cumulative = F)
  expected.noncumul <- data.frame(event.type=character(),
                            time=numeric(),
                            lineage1=character(),
                            lineage2=character(),
                            compartment1=character(),
                            compartment2=character()
  )
  expected.noncumul <- rbind(expected.noncumul,
                             list(event.type="transmission", time=4, lineage1="NA", lineage2=NA, compartment1="I_95", compartment2="I_63"),
                             list(event.type="transmission", time=6, lineage1="NA", lineage2=NA, compartment1="I_73", compartment2="I_95"),
                             list(event.type="transmission", time=3, lineage1="NA", lineage2=NA, compartment1="I_20", compartment2="I_73"),
                             list(event.type="transmission", time=2, lineage1="NA", lineage2=NA, compartment1="I_94", compartment2="I_20"),
                             list(event.type="transmission", time=1, lineage1="NA", lineage2=NA, compartment1="I_97", compartment2="I_20"),
                             stringsAsFactors=F)
  row.names(expected.noncumul) <- 1:5
  checkEquals(expected.noncumul,result.noncumul)
  
  result.cumul <- e$get.all.events()
  expected.cumul <- data.frame(event.type=character(),
                                  time=numeric(),
                                  lineage1=character(),
                                  lineage2=character(),
                                  compartment1=character(),
                                  compartment2=character()
  )
  expected.cumul <- rbind(expected.cumul,
                             list(event.type="transmission", time=11, lineage1="NA", lineage2=NA, compartment1="I_95", compartment2="I_63"),
                             list(event.type="transmission", time=5, lineage1="NA", lineage2=NA, compartment1="I_73", compartment2="I_95"),
                             list(event.type="transmission", time=2, lineage1="NA", lineage2=NA, compartment1="I_20", compartment2="I_73"),
                             list(event.type="transmission", time=0, lineage1="NA", lineage2=NA, compartment1="I_94", compartment2="I_20"),
                             list(event.type="transmission", time=1, lineage1="NA", lineage2=NA, compartment1="I_97", compartment2="I_20"),
                             stringsAsFactors=F)
  row.names(expected.cumul) <- as.character(1:5)
  checkEquals(expected.cumul,result.cumul)
}

test.get.events <- function(){
  result.noncumul <- e$get.events("transmission", cumulative=F)
  expected.noncumul <- data.frame(event.type=character(),
                         time=numeric(),
                         lineage1=character(),
                         lineage2=character(),
                         compartment1=character(),
                         compartment2=character()
  )
  expected.noncumul <- rbind(expected.noncumul,
                    list(event.type="transmission", time=4, lineage1="NA", lineage2=NA, compartment1="I_95", compartment2="I_63"),
                    list(event.type="transmission", time=6, lineage1="NA", lineage2=NA, compartment1="I_73", compartment2="I_95"),
                    list(event.type="transmission", time=3, lineage1="NA", lineage2=NA, compartment1="I_20", compartment2="I_73"),
                    list(event.type="transmission", time=2, lineage1="NA", lineage2=NA, compartment1="I_94", compartment2="I_20"),
                    list(event.type="transmission", time=1, lineage1="NA", lineage2=NA, compartment1="I_97", compartment2="I_20"),
                    stringsAsFactors=F)
  row.names(expected.noncumul) <- 1:5
  checkEquals(expected.noncumul,result.noncumul)
  
  result.cumul <- e$get.events("transmission")
  expected.cumul <- data.frame(event.type=character(),
                               time=numeric(),
                               lineage1=character(),
                               lineage2=character(),
                               compartment1=character(),
                               compartment2=character()
  )
  expected.cumul <- rbind(expected.cumul,
                          list(event.type="transmission", time=11, lineage1="NA", lineage2=NA, compartment1="I_95", compartment2="I_63"),
                          list(event.type="transmission", time=5, lineage1="NA", lineage2=NA, compartment1="I_73", compartment2="I_95"),
                          list(event.type="transmission", time=2, lineage1="NA", lineage2=NA, compartment1="I_20", compartment2="I_73"),
                          list(event.type="transmission", time=0, lineage1="NA", lineage2=NA, compartment1="I_94", compartment2="I_20"),
                          list(event.type="transmission", time=1, lineage1="NA", lineage2=NA, compartment1="I_97", compartment2="I_20"),
                          stringsAsFactors=F)
  row.names(expected.cumul) <- 1:5
  checkEquals(expected.cumul,result.cumul)
}
