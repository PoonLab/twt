# treeswithintrees

An R package (twt) for the coalescent (reverse time) simulation of pathogen trees within host transmission trees.

## Description
[treeswithintrees](http://github.com/PoonLab/twt) (twt) is an R package for the [coalescent](https://en.wikipedia.org/wiki/Coalescent_theory) (reverse time) simulation of host (transmission) trees and pathogen trees within host trees.  It is designed to be flexible so that it can accommodate a range of models at different levels of diversity, such as:
* [compartmental epidemic models](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology), *e.g.,* susceptible-infected-recovered (SIR) models
* [cospeciation](https://en.wikipedia.org/wiki/Cospeciation)/cophylogeny models
* models of [compartmentalized within-host evolution](https://veg.github.io/hyphy-site/resources/compartmentalization_detection_ppt.pdf), *e.g.,* migration between the blood and genital tract

twt attempts to accommodate these models with a common set of basic components: [compartments](Compartment) and [lineages](Lineage).  A [lineage](Lineage) is a sequence of individual pathogens descending from an ancestor in the past.  A [compartment](Compartment) represents an individual environment in which lineages are contained.  [Compartments](Compartment) are grouped into [Compartment Types](CompartmentType) to make it more convenient to specify models and for more efficient simulation.  We assume that the compartments are related through by a transmission or host tree, whose shape is determined by the transmission rates among hosts.  Lineages may also migrate between compartments that are both hosts to other lineages.  

Like many other simulation programs, twt uses a conventional [Gillespie method](https://en.wikipedia.org/wiki/Gillespie_algorithm) to sample a sequence of discrete stochastic events over reverse time.  If a host tree is not specified by the user, then the host/transmission tree is simulated while simulating the coalescence of the pathogen lineages within the hosts.  Otherwise, the timing and direction of transmission events are parsed from the user tree and set as [fixed events](Events) in the simulation.


## Funding

Development of treeswithintrees was directly supported by a grant from the Government of Canada through Genome Canada and the Ontario Genomics Institute (OGI-131) and by the Canadian Institutes of Health Research (project grant PJT-155990).
