# OccupancyModeling_Stan
Species/Site occupancy modeling with Stan in R

The code included here simulates presence/absence data for some number of species/sites/replicates, and then uses the simulated data to estimate the occupancy modeling parameters. 

Lahoz-Monfort et al (2016) applied the Royle & Link (2006) occupancy-modeling framework to environmental DNA. Here, I've taken their code for Bayesian estimation of the relevant parameters and ported it to Stan from JAGS.  I have also generalized it to apply to multiple species or multiple sites per species. 

Comments welcome, of course. 

[Lahoz‐Monfort, José J., Gurutzeta Guillera‐Arroita, and Reid Tingley. "Statistical approaches to account for false‐positive errors in environmental DNA samples." Molecular Ecology Resources 16.3 (2016): 673-685.]


I have also included a single-species demonstration (Rmd file and accompanying HTML), which walks through the logic of occupancy modeling and fits parameters for a hypothetical case. 
