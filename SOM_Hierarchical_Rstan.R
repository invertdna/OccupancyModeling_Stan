
######################
library(tidyverse)
library(here)
library(rstan)
library(shinystan)
library(bayesplot)
library(broom)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(108)
setwd(here())
######################

######################    
#create test data
######################

          Nspecies <- 2
          Ndates <- 3
          Nsites <- 3
          Ntechreplicates <- 3
          
          # mock data, given true parameters and using those to generate observations            
          psi_given <- seq(0.3, 0.9, length.out = Nspecies*Nsites) %>% 
            sample()
          p11_given <- seq(0.7, 0.9, length.out = Nspecies) %>% 
            sample()
          p10_given <- seq(0.01, 0.05, length.out = Nspecies) %>% 
            sample()
          
          #set up data with different hierarchical levels, in case you want to play with these later
          testData <- expand.grid(
            "Species" = c(1:Nspecies), # model N species
            "Date" = c(1:Ndates), # samples collected on different dates
            "Site" = c(1:Nsites) # at different sites
          ) %>%
            mutate(
              K = Ntechreplicates,
              N = NA
            ) %>%
            arrange(Species)
          
          #create unique identifier for combinations of site-date-species; for use in hierarchical modeling
          SDS <- unite(data = testData,
                       col = SDS,
                       c("Site", "Date", "Species")
          ) %>% pull(SDS)
          testData$SiteDateSpecies <- match(SDS, unique(SDS)) #index for unique site-date-species combinations
          
          #create unique identifier for combinations of site-species; for use in hierarchical modeling
          SS <- unite(data = testData,
                      col = SS,
                      c("Site", "Species")
          ) %>% pull(SS)
          testData$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations
          
          testData <- testData %>% 
            mutate(psi_given = psi_given[testData$SiteSpecies],
                   p11_given = p11_given[testData$Species],
                   p10_given = p10_given[testData$Species])
          
            #Given parameters (psi, p11, p10) assigned above, generate pattern of detections in mock data
            for (i in 1:nrow(testData)){
              if(rbinom(1,1,testData$psi_given[i]) == 0) {testData$N[i] <- 0} else {
                testData$N[i] <- rbinom(1, testData$K[i], testData$p11_given[i])
              }
            }
            
          
          
##Stan Model
       sink("Stan_SOM_hierarchical.stan")
       cat(
         "data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> K[S];   // number of replicates per site (ncol)
    int<lower=0> N[S]; // number of detections among these replicates
    int z[S];   // integer flag to help estimate psi parameter
}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0,upper=1> psi[Nloc];  //commonness parameter
    real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
    real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
}

transformed parameters{/////////////////////////////////////////////////////////////////////
}

model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
    for (i in 1:S){
			z[i] ~ bernoulli(psi[L[i]]);
			p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
			N[i] ~ binomial(K[i], p[i]);
	}; 
  
  //priors
  psi ~ beta(2,2); 
  p11 ~ beta(2,2); 
  p10 ~ beta(1,10);
}

generated quantities{
}

",
fill=TRUE)
       sink()
       
       #####################
       #run Stan model
       #note this will take a while the first time you run a particular model, because it needs to compile from C++
       #####################      
       myHierarchicalModel <- stan(file = "Stan_SOM_hierarchical.stan", 
                             data = list(
                               S = nrow(testData),
                               Species = testData$Species,
                               Nspecies = length(unique(testData$Species)),
                               L = testData$SiteSpecies,
                               Nloc = length(unique(testData$SiteSpecies)),
                               K = testData$K,
                               N = testData$N,
                               z = ifelse(testData$N > 0, 1, 0)
                             ), 
                             chains = 4,   #number of chains
                             iter = 4000   #number of iterations per chain
       )
       
       myHierarchicalStanResults <- tidy(myHierarchicalModel)   
       #launch_shinystan(myHierarchicalModel)
       
       plot(myHierarchicalModel)
