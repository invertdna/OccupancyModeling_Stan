data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of sites or samples (nrow)
    int<lower=1> K[S];   // number of replicates per site (ncol)
    int<lower=0> N[S]; // detected or not at that site/replicate
}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0,upper=1> psi;  //commonness parameter
    real<lower=0,upper=1> p11; //true positive detection rate
    real<lower=0,upper=1> p10; //false positive detection rate
}

transformed parameters{/////////////////////////////////////////////////////////////////////
}

model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
    for (i in 1:S){
			p[i] = psi*p11 + (1-psi)*p10;
			N[i] ~ binomial(K[i], p[i]);
	}; 
  
  //priors
  psi ~ beta(5,5); 
  p11 ~ beta(5,5); 
  p10 ~ beta(1,50);
}

generated quantities{
    real<lower=0,upper=1> Occupancy_prob[11];    //after inferring parameters above, now calculate occupancy probability for each possible outcome (0 to 10 detections, out of 10), and give a posterior for the probability that a sample with that outcome came from an occupied habitat. Equation from Lahoz-Monfort et al. 2015

  for (i in 0:10){
       Occupancy_prob[i+1]  = (psi*(p11^i)*(1-p11)^(10-i)) 
       / ((psi*(p11^i)*(1-p11)^(10-i)) 
          + (((1-psi)*(p10^i))*((1-p10)^(10-i)))
         );
      }
}


