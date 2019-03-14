data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  //int<lower=1> Nclusters2; // number of levels for group 2 for random intercepts
  int<lower=1> NScaleBatches1; // number of batches of scale terms 
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  //int<lower=1, upper=Nclusters2> NgroupMap2[Ntotal]; // mapping variable to map each observation to group 2 
  int<lower=1, upper=NScaleBatches1> NBatchMap1[Nclusters1]; // expanding vector to map each variance term to the relevant 
                                                    // set of coefficients from rGroupsJitter1
  int y[Ntotal]; // response variable
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
  real intercept;
  real intercept_sd;
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  real betas; // constant intercept term
  real<lower=0.01> sigmaRan1[NScaleBatches1]; // random effect standard deviations for sub-batches in group 1
  //real<lower=0.01> sigmaRan2; // random effect standard deviation for group 2
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
  //vector[Nclusters2] rGroupsJitter2; // number of random jitters for each level of cluster/group 2
  real<lower=0.1> iSize; // size parameter for the nb distribution
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = exp(betas + rGroupsJitter1[NgroupMap1]);// + rGroupsJitter2[NgroupMap2]);
}
model {
  real sigmaRan1_expanded[Nclusters1]; 
  sigmaRan1 ~ gamma(gammaShape, gammaRate);
  //sigmaRan2 ~ gamma(gammaShape, gammaRate);
  betas ~ normal(intercept, intercept_sd);
  // random effects sample
  sigmaRan1_expanded = sigmaRan1[NBatchMap1];
  rGroupsJitter1 ~ normal(0, sigmaRan1_expanded);
  //rGroupsJitter2 ~ normal(0, sigmaRan2);
  // likelihood function
  y ~ neg_binomial_2(mu, iSize);
}
