// priors for model parameters
real<lower=0> r0; // the prior expected value for r0
vector<lower=0>[M] pop;
simplex[NS] si; // fixed serial interval using empirical data



