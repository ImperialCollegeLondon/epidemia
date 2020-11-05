// priors for model parameters
real<lower=0> r0; // r0 is only used if link == 2
int<lower=1, upper=3> link; // the link function 1) log, 2) scaled_logit, 3) identity
vector<lower=0>[M] pop;
int<lower=1> si_len;
simplex[NS] si; // fixed serial interval using empirical data
int<lower=0, upper=1> pop_adjust; // whether to perform population adjustment or not




