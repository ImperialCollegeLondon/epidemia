// priors for model parameters
real<lower=0> r0; // r0 is only used if link == 2
int<lower=1, upper=3> link; // the link function 1) log, 2) scaled_logit, 3) identity
matrix<lower=0>[N2, M] susc; // susceptible population over time
int<lower=1> gen_len;
simplex[NS] gen; // fixed generation distribution using empirical data
int<lower=0, upper=1> pop_adjust; // whether to perform population adjustment or not




