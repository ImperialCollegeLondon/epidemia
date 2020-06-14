
# priors for model parameters
vector[R] prior_mean_for_phi;
vector<lower=0>[R] prior_scale_for_phi;
real<lower=0> prior_scale_for_tau;
simplex[NS] pvecs[R]; // the 'pvec' for each type of observation.
vector<lower=0>[M] pop;
simplex[NS] si; // fixed serial interval using empirical data



