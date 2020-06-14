
# priors for model parameters
vector[R] prior_mean_for_phi;
vector<lower=0>[R] prior_scale_for_phi;
vector[M] prior_mean_for_mu;
vector<lower=0>[M] prior_scale_for_mu;
real<lower=0> prior_rate_for_tau;
simplex[NS] pvecs[R]; // the 'pvec' for each type of observation.
vector<lower=0>[M] pop;
simplex[NS] si; // fixed serial interval using empirical data



