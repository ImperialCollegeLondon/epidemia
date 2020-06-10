vector<lower=(prior_dist == 8 ? 0 : negative_infinity())>[prior_dist == 7 ? sum(num_normals) : K] z_beta;
real<lower=0> global[hs];
vector<lower=0>[K] local[hs];
real<lower=0> caux[hs > 0];
vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
real<lower=0> one_over_lambda[prior_dist == 6];
vector[q] z_b;
vector[len_z_T] z_T;
vector<lower=0,upper=1>[len_rho] rho;
vector<lower=0>[len_concentration] zeta;
vector<lower=0>[t] tau;

