functions {
#include functions/reverse.stan
}

data {
#include data/data_indices.stan
real<lower=0> r0; // the prior expected value for r0
vector<lower=0,upper=1>[NS] pvecs[R]; // the 'pvec' for each type of observation.
vector<lower=0>[M] pop;
simplex[NS] si; // fixed serial interval using empirical data
int<lower=1> N; 
matrix<lower=0, upper=1>[M, R] means; // mean values for observed events / total cases (for example, IFR)
}

transformed data {
#include tdata/tdata_reverse.stan

for(r in 1:R)
      pvecs_rev[r] = reverse(pvecs[r]);
}

parameters {
  vector<lower=0>[M+1] y;
  real<lower=0> phi[R+1];
  matrix<lower=0>[M+1,R] noise;
  vector[N] eta;
}

generated quantities {
  matrix[N2, M] pred[R];
#include /tparameters/infections_rt.stan

  // initialise to 0
  for (r in 1:R)
    pred[r] = rep_matrix(0, N2, M);

#include /tparameters/gen_infections.stan
  
// simulate from posterior predictive
for (r in 1:R) {
    for (m in 1:M) {
        int n0 = starts[m];
        int n1 = n0 + NC[m] - 1;
        pred[r][n0,m] = 1e-15 * infections[n0,m];
        for (i in (n0+1):n1) {
            pred[r][i,m] = noise[m, r] * means[m, r] * dot_product(sub_col(infections, n0, m, i-n0), tail(pvecs_rev[r], i-n0));
            pred[r][i,m] = neg_binomial_2_rng(pred[r][i,m] + 1e-15, phi[r]);
        }
    }
}
}

