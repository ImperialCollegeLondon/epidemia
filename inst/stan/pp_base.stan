data {
  int <lower=1> M; // number of countries
  int <lower=1> N0; // number of days for which to impute infections
  int <lower=1> starts[M]; // the start index of each group
  int<lower=1> N; 
  int<lower=1> NC[M]; // days of observed data for each group.
  int<lower=1> N2; // total period for the simulation
  int<lower=1> NS; // maximum number of simulation days for any given group
  int<lower=0> R; // number of different observation types
  matrix[NS, R] P; // the 'pvec' for each type of observation.
  matrix[M, R] means; // mean values for observed events / total cases (for example, IFR)
  real pop[M];
  real SI[NS]; // fixed SI using empirical data
}

transformed data {
  vector[NS] SI_rev; // SI in reverse order
  vector[NS] P_rev[M]; // P in reversed order

  for(i in 1:NS)
    SI_rev[i] = SI[NS-i+1];

  for(m in 1:M){
    for(i in 1:NS) {
      P_rev[m, i] = P[N2-i+1,m];
    }
  }
}

parameters {
  real<lower=0> mu[M];
  real<lower=0> y[M];
  real<lower=0> phi;
  real<lower=0> noise[M, R];
  vector[N] eta;
}

generated quantities {
  vector[N] R0_vec;
  vector[N] Rt_vec;
  matrix[N2, M] Rt = rep_matrix(0,N2,M);
  matrix[N2, M] prediction = rep_matrix(0,N2,M);
  matrix[N2, M] E_obs[R];

  // initialise to 0
  for (r in 1:R)
    E_obs[r] = rep_matrix(0, N2, M);

  { // generate vector of R0s from group means
    int idx = NC[1]+1;
    R0_vec[1:NC[1]] = rep_vector(mu[1], NC[1]);
    for (m in 2:M) {
      R0_vec[idx:(idx+NC[m]-1)] = rep_vector(mu[m], NC[m]);
      idx += NC[m];
    }
  }

  // todo: Add branching logic for different link functions.
  Rt_vec = R0_vec * 2 .* inv_logit(-eta);

  { // predict cases over time
    int idx=1;
    matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
    for (m in 1:M){
      // time indices for group m: start date, final seed date, final date
      int n0 = starts[m];
      int n1 = n0 + N0 - 1;
      int n2 = n0 + NC[m] - 1;

      // impute unadjusted Rt from the linear predictor
      Rt[n0:n2,m] = Rt_vec[idx:(idx+NC[m]-1)];
      idx += NC[m];

      prediction[n0:n1,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
      cumm_sum[n0:n1,m] = cumulative_sum(prediction[n0:n1,m]);

      for(i in (n1+1):n2) {
        real convolution = dot_product(sub_col(prediction, n0, m, i-n0), tail(SI_rev, i-n0));
        prediction[i,m] = (pop[m] - cumm_sum[i-1,m]) * (1 - exp(-Rt[i,m] * convolution / pop[m]));
        cumm_sum[i,m] = cumm_sum[i-1,m] + prediction[i,m];
      }
    }
  }
  
    // simulate from posterior predictive
    for (r in 1:R) {
        for (m in 1:M) {
            int n0 = starts[m];
            int n1 = n0 + NC[m] - 1;
            E_obs[r][n0,m] = 1e-15 * prediction[n0,m];
            for (i in (n0+1):n1) {
                E_obs[r][i,m] = noise[m, r] * means[m, r] * dot_product(sub_col(prediction, n0, m, i-n0), tail(P_rev[r], i-n0));
            }
        }
    }
}

