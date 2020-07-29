int<lower=0> N_obs; // total size of the observation vector
int obs[N_obs]; // vector of observations
int obs_group[N_obs]; // group (1 to M) to which each observation belongs
int obs_type[N_obs]; // type of observation (1 to r). 
int obs_date[N_obs]; // observation date (1 to N2)

int<lower=0> R; // number of different observation types
int<lower=0> oN[10];  // number of each observation type
int<lower=0> oK[10];  // number of predictors for each observation type
int<lower=0> K_all; // sum of the above
int<lower=0> has_ointercept[R]; // 0 means no, otherwise gives index
int<lower=0> num_ointercepts; // total intercept params

// family for each observation type
int<lower=0,upper=2> ofamily[R]; //1:poisson 2:neg_binom
int<lower=1,upper=5> olink[R]; //1:log 2:probit 3:cauchit 4:cloglog 5:identity

// not pretty, but hopefully more efficient with algorithmic diff
vector[K_all] oxbar;
// dense matrices for each observation type (maximum of 10 types)
matrix[oN[1],oK[1]] oX1; 
matrix[oN[2],oK[2]] oX2;
matrix[oN[3],oK[3]] oX3;
matrix[oN[4],oK[4]] oX4;
matrix[oN[5],oK[5]] oX5;
matrix[oN[6],oK[6]] oX6;
matrix[oN[7],oK[7]] oX7;
matrix[oN[8],oK[8]] oX8;
matrix[oN[9],oK[9]] oX9;
matrix[oN[10],oK[10]] oX10;
