int<lower=0> N_obs; // total size of the observation vector
int obs[N_obs]; // vector of observations
int obs_group[N_obs]; // group (1 to M) to which each observation belongs
int obs_type[N_obs]; // type of observation (1 to r). 
int obs_date[N_obs]; // observation date (1 to N2)

int<lower=0> R; // number of different observation types
int<lower=0> oN[10];  // number of each observation type
int<lower=0> oK[10];  // number of predictors for each observation type
int<lower=0> K_all; // sum of the above
int<lower=0> num_ointercepts; // total intercept params
int<lower=0, upper=num_ointercepts> has_ointercept[R]; // 0 means no, otherwise gives index
vector<lower=0>[NS] pvecs[R]; // the 'lag' for each type of observation.

// family for each observation type
int<lower=1,upper=2> ofamily[R]; //1:poisson 2:neg_binom
int<lower=1,upper=5> olink[R]; //1:log 2:probit 3:cauchit 4:cloglog 5:identity

// data for auxiliary parameters
int<lower=0> num_oaux; // total number aux params
int<lower=0, upper=num_oaux> has_oaux[R];
int<lower=0, upper=3> prior_dist_for_oaux[num_oaux];
vector[num_oaux] prior_mean_for_oaux;
vector<lower=0>[num_oaux] prior_scale_for_oaux;
vector<lower=0>[num_oaux] prior_df_for_oaux;

// model matrices (maximum of 10 types)
// not pretty, but hopefully more efficient with algorithmic diff
vector[K_all] oxbar;
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
