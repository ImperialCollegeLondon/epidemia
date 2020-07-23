int<lower=0> N_obs; // total size of the observation vector
int obs[N_obs]; // vector of observations
int obs_group[N_obs]; // group (1 to M) to which each observation belongs
int obs_type[N_obs]; // type of observation (1 to r). 
int obs_date[N_obs]; // observation date (1 to N2)


// dimensions
int<lower=0> NC_obs[10];  // number of each observation type
int<lower=0> oK[10];  // number of predictors for each observation type
int<lower=0> K_all; // sum of the above
int<lower=0> has_ointercept[10]; // 0 means no, otherwise gives index
int<lower=0> num_ointercepts; // total intercept params

// not pretty, but hopefully more efficient with algorithmic diff
vector[K_all] oxbar;
// dense matrices for each observation type (maximum of 10 types)
matrix[NC_obs[1],oK[1]] oX1; 
matrix[NC_obs[2],oK[2]] oX2;
matrix[NC_obs[3],oK[3]] oX3;
matrix[NC_obs[4],oK[4]] oX4;
matrix[NC_obs[5],oK[5]] oX5;
matrix[NC_obs[6],oK[6]] oX6;
matrix[NC_obs[7],oK[7]] oX7;
matrix[NC_obs[8],oK[8]] oX8;
matrix[NC_obs[9],oK[9]] oX9;
matrix[NC_obs[10],oK[10]] oX10;
