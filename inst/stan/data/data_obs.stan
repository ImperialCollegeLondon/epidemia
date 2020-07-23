int<lower=0> N_obs; // total size of the observation vector
int obs[N_obs]; // vector of observations
int obs_group[N_obs]; // group (1 to M) to which each observation belongs
int obs_type[N_obs]; // type of observation (1 to r). 
int obs_date[N_obs]; // observation date (1 to N2)


// dimensions
int<lower=0> NC_obs[10];  // number of each observation type
int<lower=0> K_obs[10];  // number of predictors for each observation type
int<lower=0> K_all; // sum of the above
int<lower=0> has_ointercept[10]; // 0 means no, otherwise gives index
int<lower=0> num_ointercepts; // total intercept params

// not pretty, but hopefully more efficient with algorithmic diff
vector[K_all] oxbar;
// dense matrices for each observation type (maximum of 10 types)
matrix[NC_obs[1],K_obs[1]] oX1; 
matrix[NC_obs[2],K_obs[2]] oX2;
matrix[NC_obs[3],K_obs[3]] oX3;
matrix[NC_obs[4],K_obs[4]] oX4;
matrix[NC_obs[5],K_obs[5]] oX5;
matrix[NC_obs[6],K_obs[6]] oX6;
matrix[NC_obs[7],K_obs[7]] oX7;
matrix[NC_obs[8],K_obs[8]] oX8;
matrix[NC_obs[9],K_obs[9]] oX9;
matrix[NC_obs[10],K_obs[10]] oX10;
