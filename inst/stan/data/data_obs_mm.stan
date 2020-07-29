
// data relating to model matrices for each observation
int<lower=0> oK[10];  // number of predictors for each observation type
int<lower=0> K_all; // sum of the above
int<lower=0> num_ointercepts; // total intercept params
int<lower=0, upper=num_ointercepts> has_ointercept[R]; // 0 means no, otherwise gives index

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