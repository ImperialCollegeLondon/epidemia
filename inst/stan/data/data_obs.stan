int<lower=0> N_obs; // total size of the observation vector
array[N_obs] int obs_group; // group (1 to M) to which each observation belongs
array[N_obs] int obs_type; // type of observation (1 to r). 
array[N_obs] int obs_date; // observation date (1 to N2)

int<lower=0> R; // number of different observation types
array[10] int<lower=0> oN;  // number of each observation type
array[R] int<lower=1> pvecs_len; // maximum lag for each i2o distribution
array[R] vector<lower=0>[NS] pvecs; // the 'i2o' for each type of observation

array[R] int<lower=0, upper=1> has_offset;
vector[N_obs] offset_;

// family for each observation type
array[R] int<lower=1,upper=5> ofamily; //1:poisson 2:neg_binom 3:quasi-poisson 4:normal 5:log_normal
array[R] int<lower=1,upper=5> olink; //1:log 2:probit 3:cauchit 4:cloglog 5:identity
  
// data for auxiliary parameters
int<lower=0> num_oaux; // total number aux params
array[R] int<lower=0, upper=num_oaux> has_oaux;
