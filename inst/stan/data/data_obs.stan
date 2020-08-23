int<lower=0> N_obs; // total size of the observation vector
int obs_group[N_obs]; // group (1 to M) to which each observation belongs
int obs_type[N_obs]; // type of observation (1 to r). 
int obs_date[N_obs]; // observation date (1 to N2)

int<lower=0> R; // number of different observation types
int<lower=0> oN[10];  // number of each observation type
int<lower=1> pvecs_len[R]; // maximum lag for each i2o distribution
vector<lower=0>[NS] pvecs[R]; // the 'i2o' for each type of observation

int<lower=0, upper=1> has_offset[R];
vector[N_obs] offset_;

// family for each observation type
int<lower=1,upper=3> ofamily[R]; //1:poisson 2:neg_binom 3:quasi-poisson
int<lower=1,upper=5> olink[R]; //1:log 2:probit 3:cauchit 4:cloglog 5:identity
  
// data for auxiliary parameters
int<lower=0> num_oaux; // total number aux params
int<lower=0, upper=num_oaux> has_oaux[R];
