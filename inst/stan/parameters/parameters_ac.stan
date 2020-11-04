vector<lower=0>[ac_nproc] ac_scale_raw; // standard normal version of scale hyperparameter
vector[ac_q] ac_noise; // noise terms for each walk at each time period

vector<lower=0>[obs_ac_nproc] obs_ac_scale_raw; 
vector[obs_ac_q] obs_ac_noise; 