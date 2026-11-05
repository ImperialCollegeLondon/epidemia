vector[ac_nproc] ac_scale = ac_scale_raw .* ac_prior_scales;
vector[ac_q] ac_beta;
vector[ac_q] ac_steps = ac_noise; // the steps of the walks (noise times scale if non-centred)

vector[obs_ac_nproc] obs_ac_scale = obs_ac_scale_raw .* obs_ac_prior_scales; 
vector[obs_ac_q] obs_ac_beta;
vector[obs_ac_q] obs_ac_steps = obs_ac_noise;
