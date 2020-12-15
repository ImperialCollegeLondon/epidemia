// hyperparameter values are set to 0 if there is no prior
vector<lower=0>[K] prior_scale;
real<lower=0> prior_scale_for_intercept[has_intercept];
vector[K] prior_mean;
vector<lower=0>[K] prior_shape;
vector[K] prior_shift;
real prior_mean_for_intercept[has_intercept];
vector<lower=0>[K] prior_df;
real<lower=0> prior_df_for_intercept[has_intercept];
real<lower=0> global_prior_df;     // for hs priors only
real<lower=0> global_prior_scale;  // for hs priors only
real<lower=0> slab_df;     // for hs prior only
real<lower=0> slab_scale;  // for hs prior only
int<lower=2> num_normals[prior_dist == 7 ? K : 0];

// additional hyperparameters for coefficients in obs regressions
vector[K_all] prior_omean;
vector<lower=0>[K_all] prior_oscale;
vector[num_ointercepts] prior_mean_for_ointercept;
vector<lower=0>[num_ointercepts] prior_scale_for_ointercept;

// and also for auxiliary variables
int<lower=0, upper=3> prior_dist_for_oaux[num_oaux];
vector[num_oaux] prior_mean_for_oaux;
vector<lower=0>[num_oaux] prior_scale_for_oaux;
vector<lower=0>[num_oaux] prior_df_for_oaux;

// and also for inf auxiliary variables
int latent = min(inf_family, 1);
int<lower=0, upper=3> prior_dist_for_inf_aux[latent];
real prior_mean_for_inf_aux[latent];
real<lower=0> prior_scale_for_inf_aux[latent];
real<lower=0> prior_df_for_inf_aux[latent];

real<lower=0> prior_scale_for_tau;
vector<lower=0>[ac_nproc] ac_prior_scales; // prior scale for hyperparameter for each walk.
vector<lower=0>[obs_ac_nproc] obs_ac_prior_scales;
