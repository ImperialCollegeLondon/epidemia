// hyperparameter values are set to 0 if there is no prior
vector<lower=0>[K] prior_scale;
array[has_intercept] real<lower=0> prior_scale_for_intercept;
vector[K] prior_mean;
vector<lower=0>[K] prior_shape;
vector[K] prior_shift;
array[has_intercept] real prior_mean_for_intercept;
vector<lower=0>[K] prior_df;
array[has_intercept] real<lower=0> prior_df_for_intercept;
real<lower=0> global_prior_df;     // for hs priors only
real<lower=0> global_prior_scale;  // for hs priors only
real<lower=0> slab_df;     // for hs prior only
real<lower=0> slab_scale;  // for hs prior only
array[prior_dist == 7 ? K : 0] int<lower=2> num_normals;

// additional hyperparameters for coefficients in obs regressions
vector[K_all] prior_omean;
vector<lower=0>[K_all] prior_oscale;
vector[num_ointercepts] prior_mean_for_ointercept;
vector<lower=0>[num_ointercepts] prior_scale_for_ointercept;

// and also for auxiliary variables
array[num_oaux] int<lower=0, upper=3> prior_dist_for_oaux;
vector[num_oaux] prior_mean_for_oaux;
vector<lower=0>[num_oaux] prior_scale_for_oaux;
vector<lower=0>[num_oaux] prior_df_for_oaux;

// and also for inf auxiliary variables
array[latent] int<lower=0, upper=3> prior_dist_for_inf_aux;
array[latent] real prior_mean_for_inf_aux;
array[latent] real<lower=0> prior_scale_for_inf_aux;
array[latent] real<lower=0> prior_df_for_inf_aux;

// and seed auxiliary variable
array[hseeds] int<lower=0, upper=3> prior_dist_for_seeds_aux;
vector[hseeds] prior_mean_for_seeds_aux;
vector<lower=0>[hseeds] prior_scale_for_seeds_aux;
vector<lower=0>[hseeds] prior_df_for_seeds_aux;

// and the actual seeds
int<lower=0, upper=9> prior_dist_for_seeds;
vector[M] prior_mean_for_seeds;
vector<lower=0>[M] prior_scale_for_seeds;
vector<lower=0>[M] prior_df_for_seeds;

vector<lower=0>[ac_nproc] ac_prior_scales; // prior scale for hyperparameter for each walk.
vector<lower=0>[obs_ac_nproc] obs_ac_prior_scales;

array[S0_fixed ? 0 : M] real prior_mean_for_S0;
array[S0_fixed ? 0 : M] real<lower=0> prior_scale_for_S0;

array[veps_fixed ? 0 : M] real prior_mean_for_veps;
array[veps_fixed ? 0 : M] real<lower=0> prior_scale_for_veps;

