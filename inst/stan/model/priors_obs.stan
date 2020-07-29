target += normal_lpdf(oz_beta | 0, 1);

// prior for the intercepts
if (num_ointercepts > 0) 
    target += normal_lpdf(ogamma | prior_mean_for_ointercept, prior_scale_for_ointercept);

// priors for auxiliary variables
for (i in 1:num_oaux) {
if (prior_dist_for_oaux[i] == 1) 
    target += normal_lpdf(oaux_raw[i] | 0, 1);
else if (prior_dist_for_oaux[i] == 2)
    target += student_t_lpdf(oaux_raw[i] | prior_df_for_oaux[i], 0, 1);
else if (prior_dist_for_oaux[i] == 3)
    target += exponential_lpdf(oaux_raw[i] | 1);
}
