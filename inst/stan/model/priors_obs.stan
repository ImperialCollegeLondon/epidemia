target += normal_lpdf(oz_beta | 0, 1);

// prior for the intercepts
if (num_ointercepts > 0) 
    target += normal_lpdf(ogamma | oprior_mean_for_intercept, oprior_scale_for_intercept);
    