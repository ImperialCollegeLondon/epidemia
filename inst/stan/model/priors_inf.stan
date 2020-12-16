if (latent) {

    // prior for coef dispersion for inf
    if (prior_dist_for_inf_aux[1] == 1) 
        target += normal_lpdf(inf_aux_raw[1] | 0, 1);
    else if (prior_dist_for_inf_aux[1] == 2)
        target += student_t_lpdf(inf_aux_raw[1] | prior_df_for_inf_aux[1], 0, 1);
    else if (prior_dist_for_inf_aux[1] == 3)
        target += exponential_lpdf(inf_aux_raw[1] | 1);

    for (m in 1:M) {
        int n0 = starts[m];
        int n1 = n0 + N0 - 1;
        int n2 = n0 + NC[m] - 1;
        vector[n2-n1] mu = mean_inf[(n1+1):n2,m];

        if (pop_adjust) 
            mu = (pop[m] - cumm_sum[n1:(n2-1),m]) .* (1 - exp(- mu / pop[m]));
        
        target += gamma_lpdf(infections[(n1+1):n2,m] | mu / inf_aux[1], 1 / inf_aux[1]);

        // ADD JACOBIAN ADJUSTMENT HERE
    }
}

