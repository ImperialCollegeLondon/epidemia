if (latent) {

    { 
    int idx = 1;
    for (m in 1:M){
      // time indices for group m: start date, final seed date, final date
      int n0 = starts[m];
      int n1 = n0 + N0 - 1;
      int n2 = n0 + NC[m] - 1;
      vector[n2-n1] mu = E_infections[(n1+1):n2,m];
      vector[n2-n1] sigma = inf_aux[1] * E_infections[(n1+1):n2,m];
      if (fixed_vtm) sigma = sqrt(sigma);
      target += normal_lpdf(infections_raw[idx:(idx+NC[m]-N0-1)]| mu, sigma);
      idx += NC[m]-N0;
    }
    
    }

    // prior for coef dispersion for inf
    if (prior_dist_for_inf_aux[1] == 1) 
        target += normal_lpdf(inf_aux_raw[1] | 0, 1);
    else if (prior_dist_for_inf_aux[1] == 2)
        target += student_t_lpdf(inf_aux_raw[1] | prior_df_for_inf_aux[1], 0, 1);
    else if (prior_dist_for_inf_aux[1] == 3)
        target += exponential_lpdf(inf_aux_raw[1] | 1);
}


if (!I0_fixed) {
    target += normal_lpdf(I0 | prior_mean_for_I0, prior_scale_for_I0);
}
