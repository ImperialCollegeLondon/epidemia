if (latent) {

    // prior for coef dispersion for inf
    if (prior_dist_for_inf_aux[1] == 1) 
        target += normal_lpdf(inf_aux_raw[1] | 0, 1);
    else if (prior_dist_for_inf_aux[1] == 2)
        target += student_t_lpdf(inf_aux_raw[1] | prior_df_for_inf_aux[1], 0, 1);
    else if (prior_dist_for_inf_aux[1] == 3)
        target += exponential_lpdf(inf_aux_raw[1] | 1);

    target += normal_lpdf(inf_noise | 0, 1);
}

