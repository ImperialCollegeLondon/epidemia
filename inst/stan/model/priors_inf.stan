// prior for coef dispersion for inf
if (inf_family > 0) {
    if (prior_dist_for_inf_aux == 1) 
        target += normal_lpdf(inf_aux_raw | 0, 1);
    else if (prior_dist_for_inf_aux == 2)
        target += student_t_lpdf(inf_aux_raw | prior_df_for_inf_aux, 0, 1);
    else if (prior_dist_for_inf_aux == 3)
        target += exponential_lpdf(inf_aux_raw | 1);   
}
