{ // predict cases over time
int idx=1;
matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
for (m in 1:M){
    // time indices for group m: start date, final seed date, final date
    int n0 = starts[m];
    int n1 = n0 + N0 - 1;
    int n2 = n0 + NC[m] - 1;
    int len;

    if (link == 1) { // log-link
        Rt_unadj[n0:n2,m] = exp(eta[idx:(idx+NC[m]-1)]);
    } else if (link == 2) { // scaled_logit
        Rt_unadj[n0:n2,m] = r0 * 2 * inv_logit(eta[idx:(idx+NC[m]-1)]);
    } else { // identity
        Rt_unadj[n0:n2,m] = eta[idx:(idx+NC[m]-1)];
    }
    
    Rt[n0:n1,m] = Rt_unadj[n0:n1,m]; 
    idx += NC[m];

    infections[n0:n1,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
    cumm_sum[n0:n1,m] = cumulative_sum(infections[n0:n1,m]);

    for (i in (n1+1):n2) {
        int start = max(n0, i - si_len);
        load[i,m] = dot_product(sub_col(infections, start, m, i - start), tail(si_rev, i - start));
        infections[i,m] = Rt_unadj[i,m] * load[i,m];
        Rt[i,m] = Rt_unadj[i,m];

        if (latent) {
            infections[i,m] = infections[i,m] + sqrt(infections[i,m] * inf_aux[1]) * inf_noise[i,m];
        
        if (pop_adjust) {
            infections[i,m] = (pop[m] - cumm_sum[i-1,m]) * (1 - exp(- infections[i,m] / pop[m]));
            Rt[i,m] =  (pop[m] - cumm_sum[i-1,m]) * Rt[i,m] / pop[m];
        }
        infectiousness[i,m] = load[i,m] / max(si);
        cumm_sum[i,m] = cumm_sum[i-1,m] + infections[i,m];
    }
}
}

