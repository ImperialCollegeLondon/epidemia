{ // predict cases over time
int idx=1;
int idx2=1;
for (m in 1:M){
    // time indices for group m: start date, final seed date, final date
    int n0 = starts[m];
    int n1 = n0 + N0 - 1;
    int n2 = n0 + NC[m] - 1;
    int len;

    if (link == 1) { // log-link
        Rt_unadj[n0:n2,m] = exp(eta[idx:(idx+NC[m]-1)]);
    } else if (link == 2) { // scaled_logit
        Rt_unadj[n0:n2,m] = carry * inv_logit(eta[idx:(idx+NC[m]-1)]);
    } else { // identity
        Rt_unadj[n0:n2,m] = eta[idx:(idx+NC[m]-1)];
    }
    
    idx += NC[m];
    infections[n0:n1,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
    cumm_sum[n0:n1,m] = cumulative_sum(infections[n0:n1,m]);

    for (i in (n1+1):n2) {
        int start = max(n0, i - gen_len);
        load[i,m] = dot_product(sub_col(infections, start, m, i - start), tail(gen_rev, i - start));
        infections[i,m] = Rt_unadj[i,m] * load[i,m];

        if (latent) { // treat as log-normal (could extend)
            real loginf = log(infections[i,m]);
            real loginfd = log(infections[i,m] + inf_aux[1]);
            infections[i,m] = exp( (3 * loginf - loginfd) / 2 + sqrt(loginfd - loginf) * inf_noise[idx2]);
        }
        
        if (pop_adjust) 
            infections[i,m] = (susc[start,m] - cumm_sum[i-1,m]) * (1 - exp(-(susc[i,m] / susc[start,m]) * (infections[i,m] / susc[start,m])));
        
        cumm_sum[i,m] = cumm_sum[i-1,m] + infections[i,m];
        idx2 += 1;
    }
}
}

