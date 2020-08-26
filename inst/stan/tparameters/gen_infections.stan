{ // predict cases over time
    int idx=1;
    matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
    for (m in 1:M){
        // time indices for group m: start date, final seed date, final date
        int n0 = starts[m];
        int n1 = n0 + N0 - 1;
        int n2 = n0 + NC[m] - 1;
        int len;

        // impute unadjusted Rt from the linear predictor
        Rt_unadj[n0:n2,m] = r0 * 2 * inv_logit(eta[idx:(idx+NC[m]-1)]);
        Rt[n0:n1,m] = Rt_unadj[n0:n1,m]; 
        idx += NC[m];

        infections[n0:n1,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
        cumm_sum[n0:n1,m] = cumulative_sum(infections[n0:n1,m]);

        for(i in (n1+1):n2) {
            int start = max(n0, i - si_len);
            real convolution = dot_product(sub_col(infections, start, m, i - start), tail(si_rev, i - start));

            real rc = Rt_unadj[i,m] * convolution;
            real sig = sqrt(log(1/rc + 1));
            infections[i,m] = rc * exp(- sig * sig / 2 + sig * infection_noise[i,m]);
            Rt[i,m] = Rt_unadj[i,m];

            if (pop_adjust) {
                infections[i,m] = (pop[m] - cumm_sum[i-1,m]) * (1 - exp(-infections[i,m] / pop[m]));
                Rt[i,m] = (pop[m] - cumm_sum[i-1,m]) * Rt[i,m] / pop[m];
            }
            
            infectiousness[i,m] = convolution / max(si);
            cumm_sum[i,m] = cumm_sum[i-1,m] + infections[i,m];
        }
    }
}

