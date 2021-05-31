{ // predict cases over time
    int idx1=1;
    int idx2=1;
    for (m in 1:M){
        // time indices for group m: start date, final seed date, final date
        int n0 = starts[m];
        int n1 = n0 + N0 - 1;
        int n2 = n0 + NC[m] - 1;
        int len;
        real vt;

        // compute Rt over the entire period
        if (link == 1) 
            Rt_unadj[n0:n2,m] = exp(eta[idx1:(idx1+NC[m]-1)]); // log-link
        else if (link == 2) 
            Rt_unadj[n0:n2,m] = carry * inv_logit(eta[idx1:(idx1+NC[m]-1)]); // scaled_logit
        else 
            Rt_unadj[n0:n2,m] = eta[idx1:(idx1+NC[m]-1)]; // identity
        
        idx1 += NC[m];
        
        infections[n0:n1,m] = rep_vector(y[m], N0); // seeded infections

        if (pop_adjust) { // initialise susceptible population
            susc[n0,m] = pops[m];
            if (!S0_fixed) susc[n0,m] *= S0[m];
        }

        for (i in n0:n2) {

            if (i > n1) { // infections after the seeding period
                int start = max(n0, i - gen_len);
                load[i,m] = dot_product(sub_col(infections, start, m, i - start), tail(gen_rev, i - start));
                E_infections[i,m] = Rt_unadj[i,m] * load[i,m];
                if (latent) infections[i,m] = infections_raw[idx2];
                else infections[i,m] = E_infections[i,m];
                idx2 += 1;
            }
            
            if (pop_adjust) {
                infections[i,m] = susc[i,m] * (1 - exp(-infections[i,m] / pops[m]));
                vt = vacc[i,m];
                if (!veps_fixed) vt *= veps[m];
                if (i != n2) susc[i+1,m] = (1 - vt) * (susc[i,m] - infections[i,m]);
            }
        }    
    }
}
