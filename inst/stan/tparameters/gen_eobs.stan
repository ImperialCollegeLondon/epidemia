{ // apply link function
    int i = 1;
    for (r in 1:R) {
        vector[oN[r]] inv_vec = segment(oeta, i, oN[r]) + 1e-15;
        E_obs[i:(i+oN[r]-1)] = linkinv(inv_vec, olink[r]);
        i += oN[r];
    }
}

{  // compute expected values of the observations
for (i in 1:N_obs) {
    int m = obs_group[i];
    int dt = obs_date[i];
    int tp = obs_type[i];
    int n0 = starts[m];
    if (dt == 1)
        E_obs[i] *= 1e-15 * infections[1,m];
    else {
        int start = max(n0, dt - pvecs_len[tp]);
        E_obs[i] *= dot_product(sub_col(infections, start, m, dt-start), tail(pvecs_rev[tp], dt-start));
    }
}
}
