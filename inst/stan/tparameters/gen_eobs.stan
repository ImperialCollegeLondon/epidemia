// probabilities of recording an observation at that date
E_obs = inv_logit(oeta);

{  // compute expected values of the observations
for (i in 1:N_obs) {
    int m = obs_group[i];
    int dt = obs_date[i];
    int tp = obs_type[i];
    int n0 = starts[m];
    if (dt == 1)
    E_obs[i] *= 1e-15 * infections[1,m];
    else
    E_obs[i] *= dot_product(sub_col(infections, n0, m, dt-n0), tail(pvecs_rev[tp], dt-n0));
}
}
  