vector[ac_nproc] ac_scale;
vector[ac_q] ac_beta;

ac_scale = ac_scale_raw * ac_prior_scale;

i = 1;
for (l in 1:ac_nproc) { // this treats ac terms as random walks for now (to be extended to AR(p))
    slice = i:(i+ac_ntime[l]-1);
    ac_beta[slice] = cumulative_sum(ac_noise[slice]);
    i += ac_ntime[l];
}