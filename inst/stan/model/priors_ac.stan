target += normal_lpdf(ac_scale_raw | 0, 1);

i = 1;
for (l in 1:ac_nproc) { 
    target += normal_lpdf(ac_noise[i:(i+ac_nperiods[l]-1)] | 0, ac_scale[l]);
    i += ac_nperiods[l];
}

