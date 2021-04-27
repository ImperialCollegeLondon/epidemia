int<lower=0, upper=1> latent; // no latent sampling iff 0
int<lower=0, upper=1> inf_family; // 0: none, 1: gamma
int<lower=0, upper=1> I0_fixed; // 0: prior on I0, 1: I0 = 0 fixed
