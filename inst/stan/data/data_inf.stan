int<lower=0, upper=1> hseeds; // are seeds hierarchical? 1 iff yes
int<lower=0, upper=1> latent; // no latent sampling iff 0
int<lower=0, upper=1> inf_family; // 0: none, 1: gamma
int<lower=0, upper=1> S0_fixed; // 0 : prior on S0, 1: S0 = 1
int<lower=0, upper=1> veps_fixed; // 0: prior on veps, 1: veps = 1
int<lower=0, upper=1> fixed_vtm; // 0: sd linear in mean, 1: var linear in mean
