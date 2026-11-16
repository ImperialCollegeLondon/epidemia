real<lower=0> carry; // carry is only used if link == 2
int<lower=1, upper=3> link; // the link function 1) log, 2) scaled_logit, 3) identity
matrix<lower=0,upper=1> [N2,M] vacc; // proportion of susceptibles vaccinated at each point in time
array[M] real<lower=0> pops; // population of each group
int<lower=1> gen_len; // length of the longest generation distribution
int<lower=1> n_gen; // number of distinct generation distributions
array[n_gen] simplex[gen_len] gen; // generation distributions, padded with zeros to gen_len
array[M] int<lower=1, upper=n_gen> gen_idx; // generation distribution used by each group
int<lower=0, upper=1> pop_adjust; // whether to perform population adjustment or not
