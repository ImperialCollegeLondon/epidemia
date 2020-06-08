
int<lower=1> M; // number of countries
int<lower=1> N0; // number of time points for which to impute infections
int<lower=1> starts[M]; // the start index of each group
int<lower=1> NC[M]; // days of observed data for each group.
int<lower=1> N2; // total period for the simulation
int<lower=1> NS; // maximum number of simulation days for any given group
int<lower=0> R; // number of different observation types


