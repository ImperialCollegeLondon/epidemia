int<lower=0> ac_nterms; // total number of calls to rw() in the formula
int<lower=0> ac_nproc; // total number of autocorrelation processes
int<lower=0> ac_q; // total number of time periods (sum of time periods for each process)
int<lower=0> ac_nnz; // total number of non-zero entries
int<lower=0> ac_ntime[ac_nproc]; // number of time periods for each process
int<lower=0, upper=ac_q-1> ac_v[ac_nnz]; // column indices from rstan::extract_sparse_matrix

