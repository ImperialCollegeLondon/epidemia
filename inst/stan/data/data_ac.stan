int<lower=0> ac_nterms; // total number of calls to rw() in the formula
int<lower=0> ac_nproc; // total number of autocorrelation processes
int<lower=0> ac_q; // total number of time periods (sum of time periods for each process)
int<lower=0> ac_nnz; // total number of non-zero entries
int<lower=0> ac_ntime[ac_nproc]; // number of time periods for each process
int<lower=-1, upper=ac_q-1> ac_v[ac_nnz]; // column indices from rstan::extract_sparse_matrix, -1 corresponds to no ac_term


int<lower=0> obs_ac_nterms; // total number of calls to rw() over all observation models
int<lower=0> obs_ac_nproc; // total number of autocorrelation processes over all observational models
int<lower=0> obs_ac_q; // total number of time periods (sum of time periods for each process)
int<lower=0> obs_ac_nnz; // total number of non-zero entries
int<lower=0> obs_ac_ntime[obs_ac_nproc]; // number of time periods for each process
int<lower=0> obs_ac_V[obs_ac_nterms, N_obs]; // matrix giving group memberships for each observation
