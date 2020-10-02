for (i in 1:ac_nterms) {
    for (j in 1:N) {
        if (ac_V[i,j] >= 0) { # if -1 then doesn't have RW component
            eta[j] += ac_beta[ac_V[i,j]];
        }
    }
}