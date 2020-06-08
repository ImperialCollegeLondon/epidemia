vector reverse(vector vec) {
    int K = rows(vec);
    vector[K] rev;
    for (k in 1:K) {
        rev[k] = vec[K-k+1];
    }
    return rev;
}

