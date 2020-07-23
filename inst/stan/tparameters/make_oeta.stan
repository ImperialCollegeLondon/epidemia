{
    int npos = 1;
    int kpos = 1;
    int i = 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox1 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox2 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox3 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox4 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox5 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox6 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox7 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox8 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox9 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            segment(oeta, npos, oN[i]) = ox10 * segment(obeta, kpos, oK[i])
        else 
            segment(oeta, npos, oN[i]) = rep_vector(0.0, oN[i]);
        // add intercept
        if (has_ointercept[i] > 0)
            segment(oeta, npos, oN[i]) += ogamma[has_ointercept[i]];
        else
             segment(oeta, npos, oN[i]) += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]);
    }
}

