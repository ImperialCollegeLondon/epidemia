{
    int npos = 1;
    int kpos = 1;
    int i = 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX1 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX2 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX3 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX4 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX5 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX6 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX7 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX8 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX9 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX10 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
    }
}

