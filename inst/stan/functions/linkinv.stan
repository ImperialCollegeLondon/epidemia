 /** 
   * Apply inverse link function to linear predictor.
   * Adapted from rstanarm
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
*/
  vector linkinv(vector eta, int link) {
    if (link == 1)      return(inv_logit(eta)); // logit
    else if (link == 2) return(Phi(eta)); // probit
    else if (link == 3) return(atan(eta) / pi() + 0.5);  // cauchit
    else if (link == 4) return(inv_cloglog(eta)); // cloglog
    else if (link == 5) return(eta); //identity
    else reject("Invalid link");
    return eta; // never reached
  }
