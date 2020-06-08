  if (K > 0) {
     eta = X * beta;
  }
  else eta = rep_vector(0.0, N);

  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* tau;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p,
                             1.0, tau, scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }

  if (t > 0) {
#include /tparameters/eta_add_Zb.stan
  }
  if (has_intercept == 1) {
    eta += gamma[1];
  }
  else {
#include /tparameters/eta_no_intercept.stan
  }
  

