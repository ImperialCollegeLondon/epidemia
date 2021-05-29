real normal_lb_rng(real mu, real sigma, real lb) {
  real p_lb = normal_cdf(lb, mu, sigma);
  real u = uniform_rng(p_lb, 1.0);
  real y = mu + sigma * inv_Phi(u);
  return y;
}