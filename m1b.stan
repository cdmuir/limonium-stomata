data {

  int<lower = 1> n_kcatc;
  int<lower = 1> n_sco;
  int<lower = 1> n_stomata;
  int<lower = 1> n_spp;
  
  real b;
  int<lower = 0, upper = 1> drought[n_stomata];
  real m;
  vector[n_stomata] sd1;
  vector[n_stomata] sp1;
  int<lower = 1> spp_stomata[n_stomata];
  
  vector[n_kcatc] kcatc;
  vector[n_sco] sco;
  int<lower = 1> spp_kcatc[n_kcatc];
  int<lower = 1> spp_sco[n_sco];
  
}
transformed data {
  
  vector[n_stomata] gsmax;
  vector[2] stomata[n_stomata];
  vector[n_kcatc] X_kcatc = kcatc - mean(kcatc);
  vector[n_sco] X_sco = sco - mean(sco);
  
  for (n in 1:n_stomata) {
    gsmax[n] = b * m * sd1[n] * sqrt(0.5 * (2 * sp1[n]) ^ 2);
    stomata[n] = [sd1[n], sp1[n]]';
  }
  
}
parameters {

  real b0_sd;  // temporary intercept 
  real<lower=0> sigma_sd;  // residual SD 
  real b0_sp;  // temporary intercept 
  real<lower=0> sigma_sp;  // residual SD 
  // parameters for multivariate linear models 
  cholesky_factor_corr[2] stomata_rescor; 
  vector<lower=0>[2] sd_1;  // group-level standard deviations
  vector<lower=0>[2] sd_2;  // group-level standard deviations
  matrix[2, n_spp] z_1;  // unscaled group-level effects
  matrix[2, n_spp] z_2;  // unscaled group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[2] L_1;
  cholesky_factor_corr[2] L_2;
  
  real b_sd_drought;
  real b_sp_drought;

  real b0_kcatc;
  vector[n_spp] b_kcatc_spp;
  real<lower = 0> sigma_kcatc;

  real b0_sco;
  vector[n_spp] b_sco_spp;
  real<lower = 0> sigma_sco;

  corr_matrix[2] Omega_rubisco;
  vector<lower = 0>[2] sigma_rubisco;

}
transformed parameters { 
  
  cov_matrix[2] Sigma_rubisco;

  // group-level effects 
  matrix[n_spp, 2] r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)';
  matrix[n_spp, 2] r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)';
  vector[n_spp] r_1_sd_1 = r_1[, 1];
  vector[n_spp] r_1_sp_2 = r_1[, 2];
  vector[n_spp] r_2_sd_1 = r_2[, 1];
  vector[n_spp] r_2_sp_2 = r_2[, 2];
  
  Sigma_rubisco = quad_form_diag(Omega_rubisco, sigma_rubisco);

} 
model { 
  
  vector[n_stomata] mu_sd = b0_sd + rep_vector(0, n_stomata);
  vector[n_stomata] mu_sp = b0_sp + rep_vector(0, n_stomata);
  
  vector[n_kcatc] mu_kcatc; 
  vector[n_sco] mu_sco; 

  vector[2] mu_rubisco_spp[n_spp];
  vector[2] rubisco[n_spp];


  // multivariate linear predictor matrix 
  vector[2] Mu[n_stomata]; 
  vector[2] sigma = [sigma_sd, sigma_sp]';
  
  // cholesky factor of residual covariance matrix 
  matrix[2, 2] LSigma = diag_pre_multiply(sigma, stomata_rescor); 

  for (n in 1:n_kcatc) { 

    mu_kcatc[n] = b0_kcatc + b_kcatc_spp[spp_kcatc[n]];

  } 

  for (n in 1:n_sco) { 

    mu_sco[n] = b0_sco + b_sco_spp[spp_sco[n]];

  } 
  
  for (n in 1:n_stomata) { 

    mu_sd[n] += r_1_sd_1[spp_stomata[n]] + 
      drought[n] * r_2_sd_1[spp_stomata[n]];
    mu_sp[n] += r_1_sp_2[spp_stomata[n]] +
      drought[n] * r_2_sp_2[spp_stomata[n]];
    Mu[n] = [mu_sd[n] + b_sd_drought * drought[n], 
             mu_sp[n] + b_sp_drought * drought[n]]';
    
  } 
  
  for (s in 1:n_spp) {
    
    mu_rubisco_spp[s] = [0, 0]';
    rubisco[s] = [b_kcatc_spp[s], b_sco_spp[s]]';
    
  }
  
  // priors
  target += student_t_lpdf(b0_sd | 3, 108, 17); 
  target += student_t_lpdf(sigma_sd | 3, 0, 17)
    - 1 * student_t_lccdf(0 | 3, 0, 17); 
  target += student_t_lpdf(b0_sp | 3, 27, 10); 
  target += student_t_lpdf(sigma_sp | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += lkj_corr_cholesky_lpdf(stomata_rescor | 1); 
  target += student_t_lpdf(sd_1 | 3, 0, 17)
    - 2 * student_t_lccdf(0 | 3, 0, 17);
  target += student_t_lpdf(sd_2 | 3, 0, 10)
    - 2 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(to_vector(z_1) | 0, 1);
  target += normal_lpdf(to_vector(z_2) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_1 | 1); 
  target += lkj_corr_cholesky_lpdf(L_2 | 1); 
  
  target += normal_lpdf(b_sd_drought | 0, 10);
  target += normal_lpdf(b_sp_drought | 0, 10);

  target += normal_lpdf(b0_kcatc | 0, 10);
  target += cauchy_lpdf(sigma_kcatc | 0, 1);

  target += normal_lpdf(b0_sco | 0, 10);
  target += cauchy_lpdf(sigma_sco | 0, 1);

  target += cauchy_lpdf(sigma_rubisco | 0, 1);
  target += lkj_corr_lpdf(Omega_rubisco | 1);
  
  target += multi_normal_lpdf(rubisco | mu_rubisco_spp, Sigma_rubisco);

  // likelihood including all constants 
  target += multi_normal_cholesky_lpdf(stomata | Mu, LSigma);
  target += normal_lpdf(X_kcatc | mu_kcatc, sigma_kcatc);
  target += normal_lpdf(X_sco | mu_sco, sigma_sco);
  
} 
generated quantities {
  
  vector[n_spp] gsmax_wd_spp;
  vector[n_spp] gsmax_ww_spp;
  vector[n_spp] kcatc_spp;
  vector[n_spp] sco_spp;

  matrix[2, 2] Rescor = multiply_lower_tri_self_transpose(stomata_rescor); 
  vector<lower = -1, upper = 1>[1] rescor; 
  corr_matrix[2] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower = -1, upper = 1>[1] cor_1;
  corr_matrix[2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower = -1, upper = 1>[1] cor_2;
  
  for (s in 1:n_spp) {
    
    gsmax_wd_spp[s] = b * m * (b0_sd + r_1_sd_1[s] + r_2_sd_1[s]) * 
      sqrt(0.5 * (2 * (b0_sp + r_1_sp_2[s] + r_2_sp_2[s])) ^ 2);
    gsmax_ww_spp[s] = b * m * (b0_sd + r_1_sd_1[s]) * 
      sqrt(0.5 * (2 * (b0_sp + r_1_sp_2[s])) ^ 2);
    kcatc_spp[s] = b0_kcatc + b_kcatc_spp[s] + mean(kcatc);
    sco_spp[s] = b0_sco + b_sco_spp[s] + mean(sco);
   
  }
  
  // take only relevant parts of residual correlation matrix 
  rescor[1] = Rescor[1, 2]; 
  // take only relevant parts of correlation matrix
  cor_1[1] = Cor_1[1,2]; 
  cor_2[1] = Cor_2[1,2]; 
  
} 
