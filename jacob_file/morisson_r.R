pacman::p_load(SeqExpMatch, dplyr)
pacman::p_load(data.table)
pacman::p_load(ggplot2, tidyr, stringr)

rm(list = ls())

#self_data = list()
#private_data = list()

n = 1000
p = 4
mu_x = 1
sigma_x = 1
sigma_e = 1
beta_T = 1
nsim_exact_test = 51
num_cores = 20
Nsim = 1

morisson_results = data.frame(
  index                                        = 1:n,
  res_size_KK14                                = rep(NA, n),
  res_size_KK14_morisson_t                     = rep(NA, n),
  res_size_KK14_morisson_T                     = rep(NA, n),
  res_size_KK21                                = rep(NA, n),
  res_size_KK21_morisson_t                     = rep(NA, n),
  res_size_KK21_morisson_T                     = rep(NA, n),
  res_size_KK21stepwise                        = rep(NA, n),
  res_size_KK21stepwise_morisson_t             = rep(NA, n),
  res_size_KK21stepwise_morisson_T             = rep(NA, n),
  match_pair_distance_KK14                     = rep(NA, n),
  match_pair_distance_KK14_morisson_t          = rep(NA, n),
  match_pair_distance_KK14_morisson_T          = rep(NA, n),
  match_pair_distance_KK21                     = rep(NA, n),
  match_pair_distance_KK21_morisson_t          = rep(NA, n),
  match_pair_distance_KK21_morisson_T          = rep(NA, n),
  match_pair_distance_KK21stepwise             = rep(NA, n),
  match_pair_distance_KK21stepwise_morisson_t  = rep(NA, n),
  match_pair_distance_KK21stepwise_morisson_T  = rep(NA, n)
  
)

#res_size = rep(NA, n)
#match_pair_distance = rep(NA, n)

#build mvnp covariates
set.seed(1984)
errors = rnorm(n, 0, sigma_e)

all_betas_and_correlations = list()
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(1, 1, 1, 0, 0, 2, 2, 1)))) #QUAD EVEN
#all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(1, 1, 1, 0, 0)))) #QUAD MORE UNEVEN with CORR
#all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN
#all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN with CORR

# res = data.frame(
#   betas = character(),
#   rho = numeric(),
#   design = character(),
#   estimate_type = character(),
#   test_type = character(),
#   beta_hat_T = numeric(),
#   pval = numeric(),
#   lambda = numeric()
# )
for (nsim in 1 : Nsim){
  cat ("nsim:", nsim, "/", Nsim, "\n")
  for (all_betas_and_correlation in all_betas_and_correlations){
    betas = all_betas_and_correlation[["betas"]]
    rho = all_betas_and_correlation[["rho"]] 
    
    Sigma = sigma_x * (matrix(rho, nrow = p, ncol = p) + diag(1 - rho, p))
    X = data.table(MASS::mvrnorm(n, rep(mu_x, p), Sigma))
    z = betas[1] * X[, 1] +
      betas[2] * X[, 2] + 
      betas[3] * X[, 1]^2 +
      betas[4] * X[, 2]^2 +
      betas[5] * X[, 1] * X[, 2] +
      betas[6] * X[, 3] +
      betas[7] * X[, 4] +
      betas[8] * X[, 3] * X[, 1]
    y = array(NA, n)
    
    #test all designs
    for (d in c("KK14", "KK21", "KK21stepwise")){ #"CRD", "BCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"
      for(m in c("", "_morisson_t", "_morisson_T")){
        if(m == ""){
          seq_des_obj = SeqDesign$new(n = n, prob_T = 0.5, p = p, design = d, response_type = "continuous", verbose = TRUE, t_0_pct = 1e-03)
        } else if(m == "_morisson_t"){
          seq_des_obj = SeqDesign$new(n = n, prob_T = 0.5, p = p, design = d, response_type = "continuous", verbose = TRUE, morrison_little_t = TRUE, t_0_pct = 1e-03)
        } else {
          seq_des_obj = SeqDesign$new(n = n, prob_T = 0.5, p = p, design = d, response_type = "continuous", verbose = TRUE, morrison_big_T = TRUE, t_0_pct = 1e-03)
        }
        
        for (t in 1 : n){
          
          seq_des_obj$add_subject_to_experiment_and_assign(X[t, ])
          if(t %% 100 == 0){
            print(t)
            print(seq_des_obj$jacobs_fun(action = 1))
          }
          morisson_results[[paste0("res_size_", d, m)]][t] = seq_des_obj$jacobs_fun(action = 5)
          #morisson_results[[paste0("match_pair_distance_", d, m)]][t] = seq_des_obj$jacobs_fun(action = 6)
          
          #res_size[t] = seq_des_obj$jacobs_fun(action = 5)
          #match_pair_distance[t] =  seq_des_obj$jacobs_fun(action = 6)
          
          w_t = seq_des_obj$w[seq_des_obj$t]
          y[t] = beta_T * w_t + z[t] + errors[t]
          seq_des_obj$add_subject_response(t = t, y = as.numeric(y[t]))
        }
        Sys.time()
      }
      # for (test_type in c("MLE-or-KM-based")){ #, "randomization-exact"
      #   for (estimate_type in c("simple_mean_difference", "continuous_regression_with_covariates")){
      #     seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = num_cores, verbose = FALSE)
      #     
      #     beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
      #     pval = seq_des_inf_obj$compute_two_sided_pval_for_treatment_effect(nsim_exact_test = nsim_exact_test)
      #     
      #     res = rbind(res, data.frame(
      #       betas = paste0(betas, collapse=""),
      #       rho = rho,
      #       design = d,
      #       estimate_type = estimate_type,
      #       test_type = test_type,
      #       beta_hat_T = beta_hat_T,
      #       pval = pval
      #     ))
      #   }
      # }
      
    }
  }
  
}

write.csv(morisson_results, "C:\\temp\\morisson_results_comp1.csv", row.names = FALSE)


morisson_results = read.csv("C:\\temp\\morisson_results_comp1.csv")




long_df = morisson_results %>%
  pivot_longer(cols = -index, names_to = "variable", values_to = "value") %>%
  
  # Step 2: Extract base name and morisson type using regex
  mutate(
    base = str_replace(variable, "_morisson_[tT]$", ""),
    morisson_type = case_when(
      str_ends(variable, "_morisson_t") ~ "morisson_t",
      str_ends(variable, "_morisson_T") ~ "morisson_T",
      TRUE ~ "original"
    )
  )

# Step 3: Plot one graph per `base` variable
unique_bases = unique(long_df$base)

for (b in unique_bases) {
  p = long_df %>%
    filter(base == b) %>%
    ggplot(aes(x = index, y = value, color = morisson_type)) +
    xlim(0,10000) +
    scale_y_log10() +
    #geom_line() +
    geom_smooth(se = FALSE, method = "loess", span = 0.2) +
    labs(title = b, x = "Index", y = "Value", color = "Version") +
    theme_minimal()
  
  print(p)
}

  

# res_mod = res %>% 
#   mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
#   group_by(betas, rho, design, estimate_type, test_type) %>%
#   summarize(mse = mean(sq_err), power = sum(rej) / n())
# 
# res_mod %>% filter(betas == "11100" & rho == 0 & estimate_type == "simple_mean_difference" & test_type == "MLE-or-KM-based") %>% as.data.frame
# res_mod %>% filter(betas == "11100" & rho == 0 & estimate_type == "continuous_regression_with_covariates" & test_type == "MLE-or-KM-based") %>% as.data.frame
# res_mod %>% filter(betas == "11100" & rho == 0 & estimate_type == "simple_mean_difference" & test_type == "randomization-exact") %>% as.data.frame
# res_mod %>% filter(betas == "11100" & rho == 0 & estimate_type == "continuous_regression_with_covariates" & test_type == "randomization-exact") %>% as.data.frame
# 
# res_mod %>% filter(betas == "11100" & rho == 0.75 & estimate_type == "simple_mean_difference" & test_type == "MLE-or-KM-based") %>% as.data.frame
# res_mod %>% filter(betas == "11100" & rho == 0.75 & estimate_type == "continuous_regression_with_covariates" & test_type == "MLE-or-KM-based") %>% as.data.frame
# res_mod %>% filter(betas == "11100" & rho == 0.75 & estimate_type == "simple_mean_difference" & test_type == "randomization-exact") %>% as.data.frame
# res_mod %>% filter(betas == "11100" & rho == 0.75 & estimate_type == "continuous_regression_with_covariates" & test_type == "randomization-exact") %>% as.data.frame
# 
# res_mod %>% filter(betas == "61200" & rho == 0 & estimate_type == "simple_mean_difference" & test_type == "MLE-or-KM-based") %>% as.data.frame
# res_mod %>% filter(betas == "61200" & rho == 0 & estimate_type == "continuous_regression_with_covariates" & test_type == "MLE-or-KM-based") %>% as.data.frame
# res_mod %>% filter(betas == "61200" & rho == 0 & estimate_type == "simple_mean_difference" & test_type == "randomization-exact") %>% as.data.frame
# res_mod %>% filter(betas == "61200" & rho == 0 & estimate_type == "continuous_regression_with_covariates" & test_type == "randomization-exact") %>% as.data.frame
# 
# res_mod %>% filter(betas == "61200" & rho == 0.75 & estimate_type == "simple_mean_difference" & test_type == "MLE-or-KM-based") %>% as.data.frame
# res_mod %>% filter(betas == "61200" & rho == 0.75 & estimate_type == "continuous_regression_with_covariates" & test_type == "MLE-or-KM-based") %>% as.data.frame
# res_mod %>% filter(betas == "61200" & rho == 0.75 & estimate_type == "simple_mean_difference" & test_type == "randomization-exact") %>% as.data.frame
# res_mod %>% filter(betas == "61200" & rho == 0.75 & estimate_type == "continuous_regression_with_covariates" & test_type == "randomization-exact") %>% as.data.frame

