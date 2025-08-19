setwd('/ihome/jbwang/chl471/bartpp')

library(optparse)
library(doParallel)
library(tidyverse)
library(matrixStats)
library(progress)
library(gtools)
library(moderate.mediation)
library(coda)

### input options 
option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "input RDS file contains data"),
  make_option(c("-m", "--method"), type = "character", help = "method use"),
  make_option(c("-t", "--start"), type = "integer", help = "i to start run"),
  make_option(c("-n", "--n"), type = "integer", help = "number of repeats run"),
  make_option(c("-s", "--size"), type = "integer", help = "sample size run in each repeat")
)
opt = parse_args(OptionParser(option_list=option_list))

input = opt$input
method = opt$method
start = opt$start
n = opt$n
size = opt$size

dat = readRDS(input)
output = paste0(method, "_case_", gsub("[^0-9]", "", input), "_", size, ".csv")

### Basic settings 
if(method %in% c("HMBART1", "HMBART2")){
  n_burn_in = 2000; n_after_burn_in = 750;
  if(size == 1000){
    options(java.parameters = "-Xmx48g")
  } else {
      options(java.parameters = "-Xmx24g")
  }
  library(bartMachine)
  set_bart_machine_num_cores(num_cores = 10)
  output = paste0(method, "_case_", gsub("[^0-9]", "", input), "_", size, ".csv")
} else if(method %in% c("VCBCMF1", "VCBCMF2", "VCBCMF3", "VCBCMF4")){
  library(SoftBart)
  source('lib/bart_mediate2.R')
  source('lib/get_clever_cov.R')
  source('lib/get_ps.R')
  n_burn_in = 5000; n_after_burn_in = 2500;
}

if(input %in% c("dat1.rds", "dat2.rds")){
  name = c("x1", "x2")
  formula_ps = t ~ x1 + x2
  formula_m = m ~ x1 + x2 + t
  formula_y = y ~ x1 + x2 + m + t
  
} else if(input %in% c("dat3.rds", "dat4.rds")){
  name = c("x1", "x2", "x3", "x4", "x5")
  formula_ps = t ~ x1 + x2 + x3 + x4 + x5
  formula_m = m ~ x1 + x2 + x3 + x4 + x5 + t
  formula_y = y ~ x1 + x2 + x3 + x4 + x5 + m + t
} else{
  dat = subset(dat,
               y >= -4.25306 &
               y <= 1.63242)
  name = c("x1", "x2", "x3", "x4")
  formula_ps = t ~ x1 + x2 + x3 + x4
  formula_m = m ~ x1 + x2 + x3 + x4 + t
  formula_y = y ~ x1 + x2 + x3 + x4 + m + t
}


for(st in start:n){
  
  data = dat[((st-1)*size+1):(st*size), ]
  
  
  if(method %in% c("VCBCMF1", "VCBCMF2", "VCBCMF3", "VCBCMF4")){
    
    ### Clever covariates and ps score
    clever_cov = get_clever_cov(data, data, formula_m, 'm', 'y', 't')
    m0_hat = clever_cov$m0_hat
    m1_hat = clever_cov$m1_hat
    pi_hat = get_ps(data, data, formula_ps)
    
    if(method == "VCBCMF1"){
      num_tree = 20; gamma = 0.5; k = 1;
    } else if(method == "VCBCMF2"){
      num_tree = 20; gamma = 0.5; k = 0.5;
    } 
    ### Fit
    out_bart = bart_mediate2(
      data_train = data, 
      data_test = data, 
      formula_m = formula_m, 
      formula_y = formula_y, 
      pi_hat_train = pi_hat, 
      pi_hat_test = pi_hat,
      m0_hat_train = m0_hat, 
      m0_hat_test = m0_hat, 
      m1_hat_train = m1_hat, 
      m1_hat_test = m1_hat,
      mediator_name = 'm', 
      outcome_name = 'y', 
      treat_name = 't',
      n_iter = n_burn_in + n_after_burn_in, 
      burnin = n_burn_in, 
      num_tree = num_tree, gamma = gamma, k = k,
      chain = 1
    )
    
    eNIE = out_bart$d_samples * out_bart$tau_samples
    eNDE = out_bart$zeta_samples
    
    outdf = data.frame(
      true_NDE = data$true_NDE, NDE_e = colMeans(eNDE), NDE_l = colQuantiles(eNDE, probs = c(0.025)), NDE_u = colQuantiles(eNDE, probs = c(0.975)), NDE_var = colVars(eNDE),
      true_NIE = data$true_NIE, NIE_e = colMeans(eNIE), NIE_l = colQuantiles(eNIE, probs = c(0.025)), NIE_u = colQuantiles(eNIE, probs = c(0.975)), NIE_var = colVars(eNIE))
    write.table(outdf, file = output, sep = ",", row.names = FALSE,
                col.names = !file.exists(output), append = file.exists(output)) 
  
  } else if(method == "MM"){
    
    if(input  == "dat5.rds"){
      RunMM = function(data, i){
        
        resulti = modmed(data = data, 
                         treatment = "t", 
                         mediator = "m", 
                         outcome = "y", 
                         moderators.cont = name, 
                         m.model = list(intercept = name, treatment = name), 
                         y.model = list(intercept = name, treatment = NULL, 
                                        mediator = NULL, tm = NULL),
                         m.scale = "continuous", 
                         y.scale = "continuous", 
                         comp.mod.cont.values = NULL, 
                         ref.mod.cont.values = data[i, name], 
                         method = "mc",
                         nmc = 1000,
                         nboot = 1000,
                         conf.level = 0.95)
        
        return(c(resulti$effects[c('TIE.ref', 'PDE.ref'), ]))
        
      }
      
    } else {
      RunMM = function(data, i){
        
        resulti = modmed(data = data, 
                         treatment = "t", 
                         mediator = "m", 
                         outcome = "y", 
                         moderators.cont = name, 
                         m.model = list(intercept = name, treatment = name), 
                         y.model = list(intercept = name, treatment = name, 
                                        mediator = name, tm = name),
                         m.scale = "continuous", 
                         y.scale = "continuous", 
                         comp.mod.cont.values = NULL, 
                         ref.mod.cont.values = data[i, name], 
                         method = "mc",
                         nmc = 1000,
                         nboot = 1000,
                         conf.level = 0.95)
        
        return(c(resulti$effects[c('TIE.ref', 'PDE.ref'), ]))
        
      }  
    }
    
    
    ### Run the simulation
    cl = makeCluster(20)
    registerDoParallel(cl)
    sim = foreach(i = 1 : nrow(data), .combine = rbind,
                  .packages = c("moderate.mediation")) %dopar% {RunMM(data, i)}
    stopCluster(cl)
    outdf = data.frame(
      true_NDE = data$true_NDE, NDE_e = sim[, 2], NDE_l = sim[, 6], NDE_u = sim[, 8], NDE_var = sim[, 4],
      true_NIE = data$true_NIE, NIE_e = sim[, 1], NIE_l = sim[, 5], NIE_u = sim[, 7], NIE_var = sim[, 3])
    write.table(outdf, file = output, sep = ",", row.names = FALSE,
                col.names = !file.exists(output), append = file.exists(output)) 
    
  } else {

    ### HC
    hcdata = data %>% select(-c(t, m, y, true_NDE, true_NIE))
    dist_matrix = dist(scale(hcdata))
    hc = hclust(dist_matrix)
    clusters = cutree(hc, k = as.integer(size / 5))
    if(input %in% c("dat5.rds", "dat6.rds")){
      clusters = cutree(hc, k = as.integer(size / 25))
    }
    
    k_cvs = c(2, 3, 5)
    nu_q_cvs = list(c(3, 0.9), c(3, 0.99), c(10, 0.75))
    
    ### Standardization
    m_mean = mean(data$m); m_sd = sd(data$m);
    y_mean = mean(data$y); y_sd = sd(data$y);
    data$m = (data$m - m_mean) / m_sd;
    data$y = (data$y - y_mean) / y_sd;
    
    ### Fit propensity score
    fit_ps = bartMachine(X = data[, name], y = as.factor(data$t), 
                         num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in)
    data$ps = fit_ps$y_hat_train
    
    ### M model
    if(method == "HMBART1"){
      fit_m = bartMachineCV(X = data[, c(name, 't', 'ps')], y = data$m,
                            num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in,
                            num_tree_cvs = c(30), k_cvs = k_cvs, nu_q_cvs = nu_q_cvs)
    } else if(method == "HMBART2"){
      fit_m = bartMachine(X = data[, c(name, 't', 'ps')], y = data$m, num_tree = 30,
                          num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in)
    }
  
    
    #fit_m = bartMachine(X = data[, c(name, 't', 'ps')], y = data$m, 
    #                    num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in)
    residuals = fit_m$residuals
    sigsqs = sqrt(get_sigsqs(fit_m))
    
    ### Add M estimates
    data_m = data[, c(name, 't', 'ps')]
    data_m$t = 0
    m0 = predict(fit_m, data_m)
    data$m0 = m0
    data_m$t = 1
    m1 = predict(fit_m, data_m)
    data$m1 = m1
    
    ### Y model
    #fit_y = bartMachine(X = data[, c(name, 't', 'ps', 'm', 'm0', 'm1')], y = data$y, 
    #                      num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in)
    if(method == "HMBART1"){
      fit_y = bartMachineCV(X = data[, c(name, 't', 'ps', 'm', 'm0', 'm1')], y = data$y,
                            num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in,
                            num_tree_cvs = c(30), k_cvs = k_cvs, nu_q_cvs = nu_q_cvs)
    } else if(method == "HMBART2"){
      fit_y = bartMachine(X = data[, c(name, 't', 'ps', 'm', 'm0', 'm1')], y = data$y, num_tree = 30,
                            num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in)
    }
    
    ### Construct dataset
    data_m = data[, c(name, 't', 'ps')]
    data_m$t = 0
    m0 = bart_machine_get_posterior(fit_m, data_m)$y_hat_posterior_samples
    data_m$t = 1
    m1 = bart_machine_get_posterior(fit_m, data_m)$y_hat_posterior_samples
    
    ### HMBART1
    data_y_m0 = NULL; data_y_m1 = NULL;
    data_y = data[, c(name, 't', 'ps', 'm', 'm0', 'm1')]
    for(i in 1:nrow(data)){
      error = rnorm(n = n_after_burn_in, mean = 0, sd = sigsqs)
      data_i = data_y[replicate(n_after_burn_in, i), ]
      data_i$m = m0[i, ] + error
      data_y_m0 = rbind(data_y_m0, data_i)
      data_i$m = m1[i, ] + error
      data_y_m1 = rbind(data_y_m1, data_i)
    }
    
    ### Predict
    library(parallel)
    breaks = c(seq(from = 1, to = nrow(data_y_m0), by = 5 * n_after_burn_in * n_after_burn_in), nrow(data_y_m0)+1)
    get_post_samples_parallel <- function(data_set, t_value) {
      data_set$t <- t_value
      out <- mclapply(seq_along(breaks)[-length(breaks)], function(i) {
        bart_machine_get_posterior(fit_y, data_set[breaks[i]:(breaks[i+1]-1), ])$y_hat_posterior_samples
      }, mc.cores = 20)
      do.call(rbind, out)
    }
    
    post_y0m0 <- get_post_samples_parallel(data_y_m0, t_value = 0)
    post_y1m0 <- get_post_samples_parallel(data_y_m0, t_value = 1)
    post_y1m1 <- get_post_samples_parallel(data_y_m1, t_value = 1)
    
    ### Summarize
    NDE_e = c(); NDE_l = c(); NDE_u = c(); NDE_var = c();
    NIE_e = c(); NIE_l = c(); NIE_u = c(); NIE_var = c();
    for(i in 1:nrow(data)){
      start = (i-1) * n_after_burn_in + 1; end = i * n_after_burn_in;
      post_y0m0_i = c(post_y0m0[start:end, ])
      post_y1m0_i = c(post_y1m0[start:end, ])
      post_y1m1_i = c(post_y1m1[start:end, ])
      eNIE = (post_y1m1_i - post_y1m0_i) * y_sd
      eNDE = (post_y1m0_i - post_y0m0_i) * y_sd
      
      NDE_e = c(NDE_e, mean(eNDE)); NDE_l = c(NDE_l, quantile(eNDE, 0.025)); NDE_u = c(NDE_u, quantile(eNDE, 0.975)); NDE_var = c(NDE_var, var(eNDE));
      NIE_e = c(NIE_e, mean(eNIE)); NIE_l = c(NIE_l, quantile(eNIE, 0.025)); NIE_u = c(NIE_u, quantile(eNIE, 0.975)); NIE_var = c(NIE_var, var(eNIE));
    }
    
    outdf1 = data.frame(
      true_NDE = data$true_NDE, NDE_e = NDE_e, NDE_l = NDE_l, NDE_u = NDE_u, NDE_var = NDE_var,
      true_NIE = data$true_NIE, NIE_e = NIE_e, NIE_l = NIE_l, NIE_u = NIE_u, NIE_var = NIE_var)
    
    
    result_df = data.frame(
      NDE_e = outdf1$NDE_e, NDE_l = outdf1$NDE_l, NDE_u = outdf1$NDE_u, NDE_var = outdf1$NDE_var,
      NIE_e = outdf1$NIE_e, NIE_l = outdf1$NIE_l, NIE_u = outdf1$NIE_u, NIE_var = outdf1$NIE_var,
      cluster = clusters
    )
    
    # Prepare new columns for shrunk estimates
    result_df$NDE_e_shrunk <- NA
    result_df$NDE_l_shrunk <- NA
    result_df$NDE_u_shrunk <- NA
    
    result_df$NIE_e_shrunk <- NA
    result_df$NIE_l_shrunk <- NA
    result_df$NIE_u_shrunk <- NA
    
    for(cl in unique(result_df$cluster)) {
      idx <- which(result_df$cluster == cl)
      
      # If only one in the cluster, do not shrink
      if(length(idx) == 1) {
        result_df$NDE_e_shrunk[idx] <- result_df$NDE_e[idx]
        result_df$NDE_l_shrunk[idx] <- result_df$NDE_l[idx]
        result_df$NDE_u_shrunk[idx] <- result_df$NDE_u[idx]
        result_df$NIE_e_shrunk[idx] <- result_df$NIE_e[idx]
        result_df$NIE_l_shrunk[idx] <- result_df$NIE_l[idx]
        result_df$NIE_u_shrunk[idx] <- result_df$NIE_u[idx]
        next
      }
      
      # Normal shrinkage for clusters with >1
      nde_hat <- result_df$NDE_e[idx]
      nde_var <- result_df$NDE_var[idx]
      nde_lower <- result_df$NDE_l[idx]
      nde_upper <- result_df$NDE_u[idx]
      mu_hat_nde <- mean(nde_hat)
      tau2_nde <- var(nde_hat)
      w_nde <- tau2_nde / (tau2_nde + nde_var)
      
      nie_hat <- result_df$NIE_e[idx]
      nie_var <- result_df$NIE_var[idx]
      nie_lower <- result_df$NIE_l[idx]
      nie_upper <- result_df$NIE_u[idx]
      mu_hat_nie <- mean(nie_hat)
      tau2_nie <- var(nie_hat)
      w_nie <- tau2_nie / (tau2_nie + nie_var)
      
      # Shrunk NDE
      result_df$NDE_e_shrunk[idx] <- w_nde * nde_hat + (1 - w_nde) * mu_hat_nde
      result_df$NDE_l_shrunk[idx] <- result_df$NDE_e_shrunk[idx] + (nde_lower - nde_hat)
      result_df$NDE_u_shrunk[idx] <- result_df$NDE_e_shrunk[idx] + (nde_upper - nde_hat)
      
      # Shrunk NIE
      result_df$NIE_e_shrunk[idx] <- w_nie * nie_hat + (1 - w_nie) * mu_hat_nie
      result_df$NIE_l_shrunk[idx] <- result_df$NIE_e_shrunk[idx] + (nie_lower - nie_hat)
      result_df$NIE_u_shrunk[idx] <- result_df$NIE_e_shrunk[idx] + (nie_upper - nie_hat)
    }
    
    outdf2 = data.frame(
      true_NDE = data$true_NDE, NDE_e = result_df$NDE_e_shrunk, NDE_l = result_df$NDE_l_shrunk, NDE_u = result_df$NDE_u_shrunk, NDE_var = outdf1$NDE_var,
      true_NIE = data$true_NIE, NIE_e = result_df$NIE_e_shrunk, NIE_l = result_df$NIE_l_shrunk, NIE_u = result_df$NIE_u_shrunk, NIE_var = outdf1$NIE_var)
    
    write.table(outdf2, file = output, sep = ",", row.names = FALSE,
                col.names = !file.exists(output), append = file.exists(output)) 
  }

}