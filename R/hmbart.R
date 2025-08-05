#' @importFrom bartMachine bartMachine bartMachineCV bart_machine_get_posterior get_sigsqs
#' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster
#' @export
hmbart = function(data, X, t, m, y, CV = FALSE,
                  num_trees = 50, num_tree_cvs = c(50), k_cvs = c(2, 3, 5),
                  nu_q_cvs = list(c(3, 0.9), c(3, 0.99), c(10, 0.75)),
                  emb_shrink = TRUE, cluster_size = 10,
                  n_burn_in = 2000, n_after_burn_in = 500, n_process_samples = 1e5, seed = 42, fix = FALSE, 
                  argT = list(), argM = list(), argY = list()) {
  
  set.seed(seed)
  model_seed = NULL
  if(fix){
    set_bart_machine_num_cores(num_cores = 1)
    model_seed = seed
  }
  ### Extract variables
  X_name = X
  data$t = data[, t]
  data$m = data[, m]
  data$y = data[, y]
  n = nrow(data)

  ### Standardization
  m_mean = mean(data$m); m_sd = sd(data$m);
  y_mean = mean(data$y); y_sd = sd(data$y);
  data$m = (data$m - m_mean) / m_sd;
  data$y = (data$y - y_mean) / y_sd;

  ### Fit propensity score
  args_ps = c(list(X = data[, X_name], y = data$t, num_trees = num_trees, 
                   num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = model_seed), argT)
  fit_ps = do.call(bartMachine, args_ps)
  data$ps = fit_ps$y_hat_train
  rm(fit_ps)

  ### M model
  if(CV){
    args_m = c(list(X = data[, c(X_name, 't', 'ps')], y = data$m,
                    num_tree_cvs = num_tree_cvs, k_cvs = k_cvs, nu_q_cvs = nu_q_cvs,
                    num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = model_seed), argM)
    fit_m = do.call(bartMachineCV, args_m)
  }else{
    args_m = c(list(X = data[, c(X_name, 't', 'ps')], y = data$m, num_trees = num_trees,
                    num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = model_seed), argM)
    fit_m = do.call(bartMachine, args_m)
  }
  sigsqs = sqrt(get_sigsqs(fit_m))

  ### Add M estimates for outcome model
  data_m = data[, c(X_name, 't', 'ps')]
  data_m$t = 0
  m0 = predict(fit_m, data_m)
  data$m0 = m0
  data_m$t = 1
  m1 = predict(fit_m, data_m)
  data$m1 = m1

  ### Construct dataset
  data_m = data[, c(X_name, 't', 'ps')]
  data_m$t = 0
  m0 = bart_machine_get_posterior(fit_m, data_m)$y_hat_posterior_samples
  data_m$t = 1
  m1 = bart_machine_get_posterior(fit_m, data_m)$y_hat_posterior_samples

  cat('\n', 'Construct Dataset', '\n', sep = '')
  pb = txtProgressBar(min = 1, max = n - 1, style = 3)
  data_y_m0_list = vector('list', n)
  data_y_m1_list = vector('list', n)
  data_y_m0 = NULL; data_y_m1 = NULL;
  data_y = data[, c(X_name, 't', 'ps', 'm', 'm0', 'm1')]
  for(i in 1:n){
    error = rnorm(n = n_after_burn_in, mean = 0, sd = sigsqs)
    data_i = data_y[replicate(n_after_burn_in, i), ]
    data_i$m = m0[i, ] + error
    data_y_m0_list[[i]] = data_i
    data_i$m = m1[i, ] + error
    data_y_m1_list[[i]] = data_i
    setTxtProgressBar(pb, i)
  }
  close(pb)
  data_y_m0 = do.call(rbind, data_y_m0_list)
  data_y_m1 = do.call(rbind, data_y_m1_list)
  rm(fit_m); rm(data_y_m0_list); rm(data_y_m1_list);

  ### Y model
  if(CV){
    args_Y = c(list(X = data[, c(X_name, 't', 'ps', 'm', 'm0', 'm1')], y = data$y,
                    num_tree_cvs = num_tree_cvs, k_cvs = k_cvs, nu_q_cvs = nu_q_cvs,
                    num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, 
                    serialize = TRUE, seed = model_seed), argY)
    fit_y = do.call(bartMachineCV, args_Y)
  }else{
    args_Y = c(list(X = data[, c(X_name, 't', 'ps', 'm', 'm0', 'm1')], y = data$y, num_trees = num_trees,
                    num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, 
                    serialize = TRUE, seed = model_seed), argY)
    fit_y = do.call(bartMachine, args_Y)
  }

  ### Posterior sampling
  by_value = ifelse(nrow(data_y_m0) > 5 * n_process_samples, n_process_samples, as.integer(nrow(data_y_m0) / 5))
  breaks = c(seq(from = 1, to = nrow(data_y_m0), by = by_value), nrow(data_y_m0) + 1)

  get_post_samples_parallel = function(data_set, t_value) {
    data_set$t = t_value
    mc.cores = max(1, detectCores() - 1)
    cl = makeCluster(mc.cores)
    clusterExport(cl, varlist = c('fit_y', 'data_set', 'breaks', 'bart_machine_get_posterior'), envir = environment())
    out = parLapply(cl, seq_along(breaks)[-length(breaks)], function(i) {
      bart_machine_get_posterior(fit_y, data_set[breaks[i]:(breaks[i+1] - 1), ])$y_hat_posterior_samples
    })
    stopCluster(cl)
    do.call(rbind, out)
 }
  
  cat('\n', 'Posterior Sampling: y0m0', '\n', sep = '')
  post_y0m0 = get_post_samples_parallel(data_y_m0, t_value = 0)
  cat('\n', 'Posterior Sampling: y1m0', '\n', sep = '')
  post_y1m0 = get_post_samples_parallel(data_y_m0, t_value = 1)
  cat('\n', 'Posterior Sampling: y1m1', '\n', sep = '')
  post_y1m1 = get_post_samples_parallel(data_y_m1, t_value = 1)

  ### Summarize
  TE = c(); TE.l = c(); TE.u = c(); TE.var = c();
  NIE = c(); NIE.l = c(); NIE.u = c(); NIE.var = c();
  NDE = c(); NDE.l = c(); NDE.u = c(); NDE.var = c();

  TE_all = NULL; NDE_all = NULL; NIE_all = NULL; 
  for(i in 1:n){
    start = (i-1) * n_after_burn_in + 1; end = i * n_after_burn_in;
    post_y0m0_i = c(post_y0m0[start:end, ])
    post_y1m0_i = c(post_y1m0[start:end, ])
    post_y1m1_i = c(post_y1m1[start:end, ])
    eTE = (post_y1m1_i - post_y0m0_i) * y_sd
    eNIE = (post_y1m1_i - post_y1m0_i) * y_sd
    eNDE = (post_y1m0_i - post_y0m0_i) * y_sd

    TE = c(TE, mean(eTE)); NIE = c(NIE, mean(eNIE)); NDE = c(NDE, mean(eNDE));
    TE.l = c(TE.l, quantile(eTE, 0.025)); NIE.l = c(NIE.l, quantile(eNIE, 0.025)); NDE.l = c(NDE.l, quantile(eNDE, 0.025));
    TE.u = c(TE.u, quantile(eTE, 0.975)); NIE.u = c(NIE.u, quantile(eNIE, 0.975)); NDE.u = c(NDE.u, quantile(eNDE, 0.975));
    TE.var = c(TE.var, var(eTE)); NIE.var = c(NIE.var, var(eNIE)); NDE.var = c(NDE.var, var(eNDE));

  }

  ### Empirical Bayes Shrinkage
  if(emb_shrink){
      
    ### Obtain clusters
    dist_matrix = dist(scale(data[, X_name]))
    hc = hclust(dist_matrix)
    clusters = cutree(hc, k = as.integer(n / cluster_size))

    ### Prepare new columns for shrunk estimates
    TE.emb = TE; TE.l.emb = TE.l; TE.u.emb = TE.u; 
    NIE.emb = NIE; NIE.l.emb = NIE.l; NIE.u.emb = NIE.u; 
    NDE.emb = NDE; NDE.l.emb = NDE.l; NDE.u.emb = NDE.u;  

    for(cl in unique(clusters)) {
      idx = which(clusters == cl)
      
      ### If only one in the cluster, skip
      if(length(idx) == 1) {
        next
      }     
      
      ### Shrunk NDE
      nde_hat = NDE[idx]
      nde_var = NDE.var[idx]
      nde_lower = NDE.l[idx]
      nde_upper = NDE.u[idx]
      mu_hat_nde = mean(nde_hat)
      tau2_nde = var(nde_hat)
      w_nde = tau2_nde / (tau2_nde + nde_var)
      NDE.emb[idx] = w_nde * nde_hat + (1 - w_nde) * mu_hat_nde
      NDE.l.emb[idx] = NDE.emb[idx] + (nde_lower - nde_hat)
      NDE.u.emb[idx] = NDE.emb[idx] + (nde_upper - nde_hat)
      
      ### Shrunk NIE
      nie_hat = NIE[idx]
      nie_var = NIE.var[idx]
      nie_lower = NIE.l[idx]
      nie_upper = NIE.u[idx]
      mu_hat_nie = mean(nie_hat)
      tau2_nie = var(nie_hat)
      w_nie = tau2_nie / (tau2_nie + nie_var)
      NIE.emb[idx] = w_nie * nie_hat + (1 - w_nie) * mu_hat_nie
      NIE.l.emb[idx] = NIE.emb[idx] + (nie_lower - nie_hat)
      NIE.u.emb[idx] = NIE.emb[idx] + (nie_upper - nie_hat)
      
      ### Deal with TE
      TE.emb[idx] = NIE.emb[idx] + NDE.emb[idx]
      TE.l.emb[idx] = TE.emb[idx] + (TE.l[idx] - TE[idx])
      TE.u.emb[idx] = TE.emb[idx] + (TE.u[idx] - TE[idx])
            
    }

    ### Assign shrunk results
    TE = TE.emb; TE.l = TE.l.emb; TE.u = TE.u.emb
    NIE = NIE.emb; NIE.l = NIE.l.emb; NIE.u = NIE.u.emb
    NDE = NDE.emb; NDE.l = NDE.l.emb; NDE.u = NDE.u.emb
  
  }

  return(list(data = data,
              X_name = X_name,
              h_effects = data.frame(
                TE = TE, TE.l = TE.l, TE.u = TE.u,
                NDE = NDE, NDE.l = NDE.l, NDE.u = NDE.u,
                NIE = NIE, NIE.l = NIE.l, NIE.u = NIE.u),
              p_effects = lapply(list(
                TE = mean(TE), TE.l = quantile(TE, 0.025), TE.u = quantile(TE, 0.975),
                NIE = mean(NIE), NIE.l = quantile(NIE, 0.025), NIE.u = quantile(NIE, 0.975),
                NDE = mean(NDE), NDE.l = quantile(NDE, 0.025), NDE.u = quantile(NDE, 0.975)
              ), unname)))

}