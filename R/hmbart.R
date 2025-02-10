#' @importFrom bartMachine bartMachine bartMachineCV bart_machine_get_posterior
#' @export
hmbart = function(data, X, t, m, y, CV = FALSE,
                  num_trees = 100, num_tree_cvs = c(100), k_cvs = c(2, 3, 5),
                  nu_q_cvs = list(c(3, 0.9), c(3, 0.99), c(10, 0.75)),
                  n_burn_in = 2000, n_after_burn_in = 500, n_process_samples = 1e5, seed = 42, fix = FALSE) {
  
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
  fit_ps = bartMachine(X = data[, X_name], y = data$t, num_trees = num_trees,
                       num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = model_seed)
  data$ps = fit_ps$y_hat_train
  rm(fit_ps)

  ### M model
  if(CV){
    fit_m = bartMachineCV(X = data[, c(X_name, 't', 'ps')], y = data$m,
                          num_tree_cvs = num_tree_cvs, k_cvs = k_cvs, nu_q_cvs = nu_q_cvs,
                          num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = model_seed)
  }else{
    fit_m = bartMachine(X = data[, c(X_name, 't', 'ps')], y = data$m, num_trees = num_trees,
                        num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = model_seed)
  }

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
  data_y_m0 = NULL; data_y_m1 = NULL;
  data_y = data[, c(X_name, 't', 'ps', 'm', 'm0', 'm1')]
  for(i in 1:n){
    data_i = data_y[replicate(n_after_burn_in, i), ]
    data_i$m = m0[i, ]
    data_y_m0 = rbind(data_y_m0, data_i)
    data_i$m = m1[i, ]
    data_y_m1 = rbind(data_y_m1, data_i)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  rm(fit_m)

  ### Y model
  if(CV){
    fit_y = bartMachineCV(X = data[, c(X_name, 't', 'ps', 'm', 'm0', 'm1')], y = data$y,
                          num_tree_cvs = num_tree_cvs, k_cvs = k_cvs, nu_q_cvs = nu_q_cvs,
                          num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = model_seed)
  }else{
    fit_y = bartMachine(X = data[, c(X_name, 't', 'ps', 'm', 'm0', 'm1')], y = data$y, num_trees = num_trees,
                        num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = model_seed)
  }

  ### Predict
  by_value = ifelse(nrow(data_y_m0) > 5 * n_process_samples, n_process_samples, as.integer(nrow(data_y_m0) / 5))
  breaks = c(seq(from = 1, to = nrow(data_y_m0), by = by_value), nrow(data_y_m0) + 1)
  post_y0m0 = NULL; post_y1m0 = NULL; post_y1m1 = NULL;

  ### Posterior sampling: y0m0
  cat('\n', 'Posterior Sampling: y0m0', '\n', sep = '')
  pb = txtProgressBar(min = 1, max = length(breaks) - 1, style = 3)
  data_y_m0$t = 0
  for(i in 1:(length(breaks)-1)){
    post_y0m0 = rbind(post_y0m0, bart_machine_get_posterior(fit_y, data_y_m0[breaks[i]:(breaks[i+1]-1),])$y_hat_posterior_samples)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  ### Posterior sampling: y1m0
  cat('\n', 'Posterior Sampling: y1m0', '\n', sep = '')
  pb = txtProgressBar(min = 1, max = length(breaks) - 1, style = 3)
  data_y_m0$t = 1
  for(i in 1:(length(breaks)-1)){
    post_y1m0 = rbind(post_y1m0, bart_machine_get_posterior(fit_y, data_y_m0[breaks[i]:(breaks[i+1]-1),])$y_hat_posterior_samples)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  ### Posterior sampling: y1m1
  cat('\n', 'Posterior Sampling: y1m1', '\n', sep = '')
  pb = txtProgressBar(min = 1, max = length(breaks) - 1, style = 3)
  data_y_m1$t = 1
  for(i in 1:(length(breaks)-1)){
    post_y1m1 = rbind(post_y1m1, bart_machine_get_posterior(fit_y, data_y_m1[breaks[i]:(breaks[i+1]-1),])$y_hat_posterior_samples)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  ### Summarize
  TE = c(); TE.l = c(); TE.u = c();
  NDE = c(); NDE.l = c(); NDE.u = c();
  NIE = c(); NIE.l = c(); NIE.u = c();

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
