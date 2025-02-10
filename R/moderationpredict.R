#' @importFrom bartMachine bartMachine bart_machine_get_posterior
#' @export
moderationpredict = function(hmbart_obj, NDE_moderator_names = NULL, NDE_moderator_values = NULL,
                             NIE_moderator_names = NULL, NIE_moderator_values = NULL, 
                             num_trees = 100, n_burn_in = 2000, n_after_burn_in = 500, seed = 42) {
  
  set.seed(seed)
  set_bart_machine_num_cores(num_cores = 1)
  hmbart_obj$effects = hmbart_obj$h_effects
  
  ### Extract data
  NDE_data = data.frame(hmbart_obj$data[, NDE_moderator_names])
  NIE_data = data.frame(hmbart_obj$data[, NIE_moderator_names])
  colnames(NDE_data) = NDE_moderator_names
  colnames(NIE_data) = NIE_moderator_names
  
  ### Model
  NDE_model = bartMachine(X = NDE_data, y = hmbart_obj$effects$NDE, num_trees = num_trees,
                          num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = seed)
  NIE_model = bartMachine(X = NIE_data, y = hmbart_obj$effects$NIE, num_trees = num_trees,
                          num_burn_in = n_burn_in, num_iterations_after_burn_in = n_after_burn_in, seed = seed)
  
  ### Predict
  NDE_pred_raw = bart_machine_get_posterior(NDE_model, NDE_moderator_values)$y_hat_posterior_samples
  NIE_pred_raw = bart_machine_get_posterior(NIE_model, NIE_moderator_values)$y_hat_posterior_samples
  NDE_pred = data.frame(NDE = apply(NDE_pred_raw, 1, mean),
                        NDE.l = apply(NDE_pred_raw, 1, quantile, probs = c(0.025)),
                        NDE.u = apply(NDE_pred_raw, 1, quantile, probs = c(0.975)))
  NIE_pred = data.frame(NIE = apply(NIE_pred_raw, 1, mean),
                        NIE.l = apply(NIE_pred_raw, 1, quantile, probs = c(0.025)),
                        NIE.u = apply(NIE_pred_raw, 1, quantile, probs = c(0.975)))
  
  return(list(NDE_pred = NDE_pred, NIE_pred = NIE_pred))
  
}
