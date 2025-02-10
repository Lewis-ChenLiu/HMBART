#' @importFrom mgcv gam
#' @importFrom ggplot2 ggplot aes geom_linerange geom_point geom_line geom_ribbon labs theme_minimal ylim
#' @importFrom gridExtra grid.arrange
#' @export
dependentplot = function(hmbart_obj, varname, ylims_NDE = NULL, ylims_NIE = NULL) {
  
  hmbart_obj$effects = hmbart_obj$h_effects
  ### Prepare the data for NDE
  X = data.frame(
    var = hmbart_obj$data[, varname], 
    NDE = hmbart_obj$effects$NDE, 
    NIE = hmbart_obj$effects$NIE)
  not_sig_NDE = as.numeric(hmbart_obj$effects$NDE.l < 0) * as.numeric(hmbart_obj$effects$NDE.u > 0)
  dat_NDE = data.frame(
    var = X$var,
    NDE = X$NDE,
    NDE_0 = ifelse(not_sig_NDE == 1, X$NDE, NA),
    NDE_1 = ifelse(not_sig_NDE == 0, X$NDE, NA),
    NDE_0.l = ifelse(not_sig_NDE == 1, hmbart_obj$effects$NDE.l, NA),
    NDE_0.u = ifelse(not_sig_NDE == 1, hmbart_obj$effects$NDE.u, NA),
    NDE_1.l = ifelse(not_sig_NDE == 0, hmbart_obj$effects$NDE.l, NA),
    NDE_1.u = ifelse(not_sig_NDE == 0, hmbart_obj$effects$NDE.u, NA)
  )
  ### Splines fit for NDE
  gam_model = gam(NDE ~ s(var), data = X)
  pred = predict(gam_model, newdata = X, se.fit = TRUE)
  dat_NDE$NDE_gam = pred$fit
  dat_NDE$NDE_gam.l = pred$fit - 1.96 * pred$se.fit
  dat_NDE$NDE_gam.u = pred$fit + 1.96 * pred$se.fit
  
  ### Plot for NDE
  if(is.null(ylims_NDE)){
    ylims = c(min(hmbart_obj$effects$NDE.l) - 0.5, max(hmbart_obj$effects$NDE.u) + 0.5)
  }else{
    ylims = ylims_NDE
  }
  if(mean(is.na(dat_NDE$NDE_1)) == 1){
    plot_NDE_est = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NDE_0.l, ymax = NDE_0.u), color = "lightgrey", width = 0.1) +
      geom_point(aes(y = NDE_0), color = "black", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +   theme_minimal() 
  }else if(mean(is.na(dat_NDE$NDE_0)) == 1){
    plot_NDE_est = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NDE_1.l, ymax = NDE_1.u), color = "#fee090", width = 0.1) +
      geom_point(aes(y = NDE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +   theme_minimal()  
  }else{
    plot_NDE_est = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NDE_0.l, ymax = NDE_0.u), color = "lightgrey", width = 0.1) +
      geom_linerange(aes(ymin = NDE_1.l, ymax = NDE_1.u), color = "#fee090", width = 0.1) +
      geom_point(aes(y = NDE_0), color = "black", size = 0.5) +
      geom_point(aes(y = NDE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +   theme_minimal()  
  }
  plot_NDE_gam = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) + 
    geom_point(aes(y = NDE), color = "#1f77b4", size = 0.5) +
    geom_line(aes(y = NDE_gam), color = "#1f77b4") +
    geom_ribbon(aes(ymin = NDE_gam.l, ymax = NDE_gam.u), alpha = 0.2, fill = "#1f77b4") + 
    labs(
      title = "",
      x = varname,
      y = "NDE"
    ) +   theme_minimal()
  
  ### Prepare the data for NIE
  not_sig_NIE = as.numeric(hmbart_obj$effects$NIE.l < 0) * as.numeric(hmbart_obj$effects$NIE.u > 0)
  dat_NIE = data.frame(
    var = X$var,
    NIE = X$NIE,
    NIE_0 = ifelse(not_sig_NIE == 1, X$NIE, NA),
    NIE_1 = ifelse(not_sig_NIE == 0, X$NIE, NA),
    NIE_0.l = ifelse(not_sig_NIE == 1, hmbart_obj$effects$NIE.l, NA),
    NIE_0.u = ifelse(not_sig_NIE == 1, hmbart_obj$effects$NIE.u, NA),
    NIE_1.l = ifelse(not_sig_NIE == 0, hmbart_obj$effects$NIE.l, NA),
    NIE_1.u = ifelse(not_sig_NIE == 0, hmbart_obj$effects$NIE.u, NA)
  )
  ### Splines fit for NIE
  gam_model = gam(NIE ~ s(var), data = X)
  pred = predict(gam_model, newdata = X, se.fit = TRUE)
  dat_NIE$NIE_gam = pred$fit
  dat_NIE$NIE_gam.l = pred$fit - 1.96 * pred$se.fit
  dat_NIE$NIE_gam.u = pred$fit + 1.96 * pred$se.fit
  
  ### Plot for NIE
  if(is.null(ylims_NIE)){
    ylims = c(min(hmbart_obj$effects$NIE.l) - 0.5, max(hmbart_obj$effects$NIE.u) + 0.5)
  }else{
    ylims = ylims_NIE
  }
  if(mean(is.na(dat_NIE$NIE_1)) == 1){
    plot_NIE_est = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NIE_0.l, ymax = NIE_0.u), color = "lightgrey", width = 0.1) +
      geom_point(aes(y = NIE_0), color = "black", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +   theme_minimal() 
  }else if(mean(is.na(dat_NIE$NIE_0)) == 1){
    plot_NIE_est = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NIE_1.l, ymax = NIE_1.u), color = "#fee090", width = 0.1) +
      geom_point(aes(y = NIE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +   theme_minimal() 
  }else{
    plot_NIE_est = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NIE_0.l, ymax = NIE_0.u), color = "lightgrey", width = 0.1) +
      geom_linerange(aes(ymin = NIE_1.l, ymax = NIE_1.u), color = "#fee090", width = 0.1) +
      geom_point(aes(y = NIE_0), color = "black", size = 0.5) +
      geom_point(aes(y = NIE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +   theme_minimal()  
  }
  plot_NIE_gam = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) + 
    geom_point(aes(y = NIE), color = "#1f77b4", size = 0.5) +
    geom_line(aes(y = NIE_gam), color = "#1f77b4") +
    geom_ribbon(aes(ymin = NIE_gam.l, ymax = NIE_gam.u), alpha = 0.2, fill = "#1f77b4") + 
    labs(
      title = "",
      x = varname,
      y = "NIE"
    ) +   theme_minimal()
  
  ### Combine and display
  combined_plot = grid.arrange(plot_NDE_est, plot_NIE_est,
                               plot_NDE_gam, plot_NIE_gam, ncol = 2, nrow = 2)
  print(combined_plot)
  
}
