#' @importFrom mgcv gam
#' @importFrom ggplot2 ggplot aes geom_linerange geom_point geom_line geom_smooth geom_ribbon labs theme_minimal ylim
#' @importFrom gridExtra grid.arrange
#' @export
dependentplot = function(hmbart_obj, varname, categorical = FALSE, TE = FALSE, ylims_TE = NULL, ylims_NDE = NULL, ylims_NIE = NULL) {
  
  hmbart_obj$effects = hmbart_obj$h_effects
  ### For variable have unique values less than 10 (default for spline)
  k = min(max(length(unique(hmbart_obj$data[, varname])) - 2, 3), 10)
  ### Prepare the data for TE
  X = data.frame(
    var = hmbart_obj$data[, varname], 
    TE = hmbart_obj$effects$TE, 
    NDE = hmbart_obj$effects$NDE, 
    NIE = hmbart_obj$effects$NIE)
  
  ### Settings for categorical
  if(categorical){
    k = 1
    X$var = as.factor(X$var)
  }
  not_sig_TE = as.numeric(hmbart_obj$effects$TE.l < 0) * as.numeric(hmbart_obj$effects$TE.u > 0)
  dat_TE = data.frame(
    var = X$var,
    TE = X$TE,
    TE_0 = ifelse(not_sig_TE == 1, X$TE, NA),
    TE_1 = ifelse(not_sig_TE == 0, X$TE, NA),
    TE_0.l = ifelse(not_sig_TE == 1, hmbart_obj$effects$TE.l, NA),
    TE_0.u = ifelse(not_sig_TE == 1, hmbart_obj$effects$TE.u, NA),
    TE_1.l = ifelse(not_sig_TE == 0, hmbart_obj$effects$TE.l, NA),
    TE_1.u = ifelse(not_sig_TE == 0, hmbart_obj$effects$TE.u, NA)
  )
  ### Splines fit for TE
  if(k > 3){
    gam_model = gam(TE ~ s(var, k = k), data = X)
    pred = predict(gam_model, newdata = X, se.fit = TRUE)
    dat_TE$TE_gam = pred$fit
    dat_TE$TE_gam.l = pred$fit - 1.96 * pred$se.fit
    dat_TE$TE_gam.u = pred$fit + 1.96 * pred$se.fit
  }
  
  ### Plot for TE
  if(is.null(ylims_TE)){
    ylims = c(min(hmbart_obj$effects$TE.l) - 0.5, max(hmbart_obj$effects$TE.u) + 0.5)
  }else{
    ylims = ylims_TE
  }
  if(mean(is.na(dat_TE$TE_1)) == 1){
    plot_TE_est = ggplot(dat_TE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = TE_0.l, ymax = TE_0.u), color = "lightgrey") +
      geom_point(aes(y = TE_0), color = "black", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "TE"
      ) +   theme_minimal() 
    plot_TE_est1 = ggplot(dat_TE, aes(x = var, y = TE_0)) +
      ylim(ylims) +
      geom_smooth(aes(y = TE_0.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = TE_0.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_point(color = "black", size = 0.5) +
      geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = "black") +
      labs(
        title = "",
        x = varname,
        y = "TE"
      ) +
      theme_minimal()
  }else if(mean(is.na(dat_TE$TE_0)) == 1){
    plot_TE_est = ggplot(dat_TE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = TE_1.l, ymax = TE_1.u), color = "#fee090") +
      geom_point(aes(y = TE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "TE"
      ) +   theme_minimal()  
    plot_TE_est1 = ggplot(dat_TE, aes(x = var, y = TE_1)) +
      ylim(ylims) +
      geom_smooth(aes(y = TE_1.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_smooth(aes(y = TE_1.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_point(color = "#ff7f0e", size = 0.5) +
      geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = "#ff7f0e") +
      labs(
        title = "",
        x = varname,
        y = "TE"
      ) +
      theme_minimal()
  }else{
    plot_TE_est = ggplot(dat_TE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = TE_0.l, ymax = TE_0.u), color = "lightgrey") +
      geom_linerange(aes(ymin = TE_1.l, ymax = TE_1.u), color = "#fee090") +
      geom_point(aes(y = TE_0), color = "black", size = 0.5) +
      geom_point(aes(y = TE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "TE"
      ) +   theme_minimal()  
    plot_TE_est1 = ggplot(dat_TE, aes(x = var)) +
      ylim(ylims) +
      geom_smooth(aes(y = TE_0.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = TE_0.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = TE_1.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_smooth(aes(y = TE_1.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_point(aes(y = TE_0), color = "black", size = 0.5) +
      geom_point(aes(y = TE_1), color = "#ff7f0e", size = 0.5) +
      geom_smooth(aes(y = TE_0), method = "gam", formula = y ~ s(x), se = FALSE, color = "black") +
      geom_smooth(aes(y = TE_1), method = "gam", formula = y ~ s(x), se = FALSE, color = "#ff7f0e") +
      labs(
        title = "",
        x = varname,
        y = "TE"
      ) +
      theme_minimal()
  }
  if(k > 3){
    plot_TE_gam = ggplot(dat_TE, aes(x = var)) + ylim(ylims) + 
      geom_point(aes(y = TE), color = "#1f77b4", size = 0.5) +
      geom_line(aes(y = TE_gam), color = "#1f77b4") +
      geom_ribbon(aes(ymin = TE_gam.l, ymax = TE_gam.u), alpha = 0.2, fill = "#1f77b4") + 
      labs(
        title = "",
        x = varname,
        y = "TE"
      ) +   theme_minimal()
  }
  
  ### Prepare the data for NDE
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
  if(k > 3){
    gam_model = gam(NDE ~ s(var, k = k), data = X)
    pred = predict(gam_model, newdata = X, se.fit = TRUE)
    dat_NDE$NDE_gam = pred$fit
    dat_NDE$NDE_gam.l = pred$fit - 1.96 * pred$se.fit
    dat_NDE$NDE_gam.u = pred$fit + 1.96 * pred$se.fit
  }
  
  ### Plot for NDE
  if(is.null(ylims_NDE)){
    ylims = c(min(hmbart_obj$effects$NDE.l) - 0.5, max(hmbart_obj$effects$NDE.u) + 0.5)
  }else{
    ylims = ylims_NDE
  }
  if(mean(is.na(dat_NDE$NDE_1)) == 1){
    plot_NDE_est = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NDE_0.l, ymax = NDE_0.u), color = "lightgrey") +
      geom_point(aes(y = NDE_0), color = "black", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +   theme_minimal() 
    plot_NDE_est1 = ggplot(dat_NDE, aes(x = var)) +
      ylim(ylims) +
      geom_smooth(aes(y = NDE_0.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = NDE_0.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_point(aes(y = NDE_0), color = "black", size = 0.5) +
      geom_smooth(aes(y = NDE_0), method = "gam", formula = y ~ s(x), se = FALSE, color = "black") +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +
      theme_minimal()
  }else if(mean(is.na(dat_NDE$NDE_0)) == 1){
    plot_NDE_est = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NDE_1.l, ymax = NDE_1.u), color = "#fee090") +
      geom_point(aes(y = NDE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +   theme_minimal()  
    plot_NDE_est1 = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) +
      geom_smooth(aes(y = NDE_1.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_smooth(aes(y = NDE_1.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_point(aes(y = NDE_1), color = "#ff7f0e", size = 0.5) +
      geom_smooth(aes(y = NDE_1), method = "gam", formula = y ~ s(x), se = FALSE, color = "#ff7f0e") +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +
      theme_minimal()
  }else{
    plot_NDE_est = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NDE_0.l, ymax = NDE_0.u), color = "lightgrey") +
      geom_linerange(aes(ymin = NDE_1.l, ymax = NDE_1.u), color = "#fee090") +
      geom_point(aes(y = NDE_0), color = "black", size = 0.5) +
      geom_point(aes(y = NDE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +   theme_minimal()  
    plot_NDE_est1 = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) +
      geom_smooth(aes(y = NDE_0.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = NDE_0.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = NDE_1.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_smooth(aes(y = NDE_1.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_point(aes(y = NDE_0), color = "black", size = 0.5) +
      geom_point(aes(y = NDE_1), color = "#ff7f0e", size = 0.5) +
      geom_smooth(aes(y = NDE_0), method = "gam", formula = y ~ s(x), se = FALSE, color = "black") +
      geom_smooth(aes(y = NDE_1), method = "gam", formula = y ~ s(x), se = FALSE, color = "#ff7f0e") +
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +
      theme_minimal()
  }
  if(k > 3){
    plot_NDE_gam = ggplot(dat_NDE, aes(x = var)) + ylim(ylims) + 
      geom_point(aes(y = NDE), color = "#1f77b4", size = 0.5) +
      geom_line(aes(y = NDE_gam), color = "#1f77b4") +
      geom_ribbon(aes(ymin = NDE_gam.l, ymax = NDE_gam.u), alpha = 0.2, fill = "#1f77b4") + 
      labs(
        title = "",
        x = varname,
        y = "NDE"
      ) +   theme_minimal()
  }
  
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
  if(k > 3){
    gam_model = gam(NIE ~ s(var, k = k), data = X)
    pred = predict(gam_model, newdata = X, se.fit = TRUE)
    dat_NIE$NIE_gam = pred$fit
    dat_NIE$NIE_gam.l = pred$fit - 1.96 * pred$se.fit
    dat_NIE$NIE_gam.u = pred$fit + 1.96 * pred$se.fit
  }
  
  ### Plot for NIE
  if(is.null(ylims_NIE)){
    ylims = c(min(hmbart_obj$effects$NIE.l) - 0.5, max(hmbart_obj$effects$NIE.u) + 0.5)
  }else{
    ylims = ylims_NIE
  }
  if(mean(is.na(dat_NIE$NIE_1)) == 1){
    plot_NIE_est = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NIE_0.l, ymax = NIE_0.u), color = "lightgrey") +
      geom_point(aes(y = NIE_0), color = "black", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +   theme_minimal() 
    plot_NIE_est1 = ggplot(dat_NIE, aes(x = var)) +
      ylim(ylims) +
      geom_smooth(aes(y = NIE_0.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = NIE_0.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_point(aes(y = NIE_0), color = "black", size = 0.5) +
      geom_smooth(aes(y = NIE_0), method = "gam", formula = y ~ s(x), se = FALSE, color = "black") +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +
      theme_minimal()
  }else if(mean(is.na(dat_NIE$NIE_0)) == 1){
    plot_NIE_est = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NIE_1.l, ymax = NIE_1.u), color = "#fee090") +
      geom_point(aes(y = NIE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +   theme_minimal() 
    plot_NIE_est1 = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) +
      geom_smooth(aes(y = NIE_1.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_smooth(aes(y = NIE_1.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_point(color = "#ff7f0e", size = 0.5) +
      geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = "#ff7f0e") +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +
      theme_minimal()
  }else{
    plot_NIE_est = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) + 
      geom_linerange(aes(ymin = NIE_0.l, ymax = NIE_0.u), color = "lightgrey") +
      geom_linerange(aes(ymin = NIE_1.l, ymax = NIE_1.u), color = "#fee090") +
      geom_point(aes(y = NIE_0), color = "black", size = 0.5) +
      geom_point(aes(y = NIE_1), color = "#ff7f0e", size = 0.5) +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +   theme_minimal()  
    plot_NIE_est1 = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) +
      geom_smooth(aes(y = NIE_0.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = NIE_0.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "lightgrey", linetype = "dashed") +
      geom_smooth(aes(y = NIE_1.l), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_smooth(aes(y = NIE_1.u), method = "gam", formula = y ~ s(x), se = FALSE, color = "#fee090", linetype = "dashed") +
      geom_point(aes(y = NIE_0), color = "black", size = 0.5) +
      geom_point(aes(y = NIE_1), color = "#ff7f0e", size = 0.5) +
      geom_smooth(aes(y = NIE_0), method = "gam", formula = y ~ s(x), se = FALSE, color = "black") +
      geom_smooth(aes(y = NIE_1), method = "gam", formula = y ~ s(x), se = FALSE, color = "#ff7f0e") +
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +
      theme_minimal()
  }
  if(k > 3){
    plot_NIE_gam = ggplot(dat_NIE, aes(x = var)) + ylim(ylims) + 
      geom_point(aes(y = NIE), color = "#1f77b4", size = 0.5) +
      geom_line(aes(y = NIE_gam), color = "#1f77b4") +
      geom_ribbon(aes(ymin = NIE_gam.l, ymax = NIE_gam.u), alpha = 0.2, fill = "#1f77b4") + 
      labs(
        title = "",
        x = varname,
        y = "NIE"
      ) +   theme_minimal()
  }
  
  ### Combine and display
  if(TE){
    if(k > 3){
      combined_plot = grid.arrange(plot_TE_est, plot_NDE_est, plot_NIE_est,
                                   plot_TE_est1, plot_NDE_est1, plot_NIE_est1,
                                   plot_TE_gam, plot_NDE_gam, plot_NIE_gam, ncol = 3, nrow = 3)
    } else {
      combined_plot = grid.arrange(plot_TE_est, plot_NDE_est, plot_NIE_est, 
                                   plot_TE_est1, plot_NDE_est1, plot_NIE_est1,ncol = 3, nrow = 2)
    }
  }else{
    if(k > 3){
      combined_plot = grid.arrange(plot_NDE_est, plot_NIE_est,
                                   plot_NDE_est1, plot_NIE_est1,
                                   plot_NDE_gam, plot_NIE_gam, ncol = 2, nrow = 3)
    } else {
      combined_plot = grid.arrange(plot_NDE_est, plot_NIE_est, 
                                   plot_NDE_est1, plot_NIE_est1, ncol = 2, nrow = 2)
    }
  }
  invisible(combined_plot)
  
}
