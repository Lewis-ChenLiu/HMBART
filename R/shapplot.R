#' @importFrom xgboost xgboost
#' @importFrom ggplot2 ggtitle theme element_text
#' @importFrom gridExtra grid.arrange
#' @export
shapplot = function(hmbart_obj, moderator_names = NULL, TE = FALSE, plot_each = FALSE, seed = 42) {
  
  library(SHAPforxgboost)
  ### Model
  if(is.null(moderator_names)){moderator_names = hmbart_obj$X_name}
  X = as.matrix(hmbart_obj$data[, moderator_names])
  colnames(X) = moderator_names
  hmbart_obj$effects = hmbart_obj$h_effects
  
  set.seed(seed)
  TE_model = xgboost(data = X, label = hmbart_obj$effects$TE, nrounds = 50, verbose = F)
  NDE_model = xgboost(data = X, label = hmbart_obj$effects$NDE, nrounds = 50, verbose = F)
  NIE_model = xgboost(data = X, label = hmbart_obj$effects$NIE, nrounds = 50, verbose = F)
  
  ### Plot
  shap_TE = shap.prep(TE_model, X_train = X)
  shap_NDE = shap.prep(NDE_model, X_train = X)
  shap_NIE = shap.prep(NIE_model, X_train = X)
  shap_TE = shap.plot.summary(shap_TE, min_color_bound = "#fee090", max_color_bound = "#a50026") + 
    ggtitle("SHAP Summary - TE") + theme(plot.title = element_text(hjust = 0.5, ))
  shap_NDE = shap.plot.summary(shap_NDE, min_color_bound = "#fee090", max_color_bound = "#a50026") + 
    ggtitle("SHAP Summary - NDE") + theme(plot.title = element_text(hjust = 0.5, ))
  shap_NIE = shap.plot.summary(shap_NIE, min_color_bound = "#fee090", max_color_bound = "#a50026") + 
    ggtitle("SHAP Summary - NIE") + theme(plot.title = element_text(hjust = 0.5, ))
  
  ### Display
  if(TE){
    if(plot_each){
      print(shap_TE); print(shap_NDE); print(shap_NIE);
    }else{
      invisible(grid.arrange(shap_TE, shap_NDE, shap_NIE, ncol = 1, nrow = 3))
    }
  }else{
    if(plot_each){
      print(shap_NDE); print(shap_NIE);
    }else{
      invisible(grid.arrange(shap_NDE, shap_NIE, ncol = 1, nrow = 2))
    }
  }
  
}
