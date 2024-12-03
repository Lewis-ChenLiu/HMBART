#' @importFrom xgboost xgboost
#' @importFrom ggplot2 ggtitle theme element_text
#' @importFrom gridExtra grid.arrange
#' @export
shapplot = function(hmbart_obj) {
  
  library(SHAPforxgboost)
  ### Model
  X = as.matrix(hmbart_obj$data[, hmbart_obj$X_name])
  colnames(X) = hmbart_obj$X_name
  NDE_model = xgboost(data = X, label = hmbart_obj$effects$NDE, nrounds = 50, verbose = F)
  NIE_model = xgboost(data = X, label = hmbart_obj$effects$NIE, nrounds = 50, verbose = F)
  
  ### Plot
  shap_NDE = shap.prep(NDE_model, X_train = X)
  shap_NIE = shap.prep(NIE_model, X_train = X)
  shap_NDE = shap.plot.summary(shap_NDE, min_color_bound = "#fee090", max_color_bound = "#a50026") + 
    ggtitle("SHAP Summary - NDE") + theme(plot.title = element_text(hjust = 0.5, ))
  shap_NIE = shap.plot.summary(shap_NIE, min_color_bound = "#fee090", max_color_bound = "#a50026") + 
    ggtitle("SHAP Summary - NIE") + theme(plot.title = element_text(hjust = 0.5, ))
  
  ### Display
  combined_plot = grid.arrange(shap_NDE, shap_NIE, ncol = 2, nrow = 1)
  print(combined_plot)

}
