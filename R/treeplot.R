#' @importFrom rpart rpart rpart.control
#' @importFrom rpart.plot prp
#' @export
treeplot = function(hmbart_obj, moderator_names = NULL, cp = 0.01, seed = 42) {
  
  set.seed(seed)
  hmbart_obj$effects = hmbart_obj$h_effects
  ### Model
  if(is.null(moderator_names)){moderator_names = hmbart_obj$X_name}
  X = data.frame(hmbart_obj$data[, moderator_names])
  colnames(X) = moderator_names

  tree_NDE = rpart(hmbart_obj$effects$NDE ~ ., data = X, control = rpart.control(cp = cp))
  tree_NIE = rpart(hmbart_obj$effects$NIE ~ ., data = X, control = rpart.control(cp = cp))
  
  ### Plot and display
  layout(matrix(1:2, ncol = 2), widths = c(1, 1))
  prp(tree_NDE, main = "", yesno = 1, extra = "auto")
  title(main = "Tree Summary - NDE", adj = 0.5)
  prp(tree_NIE, main = "", yesno = 1, extra = "auto")
  title(main = "Tree Summary - NIE", adj = 0.5)
  layout(1)
  
}
