#' @importFrom rpart rpart rpart.control path.rpart
#' @importFrom ggplot2 ggplot aes geom_density facet_wrap theme_minimal labs
#' @importFrom gridExtra grid.arrange
#' @export
subgroupplot = function(hmbart_obj, moderator_names = NULL, cp = 0.01, TE = FALSE, seed = 42) {
  
  set.seed(seed)
  hmbart_obj$effects = hmbart_obj$h_effects
  ### Model
  if(is.null(moderator_names)){moderator_names = hmbart_obj$X_name}
  X = data.frame(hmbart_obj$data[, moderator_names])
  colnames(X) = moderator_names
  
  tree_TE = rpart(hmbart_obj$effects$TE ~ ., data = X, control = rpart.control(cp = cp))
  tree_NDE = rpart(hmbart_obj$effects$NDE ~ ., data = X, control = rpart.control(cp = cp))
  tree_NIE = rpart(hmbart_obj$effects$NIE ~ ., data = X, control = rpart.control(cp = cp))
  
  ### Assign leaf for each data points
  getsubgroup = function(tree, effect){
    
    node_assignments = tree$where
    existing_nodes = as.numeric(row.names(tree$frame))
    valid_nodes = intersect(unique(node_assignments), existing_nodes)
    rules = path.rpart(tree, nodes = valid_nodes, print.it = FALSE)
    rule_strings = lapply(rules, function(x) paste(x[x != "root"], collapse = ", "))
    names = names(rule_strings)
    value = c()
    group = c()
    for(i in 1:length(names)){
      name = names[i]
      subgroup = effect[node_assignments == name]
      value = c(value, subgroup)
      group = c(group, rep(rule_strings[[i]], length(subgroup)))
    }
    df = data.frame(value = value, group = group)
    return(df)
    
  }
  
  ### Plot and display
  if(TE){
    
    df = getsubgroup(tree_TE, hmbart_obj$effects$TE)
    p1 = ggplot(df, aes(x = value)) +
      geom_density(fill = "lightblue", alpha = 0.5) +
      facet_wrap(~ group, ncol = 1, scales = "free_y") +
      theme_minimal() +
      labs(x = "TE", y = "Density", title = "Subgroup Summary - TE")
  }
  
  ### NDE
  df = getsubgroup(tree_NDE, hmbart_obj$effects$NDE)
  p2 = ggplot(df, aes(x = value)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    facet_wrap(~ group, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(x = "NDE", y = "Density", title = "Subgroup Summary - NDE")
  ### NIE
  df = getsubgroup(tree_NIE, hmbart_obj$effects$NIE)
  p3 = ggplot(df, aes(x = value)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    facet_wrap(~ group, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(x = "NIE", y = "Density", title = "Subgroup Summary - NIE")
  
  if(TE){
    combined_plot = grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)
  } else {
    combined_plot = grid.arrange(p2, p3, ncol = 2, nrow = 1)
  }
  invisible(combined_plot)
  
}