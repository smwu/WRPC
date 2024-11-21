#===============================
# WRPC plotting functions
# ==============================

library(tidyverse)

plot_wrpc_global_pattern_profiles <- function(res, item_labels = NULL, 
                                         item_title = "Item",
                                  categ_labels = NULL, 
                                  categ_title = "Consumption Level",
                                  class_labels = NULL, 
                                  class_title = "Dietary Pattern", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca", "wrpc"))) {
    stop("res must be an object of class `swolca`, `wolca` or `wrpc`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  est_item_probs <- res$estimates$theta_global_med
  
  # Obtain theta modes as the category with highest probability for each item
  mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
  
  # Define item, latent class, and category names if not specified
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  } else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ",
                res$data_vars$J))
  }
  K <- dim(mode_item_probs)[2]
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", K))
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  } else if (length(categ_labels) != res$data_vars$R) {
    stop(paste0("length of categ_labels must equal the number of exposure categories, R = ", 
                res$data_vars$R))
  }
  rownames(mode_item_probs) <- item_labels
  colnames(mode_item_probs) <- class_labels
  mode_item_probs$Item <- rownames(mode_item_probs)
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Class <- Item <- Level <- NULL
  
  # Create plot
  mode_plot <- mode_item_probs %>% 
    tidyr::gather("Class", "Level", -Item) 
  mode_plot %>% ggplot2::ggplot(ggplot2::aes(x=Class, 
                                             y=factor(Item, 
                                                      levels = rev(item_labels)), 
                                             fill = factor(Level))) + 
    ggplot2::geom_tile(color="black", linewidth = 0.3) + 
    ggplot2::scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                               name = categ_title, labels = categ_labels) +
    ggplot2::labs(x = class_title, y = item_title) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete() + 
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "top")
}



plot_wrpc_local_pattern_profiles <- function(res, item_labels = NULL, 
                                              item_title = "Item",
                                              categ_labels = NULL, 
                                              categ_title = "Consumption Level",
                                              subgroup_labels = NULL, 
                                              subgroup_title = "Subgroup", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca", "wrpc"))) {
    stop("res must be an object of class `swolca`, `wolca` or `wrpc`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  est_item_probs <- res$estimates$theta_local_med
  
  # Obtain theta modes as the category with highest probability for each item
  mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
  
  # Define item, subgroup, and category names if not specified
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  } else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ",
                res$data_vars$J))
  }
  H <- dim(mode_item_probs)[2]
  if (is.null(subgroup_labels)) {
    subgroup_labels <- 1:H
  } else if (length(subgroup_labels) != H) {
    stop(paste0("length of subgroup_labels must equal the number of subgroups, H = ", H))
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  } else if (length(categ_labels) != res$data_vars$R) {
    stop(paste0("length of categ_labels must equal the number of exposure categories, R = ", 
                res$data_vars$R))
  }
  rownames(mode_item_probs) <- item_labels
  colnames(mode_item_probs) <- subgroup_labels
  mode_item_probs$Item <- rownames(mode_item_probs)
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Subgroup <- Item <- Level <- NULL
  
  # Create plot
  mode_plot <- mode_item_probs %>% 
    tidyr::gather("Subgroup", "Level", -Item) 
  mode_plot %>% ggplot2::ggplot(ggplot2::aes(x=factor(Subgroup, 
                                                      levels = subgroup_labels),
                                             y=factor(Item, 
                                                      levels = rev(item_labels)), 
                                             fill = factor(Level))) + 
    ggplot2::geom_tile(color="black", linewidth = 0.3) + 
    ggplot2::scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                               name = categ_title, labels = categ_labels) +
    ggplot2::labs(x = subgroup_title, y = item_title) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete() + 
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "top")
}


plot_wrpc_allocation <- function(res, item_labels = NULL, 
                                 item_title = "Item",
                                 subgroup_labels = NULL, 
                                 subgroup_title = "Subgroup", ...) {

  # Check object class
  if (!(class(res) %in% c("swolca", "wolca", "wrpc"))) {
    stop("res must be an object of class `swolca`, `wolca` or `wrpc`, resulting 
       from a call to one of these functions")
  }
  
  # Obtain nu estimates
  est_alloc_prob <- as.data.frame(t(res$estimates$nu_med))
  
  # Define item, subgroup, and category names if not specified
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  } else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ",
                res$data_vars$J))
  }
  H <- dim(est_alloc_prob)[2]
  if (is.null(subgroup_labels)) {
    subgroup_labels <- 1:H
  } else if (length(subgroup_labels) != H) {
    stop(paste0("length of subgroup_labels must equal the number of subgroups, H = ", H))
  }

  rownames(est_alloc_prob) <- item_labels
  colnames(est_alloc_prob) <- subgroup_labels
  est_alloc_prob$Item <- rownames(est_alloc_prob)
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Subgroup <- Item <- NULL
  
  # Create plot
      # mode_plot <- est_alloc_prob %>% 
      #   tidyr::gather("Subgroup", "Level", -Item) 
  mode_plot <- est_alloc_prob %>%
    pivot_longer(cols = -Item, names_to = "Subgroup", values_to = "Value")
  mode_plot %>% ggplot2::ggplot(ggplot2::aes(x=factor(Subgroup, 
                                                      levels = subgroup_labels),
                                             y=factor(Item, 
                                                      levels = rev(item_labels)),
                                             fill = Value)) + 
    ggplot2::geom_tile(color="black", linewidth = 0.1) + 
    scale_fill_gradient2(low = "#D95F02", mid = "white", high = "#1B9E77", 
                         midpoint = 0.5) + 
    geom_text(aes(label = round(Value, 2)), color = "black", size = 2.5) +  # Adds values
    ggplot2::labs(x = subgroup_title, y = item_title) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete() + 
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "right")
}

# Plot local profiles but gray out tiles that have global allocation probability
# nu_med >= 0.5
plot_wrpc_local_profiles_allocation <- function(res, item_labels = NULL, 
                                             item_title = "Item",
                                             categ_labels = NULL, 
                                             categ_title = "Consumption Level",
                                             subgroup_labels = NULL, 
                                             subgroup_title = "Subgroup", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca", "wrpc"))) {
    stop("res must be an object of class `swolca`, `wolca` or `wrpc`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  est_item_probs <- res$estimates$theta_local_med
  
  # Obtain theta modes as the category with highest probability for each item
  mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
  
  # Define item, subgroup, and category names if not specified
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  } else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ",
                res$data_vars$J))
  }
  H <- dim(mode_item_probs)[2]
  if (is.null(subgroup_labels)) {
    subgroup_labels <- 1:H
  } else if (length(subgroup_labels) != H) {
    stop(paste0("length of subgroup_labels must equal the number of subgroups, H = ", H))
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  } else if (length(categ_labels) != res$data_vars$R) {
    stop(paste0("length of categ_labels must equal the number of exposure categories, R = ", 
                res$data_vars$R))
  }
  rownames(mode_item_probs) <- item_labels
  colnames(mode_item_probs) <- subgroup_labels
  mode_item_probs$Item <- rownames(mode_item_probs)
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Subgroup <- Item <- Level <- NULL 
  # Create plot data
  mode_plot <- mode_item_probs %>% 
    tidyr::gather("Subgroup", "Level", -Item)
  
  # Obtain nu allocation estimates
  est_alloc_prob <- as.data.frame(t(res$estimates$nu_med))
  # nu plot values
  colnames(est_alloc_prob) <- subgroup_labels
  nu_vals <- est_alloc_prob %>%
    mutate(Item = item_labels) %>%
    pivot_longer(cols = -Item, names_to = "Subgroup", values_to = "Value")
  
  # Add in nu values
  mode_plot <- mode_plot %>%
    left_join(nu_vals, by = join_by(Subgroup, Item))
  # Set global values to NA
  mode_plot <- mode_plot %>%
    mutate(Global_bool = (Value > 0.5))
  # mode_plot <- mode_plot %>%
  #   mutate(Value = ifelse(Value > 0.5, NA, Value),
  #          Level = ifelse(Value > 0.5, NA, Level))
  # alloc_bool <- as.matrix(est_alloc_prob > 0.5)
  # mode_item_probs[alloc_bool] <- NA
  
  mode_plot %>% ggplot2::ggplot(ggplot2::aes(x=factor(Subgroup, 
                                                      levels = subgroup_labels),
                                             y=factor(Item, 
                                                      levels = rev(item_labels)), 
                                             fill = factor(Level))) + 
    # ggplot2::geom_tile(color="black", linewidth = 0.1, 
    #                    aes(fill = factor(ifelse(Global_bool, NA, Level),
    #                                      levels = categ_labels, 
    #                                      labels = categ_labels))) + 
    ggplot2::geom_tile(color="black", linewidth = 0.1) + 
    ggplot2::scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                               name = categ_title, labels = categ_labels,
                               na.value = "gray80") +
    ggplot2::geom_tile(data = subset(mode_plot, Global_bool), fill = "gray80",
                       color="black", linewidth = 0.1) +
    geom_text(data = subset(mode_plot, !Global_bool),
              aes(label = round(Value, 2)), color = "black", size = 2.5) + 
    # geom_text(aes(label = round(Value, 2)), color = "black", size = 2.5,
    #           na.rm = TRUE) + 
    ggplot2::labs(x = subgroup_title, y = item_title) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete() + 
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "top")
}


# normalize: normalize within each subgroup
plot_class_subgroup_dist <- function(res,  
                                     subgroup_labels = NULL, 
                                     subgroup_title = "Subgroup",
                                     class_labels = NULL,
                                     class_title = "Class", 
                                     normalize = TRUE) {
  
  grid <- table(res$estimates$c_all, res$data_vars$h_all)
  if (normalize) {
    grid <- prop.table(grid, margin = 2)
  }
  if (is.null(class_labels)) {
    class_labels <- 1:res$estimates$K_red
  }
  if (is.null(subgroup_labels)) {
    subgroup_labels <- 1:res$data_vars$H
  }
  grid_plot <- as.data.frame(grid)
  colnames(grid_plot) <- c(class_title, subgroup_title, "Frequency")
  
  
  grid_plot %>% ggplot(aes(x = Subgroup, y = Frequency, 
                           group = factor(Class, levels = rev(class_labels)),
                           fill = factor(Class, levels = class_labels))) + 
    geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(x = subgroup_title, y = "Frequency", fill = class_title) +
    coord_flip() + 
    ggplot2::theme_classic() +
    # guides(fill = guide_legend(
    #   title = class_title,
    #   labels = class_labels
    # )) + 
    # ggplot2::scale_x_discrete() + 
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "right")
}


plot_conf_cmd_dist <- function(data, normalize = TRUE) {
  
  grid <- table(data$conf_risk, data$cmd_risk)
  if (normalize) {
    grid <- prop.table(grid, margin = 2)
  }
  conf_labels <- 1:length(unique(data$conf_risk))
  cmd_labels <- 1:length(unique(data$cmd_risk))
  grid_plot <- as.data.frame(grid)
  colnames(grid_plot) <- c("Non-CMD Group", "CMD Group", "Frequency")
  grid_plot %>% ggplot(aes(x = `CMD Group`, y = Frequency, 
                           group = factor(`Non-CMD Group`, levels = rev(conf_labels)),
                           fill = factor(`Non-CMD Group`, levels = rev(conf_labels)))) + 
    geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(x = "CMD Group", y = "Frequency", fill = "Non-CMD Group") +
    coord_flip() + 
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "right")
}



plot_global_pattern_probs <- function(res, item_labels = NULL, categ_labels = NULL, 
                               categ_title = "Consumption Level",
                               class_labels = NULL, num_rows = 4,
                               class_title = "Dietary Pattern", 
                               y_title = "Consumption Level Probability", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca", "wrpc"))) {
    stop("res must be an object of class `swolca`, `wolca` or `wrpc`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  est_item_probs <- res$estimates$theta_global_med

  # Get number of latent classes
  K <- dim(est_item_probs)[2]
  
  # Define item, latent class, and category names if not specified
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  } else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ",
                res$data_vars$J))
  }
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", K))
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  } else if (length(categ_labels) != res$data_vars$R) {
    stop(paste0("length of categ_labels must equal the number of exposure categories, R = ", 
                res$data_vars$R))
  }
  
  dimnames(est_item_probs)[[1]] <- item_labels
  dimnames(est_item_probs)[[2]] <- class_labels
  
  # Convert to dataframe with each row corresponding to a value in the array
  # Use base R instead of reshape2 to reduce number of package imports
  # theta_plot <- reshape2::melt(est_item_probs, level = 2)
  theta_plot <- data.frame(expand.grid(lapply(dim(est_item_probs), seq_len)), 
                           value = as.vector(est_item_probs))
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Item <- Class <- Probability <- Level <- NULL
  
  # Create plot
  colnames(theta_plot) <- c("Item", "Class", "Level", "Probability")
  theta_plot %>%
    ggplot2::ggplot(ggplot2::aes(x = factor(Class, labels = class_labels), 
                                 y = Probability,
                                 fill = factor(Level))) + 
    ggplot2::geom_bar(stat = "identity", position = "stack") + 
    ggplot2::facet_wrap(factor(Item, labels = item_labels) ~ ., nrow = num_rows) + 
    ggplot2::scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                               name = categ_title, labels = categ_labels) +
    ggplot2::theme_bw() + 
    ggplot2::labs(x = class_title, y = y_title) + 
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 11, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "top",
                   strip.text = ggplot2::element_text(size = 9),
                   strip.background = ggplot2::element_rect(fill = "gray90"))
}

