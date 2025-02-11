#===============================
# WRPC plotting functions
# ==============================

library(tidyverse)

#' Plot theta modal exposure categories for each global latent class
#'
#' @description
#' `plot_wrpc_global_pattern_profiles` plots a heatmap of the global latent 
#' class patterns, where the patterns are defined by the category with the 
#' highest probability (i.e., the model category) for each exposure item.
#'
#' @param res An object of class `"wrpc"`, resulting from a call to [wrpc()]. 
#' @param item_labels String vector of names for the exposure items. Jx1. If
#' `NULL` (default), numbers from 1 to J are used.
#' @param item_title String specifying the title for the exposure items. 
#' Default is `"Item"`.
#' @param categ_labels String vector of names for the item categories. Rx1. If
#' `NULL` (default), numbers from 1 to R are used.
#' @param categ_title String specifying the title for the item categories. 
#' Default is `"Consumption Level"`.
#' @param class_labels String vector of names for the latent classes. Kx1. 
#' If `NULL` (default), numbers from 1 to K are used, where K is the final 
#' determined number of latent classes.
#' @param class_title String specifying the title for the latent classes. 
#' Default is `"Dietary Pattern"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a heatmap of the latent class patterns.
#' 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_brewer labs 
#' theme_classic scale_x_discrete theme element_text
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
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


#' Plot theta modal exposure categories for each local latent class
#'
#' @description
#' `plot_wrpc_local_pattern_profiles` plots a heatmap of the local latent 
#' class patterns, where the patterns are defined by the category with the 
#' highest probability (i.e., the model category) for each exposure item.
#'
#' @param res An object of class `"wrpc"`, resulting from a call to [wrpc()]. 
#' @param item_labels String vector of names for the exposure items. Jx1. If
#' `NULL` (default), numbers from 1 to J are used.
#' @param item_title String specifying the title for the exposure items. 
#' Default is `"Item"`.
#' @param categ_labels String vector of names for the item categories. Rx1. If
#' `NULL` (default), numbers from 1 to R are used.
#' @param categ_title String specifying the title for the item categories. 
#' Default is `"Consumption Level"`.
#' @param subgroup_labels String vector of names for the subpopulations. Hx1. 
#' If `NULL` (default), numbers from 1 to H are used, where H is the final 
#' determined number of subpopulations.
#' @param subgroup_title String specifying the title for the subpopulations. 
#' Default is `"Subgroup"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a heatmap of the latent class patterns.
#' 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_brewer labs 
#' theme_classic scale_x_discrete theme element_text
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
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


#' Plot nu global-local allocation probabilities for each subpopulation and item
#'
#' @description
#' `plot_wrpc_allocation` plots a heatmap of the global-local allocation 
#' probabilities for each subpopulation and item.
#'
#' @param res An object of class `"wrpc"`, resulting from a call to [wrpc()]. 
#' @param item_labels String vector of names for the exposure items. Jx1. If
#' `NULL` (default), numbers from 1 to J are used.
#' @param item_title String specifying the title for the exposure items. 
#' Default is `"Item"`.
#' @param subgroup_labels String vector of names for the subpopulations. Hx1. 
#' If `NULL` (default), numbers from 1 to H are used, where H is the final 
#' determined number of subpopulations.
#' @param subgroup_title String specifying the title for the subpopulations. 
#' Default is `"Subgroup"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a heatmap of the latent class patterns.
#' 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_brewer labs 
#' theme_classic scale_x_discrete theme element_text
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
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


#' Plot theta modal exposure categories for local deviations.
#'
#' @description
#' `plot_wrpc_local_pattern_profiles` plots a heatmap of the local latent 
#' class patterns, specifically for subpopulation-item combinations that have 
#' local allocation (i.e., \eqn{\nu} < 0.5). Items that have global allocation 
#' are grayed out. Pattern values are defined by the category with the highest 
#' probability (i.e., the model category) for each exposure item.
#'
#' @param res An object of class `"wrpc"`, resulting from a call to [wrpc()]. 
#' @param item_labels String vector of names for the exposure items. Jx1. If
#' `NULL` (default), numbers from 1 to J are used.
#' @param item_title String specifying the title for the exposure items. 
#' Default is `"Item"`.
#' @param categ_labels String vector of names for the item categories. Rx1. If
#' `NULL` (default), numbers from 1 to R are used.
#' @param categ_title String specifying the title for the item categories. 
#' Default is `"Consumption Level"`.
#' @param subgroup_labels String vector of names for the subpopulations. Hx1. 
#' If `NULL` (default), numbers from 1 to H are used, where H is the final 
#' determined number of subpopulations.
#' @param subgroup_title String specifying the title for the subpopulations. 
#' Default is `"Subgroup"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a heatmap of the latent class patterns.
#' 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_brewer labs 
#' theme_classic scale_x_discrete theme element_text
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
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
    mutate(Item = rownames(mode_item_probs)) %>%
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


#' Plot probabilities of exposure categories, theta, for each global latent class
#'
#' @description
#' `plot_wrpc_global_pattern_probs` plots a grouped barplot of the probability 
#' of the exposure categories, for each exposure item and each global latent class.
#' 
#' @inheritParams plot_global_pattern_profiles
#' @param y_title String specifying the title for the y-axis. Default is 
#' `"Consumption Level Probability"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a grouped barplot of the probability of 
#' the exposure categories, for each exposure item and each global latent class
#' 
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap scale_fill_brewer labs 
#' theme_bw theme element_text element_blank element_rect
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
plot_wrpc_global_pattern_probs <- function(res, item_labels = NULL, categ_labels = NULL, 
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



#' Plot the distribution of global classes by subpopulations
#' 
#' @param normalize Boolean indicating whether to normalize to sum to 1 within 
#' each subpopulation. Default is `TRUE`.
#' @param weights Boolean indicating if the proportions should incorporate 
#' survey weights. Default is `TRUE`.
#' 
#' @return
#' Returns a `ggplot2` object displaying a horizontal stacked bar plot of the 
#' distribution of global classes for each subpopulation.
#' 
plot_wrpc_class_subgroup_dist <- function(res, weights = TRUE, 
                                         subgroup_labels = NULL, 
                                         subgroup_title = "Subgroup",
                                         class_labels = NULL,
                                         class_title = "Class", 
                                         normalize = TRUE) {
  if (weights) {
    svy_data <- data.frame(sampling_wt = res$data_vars$sampling_wt, 
                           subgroup = as.factor(res$data_vars$h_all),
                           class = as.factor(res$estimates$c_all))
    svy_design_wts_only <- survey::svydesign(id = ~1, 
                                             weights = ~sampling_wt, 
                                             data = svy_data)
    if (normalize) {
      # row proportions
      temp <- as.data.frame(survey::svyby(~class, ~subgroup, svy_design_wts_only, 
                                          svymean, na.rm = TRUE))
    } else {
      # totals
      temp <- as.data.frame(survey::svyby(~class, ~subgroup, svy_design_wts_only, 
                                          svytotal, na.rm = TRUE))
    }
    
    grid <- temp %>% 
      select(subgroup:class4) %>%
      pivot_longer(cols = class1:class4, names_to = "class_v1", 
                                         values_to = "Freq") %>%
      select(class_v1, subgroup, Freq)
    grid$class <- stringr::str_remove(grid$class_v1, "class")
    grid <- grid %>% select(class, subgroup, Freq)
    
    # knitr::kable(grid, digits = 3, booktabs = TRUE)
    # colSums(grid)
  } else {
    grid <- table(res$estimates$c_all, res$data_vars$h_all)
    if (normalize) {
      grid <- prop.table(grid, margin = 2)
    }
  }
  
  if (is.null(class_labels)) {
    class_labels <- 1:res$estimates$K_red
  }
  if (is.null(subgroup_labels)) {
    subgroup_labels <- 1:res$data_vars$H
  }
  grid_plot <- as.data.frame(grid)
  colnames(grid_plot) <- c(class_title, subgroup_title, "Frequency")
  
  
  grid_plot %>% ggplot(aes(x = factor(Subgroup, levels = rev(1:res$data_vars$H), 
                                      labels = rev(subgroup_labels)), y = Frequency, 
                           group = factor(Class, levels = rev(class_labels)),
                           fill = factor(Class, levels = class_labels))) + 
    geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(x = subgroup_title, y = "Frequency", fill = class_title) +
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

