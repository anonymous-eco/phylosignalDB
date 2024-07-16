#' Calculate Gower distance
#'
#' `gower_dist()` calculates Gower distance among observations or species.
#'
#' @param x A data frame. The columns usually represent trait data, and the row names are species names.
#' @param type A list for specifying the variable types of the columns in `x`.
#' Default is numeric type. More details in [cluster::daisy()].
#' @param dist_format The class of the return value. Default is "matrix".
#' @returns A matrix or dist object containing the Gower distance among the rows of `x`.
#' @references
#'    Gower, J.C. (1971) A general coefficient of similarity and some of its properties. Biometrics: 857-871.
#'
#'    Kaufman, L. & Rousseeuw, P.J. (1990) Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.
#'
#' @seealso [cluster::daisy()] which this function wraps.
#' @examples
#' data("turtles")
#' # Continuous trait
#' trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df)
#'
#' # Nominal discrete trait
#' trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
#'
#' # Ordinal discrete trait
#' trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
#'
#' # Multi-trait Combinations
#' trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
#'                        row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
#'
#' @export
gower_dist <- function(x, type = list(),
                       dist_format = c("matrix", "dist")){
  if (dist_format[1] %in% c("matrix", "dist")) {
    dist_format <- dist_format[1]
  } else {
    stop("The 'value_format' must be one in c('matrix', 'dist').")
  }
  trait_dist <- cluster::daisy(x=x, metric="gower", type=type)
  if (dist_format == "matrix") {
    trait_dist <- as.matrix(trait_dist)
  }
  return(trait_dist)
}

#' Calculate M statistic
#'
#' `M_stat` calculates the value of M statistic as a measurement of the strength of
#' the phylogenetic signal for the trait(s). The trait(s) could be continuous, discrete, or multi-variable.
#' Blomberg and Garland (2002) provided a widely accepted statistical definition of
#' the phylogenetic signal, which is the "tendency for related species to resemble
#' each other more than they resemble species drawn at random from the tree".
#' The M statistic strictly adheres to the definition of phylogenetic signal,
#' formulating an index and developing a method of testing in strict accordance
#' with the definition, instead of relying on correlation analysis or evolutionary models.
#' The novel method equivalently expressed the textual definition of the phylogenetic signal
#' as an inequality equation of the phylogenetic and trait distances and constructed the M statistic.
#' For more details, please refer to Yao & Yuan (2024).
#'
#' @param trait_dist A distance object of class `matrix` or [dist].
#'    Its row and column names should match the tip labels of the phylogenetic tree (`phy`).
#'    The functions [gower_dist()] and [cluster::daisy()] can be used to calculate distances using trait data.
#' @param phy A phylogenetic tree of class [phylo], [phylo4] or [phylo4d].
#' @param auto_multi2di A logical switch, `TRUE` or `FALSE`. Default is `TRUE`,
#'    then function [multi2di()] in `ape` package will be called to make the phylogeney (tree)
#'    be dichotomous if the tree (`phy`) contains some polytomies.
#' @returns A value that lies between 0 and 1, inclusive.
#' @references
#'    Blomberg, S.P. & Garland, T., Jr (2002) Tempo and mode in evolution: phylogenetic inertia,
#'    adaptation and comparative methods. Journal of Evolutionary Biology, 15(6): 899-910.
#'    Yao, L. & Yuan, Y. (2024) A unified method for detecting phylogenetic signals in continuous traits,
#'    discrete traits, and multi-trait combinations. (manuscript)
#'
#' @seealso [M_rand_perm()] [phylosignal_M()]
#' @examples
#' data("turtles")
#' # Continuous trait
#' trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df)
#' M_stat(trait_dist, turtles$phylo4)
#'
#' # Nominal discrete trait
#' trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
#' M_stat(trait_dist, turtles$phylo4)
#'
#' # Ordinal discrete trait
#' trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
#' M_stat(trait_dist, turtles$phylo4)
#'
#' # Multi-trait Combinations
#' trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
#'                        row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
#' M_stat(trait_dist, turtles$phylo4)
#'
#' @export
M_stat <- function(trait_dist = NULL, phy = NULL, auto_multi2di = TRUE){
  if (length(phy)*length(trait_dist) == 0) {
    stop("The 'phy' and 'trait_dist' cannot be NULL.")
  }
  data_dist <- trait_dist # Create a copy of trait_dist named data_dist.
  if (TRUE %in% (c("phylo", "phylo4", "phylo4d") %in% class(phy))) {
    phy <- as(phy, "phylo4")
  } else {
    stop("The 'phy' should be a phylogenetic tree of class 'phylo', 'phylo4' or 'phylo4d'.")
  }
  if (!phylobase::isRooted(phy)) {
    stop("The 'phy' tree must be rooted.")
  }
  if (!phylobase::isUltrametric(phy)) {
    stop("The 'phy' tree is not ultrametric such that all tips are at equal distance from the root node.")
  }
  if (phylobase::hasPoly(phy)) {
    if (auto_multi2di) {
      warning("The 'phy' tree contains some polytomies. Function multi2di() in ape package have been called to make phylogeney(tree) be dichotomous.", call. = FALSE)
      phy <- as(ape::multi2di(as(phy, "phylo")), "phylo4")
    } else {
      stop("The 'phy' tree contains some polytomies. Function multi2di() in ape package maybe helpful. You also can set 'auto_multi2di' to be TRUE.")
    }
  }
  S <- phylobase::nTips(phy)
  label_vec <- phylobase::tipLabels(phy)
  # Sort data_dist/trait_dist by the naming and numbering according to phy's tips.
  data_dist <- as.matrix(data_dist)
  data_dist <- data_dist[label_vec, label_vec]
  # Create an index list of internal nodes.
  internal_nodes <- phylobase::getNode(phy, type = "internal")
  k <- NULL
  nodes_info <- foreach::foreach(k=1:length(internal_nodes)) %do% {
    temp_tips_id <- phylobase::descendants(phy, internal_nodes[k], type = "tips")
    temp_tips_num <- length(temp_tips_id)
    temp_tips_label <- names(temp_tips_id)
    temp_child_nodes <- phylobase::children(phy, internal_nodes[k])
    temp_child_tips_1left <- phylobase::descendants(phy, temp_child_nodes[1], type = "tips")
    temp_child_tips_2right <- phylobase::descendants(phy, temp_child_nodes[2], type = "tips")
    list(node_id = internal_nodes[k],
         desc_tips_id = temp_tips_id,
         desc_tips_num = temp_tips_num,
         desc_tips_label = temp_tips_label,
         child_nodes_id = temp_child_nodes,
         child_tipsid_1left = temp_child_tips_1left,
         child_tipsid_2right = temp_child_tips_2right)
  }
  names(nodes_info) <- as.character(internal_nodes)
  # Filter out the internal nodes that have at least 3 tips as descendants.
  kk <- NULL
  rootnodes_subtree <- foreach::foreach(kk=1:length(internal_nodes), .combine = "c") %do% {
    temp_node_info <- nodes_info[[as.character(internal_nodes[kk])]]
    if (temp_node_info$desc_tips_num > 2) {
      temp_node <- temp_node_info$node_id
    } else {
      temp_node <- NULL
    }
    temp_node
  }
  N_subtree <- length(rootnodes_subtree)
  # Assign equal weight to each subtree.
  w_vec <- rep(1/N_subtree, N_subtree)
  subtree_scores_0 <- rootnodes_subtree - rootnodes_subtree
  # Calculate the M value of the observed sample data.
  i <- NULL
  M_obs <- foreach::foreach(i=1:N_subtree, .combine="+") %do% {
    temp_node_info <- nodes_info[[as.character(rootnodes_subtree[i])]]
    dist_child1 <- as.vector(as.dist(data_dist[temp_node_info$child_tipsid_1left,
                                               temp_node_info$child_tipsid_1left]))
    dist_child2 <- as.vector(as.dist(data_dist[temp_node_info$child_tipsid_2right,
                                               temp_node_info$child_tipsid_2right]))
    dist_mix <- as.vector(data_dist[temp_node_info$child_tipsid_1left,
                                    temp_node_info$child_tipsid_2right])
    average_child1 <- ifelse(length(dist_child1) > 0, mean(dist_child1), 0)
    average_child2 <- ifelse(length(dist_child2) > 0, mean(dist_child2), 0)
    average_mix <- mean(dist_mix)
    subtree_score <- ifelse(average_mix >= average_child1, 0.5, 0) +
      ifelse(average_mix >= average_child2, 0.5, 0)
    subtree_score*w_vec[i]
  }
  names(M_obs) <- "M_stat"
  return(M_obs)
}

#' Calculate M statistics after random permutations
#'
#' `M_rand_perm` calculates M statistic for trait(s) after randomly permuting the species names or tip labels in phylogeny.
#' The M statistic is a unified method for detecting phylogenetic signals in continuous traits,
#' discrete traits, and multi-trait combinations.
#' Blomberg and Garland (2002) provided a widely accepted statistical definition of
#' the phylogenetic signal, which is the "tendency for related species to resemble
#' each other more than they resemble species drawn at random from the tree".
#' The M statistic strictly adheres to the definition of phylogenetic signal,
#' formulating an index and developing a method of testing in strict accordance
#' with the definition, instead of relying on correlation analysis or evolutionary models.
#' The novel method equivalently expressed the textual definition of the phylogenetic signal
#' as an inequality equation of the phylogenetic and trait distances and constructed the M statistic.
#' For more details, please refer to Yao & Yuan (2024).
#'
#' @param trait_dist A distance object of class `matrix` or [dist].
#'    Its row and column names should match the tip labels of the phylogenetic tree (`phy`).
#'    The functions [gower_dist()] and [cluster::daisy()] can be used to calculate distances using trait data.
#' @param phy A phylogenetic tree of class [phylo], [phylo4] or [phylo4d].
#' @param reps An integer. The number of random permutations.
#' @param auto_multi2di A logical switch, `TRUE` or `FALSE`. Default is `TRUE`,
#'    then function [multi2di()] in `ape` package will be called to make the phylogeney (tree)
#'    be dichotomous if the tree (`phy`) contains some polytomies.
#' @returns A list object containing two components.
#'    Component `$permuted` is the vector of M values obtained after random permutation for `reps` times;
#'    component `$observed` is the value of M statistic obtained from the original input data.
#' @references
#'    Blomberg, S.P. & Garland, T., Jr (2002) Tempo and mode in evolution: phylogenetic inertia,
#'    adaptation and comparative methods. Journal of Evolutionary Biology, 15(6): 899-910.
#'    Yao, L. & Yuan, Y. (2024) A unified method for detecting phylogenetic signals in continuous traits,
#'    discrete traits, and multi-trait combinations. (manuscript)
#'
#' @seealso [M_stat()] [phylosignal_M()]
#' @examples
#' data("turtles")
#' # Continuous trait
#' trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df)
#' M_rand_perm(trait_dist, turtles$phylo4, reps = 99) # reps=999 better
#'
#' # Nominal discrete trait
#' trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
#' M_rand_perm(trait_dist, turtles$phylo4, reps = 99) # reps=999 better
#'
#' # Ordinal discrete trait
#' trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
#' M_rand_perm(trait_dist, turtles$phylo4, reps = 99) # reps=999 better
#'
#' # Multi-trait Combinations
#' trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
#'                        row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
#' M_rand_perm(trait_dist, turtles$phylo4, reps = 99) # reps=999 better
#'
#' @export
M_rand_perm <- function(trait_dist = NULL, phy = NULL, reps = 999, auto_multi2di = TRUE){
  if (length(phy)*length(trait_dist) == 0) {
    stop("The 'phy' and 'trait_dist' cannot be NULL.")
  }
  data_dist <- trait_dist # Create a copy of trait_dist named data_dist.
  if ((!is.numeric(reps))) {
    stop("The 'reps' must be a positive integer.")
  } else if (reps < 1) {
    stop("The 'reps' must be a positive integer.")
  } else {
    reps <- round(reps) # coerce to a positive integer
  }
  if (TRUE %in% (c("phylo", "phylo4", "phylo4d") %in% class(phy))) {
    phy <- as(phy, "phylo4")
  } else {
    stop("The 'phy' should be a phylogenetic tree of class 'phylo', 'phylo4' or 'phylo4d'.")
  }
  if (!phylobase::isRooted(phy)) {
    stop("The 'phy' tree must be rooted.")
  }
  if (!phylobase::isUltrametric(phy)) {
    stop("The 'phy' tree is not ultrametric such that all tips are at equal distance from the root node.")
  }
  if (phylobase::hasPoly(phy)) {
    if (auto_multi2di) {
      warning("The 'phy' tree contains some polytomies. Function multi2di() in ape package have been called to make phylogeney(tree) be dichotomous.", call. = FALSE)
      phy <- as(ape::multi2di(as(phy, "phylo")), "phylo4")
    } else {
      stop("The 'phy' tree contains some polytomies. Function multi2di() in ape package maybe helpful. You also can set 'auto_multi2di' to be TRUE.")
    }
  }
  S <- phylobase::nTips(phy)
  label_vec <- phylobase::tipLabels(phy)
  # Sort data_dist/trait_dist by the naming and numbering according to phy's tips
  data_dist <- as.matrix(data_dist)
  data_dist <- data_dist[label_vec, label_vec]
  # Create an index list of internal nodes
  internal_nodes <- phylobase::getNode(phy, type = "internal")
  k <- NULL
  nodes_info <- foreach::foreach(k=1:length(internal_nodes)) %do% {
    temp_tips_id <- phylobase::descendants(phy, internal_nodes[k], type = "tips")
    temp_tips_num <- length(temp_tips_id)
    temp_tips_label <- names(temp_tips_id)
    temp_child_nodes <- phylobase::children(phy, internal_nodes[k])
    temp_child_tips_1left <- phylobase::descendants(phy, temp_child_nodes[1], type = "tips")
    temp_child_tips_2right <- phylobase::descendants(phy, temp_child_nodes[2], type = "tips")
    list(node_id = internal_nodes[k],
         desc_tips_id = temp_tips_id,
         desc_tips_num = temp_tips_num,
         desc_tips_label = temp_tips_label,
         child_nodes_id = temp_child_nodes,
         child_tipsid_1left = temp_child_tips_1left,
         child_tipsid_2right = temp_child_tips_2right)
  }
  names(nodes_info) <- as.character(internal_nodes)
  # Filter out the internal nodes that have at least 3 tips as descendants
  kk <- NULL
  rootnodes_subtree <- foreach::foreach(kk=1:length(internal_nodes), .combine = "c") %do% {
    temp_node_info <- nodes_info[[as.character(internal_nodes[kk])]]
    if (temp_node_info$desc_tips_num > 2) {
      temp_node <- temp_node_info$node_id
    } else {
      temp_node <- NULL
    }
    temp_node
  }
  N_subtree <- length(rootnodes_subtree)
  # Assign equal weight to each subtree.
  w_vec <- rep(1/N_subtree, N_subtree)
  subtree_scores_0 <- rootnodes_subtree - rootnodes_subtree
  # Calculate the M value of the observed sample data.
  i <- NULL
  M_obs <- foreach::foreach(i=1:N_subtree, .combine="+") %do% {
    temp_node_info <- nodes_info[[as.character(rootnodes_subtree[i])]]
    dist_child1 <- as.vector(as.dist(data_dist[temp_node_info$child_tipsid_1left,
                                               temp_node_info$child_tipsid_1left]))
    dist_child2 <- as.vector(as.dist(data_dist[temp_node_info$child_tipsid_2right,
                                               temp_node_info$child_tipsid_2right]))
    dist_mix <- as.vector(data_dist[temp_node_info$child_tipsid_1left,
                                    temp_node_info$child_tipsid_2right])
    average_child1 <- ifelse(length(dist_child1) > 0, mean(dist_child1), 0)
    average_child2 <- ifelse(length(dist_child2) > 0, mean(dist_child2), 0)
    average_mix <- mean(dist_mix)
    subtree_score <- ifelse(average_mix >= average_child1, 0.5, 0) +
      ifelse(average_mix >= average_child2, 0.5, 0)
    subtree_score*w_vec[i]
  }
  # Calculate the simulated M values for random permutations, reps times.
  j <- NULL
  M_reps <- foreach::foreach(j = 1:reps, .combine = "c") %do% {
    each_perm <- sample(1:S, size = S)
    data_dist_perm <- data_dist[each_perm, each_perm]
    kkk <- NULL
    M_perm <- foreach::foreach(kkk=1:N_subtree, .combine="+") %do% {
      temp_node_info <- nodes_info[[as.character(rootnodes_subtree[kkk])]]
      dist_child1 <- as.vector(as.dist(data_dist_perm[temp_node_info$child_tipsid_1left,
                                                      temp_node_info$child_tipsid_1left]))
      dist_child2 <- as.vector(as.dist(data_dist_perm[temp_node_info$child_tipsid_2right,
                                                      temp_node_info$child_tipsid_2right]))
      dist_mix <- as.vector(data_dist_perm[temp_node_info$child_tipsid_1left,
                                           temp_node_info$child_tipsid_2right])
      average_child1 <- ifelse(length(dist_child1) > 0, mean(dist_child1), 0)
      average_child2 <- ifelse(length(dist_child2) > 0, mean(dist_child2), 0)
      average_mix <- mean(dist_mix)
      subtree_score <- ifelse(average_mix >= average_child1, 0.5, 0) +
        ifelse(average_mix >= average_child2, 0.5, 0)
      subtree_score*w_vec[kkk]
    }
    M_perm
  }
  return(list(M_permuted = M_reps,
              M_observed = M_obs))
}

#' Measure and test phylogenetic signal with M statistic
#'
#' `phylosignal_M` computes the M statistic for trait(s) and evaluates
#' its statistical significance through a random permutation test.
#' The M statistic is a unified method for detecting phylogenetic signals in continuous traits,
#' discrete traits, and multi-trait combinations.
#' Blomberg and Garland (2002) provided a widely accepted statistical definition of
#' the phylogenetic signal, which is the "tendency for related species to resemble
#' each other more than they resemble species drawn at random from the tree".
#' The M statistic strictly adheres to the definition of phylogenetic signal,
#' formulating an index and developing a method of testing in strict accordance
#' with the definition, instead of relying on correlation analysis or evolutionary models.
#' The novel method equivalently expressed the textual definition of the phylogenetic signal
#' as an inequality equation of the phylogenetic and trait distances and constructed the M statistic.
#' For more details, please refer to Yao & Yuan (2024).
#'
#' @param trait_dist A distance object of class `matrix` or [dist].
#'    Its row and column names should match the tip labels of the phylogenetic tree (`phy`).
#'    The functions [gower_dist()] and [cluster::daisy()] can be used to calculate distances using trait data.
#' @param phy A phylogenetic tree of class [phylo], [phylo4] or [phylo4d].
#' @param reps An integer. The number of random permutations.
#' @param auto_multi2di A logical switch, `TRUE` or `FALSE`. Default is `TRUE`,
#'    then function [multi2di()] in `ape` package will be called to make the phylogeney (tree)
#'    be dichotomous if the tree (`phy`) contains some polytomies.
#' @param output_M_permuted A logical switch, `TRUE` or `FALSE`. Default is `FALSE`.
#'    If this logical switch is set to `TRUE`, the returned list will include the vector
#'    of M values obtained after random permutations.
#' @returns A list object containing two components.
#'    Component `$permuted` is the vector of M values obtained after random permutation for `reps` times;
#'    component `$observed` is the value of M statistic obtained from the original input data.
#' @references
#'    Blomberg, S.P. & Garland, T., Jr (2002) Tempo and mode in evolution: phylogenetic inertia,
#'    adaptation and comparative methods. Journal of Evolutionary Biology, 15(6): 899-910.
#'    Yao, L. & Yuan, Y. (2024) A unified method for detecting phylogenetic signals in continuous traits,
#'    discrete traits, and multi-trait combinations. (manuscript)
#'
#' @seealso [M_stat()] [M_rand_perm()]
#' @examples
#' data("turtles")
#' # Continuous trait
#' trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df)
#' phylosignal_M(trait_dist, turtles$phylo4, reps = 99) # reps=999 better
#'
#' # Nominal discrete trait
#' trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
#' phylosignal_M(trait_dist, turtles$phylo4, reps = 99) # reps=999 better
#'
#' # Ordinal discrete trait
#' trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
#' phylosignal_M(trait_dist, turtles$phylo4, reps = 99) # reps=999 better
#'
#' # Multi-trait Combinations
#' trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
#'                        row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
#' phylosignal_M(trait_dist, turtles$phylo4, reps = 99) # reps=999 better
#'
#' @export
phylosignal_M <- function(trait_dist = NULL, phy = NULL, reps = 999,
                          auto_multi2di = TRUE,
                          output_M_permuted = FALSE) {
  if (length(phy)*length(trait_dist) == 0) {
    stop("The 'phy' and 'trait_dist' cannot be NULL.")
  }
  if ((!is.numeric(reps))) {
    stop("The 'reps' must be a positive integer.")
  } else if (reps < 1) {
    stop("The 'reps' must be a positive integer.")
  } else {
    reps <- round(reps) # coerce to a positive integer
  }
  if (!is.logical(output_M_permuted)) {
    stop("The 'output_M_permuted' must be a logical value (TRUE or FALSE).")
  }
  if (TRUE %in% (c("phylo", "phylo4", "phylo4d") %in% class(phy))) {
    phy <- as(phy, "phylo4")
  } else {
    stop("The 'phy' should be a phylogenetic tree of class 'phylo', 'phylo4' or 'phylo4d'.")
  }
  if (!phylobase::isRooted(phy)) {
    stop("The 'phy' tree must be rooted.")
  }
  if (!phylobase::isUltrametric(phy)) {
    stop("The 'phy' tree is not ultrametric such that all tips are at equal distance from the root node.")
  }
  if (phylobase::hasPoly(phy)) {
    if (auto_multi2di) {
      warning("The 'phy' tree contains some polytomies. Function multi2di() in ape package have been called to make phylogeney(tree) be dichotomous.", call. = FALSE)
      phy <- as(ape::multi2di(as(phy, "phylo")), "phylo4")
    } else {
      stop("The 'phy' tree contains some polytomies. Function multi2di() in ape package maybe helpful. You also can set 'auto_multi2di' to be TRUE.")
    }
  }
  M_perm <- M_rand_perm(phy = phy, trait_dist = trait_dist,
                        reps = reps)
  M_observed <- M_perm$M_observed
  M_permuted <- M_perm$M_permuted
  # The rank() in ascending order handles ties more freely than the order().
  # "average": There is a tendency to exaggerate p-values when the number of tips/species is small.
  # rank_order <- rank(c(M_observed, M_permuted),
  #                    ties.method = "average", na.last = TRUE)
  # "max": Due to the right-tailed test, there is a tendency to reduce p-values when the number of tips/species is small.
  rank_order <- rank(c(M_observed, M_permuted),
                     ties.method = "max", na.last = NA)
  N <- length(rank_order)
  # Reverse order.
  pvalue <- (N+1-rank_order[1])/N
  if(output_M_permuted){
    result <- list(stat = M_observed, pvalue = pvalue, M_permuted = M_permuted)
  } else {
    result <- list(stat = M_observed, pvalue = pvalue)
  }
  return(result)
}






