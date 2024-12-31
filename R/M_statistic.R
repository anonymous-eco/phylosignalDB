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
#'
#' @param trait_dist A distance object of class `matrix` or [dist].
#'    Its row and column names should match the tip labels of the phylogenetic tree (`phy`).
#'    The functions [gower_dist()] and [cluster::daisy()] can be used to calculate distances using trait data.
#' @param phy A phylogenetic tree of class [phylo].
#' @param auto_multi2di A logical switch, `TRUE` or `FALSE`. Default is `TRUE`,
#'    then function [multi2di()] in `ape` package will be called to make the phylogeney (tree)
#'    be dichotomous if the tree (`phy`) contains some polytomies.
#' @returns A value that lies between 0 and 1, inclusive.
#' @references
#'    Blomberg, S.P. & Garland, T., Jr (2002) Tempo and mode in evolution: phylogenetic inertia, adaptation and comparative methods. Journal of Evolutionary Biology, 15(6): 899-910.
#'
#'
#' @seealso [M_rand_perm()] [phylosignal_M()]
#' @examples
#' data("turtles")
#' # Continuous trait
#' trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df)
#' M_stat(trait_dist, turtles$phylo)
#'
#' # Nominal discrete trait
#' trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
#' M_stat(trait_dist, turtles$phylo)
#'
#' # Ordinal discrete trait
#' trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
#' M_stat(trait_dist, turtles$phylo)
#'
#' # Multi-trait Combinations
#' trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
#'                        row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
#' M_stat(trait_dist, turtles$phylo)
#'
#' @export
M_stat <- function(trait_dist = NULL, phy = NULL, auto_multi2di = TRUE){
  if (length(phy)*length(trait_dist) == 0) {
    stop("The 'phy' and 'trait_dist' cannot be NULL.")
  }
  if (!("phylo" %in% class(phy))) {
    stop("The 'phy' should be a phylogenetic tree of class 'phylo'.")
  }
  if (!ape::is.rooted(phy)) {
    stop("The 'phy' tree must be rooted.")
  }
  if (!ape::is.ultrametric(phy)) {
    stop("The 'phy' tree is not ultrametric such that all tips are at equal distance from the root node.")
  }
  if (!ape::is.binary(phy)) {
    if (auto_multi2di) {
      warning("The 'phy' tree contains some polytomies. Function multi2di() in ape package have been called to make phylogeney(tree) be dichotomous.", call. = FALSE)
      phy <- ape::multi2di(phy)
    } else {
      stop("The 'phy' tree contains some polytomies. Function multi2di() in ape package maybe helpful. You also can set 'auto_multi2di' to be TRUE.")
    }
  }
  # Deal with the NA values in trait_dist/data_dist
  data_dist <- as.matrix(trait_dist) # Create a copy of trait_dist named data_dist.
  diag(data_dist) <- 0 # The diagonal values are equal to zero.
  row_col_narm <- (colSums(data_dist, na.rm = TRUE) > 0)
  data_dist <- data_dist[row_col_narm, row_col_narm]
  # Prune the phylogenetic tree
  phy <- castor::get_subtree_with_tips(tree = phy, only_tips = rownames(data_dist),
                                       collapse_monofurcations = TRUE,
                                       force_keep_root = FALSE)[["subtree"]]
  phy$node.label <- as.character(sort(unique(phy$edge[,1])))
  label_vec <- phy$tip.label
  S <- length(phy$tip.label)
  # Sort data_dist/trait_dist by the naming and numbering according to phy's tips.
  data_dist <- data_dist[label_vec, label_vec]
  # Get lists of internal nodes and tips about child nodes and descendant tips.
  edge_vec <- as.character(phy$edge[,2])
  names(edge_vec) <- as.character(phy$edge[,1])
  child_nodes <- split(edge_vec, names(edge_vec))
  descendant_tips_internal_nodes <- castor::get_subtrees_at_nodes(tree = phy,
                                                                  nodes = phy$node.label)[["new2old_tips"]]
  names(descendant_tips_internal_nodes) <- phy$node.label
  tips_vec <- 1:(ape::Ntip(phy))
  names(tips_vec) <- as.character(tips_vec)
  descendant_tips <- c(descendant_tips_internal_nodes, as.list(tips_vec))

  # Filter out the internal nodes that have at least 3 tips as descendants.
  rootnodes_subtree <- names(descendant_tips)[unlist(lapply(descendant_tips, function(x) ifelse(length(x) > 2, TRUE, FALSE)))]
  # subtree score
  N_subtree <- length(rootnodes_subtree)
  score_subtree <- rep(NA, N_subtree)
  names(score_subtree) <- rootnodes_subtree
  for (i in 1:N_subtree) {
    node_two <- child_nodes[[rootnodes_subtree[i]]]
    if (length(descendant_tips[[node_two[1]]]) > 1) {
      dist_child1 <- data_dist[descendant_tips[[node_two[1]]],descendant_tips[[node_two[1]]]]
      dist_child1 <- dist_child1[upper.tri(dist_child1)]
      average_child1 <- mean(dist_child1, na.rm = TRUE)
    } else {average_child1 <- 0}
    if (length(descendant_tips[[node_two[2]]]) > 1) {
      dist_child2 <- data_dist[descendant_tips[[node_two[2]]],descendant_tips[[node_two[2]]]]
      dist_child2 <- dist_child2[upper.tri(dist_child2)]
      average_child2 <- mean(dist_child2, na.rm = TRUE)
    } else {average_child2 <- 0}
    dist_mix <- data_dist[descendant_tips[[node_two[1]]], descendant_tips[[node_two[2]]]]
    average_mix <- mean(dist_mix, na.rm = TRUE)

    score_subtree[i] <- ifelse(average_mix >= average_child1, 0.5, 0) +
      ifelse(average_mix >= average_child2, 0.5, 0)
  }
  M_obs <- mean(score_subtree, na.rm = TRUE)
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
#'
#' @param trait_dist A distance object of class `matrix` or [dist].
#'    Its row and column names should match the tip labels of the phylogenetic tree (`phy`).
#'    The functions [gower_dist()] and [cluster::daisy()] can be used to calculate distances using trait data.
#' @param phy A phylogenetic tree of class [phylo].
#' @param reps An integer. The number of random permutations.
#' @param auto_multi2di A logical switch, `TRUE` or `FALSE`. Default is `TRUE`,
#'    then function [multi2di()] in `ape` package will be called to make the phylogeney (tree)
#'    be dichotomous if the tree (`phy`) contains some polytomies.
#' @param cores Number of cores to be used in parallel processing.
#' Default is 1, indicating no parallel computation is performed.
#' If set to 0, parallel computation is executed using `parallel::detectCores() - 1` number of cores.
#' @returns A list object containing two components.
#'    Component `$permuted` is the vector of M values obtained after random permutation for `reps` times;
#'    component `$observed` is the value of M statistic obtained from the original input data.
#' @references
#'    Blomberg, S.P. & Garland, T., Jr (2002) Tempo and mode in evolution: phylogenetic inertia, adaptation and comparative methods. Journal of Evolutionary Biology, 15(6): 899-910.
#'
#'
#' @seealso [M_stat()] [phylosignal_M()]
#' @examples
#' data("turtles")
#' # Continuous trait
#' trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df)
#' M_rand_perm(trait_dist, turtles$phylo, reps = 99) # reps=999 better
#'
#' # Nominal discrete trait
#' trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
#' M_rand_perm(trait_dist, turtles$phylo, reps = 99) # reps=999 better
#'
#' # Ordinal discrete trait
#' trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
#' M_rand_perm(trait_dist, turtles$phylo, reps = 99) # reps=999 better
#'
#' # Multi-trait Combinations
#' trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
#'                        row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
#' M_rand_perm(trait_dist, turtles$phylo, reps = 99) # reps=999 better
#'
#' @export
M_rand_perm <- function(trait_dist = NULL, phy = NULL, reps = 999, auto_multi2di = TRUE,
                        cores = 1){
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
  if (FALSE %in% (cores %in% c(1, 0))) {
    stop("The 'cores' must be 1 (no parallel computation) or 0 (parallel computation).")
  }
  if (!("phylo" %in% class(phy))) {
    stop("The 'phy' should be a phylogenetic tree of class 'phylo'.")
  }
  if (!ape::is.rooted(phy)) {
    stop("The 'phy' tree must be rooted.")
  }
  if (!ape::is.ultrametric(phy)) {
    stop("The 'phy' tree is not ultrametric such that all tips are at equal distance from the root node.")
  }
  if (!ape::is.binary(phy)) {
    if (auto_multi2di) {
      warning("The 'phy' tree contains some polytomies. Function multi2di() in ape package have been called to make phylogeney(tree) be dichotomous.", call. = FALSE)
      phy <- ape::multi2di(phy)
    } else {
      stop("The 'phy' tree contains some polytomies. Function multi2di() in ape package maybe helpful. You also can set 'auto_multi2di' to be TRUE.")
    }
  }
  # Deal with the NA values in trait_dist/data_dist
  data_dist <- as.matrix(trait_dist) # Create a copy of trait_dist named data_dist.
  diag(data_dist) <- 0 # The diagonal values are equal to zero.
  row_col_narm <- (colSums(data_dist, na.rm = TRUE) > 0)
  data_dist <- data_dist[row_col_narm, row_col_narm]
  # Prune the phylogenetic tree
  phy <- castor::get_subtree_with_tips(tree = phy, only_tips = rownames(data_dist),
                                       collapse_monofurcations = TRUE,
                                       force_keep_root = FALSE)[["subtree"]]
  phy$node.label <- as.character(sort(unique(phy$edge[,1])))
  label_vec <- phy$tip.label
  S <- ape::Ntip(phy)
  # Sort data_dist/trait_dist by the naming and numbering according to phy's tips.
  data_dist <- data_dist[label_vec, label_vec]

  # Get lists of internal nodes and tips about child nodes and descendant tips.
  edge_vec <- as.character(phy$edge[,2])
  names(edge_vec) <- as.character(phy$edge[,1])
  child_nodes <- split(edge_vec, names(edge_vec))
  descendant_tips_internal_nodes <- castor::get_subtrees_at_nodes(tree = phy,
                                                                  nodes = phy$node.label)[["new2old_tips"]]
  names(descendant_tips_internal_nodes) <- phy$node.label
  tips_vec <- 1:S
  names(tips_vec) <- as.character(tips_vec)
  descendant_tips <- c(descendant_tips_internal_nodes, as.list(tips_vec))

  num_descendant_tips <- unlist(lapply(descendant_tips, function(x) length(x)))
  # Filter out the internal nodes that have at least 3 tips as descendants.
  rootnodes_subtree <- names(descendant_tips)[unlist(lapply(descendant_tips, function(x) ifelse(length(x) > 2, TRUE, FALSE)))]

  num_tips_sampled <- round(mean(num_descendant_tips[rootnodes_subtree]))

  # subtree score
  N_subtree <- length(rootnodes_subtree)
  score_subtree <- rep(NA, N_subtree)
  names(score_subtree) <- rootnodes_subtree
  for (i in 1:N_subtree) {
    node_two <- child_nodes[[rootnodes_subtree[i]]]
    if (length(descendant_tips[[node_two[1]]]) > 1) {
      dist_child1 <- data_dist[descendant_tips[[node_two[1]]],descendant_tips[[node_two[1]]]]
      dist_child1 <- dist_child1[upper.tri(dist_child1)]
      average_child1 <- mean(dist_child1, na.rm = TRUE)
    } else {average_child1 <- 0}
    if (length(descendant_tips[[node_two[2]]]) > 1) {
      dist_child2 <- data_dist[descendant_tips[[node_two[2]]],descendant_tips[[node_two[2]]]]
      dist_child2 <- dist_child2[upper.tri(dist_child2)]
      average_child2 <- mean(dist_child2, na.rm = TRUE)
    } else {average_child2 <- 0}
    dist_mix <- data_dist[descendant_tips[[node_two[1]]], descendant_tips[[node_two[2]]]]
    average_mix <- mean(dist_mix, na.rm = TRUE)

    score_subtree[i] <- ifelse(average_mix >= average_child1, 0.5, 0) +
      ifelse(average_mix >= average_child2, 0.5, 0)
  }
  M_obs <- mean(score_subtree, na.rm = TRUE)
  # Calculate the simulated M values for random permutations, reps times.
  if (cores == 1) { # no parallel computation
    M_perm <- rep(NA, reps)
    for (j in 1:reps) {
      each_perm <- sample(1:S, size = S)
      data_dist_perm <- data_dist[each_perm, each_perm]

      score_subtree <- rep(NA, N_subtree)
      names(score_subtree) <- rootnodes_subtree
      for (i in 1:N_subtree) {
        node_two <- child_nodes[[rootnodes_subtree[i]]]
        if (length(descendant_tips[[node_two[1]]]) > 1) {
          dist_child1 <- data_dist_perm[descendant_tips[[node_two[1]]],descendant_tips[[node_two[1]]]]
          dist_child1 <- dist_child1[upper.tri(dist_child1)]
          average_child1 <- mean(dist_child1, na.rm = TRUE)
        } else {average_child1 <- 0}
        if (length(descendant_tips[[node_two[2]]]) > 1) {
          dist_child2 <- data_dist_perm[descendant_tips[[node_two[2]]],descendant_tips[[node_two[2]]]]
          dist_child2 <- dist_child2[upper.tri(dist_child2)]
          average_child2 <- mean(dist_child2, na.rm = TRUE)
        } else {average_child2 <- 0}
        dist_mix <- data_dist_perm[descendant_tips[[node_two[1]]], descendant_tips[[node_two[2]]]]
        average_mix <- mean(dist_mix, na.rm = TRUE)
        score_subtree[i] <- ifelse(average_mix >= average_child1, 0.5, 0) +
          ifelse(average_mix >= average_child2, 0.5, 0)
      }
      M_perm[j] <- mean(score_subtree, na.rm = TRUE)
    }
  }
  if (cores == 0) { # parallel computation
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    k <- NULL
    M_perm <- foreach::foreach(k = 1:reps, .combine = "c") %dopar% {
      each_perm <- sample(1:S, size = S)
      data_dist_perm <- data_dist[each_perm, each_perm]
      score_subtree <- rep(NA, N_subtree)
      names(score_subtree) <- rootnodes_subtree
      i <- NULL
      for (i in 1:N_subtree) {
        node_two <- child_nodes[[rootnodes_subtree[i]]]
        if (length(descendant_tips[[node_two[1]]]) > 1) {
          dist_child1 <- data_dist_perm[descendant_tips[[node_two[1]]],descendant_tips[[node_two[1]]]]
          dist_child1 <- dist_child1[upper.tri(dist_child1)]
          average_child1 <- mean(dist_child1, na.rm = TRUE)
        } else {average_child1 <- 0}
        if (length(descendant_tips[[node_two[2]]]) > 1) {
          dist_child2 <- data_dist_perm[descendant_tips[[node_two[2]]],descendant_tips[[node_two[2]]]]
          dist_child2 <- dist_child2[upper.tri(dist_child2)]
          average_child2 <- mean(dist_child2, na.rm = TRUE)
        } else {average_child2 <- 0}
        dist_mix <- data_dist_perm[descendant_tips[[node_two[1]]], descendant_tips[[node_two[2]]]]
        average_mix <- mean(dist_mix, na.rm = TRUE)
        score_subtree[i] <- ifelse(average_mix >= average_child1, 0.5, 0) +
          ifelse(average_mix >= average_child2, 0.5, 0)
      }
      mean(score_subtree, na.rm = TRUE)
    }
    parallel::stopCluster(cl)
  }
  return(list(M_permuted = M_perm,
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
#'
#' @param trait_dist A distance object of class `matrix` or [dist].
#'    Its row and column names should match the tip labels of the phylogenetic tree (`phy`).
#'    The functions [gower_dist()] and [cluster::daisy()] can be used to calculate distances using trait data.
#' @param phy A phylogenetic tree of class [phylo].
#' @param reps An integer. The number of random permutations.
#' @param auto_multi2di A logical switch, `TRUE` or `FALSE`. Default is `TRUE`,
#'    then function [multi2di()] in `ape` package will be called to make the phylogeney (tree)
#'    be dichotomous if the tree (`phy`) contains some polytomies.
#' @param output_M_permuted A logical switch, `TRUE` or `FALSE`. Default is `FALSE`.
#'    If this logical switch is set to `TRUE`, the returned list will include the vector
#'    of M values obtained after random permutations.
#' @param cores Number of cores to be used in parallel processing.
#'    Default is 1, indicating no parallel computation is performed.
#'    If set to 0, parallel computation is executed using `parallel::detectCores() - 1` number of cores.
#' @returns A list object containing two components.
#'    Component `$permuted` is the vector of M values obtained after random permutation for `reps` times;
#'    component `$observed` is the value of M statistic obtained from the original input data.
#' @references
#'    Blomberg, S.P. & Garland, T., Jr (2002) Tempo and mode in evolution: phylogenetic inertia, adaptation and comparative methods. Journal of Evolutionary Biology, 15(6): 899-910.
#'
#'
#' @seealso [M_stat()] [M_rand_perm()]
#' @examples
#' data("turtles")
#' # Continuous trait
#' trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df)
#' phylosignal_M(trait_dist, turtles$phylo, reps = 99) # reps=999 better
#'
#' # Nominal discrete trait
#' trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
#' phylosignal_M(trait_dist, turtles$phylo, reps = 99) # reps=999 better
#'
#' # Ordinal discrete trait
#' trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
#' phylosignal_M(trait_dist, turtles$phylo, reps = 99) # reps=999 better
#'
#' # Multi-trait Combinations
#' trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
#'                        row.names = turtles$traits$specie)
#' trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
#' phylosignal_M(trait_dist, turtles$phylo, reps = 99) # reps=999 better
#'
#' @export
phylosignal_M <- function(trait_dist = NULL, phy = NULL, reps = 999,
                              auto_multi2di = TRUE, output_M_permuted = FALSE,
                              cores = 1) {
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
  if (FALSE %in% (cores %in% c(1, 0))) {
    stop("The 'cores' must be 1 (no parallel computation) or 0 (parallel computation).")
  }
  if (!("phylo" %in% class(phy))) {
    stop("The 'phy' should be a phylogenetic tree of class 'phylo'.")
  }
  if (!ape::is.rooted(phy)) {
    stop("The 'phy' tree must be rooted.")
  }
  if (!ape::is.ultrametric(phy)) {
    stop("The 'phy' tree is not ultrametric such that all tips are at equal distance from the root node.")
  }
  if (!ape::is.binary(phy)) {
    if (auto_multi2di) {
      warning("The 'phy' tree contains some polytomies. Function multi2di() in ape package have been called to make phylogeney(tree) be dichotomous.", call. = FALSE)
      phy <- ape::multi2di(phy)
    } else {
      stop("The 'phy' tree contains some polytomies. Function multi2di() in ape package maybe helpful. You also can set 'auto_multi2di' to be TRUE.")
    }
  }
  if (!is.logical(output_M_permuted)) {
    stop("The 'output_M_permuted' must be FALSE or TRUE.")
  }
  M_perm <- M_rand_perm(trait_dist = trait_dist, phy = phy,
                        reps = reps, auto_multi2di = auto_multi2di,
                        cores = cores)
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
