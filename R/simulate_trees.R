library(TreeSimGM)
library(castor)
library(ape)
library(dplyr)

simulate_tree <- function(a, b, deathprob, rho, origin, seed = 1) {
  set.seed(seed)
  # generate phylogeny for a fixed time
  birth_dist = paste0("rgamma(shape=", b, ", scale=", a, ")") # birth events (branching times) are Gamma distributed
  death_dist = paste0("rexp(0)") # no death
  tree_complete = TreeSimGM::sim.age(age = origin, numbsim = 1, 
                                     waitsp = birth_dist, waitext = death_dist, 
                                     tiplabel = c("","","",""), 
                                     complete = FALSE, symmetric = TRUE)[[1]]
  # add death on branching times
  tree = prune_tree(tree_complete, deathprob)
  tree = ape::drop.fossil(tree)
  
  # sample tips
  tips = tree$tip.label
  sample_tip = sample(c(TRUE, FALSE), length(tips), prob = c(rho, 1 - rho), replace = TRUE)
  tree = ape::keep.tip(tree, tip = tips[sample_tip])
  # relabel
  tree$tip.label = as.character( c(0:(length(tree$tip.label) - 1)) ) 
 
  # add tree parameters
  tree$length = sum(tree$edge.length) + tree$root.edge
  tree$height = max(node.depth.edgelength(tree))
  return(tree)
}


# add death at branching times
prune_tree <- function(tree, deathprob) {
  
  # extract edge matrix, root and inner nodes
  edge = tree$edge
  edge.length = tree$edge.length
  
  # go along branches, let cells die with probability deathprob
  i = 1
  while (i <= nrow(edge)) {
    branch = edge[i, ]
    die = sample(c(TRUE, FALSE), 1, prob = c(deathprob, 1 - deathprob))
    if (die) {
      # shorten branch to time of cell death 
      edge.length[i] = runif(1, min = 0, max = edge.length[i])
      dying_node = branch[2]
      if (!dying_node %in% tree$tip.label) {
        # remove all offspring of dying cell 
        ix_prev = i + which(edge[-c(1:i), 1] < dying_node)[1]
        if (!is.na(ix_prev)) {
          edge = edge[-c((i + 1) : (ix_prev - 1)), ] 
          edge.length = edge.length[-c((i + 1) : (ix_prev - 1))] 
        } else {
          edge = edge[-c((i + 1) : nrow(edge)), ] 
          edge.length = edge.length[-c((i + 1) : length(edge.length))] 
        }
      }
    }
    i = i + 1
  }
  
  # # plot surviving nodes in color over whole tree
  # ix_in = which(apply(tree$edge, 1, function(x) all(x %in% edge)))
  # edge_color = rep("black", nrow(tree$edge))
  # edge_color[ix_in] = "blue"
  # plot(tree, edge.color = edge_color, show.tip.label = F)
  # axisPhylo()
  
  # restore phylo object 
  tree_pruned = restore_phylo(edge, edge.length)
  # add age and root age
  tree_pruned$age = tree$age
  tree_pruned$root.edge = tree$root.edge
  
  return(tree_pruned)
}


# helper function for restoring phylo object after cell death
restore_phylo <- function(edge, edge.length) {
  # find and rename tips
  ix_tips = which(!(edge[, 2] %in% edge[, 1]))
  tip.label = c(1:length(ix_tips)) 
  edge[ix_tips, 2] = tip.label
  
  # find and rename inner nodes
  innodes = unique(edge[, 1])
  Nnode = length(innodes)
  new_innodes = c((length(ix_tips) + 1):(length(ix_tips) + length(innodes)))
  innodes = setNames(innodes, new_innodes)
  edge[, 1] = sapply(edge[, 1], function(x) {x = as.integer(names(which(innodes == x))) }, simplify = T)
  if (!is.null(edge[-ix_tips, ]) & length(edge[-ix_tips, ]) != 0) {
    edge[-ix_tips, 2] = sapply(edge[-ix_tips, 2], function(x) {x = as.integer(names(which(innodes == x))) }, simplify = T)
  }
  
  # define updated phylo object
  phy = list(edge = edge, edge.length = edge.length, tip.label = as.character(tip.label), Nnode = Nnode)
  class(phy) = "phylo"
  attr(phy, "order") = "cladewise"
  return(phy)
}


# get all branching times 
get_tree_times <- function(tree) {
  # get branches
  df = data.frame(tree$edge)
  names(df) <- c("end.node", "start.node")
  
  # calculate distances to root
  dists = castor::get_all_distances_to_root(tree)
  dists = tree$age - (dists + tree$root.edge)
  dists = data.frame(end.time = zapsmall(dists), end.node = seq(1:length(dists)))
  
  # merge tables
  df = merge(dists, df, by = "end.node")
  dists = dists %>% dplyr::rename(start.node = end.node, start.time = end.time)
  df = merge(dists, df, by = "start.node")
  
  # add type
  ntips = length(tree$tip.label)
  df = df %>% dplyr::mutate(edge.type = case_when(
    start.node <= ntips ~ "External", TRUE ~ "Internal"))
  
  # add stem edge
  stem = data.frame(start.node = max(df$start.node) + 1, start.time = max(df$end.time),
                    end.node = max(df$start.node), end.time = tree$age, 
                    edge.type = "Internal")
  df =  rbind(df, stem)
  return(df)
}

