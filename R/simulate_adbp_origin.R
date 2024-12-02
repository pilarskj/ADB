library(ape)
library(treeio)
library(dplyr)
library(assertthat)
library(ggtree) # for plotting


#' Simulator of a phylogeny from an Age-Dependent Branching Process for a fixed time interval (since origin)
#' @param origin_time time of birth of the initial particle
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param Xsi_as matrix of asymetric type transition probabilities
#' @param Xsi_s matrix of symetric type transition probabilities
#' @param rho sampling probability 
#' @param rho minimum number of tips in the phylogeny
simulate_phylogeny <- function(origin_time, a, b, d = 0, rho = 1, origin_type = 0, Xsi_as = matrix(0), Xsi_s = matrix(1), min_tips = 2) {
  # assert that all inputs are correct
  ntypes = length(a)
  assert_that(all(c(length(b) == ntypes, length(d) == ntypes, 
                    dim(Xsi_as) == c(ntypes, ntypes), dim(Xsi_s) == c(ntypes, ntypes),
                    origin_time > 0, origin_type %in% c(0:ntypes),
                    d >= 0, d < 1, rho > 0, rho <= 1)),
              msg = 'The inputs do not have proper dimensions or values. Please check all parameters!')
  assert_that(all(sapply(seq(1, ntypes), function(i) {sum(2*Xsi_as[i, ]) + sum(Xsi_s[i, ]) == 1})),
              msg = 'The transition probabilities do not some to 1. Please check!')
  
  # simulate full tree
  tree = simulate_complete_tree(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s, min_tips)
  if (is.null(tree)) {
    return(NULL)
  }
  
  # prune dead particles
  dead_nodes = tree %>% as_tibble %>% filter(status == 0) %>% pull(label)
  phylogeny = drop.tip(tree, dead_nodes)
  if (is.null(phylogeny)) {
    message('All particle died! Try another seed.')
    return(NULL)
  }

  # prune unsampled particles
  tips = phylogeny@phylo$tip.label
  sampled_tips = sample(c(TRUE, FALSE), length(tips), prob = c(rho, 1 - rho), replace = TRUE)
  if (sum(sampled_tips) < min_tips) {
    message('Not enough tips are sampled! Try another seed.')
    return(NULL)
  }
  phylogeny = drop.tip(phylogeny, tips[!sampled_tips])
  phylogeny@phylo$tip.label = as.character(c(1:length(phylogeny@phylo$tip.label)))
  
  # remove height from data for exporting tree
  phylogeny@data = phylogeny@data %>% select(-status)
                                             
  # keep origin in final object
  phylogeny@phylo$origin = tree@phylo$origin
  
  return(phylogeny)
}


# Simulator of the complete Age-Dependent Branching Process for a fixed time interval (since origin)
simulate_complete_tree <- function(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s, min_tips = 2) {
  
  # initialize
  edges = matrix(nrow = 0, ncol = 2)
  edge_lengths = c()
  nodes = data.frame(
    id = integer(0), 
    height = numeric(0), 
    type = integer(0), 
    parent = integer(0), 
    leftchild = integer(0), 
    rightchild = integer(0), 
    status = integer(0) # 0 = dead, 1 = alive, 2 = divided
  )
  #nodes = matrix(nrow = 0, ncol = 7, 
  #               dimnames = list(NULL, c("id", "height", "type", "parent", "leftchild", "rightchild", "isTip")))
  
  # sample the lifetime of the first particle
  root_edge = rgamma(1, shape = b[origin_type + 1], scale = a[origin_type + 1])
  nodes = bind_rows(nodes, c(id = 1, height = origin_time - root_edge, type = origin_type, 
                             parent = NA, leftchild = NA, rightchild = NA, status = 1))
  events = nodes
  event_counter = 1
  
  while (nrow(events) > 0) {
    event = as.list(events[1, ]); events = events[-1, ] # look at one event

    if (runif(1) < d[event$type + 1]) {
      # particle dies
      nodes[event$id, "status"] = 0
    } else {
      # particle divides, create two new children
      left_id = event_counter + 1
      right_id = event_counter + 2
      event_counter = event_counter + 2
      nodes[event$id, "status"] = 2

      # sample types
      children_types = sample_types(event$type, Xsi_as, Xsi_s)

      # sample lifetimes and add new nodes
      left_lifetime = rgamma(1, shape = b[event$type + 1], scale = a[event$type + 1])
      if (event$height - left_lifetime < 0) {
        # censor lifetime
        left_lifetime = event$height
        left_node = c(id = left_id, height = event$height - left_lifetime, type = children_types[1], 
                      parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      } else {
        left_node = c(id = left_id, height = event$height - left_lifetime, type = children_types[1], 
                      parent = event$id, leftchild = NA, rightchild = NA, status = 1)
        events = bind_rows(events, left_node) # to be processed
      }

      right_lifetime = rgamma(1, shape = b[event$type + 1], scale = a[event$type + 1])
      if (event$height - right_lifetime < 0) {
        # censor lifetime, make tip
        right_lifetime = event$height
        right_node = c(id = right_id, height = event$height - right_lifetime, type = children_types[2], 
                      parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      } else {
        right_node = c(id = right_id, height = event$height - right_lifetime, type = children_types[2], 
                       parent = event$id, leftchild = NA, rightchild = NA, status = 1)
        events = bind_rows(events, right_node) # to be processed
      }
    
      # add child relationships and append the current edges and edge lengths
      nodes = bind_rows(nodes, left_node, right_node)
      nodes[event$id, c("leftchild", "rightchild")] = c(left_id, right_id)
      edges = rbind(edges, c(event$id, left_id), c(event$id, right_id))
      edge_lengths = c(edge_lengths, left_lifetime, right_lifetime)
    }
  }
  
  # check number of tips
  Nnode = sum(nodes$status == 2)
  Ntip = sum(nodes$status %in% c(0,1))
  if (Ntip < min_tips) {
    message("The simulated tree has too few tips. Try another seed.")
    return(NULL)
  }
  
  # assign labels
  nodes[which(nodes$status %in% c(0,1)), 'label'] = c(1:Ntip) 
  nodes[which(nodes$status == 2), 'label'] = c((Ntip + 1):(Ntip + Nnode)) 
  
  # use labels in edge matrix
  edges_recoded = matrix(nodes[as.double(edges), 'label'], nrow = nrow(edges), ncol = ncol(edges))
  
  # create phylogenetic tree
  phylo_tree = list(edge = edges_recoded, edge.length = edge_lengths, Nnode = Nnode, 
                    tip.label = as.character(nodes[which(nodes$status %in% c(0,1)), 'label']))
  class(phylo_tree) = "phylo"
  
  # create treedata object
  tree = as.treedata(phylo_tree)
  tree@phylo$root.edge = root_edge
  tree@phylo$origin = origin_time
  data = as_tibble(tree)
  types = nodes %>% 
    select(node = label, status, type) %>% 
    mutate(type = as.factor(type)) %>%
    arrange(node) %>%
    as_tibble
  tree@data = types
  
  return(tree)
}


# Helper function for simulators of ADBP for sampling types of offspring upon division
sample_types <- function(parent_type, Xsi_as, Xsi_s) {
  ntypes = ncol(Xsi_s)
  
  cum_prob = 0
  for (i in seq(0, ntypes - 1)) {
    # symmetric case
    cum_prob = cum_prob + Xsi_s[parent_type + 1, i + 1];
    if (runif(1) < cum_prob) {
      return(c(i, i))
    }
    # asymmetric case
    cum_prob = cum_prob + 2*Xsi_as[parent_type + 1, i + 1];
    if (runif(1) < cum_prob) {
      if (runif(1) < 0.5) {
        return(c(parent_type, i))
      } else {
        return(c(i, parent_type))
      }
    }
  }
}


# # Example
# origin_time = 10
# origin_type = 0
# a = c(2, 3)
# b = c(1, 2)
# d = c(0.2, 0.1)
# Xsi_as = rbind(c(0, 0.1), c(0.2, 0))
# Xsi_s = rbind(c(0.5, 0.3), c(0.1, 0.5))
# rho = 0.8
# set.seed(1)

# tree = simulate_complete_tree(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s)
# ggtree(tree) + geom_point(aes(color = type)) + geom_rootedge() + geom_tiplab() + theme_tree2()
# ggtree(tree) + geom_rootedge() + geom_point(aes(x = x - branch.length, color = type), size = 2)
# phylogeny = simulate_phylogeny(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s, rho)
# ggtree(phylogeny) + geom_point(aes(color = type))
# write.beast.newick(phylogeny)


# other:
# beast = read.beast.newick(textConnection("<newick_string>;")) # to plot tree from Java simulator
# beast = read.beast.newick("~/intellij/ADBP/examples/tree_typed.newick") # example
# beast@data = beast@data %>% mutate(type = as.factor(type))
# ggtree(beast) + geom_point(aes(color = type))

