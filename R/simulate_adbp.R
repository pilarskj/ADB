library(ape)
library(treeio)
library(dplyr)
library(assertthat)
library(ggtree) # for plotting


#' Simulator of a phylogeny from an Age-Dependent Branching Process
#' @param origin_time time of birth of the initial particle
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param Xsi_as matrix of asymetric type transition probabilities
#' @param Xsi_s matrix of symetric type transition probabilities
simulate_phylogeny <- function(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s, rho, min_tips = 2) {
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
  dead_nodes = tree %>% as_tibble %>% filter(!is.na(label) & height > 0) %>% pull(label)
  phylogeny = drop.tip(tree, dead_nodes)
  if (is.null(phylogeny)) {
    print('All particle died! Try another seed.')
    return(NULL)
  }

  # prune unsampled particles
  tips = phylogeny@phylo$tip.label
  sampled_tips = sample(c(TRUE, FALSE), length(tips), prob = c(rho, 1 - rho), replace = TRUE)
  if (sum(sampled_tips) < min_tips) {
    print('Not enough tips are sampled! Try another seed.')
    return(NULL)
  }
  phylogeny = drop.tip(phylogeny, tips[!sampled_tips])
  
  # remove height from data for exporting tree
  phylogeny@data = phylogeny@data %>% select(-height)
  
  return(phylogeny)
}


# Simulator of the complete Age-Dependent Branching Process
simulate_complete_tree <- function(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s, min_tips = 2) {
  
  # initialize
  edges = matrix(nrow = 0, ncol = 2)
  edge_lengths = c()
  nodes = matrix(nrow = 0, ncol = 7, 
                 dimnames = list(NULL, c("id", "height", "type", "parent", "leftchild", "rightchild", "isTip")))
  
  # sample the lifetime of the first particle
  root_edge = rgamma(1, shape = b[origin_type + 1], scale = a[origin_type + 1])
  nodes = rbind(nodes, c(1, origin_time - root_edge, origin_type, NA, NA, NA, F))
  events = rbind(matrix(nrow = 2, ncol = 7), nodes) # add dummy rows due to type conversions
  event_counter = 1
  
  while (nrow(events) > 2) {
    event = as.list(events[3, ]); events = events[-3, ] # look at one event

    if (runif(1) < d[event$type + 1]) {
      # particle dies
      nodes[event$id, "isTip"] = T
    } else {
      # particle divides, create two new children
      left_id = event_counter + 1
      right_id = event_counter + 2
      event_counter = event_counter + 2

      # sample types
      children_types = sample_types(event$type, Xsi_as, Xsi_s)
      left_type = children_types[1]
      right_type = children_types[2]

      # sample lifetimes and add new nodes
      left_lifetime = rgamma(1, shape = b[event$type + 1], scale = a[event$type + 1])
      if (event$height - left_lifetime < 0) {
        # censor lifetime, make tip
        left_lifetime = event$height
        nodes = rbind(nodes, c(left_id, 0, left_type, event$id, NA, NA, T), deparse.level = 0)
      } else {
        left_node = c(left_id, event$height - left_lifetime, left_type, event$id, NA, NA, F)
        nodes = rbind(nodes, left_node, deparse.level = 0)
        events = rbind(events, left_node, deparse.level = 0) # to be processed
      }

      right_lifetime = rgamma(1, shape = b[event$type + 1], scale = a[event$type + 1])
      if (event$height - right_lifetime < 0) {
        # censor lifetime, make tip
        right_lifetime = event$height
        nodes = rbind(nodes, c(right_id, 0, right_type, event$id, NA, NA, T), deparse.level = 0)
      } else {
        right_node = c(right_id, event$height - right_lifetime, right_type, event$id, NA, NA, F)
        nodes = rbind(nodes, right_node, deparse.level = 0)
        events = rbind(events, right_node, deparse.level = 0) # to be processed
      }
    
      # add child relationships and append the current edges and edge lengths
      nodes[event$id, c("leftchild", "rightchild")] = c(left_id, right_id)
      edges = rbind(edges, c(event$id, left_id), c(event$id, right_id), deparse.level = 0)
      edge_lengths = c(edge_lengths, left_lifetime, right_lifetime)
    }
  }
  
  # check number of tips
  nodes = data.frame(nodes)
  Nnode = sum(!nodes$isTip)
  Ntip = sum(nodes$isTip)
  if (Ntip < min_tips) {
    print("The simulated tree has too few tips. Try another seed.")
    return(NULL)
  }
  
  # assign labels
  nodes[which(nodes$isTip == T), 'label'] = c(1:Ntip) 
  nodes[which(nodes$isTip == F), 'label'] = c((Ntip + 1):(Ntip + Nnode)) 
  
  # use labels in edge matrix
  edges_recoded = matrix(nodes[as.double(edges), 'label'], nrow = nrow(edges), ncol = ncol(edges))
  
  # create phylogenetic tree
  phylo_tree = list(edge = edges_recoded, edge.length = edge_lengths, Nnode = Nnode, 
                    tip.label = as.character(nodes[which(nodes$isTip == T), 'label']))
  class(phylo_tree) = "phylo"
  
  # create treedata object
  tree = as.treedata(phylo_tree)
  tree@phylo$root.edge = root_edge
  data = as_tibble(tree)
  types = nodes %>% 
    select(node = label, height, type) %>% 
    mutate(type = as.factor(type)) %>%
    arrange(node) %>%
    as_tibble
  tree@data = types
  
  return(tree)
}

  
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


# Example
origin_time = 10
origin_type = 0
a = c(2, 3)
b = c(1, 2)
d = c(0.2, 0.1)
Xsi_as = rbind(c(0, 0.1), c(0.2, 0))
Xsi_s = rbind(c(0.5, 0.3), c(0.1, 0.5))
rho = 0.8
set.seed(1)

# tree = simulate_complete_tree(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s)
# ggtree(tree) + geom_point(aes(color = type)) + geom_rootedge() + geom_tiplab() 
# ggtree(tree) + geom_rootedge() + geom_point(aes(x = x - branch.length, color = type), size = 2)
phylogeny = simulate_phylogeny(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s, rho)
ggtree(phylogeny) + geom_point(aes(color = type)) 
write.beast.newick(phylogeny)
#for (i in seq(1, 1000)) {
#  simulate_phylogeny(origin_time, origin_type, a, b, d, Xsi_as, Xsi_s, rho)
#}


# other:
# beast = read.beast.newick(textConnection("<newick_string>;")) # to plot tree from Java simulator
# beast = read.beast.newick("~/intellij/ADBP/examples/tree_typed.newick") # example
# beast@data = beast@data %>% mutate(type = as.factor(type))
# ggtree(beast) + geom_point(aes(color = type))

