# Functions for merging level 1 and level 2

edges_to_matrix_tree1 <- function(edge_list) {
  # edge_list: list of edges (a, b), each element a vector of two integers
  
  all_vars <- sort(unique(unlist(edge_list)))
  d <- length(all_vars)
  
  # Build undirected graph from edge list
  g <- make_empty_graph(n = d) %>%
    add_edges(unlist(edge_list)) %>%
    as_undirected()
  
  if (!is_connected(g) || ecount(g) != d - 1) {
    stop("Edge list must form a connected spanning tree on d nodes.")
  }
  
  # Assign node labels to 1..d if needed
  V(g)$name <- as.character(all_vars)
  
  # Get vine variable order from BFS traversal
  root <- as.character(all_vars[1])
  bfs_res <- bfs(g, root = root, order = TRUE)
  vine_order <- as.integer(bfs_res$order[!is.na(bfs_res$order)])
  
  # Initialize vine matrix
  M <- matrix(0, d, d)
  M[1, d] <- vine_order[1]
  for (i in 2:d) {
    M[i, d - i + 1] <- vine_order[i]
  }
  
  # Fill Tree 1 entries
  # Match edge direction to vine_order (a conditioned on b)
  for (e in 1:(d - 1)) {
    # Get the conditioned variable (counter-diagonal entry)
    a <- M[d - e + 1, e]
    
    # Find the neighbor b such that (a, b) is in the edge list
    edge_candidates <- Filter(
      function(x) a %in% x,
      edge_list
    )
    # Among them, select the b that is earlier in the vine_order
    b <- NULL
    for (edge in edge_candidates) {
      other <- setdiff(edge, a)
      if (which(vine_order == other) < which(vine_order == a)) {
        b <- other
        break
      }
    }
    if (is.null(b)) stop(paste("No valid parent found for node", a))
    
    M[1, e] <- b  # fill in conditioned-on variable
  }
  
  return(M)
}

# Function to get all edges in tree 1 in the structure
generate_edge_list_tree1 <- function(vine_struct) {
  d <- vine_struct$d
  
  # Preallocate edge list matrix
  edge_list <- matrix(NA, nrow = d-1, ncol = 2)
  colnames(edge_list) <- paste0("e", 1:2)
  
  edge_list[, "e1"] <- vine_struct$order[1:(d-1)]
  edge_list[, "e2"] <- vine_struct$order[vine_struct$struct_array[[1]]][1:(d-1)]
  
  return(edge_list)
}


get_bridging_edges = function(rvs, bridge_var, nvars){
  
  d = rvs$d
  edges = generate_edge_list_tree1(rvs)
  
  ## Map the edges to the full data set
  mapping = cbind("original" = 1:d, "remapped" = bridge_var + nvars * seq(0, d-1))
  
  # Remap function
  remap_values <- function(edges, mapping) {
    apply(edges, c(1,2), function(x) mapping[mapping[, "original"] == x, "remapped"])
  }
  
  edges_remapped <- remap_values(edges, mapping)
  
  # return list with bridging edges
  return(edges_remapped) 
}


merge_edges_fixed = function(rvs_level1, rvs_level2, bridge_var){
  
  nlocs = rvs_level1$d
  nvars = rvs_level2$d
  
  # Get the edges of both levels
  edges_level1 = generate_edge_list_tree1(rvs_level1)
  edges_level2 = generate_edge_list_tree1(rvs_level2)
  
  # Create the edges of level 2 for all variables
  edges_level_2_all = matrix(NA, nrow = 0, ncol = 2)
  colnames(edges_level_2_all) <- paste0("e", 1:2)
  
  for (iLoc in 1:nlocs){
    edges_level_2_all = rbind(edges_level_2_all, edges_level2 + (iLoc - 1) * nvars)
  }
  
  # Get bridging edges
  bridging_edges = get_bridging_edges(rvs_level1, bridge_var, nvars)
  
  # Combine everything
  edges_all = rbind(edges_level_2_all, bridging_edges)
  edges_list = lapply(seq_len(nrow(edges_all)), function(i) as.integer(edges_all[i, ]))
  
  # Construct R vine matrix
  rvine_mat_merged = edges_to_matrix_tree1(edges_list)
  
  # Transfrom it to R vine structure
  order_level3 = rev(rvine_mat_merged[cbind(1:(nlocs*nvars), (nlocs*nvars):1)])
  tree1 = list(rvine_mat_merged[1, -(nlocs*nvars)])
  
  return(rvine_structure(order = order_level3, struct_array = tree1))
}


merge_edges_individual = function(rvs_level1, rvs_level2, bridge_var){
  
  nlocs = nrow(rvs_level1[[1]])
  nvars = rvs_level2$d
  
  # Get the edges of both levels
  edges_level1 = lapply(rvs_level1, FUN = generate_edge_list_tree1)
  edges_level2 = generate_edge_list_tree1(rvs_level2)
  
  # Create the edges of level 2 for all variables
  edges_level_1_all = matrix(NA, nrow = 0, ncol = 2)
  colnames(edges_level_1_all) <- paste0("e", 1:2)
  
  for (iLoc in 1:nvars){
    edges_level_1_all = rbind(edges_level_1_all, edges_level1[[iLoc]] + (iLoc - 1) * nlocs)
  }
  
  # Get bridging edges
  bridging_edges = get_bridging_edges(rvs_level2, bridge_var, nlocs)
  
  # Combine everything
  edges_all = rbind(edges_level_1_all, bridging_edges)
  edges_list = lapply(seq_len(nrow(edges_all)), function(i) as.integer(edges_all[i, ]))
  
  # Construct R vine matrix
  rvine_mat_merged = edges_to_matrix_tree1(edges_list)
  
  # Transfrom it to R vine structure
  order_level3 = rev(rvine_mat_merged[cbind(1:(nlocs*nvars), (nlocs*nvars):1)])
  tree1 = list(rvine_mat_merged[1, -(nlocs*nvars)])
  
  return(rvine_structure(order = order_level3, struct_array = tree1))
}
# 
# merge_structures_inverse = function(rvs_level1, rvs_level2, bridge_var){
#   nvars <- length(rvs_level1)
#   nlocs <- dim(rvs_level1[[1]])[1]
# 
#   rvs_level1_mapped <- mapply(function(matrix, integer) {
#     matrix[matrix != 0] <- matrix[matrix != 0] + integer
#     return(matrix)
#   }, rvs_level1, (c(0:nlocs) * nlocs)[order(rev(rvs_level2$order))], SIMPLIFY = FALSE)
# 
#   # Stack the matrices all together
#   # Function to extract antidiagonals from a matrix
#   extract_antidiagonals = function(mat) {
#     n = nrow(mat)
#     anti_diag = mat[cbind(1:n, n:1)]
#     return(rev(anti_diag))
#   }
#   order_level3 = unlist(lapply(rvs_level1_mapped[rvs_level2$order], extract_antidiagonals))
# 
#   ## Reuse the
#   rvs_level1_mapped_tmp = lapply(rvs_level1_mapped, function(x) {
#     x[cbind(1:nrow(x), nrow(x):1)] = 0
#     return(x)
#   } )
# 
#   rvs_level1_mapped_stacked = do.call(cbind, rvs_level1_mapped_tmp[rvs_level2$order])
# 
#   # Set up the matrix
#   matrix_level3 = matrix(0, nrow = nlocs*nvars, ncol = nlocs * nvars)
#   matrix_level3[1:nvars,] = rvs_level1_mapped_stacked[1:nvars,]
#   matrix_level3[cbind(1:(nlocs*nvars), (nlocs*nvars):1)] = rev(order_level3)
# 
#   ## Create mapping for level 2
#   mat_level2 = struct_to_matrix(rvs_level2)
#   mat_level2_mapped = mat_level2
#   replacement_values <- (c(0:nlocs) * nlocs)[order(rev(rvs_level2$order))]+bridge_var
# 
#   mapping <- setNames(replacement_values, 1:nvars)
# 
#   mat_level2_mapped[mat_level2_mapped %in% 1:nvars] <- mapping[mat_level2_mapped[mat_level2_mapped %in% 1:nvars]]
#   mat_level2_mapped[cbind(1:nvars, nvars:1)] = 0
# 
#   for (iTree in (1:nvars)){
#     idx = which(matrix_level3[iTree,] == 0) #1:(nvars*nlocs - iTree + 1)
#     matrix_level3[iTree, idx] = c(mat_level2_mapped[1:iTree,]) [1:length(idx)]#[1:iTree,-(nvars:(nvars-iTree+1))]
#   }
# 
#   # Fill in missing values: Reuse the order of level 1 structures
#   orders_level1 = rev(lapply(rvs_level1_mapped[rvs_level2$order], extract_antidiagonals))
#   mapping_missings = cbind(c(1:nvars), rev(rvs_level2$order))
# 
#   for (iTree in (2:nvars)){
#     idx = rev(which(matrix_level3[iTree, 1:(nvars*nlocs - iTree + 1)] == 0))
#     id_var = sapply(idx, function(val) {
#       which(sapply(orders_level1, function(x) val %in% x))
#     })
# 
#     for (j in 1:length(idx)){
#       order_tmp = orders_level1[[j]]
#       matrix_level3[iTree:min(nlocs,nvars), idx[j]] = rev(order_tmp)[2:(nvars-iTree+2)]
#     }
#   }
# 
#   # Transform level 3 matrix into structure
# 
#   struct_array_merged = list()
# 
#   for (iTree in 1:nvars){
#     struct_array_merged = append(struct_array_merged, list(matrix_level3[iTree,1:(nvars*nlocs - iTree) ]))
#   }
# 
#   rvs_level3_merged = rvine_structure(order = order_level3, struct_array = struct_array_merged)
# 
#   return(rvs_level3_merged)
# }
# 
# 
# merge_structures = function(rvs_level1, rvs_level2){
#   nlocs <- rvs_level1$d
#   nvars <- rvs_level2$d
# 
#   # Create mapping
#   mapping = matrix(c(1:nlocs, rev(rvs_level1$order), (rev(rvs_level1$order)-1) * nvars + 1), ncol = 3)
# 
#   ### Order level 3
#   order_level3 <- rep(NA, nlocs * nvars)
# 
#   j <- 1
#   for (i in seq(1, length(order_level3), by = nvars)) {
#     order_level3[i:(i + nvars - 1)] <- rvs_level2$order + nvars * (rvs_level1$order[j]-1)
#     j <- j + 1
#   }
# 
#   ### First nvars-1 rows
# 
#   # Create mappings
#   mapping_level2 = matrix(c(1:(nvars*nlocs)), ncol = nlocs)[, rev(rvs_level1$order)]
#   tmp = c()
#   for (i in rev(rvs_level1$order)){
#     tmp = c(tmp, rev(rvs_level2$order) + nvars * (i-1))
#   }
#   mapping_2 = matrix(c(1:(nvars*nlocs), tmp), ncol = 2)
# 
#   # Get connections between the trees
#   connections_treei = c()
#   for (iRow in 1:(nvars-1)){
#     rowi_level1 = unlist(rvs_level1$order[rvs_level1$struct_array[[iRow]]])
# 
#     fill_vars = c()
#     if (iRow > 1){
#       for (iVar in iRow:2){
#         fill_vars = c(fill_vars, order_level3[length(order_level3)-iVar+1])
#       }
#     }
# 
#     connections_treei =  cbind(connections_treei, c(fill_vars, rev(nvars * (rowi_level1-1) + (rvs_level2$order)[nvars])))
#   }
#   # Initialize the list for the structure
#   struct_array = list()
# 
#   # Create the trees
#   for (iRow in 1:(nvars-1)){
#     rowi_level2 <- unlist(rvs_level2$order[rvs_level2$struct_array[[iRow]]])
#     rowi_level3 <- c()
#     for (iLoc in nlocs:1){
#       rowi_level3 = c(rowi_level3, mapping_level2[match(rowi_level2, mapping_2[,1]), iLoc], connections_treei[iLoc-1, 1:iRow])
#     }
# 
#     struct_array = append(struct_array,list(rowi_level3))
#   }
#   # Create the R vine
#   rvs_level3 <- rvine_structure(
#     order = order_level3,
#     struct_array = struct_array,
#     is_natural_order = FALSE)
# 
#   return(rvs_level3)
# }
# 
# rotate_until_last <- function(perm, i) {
#   stopifnot(i %in% perm)
#   pos <- which(perm == i)
#   n <- length(perm)
#   # Perform cyclic rotation: move everything from pos+1 to end, then from start to pos
#   rotated <- c(perm[(pos + 1):n], perm[1:pos])
#   return(rotated)
# }
# 
# rotate_vine <- function(rvs, var_index) {
#   rvs_rotated <- rvs
#   rvs_rotated$order <- rotate_until_last(rvs$order, var_index)
#   return(rvs_rotated)
# }