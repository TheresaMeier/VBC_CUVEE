# Functions for tree manipulations and generation

get_spatial_mask <- function(rows, cols) {
  # Grid dimensions
  n_cells <- rows * cols
  
  # Initialize sparse adjacency matrix
  adj_matrix <- Matrix(0, nrow = n_cells, ncol = n_cells, sparse = TRUE)
  
  # Function to convert (row, col) to index in vectorized form
  cell_index <- function(r, c) {
    if (r >= 1 && r <= rows && c >= 1 && c <= cols) {
      return((r - 1) * cols + c)
    } else {
      return(NA) # Out of bounds
    }
  }
  
  # Neighbor offsets (left, right, up, down, diagonals)
  neighbors <- list(
    c(-1, 0), c(1, 0), c(0, -1), c(0, 1),
    c(-1, -1), c(-1, 1), c(1, -1), c(1, 1)
  )
  
  # Loop through each grid cell
  for (r in 1:rows) {
    for (c in 1:cols) {
      current_index <- cell_index(r, c)
      for (offset in neighbors) {
        neighbor_index <- cell_index(r + offset[1], c + offset[2])
        if (!is.na(neighbor_index)) {
          adj_matrix[current_index, neighbor_index] <- 1
        }
      }
    }
  }
  
  # Convert to symmetric matrix (undirected adjacency)
  adj_matrix <- forceSymmetric(adj_matrix, uplo = "L")
  return(adj_matrix)
}

sample_tree_wilson <- function(g, root = 1, seed = 123) {
  stopifnot(is_connected(g), is_simple(g), !is_directed(g))
  set.seed(seed)
  root <- sample(V(g), 1)
  n <- vcount(g)
  visited <- rep(FALSE, n)
  visited[root] <- TRUE
  tree_edges <- c()
  
  while (any(!visited)) {
    start <- sample(V(g)[!visited], 1)
    path <- integer(0)
    current <- as.integer(start)
    prev <- rep(NA_integer_, n)
    
    while (!visited[current]) {
      neighbors <- as.integer(neighbors(g, current))
      next_node <- sample(neighbors, 1)
      prev[current] <- next_node
      current <- next_node
    }
    
    current <- as.integer(start)
    seen <- rep(NA_integer_, n)
    while (!visited[current]) {
      seen[current] <- prev[current]
      visited[current] <- TRUE
      current <- prev[current]
    }
    
    for (i in which(!is.na(seen))) {
      edge <- get_edge_ids(g, c(i, seen[i]), directed = FALSE)
      tree_edges <- c(tree_edges, edge)
    }
  }
  
  subgraph_from_edges(g, E(g)[tree_edges])
}

spanning_tree_to_rvine_structure <- function(spanning_tree) {
  # Turn a spanning tree into an R-vine structure
  order <- bfs(spanning_tree, root = 1, order = TRUE)$order
  edges <- as_edgelist(spanning_tree, names = FALSE)
  
  row1 <- rep(NA, nrow(edges))
  j <- 1
  
  for (i in order) {
    ind_i <- which(edges[, 1] == i)
    ind_i_2 <- which(edges[, 2] == i)
    
    if (length(ind_i) != 0) {
      for (k in 1:length(ind_i)) {
        row1[j] <- edges[ind_i[k], 1]
        edges[ind_i[k], ] <- NA
        j <- j + 1
      }
    }
    if (length(ind_i_2) != 0) {
      for (k in 1:length(ind_i_2)) {
        row1[j] <- edges[ind_i_2[k], 2]
        edges[ind_i_2[k], ] <- NA
        j <- j + 1
      }
    }
    ind_i <- c()
    ind_i_2 <- c()
  }
  
  rvine_struct <- rvine_structure(order = rev(as.vector(order)), struct_array = list(rev(row1)))
  return(rvine_struct)
}

# Function to get all the edges in the structure
generate_edge_list <- function(vine_struct) {
  d <- vine_struct$d
  
  # Preallocate edge list matrix
  edge_list <- matrix(NA, nrow = sum(1:(d-1)), ncol = d)
  colnames(edge_list) <- paste0("e", 1:d)
  
  start_idx <- 1  # Row index tracker
  
  for (tree in 1:(d-1)) {
    if (length(vine_struct$struct_array) >= tree){
      num_edges <- d - tree  # Number of edges in this tree
      
      edge_list[start_idx:(start_idx + num_edges - 1), "e1"] <- vine_struct$order[1:num_edges]
      edge_list[start_idx:(start_idx + num_edges - 1), "e2"] <- vine_struct$order[vine_struct$struct_array[[1]]][1:num_edges]
      
      if (tree > 1) {
        for (j in 2:tree) {
          edge_list[start_idx:(start_idx + num_edges - 1), paste0("e", j + 1)] <- 
            vine_struct$order[vine_struct$struct_array[[j]]][1:num_edges]
        }
      }
      
      start_idx <- start_idx + num_edges
    }
  } 
  
  #edge_list = edge_list[!apply(is.na(edge_list), 1, all), ]
  #edge_list = edge_list[, !apply(is.na(edge_list), 2, all)]
  
  return(edge_list)
}