# Functions for data set manipulation and data sampling

reorder_dataset = function(data, direction = c("variable-major", "location-major"), order = NULL){
  # Extract variable names and location numbers
  var_names <- unique(sub("\\.\\d+$", "", colnames(data)))  # Extract Variables
  
  if (!(is.null(order))) var_names = var_names[order]
  
  loc_numbers <- unique(sub("^.*\\.", "", colnames(data)))  # Extract Location
  
  if (direction == "variable-major"){
    new_order = as.vector(sapply(var_names, function(v) paste0(v, ".", loc_numbers)))
  } else if (direction == "location-major"){
    new_order = as.vector(sapply(loc_numbers, function(i) paste0(var_names, ".", i)))
  } else {
    stop("Please select the direction of the reordering.")
  }
  
  # Reorder dataset
  return(data.frame(data)[, new_order])
}

create_simulation_data = function(nlocs, nvars, nsample, seed = 1, method = c("indep", "random", "given"), cov_mat = NULL){
  set.seed(seed)
  d = nlocs * nvars
  
  if (method == "indep"){
    data <- matrix(rnorm(nsample * d), nsample, d)
  } else if (method == "random"){
    # Create a random orthogonal matrix Q
    Q <- qr.Q(qr(matrix(rnorm(d^2), d, d)))
    
    # Create random eigenvalues (positive)
    eigenvalues <- runif(d, min = 0.5, max = 2)
    
    # Construct covariance matrix: Σ = Q Λ Qᵗ
    Sigma <- Q %*% diag(eigenvalues) %*% t(Q)
    
    data <- mvrnorm(n = nsample, mu = rep(0, d), Sigma = Sigma)
  } else if (method == "given"){
    if (is.null(cov_mat)) stop("Please specify a covariance matrix.") 
    data <- mvrnorm(n = nsample, mu = rep(0, d), Sigma = cov_mat)
  } else {
    stop("Please select a sampling method.")
  }
  
  colnames(data) <- make_names(nlocs, nvars)
  return(data)
}

make_names <- function(nlocs, nvars) {
  vars <- paste0("v_", letters[1:nvars])
  locs <- paste0(".", 1:nlocs)
  return(as.vector(outer(vars, locs, FUN = paste0)))
}

sample_from_level3 = function(rvs_level1, rvs_level2, nsample, seed = 1, cores = 4){
  set.seed(seed)
  
  # Merge the vines
  rvs_level3 = merge_structures(rvs_level1, rvs_level2)
  
  # Specify non-parametric pair copulas for tree_level 1: nvars - 1
  d = rvs_level3$d
  pair_copulas <- list()
  
  for (tree_level in 1:(nvars-1)) {
    n_copulas_in_tree <- d - tree_level
    tree_list <- vector("list", n_copulas_in_tree)
    
    for (i in 1:n_copulas_in_tree) {
      tree_list[[i]] <- bicop_dist(family = "tll")  # or customize kernel parameters
    }
    
    pair_copulas[[tree_level]] <- tree_list
  }
  
  # Create the vine copula model
  vine_fit <- vinecop_dist(structure = rvs_level3, pair_copulas = pair_copulas)
  
  # Sample from that vine copula
  vine_sample = rvinecop(n = nsample, vine_fit)
  
  return(vine_sample)
}

struct_to_matrix = function(vine_struct){
  d = vine_struct$d
  mat = set_antidiagonal(matrix(0, nrow = d, ncol = d), rev(vine_struct$order))
  
  for (iRow in 1:(d-1)){
    if (length(vine_struct$struct_array) >= iRow){
      mat[iRow, 1:(d-iRow)] = vine_struct$order[vine_struct$struct_array[[iRow]]]
    }
  }
  return(mat)
}

set_antidiagonal <- function(mat, values) {
  n <- nrow(mat)
  mat[cbind(1:n, n:1)] <- values
  return(mat)
}

sample_rvs_locations = function(nrows, ncols, mask, seed = 123){
  set.seed(seed)
  if (mask) {
    global_mask <- get_spatial_mask(nrows, ncols)
    
    g <- graph_from_adjacency_matrix(
      as.matrix(global_mask),
      mode = "undirected",
      weighted = TRUE, diag = FALSE)
    
    tree <- mst(g)
    rvs_tmp <- spanning_tree_to_rvine_structure(tree)
    
    data_sample = create_simulation_data(nrows*ncols, 1, 100, seed = seed, method = "random")
    
    vine_fit = vinecop(pseudo_obs(data_sample), family_set = "tll", structure = rvs_tmp)
    
    rvs = vine_fit$structure
    
  } else {
    rvs = rvine_structure_sim(nrows * ncols)
  }
  
  return(rvs)
}