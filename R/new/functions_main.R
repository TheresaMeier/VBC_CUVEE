## Main functions: compute hierarchical vines

get_nested_vine <- function(
    data, nrows, ncols, nvars,
    direction = c("loc-var", "var-loc"),
    fixed = TRUE,             # Indicator if fixed local level is assumed or not
    mask = TRUE,
    bridge_var = NULL,
    fit_levels = FALSE,
    fit_final = TRUE,
    location_vine = c("random", "dissman"),
    variable_vine = c("random", "dissman"),
    seed = 123, cores = 1) {
  
  stopifnot(ncol(data) == nrows * ncols * nvars)
  stopifnot(location_vine %in% c("random", "dissman"))
  stopifnot(variable_vine %in% c("random", "dissman"))
  
  location_vine <- match.arg(location_vine)
  variable_vine <- match.arg(variable_vine)
  
  # Set the seed
  set.seed(seed)
  
  if (direction == "var-loc"){          # original
    
    # # Ensure the correct ordering of the data set
    # colnames_tmp = colnames(data)
    # data = reorder_dataset(data, direction = "location-major")
    # 
    # if (any(colnames_tmp != colnames(data))) print("Warning: Dataset is reordered to 'location-major' format.")

    output = nested_vine_var_loc(data, nrows, ncols, nvars,
                                 mask, bridge_var, fixed,
                                 fit_levels, fit_final,
                                 location_vine, variable_vine,
                                 seed = 123, cores = 1)
    
  } else if (direction == "loc-var"){       # inverse
    
    # # Ensure the correct ordering of the data set
    # colnames_tmp = colnames(data)
    # data = reorder_dataset(data, direction = "variable-major")
    # 
    # if (any(colnames_tmp != colnames(data))) print("Warning: Dataset is reordered to 'variable-major' format.")
    # 
    output = nested_vine_loc_var(data, nrows, ncols, nvars,
                                 mask, bridge_var, fixed,
                                 fit_levels, fit_final,
                                 location_vine, variable_vine,
                                 seed = 123, cores = 1)
  } else{
    stop("Please specify a direction for the merge.")
  }
  
  return(output)
}


nested_vine_var_loc = function(data, nrows, ncols, nvars,
                               mask = TRUE,
                               bridge_var = NULL,
                               fixed = TRUE,
                               fit_levels = FALSE,
                               fit_final = TRUE,
                               location_vine = c("random", "dissman"),
                               variable_vine = c("random", "dissman"),
                               seed = 123, cores = 1){
  
  d <- ncol(data)
  nlocs = nrows * ncols

  # Sample a bridging variable
  if (is.null(bridge_var)) {
    bridge_var = sample(1:nvars, 1)
  } else if (!(bridge_var %in% c(1:nvars))){
    stop(paste("Bridging variable must be between 1 and", nvars))
  }
  
  output <- list()
  output[["bridge_var"]] <- bridge_var
  
  ### Level 1: Get global (location) vine
  ## Get mask to define possible spatial connections
  if (location_vine == "dissman") {
    if (mask) {
      global_mask <- get_spatial_mask(nrows, ncols)
      
      ktau_matrix <- wdm(data[, seq(bridge_var, d, nvars)], method = "kendall")
      weight_matrix <- global_mask * (1 - abs(ktau_matrix)) # to get maximum spanning tree
      g <- graph_from_adjacency_matrix(
        as.matrix(weight_matrix),
        mode = "undirected",
        weighted = TRUE, diag = FALSE)
      
      tree <- mst(g)
      rvs_level1 <- spanning_tree_to_rvine_structure(tree)
      
      vine_level1_tmp = vinecop(pseudo_obs(data[, ((bridge_var-1)*nlocs+1):((bridge_var-1)*nlocs+nlocs)]), trunc_lvl = 1,
                                family_set = "tll", structure = rvs_level1, cores = cores)
      
    } else {
      vine_level1_tmp = vinecop(pseudo_obs(data[, seq(bridge_var, ncol(data), nvars)]), trunc_lvl = 1,
                                family_set = "tll", cores = cores)
    }
  } else {
    if (mask) {
      global_mask <- get_spatial_mask(nrows, ncols)
      
      g <- graph_from_adjacency_matrix(
        global_mask,
        mode = "undirected",
        diag = FALSE
      )
      tree <- sample_tree_wilson(g, seed)
      rvs_level1 <- spanning_tree_to_rvine_structure(tree)
      
      vine_level1_tmp = vinecop(pseudo_obs(data[, ((bridge_var-1)*nlocs+1):((bridge_var-1)*nlocs+nlocs)]), trunc_lvl = 1,
                                family_set = "tll", structure = rvs_level1, cores = cores, tree_algorithm = "random_weighted")
      
    } else {
      vine_level1_tmp = vinecop(pseudo_obs(data[, ((bridge_var-1)*nlocs+1):((bridge_var-1)*nlocs+nlocs)]), trunc_lvl = 1,
                                family_set = "tll", cores = cores, tree_algorithm = "random_weighted")
    }
  }
  
  
  rvs_level1 = vine_level1_tmp$structure
  
  output[["rvs_level1"]] <- rvs_level1
  
  if (fit_levels) {
    # Find the vine based on the matching variable
    vine_level1 = vinecop(pseudo_obs(data[, ((bridge_var-1)*nlocs+1):((bridge_var-1)*nlocs+nlocs)]), 
                          family_set = "tll", trunc_lvl = 1,
                          structure = rvs_level1, cores = cores)  
    
    output[["vine_level1"]] = vine_level1
  }
  
  ### Level 2: Get local (inter-variable) vine
  
  if (variable_vine == "dissman") { 
    tree_alg = "mst_prim"
  } else{
    tree_alg = "random_weighted"
  }
  
  # Sample a location
  
  if (fixed){
    loc = sample(1:nlocs,1)
    
    # Fit the vine to that location
    vine_level2 =  vinecop(pseudo_obs(data[, seq(loc, nlocs*nvars, by = nlocs)]), 
                           family_set = "tll", cores = cores, tree_algorithm = tree_alg, trunc_lvl = 1) 
    
    rvs_level2 = vine_level2$structure
    
    output[["rvs_level2"]] <- rvs_level2
  } else {
    rvs_level2 = list()
    
    for (iLoc in 1:nlocs){
      # Fit the vine to each location separately
      vine_level2_iLoc =  vinecop(pseudo_obs(data[, seq(iLoc, nlocs*nvars, by = nlocs)]), 
                             family_set = "tll", cores = cores, tree_algorithm = tree_alg, trunc_lvl = 1) 
      
      output[[paste0("vine_level2_loc", iLoc)]] = vine_level2_iLoc
      output[[paste0("rvs_level2_loc", iLoc)]] = vine_level2_iLoc$structure
      rvs_level2[[iLoc]] = vine_level2_iLoc$structure
    }
  }
  
  
  if (fit_levels & fixed) {
    # Fit the vine to each location separately
    for (loc in 1:nlocs) {
      vine_level2_loc =  vinecop(pseudo_obs(data[, seq(loc + (loc-1) * (nvars-1), loc + loc * (nvars-1))]), 
                                 family_set = "tll", trunc_lvl = 1,
                                 structure = rvs_level2, cores = cores)  
      
      output[[paste0("vine_level2_loc", loc)]] = vine_level2_loc
    }
  }
  
  ### Level 3: Merge the vine structures
  if (fixed){
    rvs_level3 <- merge_edges_fixed(rvs_level1, rvs_level2, bridge_var)
  } else {
    rvs_level3 <- merge_edges_individual(rvs_level2, rvs_level1, bridge_var)
  }
  
  output[["rvs_level3"]] <- rvs_level3
  
  # Reorder data accordingly
  if (fit_final) {
    
    # Fit the corresponding R vine
    vine_level3 <- vinecop(pseudo_obs(reorder_dataset(data, direction = "location-major")),
                           family_set = "tll", trunc_lvl = 1,
                           structure = rvs_level3, cores = cores)
    
    output[["vine_level3"]] <- vine_level3
  }
  
  return(output)
}


nested_vine_loc_var = function(data, nrows, ncols, nvars,
                               mask = TRUE,
                               bridge_var = NULL,
                               fixed = TRUE,
                               fit_levels = FALSE,
                               fit_final = TRUE,
                               location_vine = c("random", "dissman"),
                               variable_vine = c("random", "dissman"),
                               seed = 123, cores = 1){
  
  d <- ncol(data)
  nlocs = nrows * ncols
  
  # Sample a bridging variable
  if (is.null(bridge_var)) {
    bridge_var = sample(1:nlocs, 1)
  } else if (!(bridge_var %in% c(1:nlocs))){
    stop(paste("Bridging variable must be between 1 and", nlocs))
  }
  
  output <- list()
  output[["bridge_var"]] <- bridge_var
  
  rvs_level1 = list()
  
  ### Level 1: Get spatial structure for each variable
  
  if (location_vine == "dissman") { 
    tree_alg_level1 = "mst_prim"
  } else{
    tree_alg_level1 = "random_weighted"
  }
  
  if (mask) {
    global_mask <- get_spatial_mask(nrows, ncols)
    
    # If fixed = TRUE, sample a variable to determine one structure for all variables/locations
    if (fixed){
      iVar = sample(1:nvars, 1)
      
      ktau_matrix <- wdm(data[, seq(iVar + (iVar-1) * (nlocs-1), iVar + iVar * (nlocs-1))], method = "kendall")
      weight_matrix <- global_mask * (1 - abs(ktau_matrix)) # to get maximum spanning tree
      g <- graph_from_adjacency_matrix(
        as.matrix(weight_matrix),
        mode = "undirected",
        weighted = TRUE, diag = FALSE)
      
      tree <- mst(g)
      rvs_level1_tmp <- spanning_tree_to_rvine_structure(tree)
      
      vine_level1_iVar = vinecop(pseudo_obs(data[, seq(iVar + (iVar-1) * (nlocs-1), iVar + iVar * (nlocs-1))]), trunc_lvl = 1,
                                 family_set = "tll", structure = rvs_level1_tmp, cores = 12, tree_algorithm = tree_alg_level1)
      
      # rvs_level1_iVar = vine_level1_iVar$structure
      # 
      # vine_level1_iVar = vinecop(pseudo_obs(data[, seq(iVar + (iVar-1) * (nlocs-1), iVar + iVar * (nlocs-1))]), trunc_lvl = 1,
      #                            family_set = "tll", structure = rvs_level1_iVar, cores = 12, tree_algorithm = tree_alg_level1)
      
      output[["vine_level1"]] = vine_level1_iVar
      output[["rvs_level1"]] = vine_level1_iVar$structure
      rvs_level1 = vine_level1_iVar$structure
      
    } else {
      for (iVar in 1:nvars){
        ktau_matrix <- wdm(data[, seq(iVar + (iVar-1) * (nlocs-1), iVar + iVar * (nlocs-1))], method = "kendall")
        weight_matrix <- global_mask * (1 - abs(ktau_matrix)) # to get maximum spanning tree
        g <- graph_from_adjacency_matrix(
          as.matrix(weight_matrix),
          mode = "undirected",
          weighted = TRUE, diag = FALSE)
        
        tree <- mst(g)
        rvs_level1_tmp <- spanning_tree_to_rvine_structure(tree)
        
        vine_level1_iVar = vinecop(pseudo_obs(data[, seq(iVar + (iVar-1) * (nlocs-1), iVar + iVar * (nlocs-1))]), trunc_lvl = 1,
                                   family_set = "tll", structure = rvs_level1_tmp, cores = 12, tree_algorithm = tree_alg_level1)
        
        rvs_level1_iVar = vine_level1_iVar$structure
        
        vine_level1_iVar = vinecop(pseudo_obs(data[, seq(iVar + (iVar-1) * (nlocs-1), iVar + iVar * (nlocs-1))]), trunc_lvl = 1,
                                   family_set = "tll", structure = rvs_level1_iVar, cores = 12, tree_algorithm = tree_alg_level1)
        
        output[[paste0("vine_level1_var", iVar)]] = vine_level1_iVar
        output[[paste0("rvs_level1_var", iVar)]] = rvs_level1_iVar
        rvs_level1[[iVar]] = rvs_level1_iVar
      }
    }
    
  } else {
    if (fixed){
      iVar = sample(1:nvars, 1)
      
      vine_level1_iVar = vinecop(pseudo_obs(data[,seq(iVar + (iVar-1) * (nlocs-1), iVar + iVar * (nlocs-1))]), trunc_lvl = 1,
                                 family_set = "tll", cores = cores, tree_algorithm = tree_alg_level1)
      
      rvs_level1_iVar = vine_level1_iVar$structure
      
      output[["vine_level1"]] = vine_level1_iVar
      output[["rvs_level1"]] = rvs_level1_iVar
      rvs_level1 = rvs_level1_iVar
      
    } else {
      for (iVar in 1:nvars){
        vine_level1_iVar = vinecop(pseudo_obs(data[,seq(iVar + (iVar-1) * (nlocs-1), iVar + iVar * (nlocs-1))]), trunc_lvl = 1,
                                   family_set = "tll", cores = cores, tree_algorithm = tree_alg_level1)
        
        rvs_level1_iVar = vine_level1_iVar$structure
        
        output[[paste0("vine_level1_var", iVar)]] = vine_level1_iVar
        output[[paste0("rvs_level1_var", iVar)]] = rvs_level1_iVar
        rvs_level1[[iVar]] = rvs_level1_iVar
      }
    }
  }
  
  ### Level 2: Get local (inter-variable) vine
  if (variable_vine == "dissman") { 
    tree_alg_level2 = "mst_prim"
  } else{
    tree_alg_level2 = "random_weighted"
  }
  
  # Determine inter-variable structure based on the bridging location
  vine_level2 = vinecop(pseudo_obs(data[, seq(bridge_var, d, nlocs)]), trunc_lvl = 1,
                        family_set = "tll", cores = cores, tree_algorithm = tree_alg_level2) 
  
  rvs_level2 = vine_level2$structure
  output[["rvs_level2"]] <- rvs_level2
  
  if (fit_levels) {
    output[["vine_level2"]] = vine_level2
  }
  
  ### Level 3: Merge the vine structures
  
  if (fixed){
    rvs_level3 <- merge_edges_fixed(rvs_level2, rvs_level1, bridge_var)
  } else {
    rvs_level3 <- merge_edges_individual(rvs_level1, rvs_level2, bridge_var)
  }
  
  output[["rvs_level3"]] <- rvs_level3
  
  # Reorder data accordingly
  if (fit_final) {
    
    # Fit the corresponding R vine
    vine_level3 <- vinecop(pseudo_obs(data),
                           family_set = "tll", trunc_lvl = 1,
                           structure = rvs_level3, cores = cores)
    
    output[["vine_level3"]] <- vine_level3
  }
  
  return(output)
}
