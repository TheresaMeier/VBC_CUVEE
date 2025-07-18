


get_GAM_decomposition = function(mp, mc, rc, 
                                 locs, 
                                 time_c, time_p,
                                 var_names,
                                 families = NULL) {
  
  # Step 1: Transform inputs to wide format tibbles for mc, mp, rc
  dfs <- list(
    "mc" = transform_to_wide_format(mc, locs, var_names, time_c),
    "mp" = transform_to_wide_format(mp, locs, var_names, time_p),
    "rc" = transform_to_wide_format(rc, locs, var_names, time_c)
  )
  
  # Step 2: Fit GAM models for each variable on each dataset
  gam_list <- map(dfs, function(df) {
    map(set_names(var_names), function(var) {
      select_best_gam(var, df, families)
    })
  })
  
  # Step 3: Add seasonality and remainder columns to each dataset
  dfs <- map2(dfs, gam_list, function(df, fit_list) {
    for (var in names(fit_list)) {
      preds <- predict(fit_list[[var]], type = "response")
      df[[paste0("seasonality_", var)]] <- preds
      df[[paste0("remainder_", var)]] <- df[[var]] - preds
    }
    df
  })
  
  # Step 4: Calculate seasonality predictions on mp data using rc and mc fits, plus delta seasonality
  # Initialize empty lists to store these seasonality predictions
  seasonality_rc <- vector("list", length(var_names))
  names(seasonality_rc) <- var_names
  
  seasonality_mc <- vector("list", length(var_names))
  names(seasonality_mc) <- var_names
  
  seasonality_delta <- vector("list", length(var_names))
  names(seasonality_delta) <- var_names
  
  for (var in var_names) {
    # Predict seasonality on mp data using rc fit
    seasonality_rc[[var]] <- predict(gam_list$rc[[var]], newdata = dfs$mp, type = "response")
    
    # Predict seasonality on mp data using mc fit
    seasonality_mc[[var]] <- predict(gam_list$mc[[var]], newdata = dfs$mp, type = "response")
    
    # Calculate delta = mp seasonality (original) - predicted mc seasonality
    seasonality_delta[[var]] <- dfs$mp[[paste0("seasonality_", var)]] - seasonality_mc[[var]]
  }
  
  # Step 5: Combine the predictions for seasonality_rc and seasonality_delta into tibbles
  seasonality_rc_df <- dfs$mp %>% dplyr::select(Id, time)
  seasonality_delta_df <- dfs$mp %>% dplyr::select(Id, time)
  
  for (var in var_names) {
    seasonality_rc_df[[var]] <- seasonality_rc[[var]]
    seasonality_delta_df[[var]] <- seasonality_delta[[var]]
  }
  
  # Step 6: Add these new seasonality dataframes to dfs list before pivoting wider
  dfs[["seasonality_rc"]] <- seasonality_rc_df
  dfs[["seasonality_delta"]] <- seasonality_delta_df
  
  # Step 7: Transform all dfs back to wide format (including new seasonality ones)
  dfs_wide <- map(dfs, function(df) {
    # Remove Lat/Lon if present
    df2 <- df %>%  dplyr::select(-any_of(c("Lat", "Lon")))
    
    # Variables to spread (exclude Id and time)
    value_vars <- setdiff(names(df2), c("time", "Id"))
    
    df2 %>%
      pivot_wider(
        id_cols = time,
        names_from = Id,
        values_from = all_of(value_vars),
        names_sep = "."
      ) %>%
      dplyr::select(-time)
  })
  
  return(dfs_wide)
}



# Transform data into wide format (with locs being a df with columns Lat and Lon)
transform_to_wide_format <- function(data, locs, vars, time) {
  
  # Build list of relevant column names: time + all variable columns matching var patterns
  selected_cols <- unlist(lapply(vars, function(v) {
    grep(paste0("^", v, "\\."), names(data), value = TRUE)
  }))
  
  # Select only relevant columns
  data <- data %>%
    dplyr::select(all_of(selected_cols)) %>%
    dplyr::mutate(time = as.Date(time))
  
  # Pivot to long format: variable values go into columns named by `.value`, and location into 'Id'
  data_long <- data %>%
    pivot_longer(
      cols = -time,
      names_to = c(".value", "Id"),
      names_sep = "\\."
    ) %>%
    mutate(
      Id = as.integer(Id),
      t = lubridate::yday(time)  # Extract day-of-year
    ) %>%
    left_join(locs, by = "Id")
  
  return(data_long)
}


# Helper to suggest families based on variable properties
suggest_families <- function(var, df) {
  y <- df[[var]]
  
  # Get indicators
  is_positive <- all(y > 0, na.rm = TRUE)
  is_non_negative <- all(y >= 0, na.rm = TRUE)
  has_zeros <- (sum(y < 0.0001, na.rm = TRUE) / length(y)) > 0.2
  skewness_val <- e1071::skewness(y, na.rm = TRUE)
  is_count <- all(y >= 0, na.rm = TRUE) && all(y == floor(y), na.rm = TRUE)
  is_binary <- length(unique(y[!is.na(y)])) == 2
  is_bounded_01 <- all(y > 0 & y < 1, na.rm = TRUE)
  
  if (is_binary) {
    return(list(binomial = binomial()))
  } else if (is_count) {
    return(list(poisson = poisson()))
  } else if (!is_non_negative) {
    return(list(gaussian = gaussian()))
  } else if (has_zeros) {
    return(list(tweedie = tw(link = "log")))
  } else if (is_bounded_01) {
    return(list(betar = mgcv::betar()))
  } else if (is_positive && skewness_val > 1) {
    return(list(
      Gamma = Gamma(link = "log"),
      gaussian = gaussian()
    ))
  } else {
    return(list(gaussian = gaussian()))
  }
}


# Function to fit the best GAM model based on AIC
select_best_gam <- function(var, df, families) {
  
  if (is.null(families)){
    families_to_try <- suggest_families(var, df)
    message(sprintf(
      "Families to try for variable %s: %s", 
      var, 
      paste(names(families_to_try), collapse = ", ")
    ))
  } else {
    families_to_try = families[[var]]
  }
  
  models <- purrr::imap(families_to_try, function(fam, fam_name) {
    
    fit <- tryCatch(
      bam(
        formula = as.formula(paste0(var, " ~ te(t, Lat, Lon, bs = c('cc', 'tp', 'tp'))")),
        data = df,
        family = fam
      ),
      error = function(e) {
        message(sprintf("Failed to fit %s with family %s", var, fam_name))
        NULL
      }
    )
    
    if (!is.null(fit)) {
      aic_val <- tryCatch(AIC(fit), error = function(e) NA)
      list(model = fit, aic = aic_val, family = fam_name)
    } else {
      list(model = fit, family = fam_name)
      message(sprintf("Failed to fit %s with family %s", var, fam_name))
    }
  })
  
  models <- purrr::compact(models)
  
  if (length(families_to_try) > 1){
    if (length(models) == 0) {
      message(sprintf("No valid models fitted for variable '%s'", var))
      return(NULL)
    }
    
    best_model <- models[[which.min(map_dbl(models, "aic"))]]
    
    message(sprintf("Best family for '%s': %s (AIC = %.2f)", var, best_model$family, best_model$aic))
    
    return(best_model$model)
  } else {
    message(sprintf("Family for '%s': %s ", var, models[[1]]$family))
    
    return(models[[1]]$model)
  }
  
}
