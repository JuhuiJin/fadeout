######################################
## Solve NA's 
######################################

process_nas <- function(df, vars) {
  df <- df %>%
    dplyr::mutate(
      across(
        all_of(vars), 
        ~ if_else(.x <= 0, NA_real_, .x)
      )
    )
  
  return(df)
}

factorize <- function(df, vars) {
  df <- df %>%
    dplyr::mutate(
      across(
        all_of(vars), 
        ~ as.factor(.x)
      )
    )
}

unlabel <- function(df, vars) {
  df <- df %>%
    dplyr::mutate(
      across(
        all_of(vars), 
        ~ zap_labels(.x)
      )
    )
}

######################################
## Standardize Scores 
######################################

standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}


######################################
## Instrumental Variables 
######################################

# Feed in df, outcome var, w, aux_vars 

run_iv <- function(df, grades, w, aux_vars) {
  iv_coefs <- data.frame(
    grade = grades, 
    exp_estimates = NA
  )
  
  for (g in grades) {
    outcome_var <- paste0("z_score", g)
    
    formula_str <- paste0(
      outcome_var, 
      "~", paste(c(w, aux_vars), collapse = "+"), "| interval +" , paste(aux_vars, collapse = "+")
    )
    
    model_iv <- ivreg(as.formula(formula_str), data = df)
    
    iv_coefs$exp_estimates[iv_coefs$grade == g] <- coef(model_iv)[[w]]
    
    for (var in c(aux_vars)) {
      iv_coefs[[paste0("coef_", var)]][iv_coefs$grade == g] <- coef(model_iv)[[paste(var)]]
    }
  }
  
  return(iv_coefs)
}

######################################
## Ordinary Least Squares 
######################################

run_ols_estimates <- function(df, grades, w, aux_vars = NULL){
  df <- df %>% 
    filter(if_all(all_of(c(w, aux_vars)), ~ !is.na(.)))
  
  ols_coefs <- data.frame(
    grade = grades
  )
  
  for (g in grades) {
    formula_str <- paste0("z_score", g, "~", w)
    if (length(aux_vars) > 0) {
      formula_str <- paste0("z_score", g, "~", w, "+", paste0(aux_vars, collapse = "+"))
    }
    formula <- as.formula(formula_str)
    
    mod <- feols(formula, data = df)
    
    ols_coefs$ols_estimates[ols_coefs$grade == g] <- coef(mod)[[w]]
    
    for (var in aux_vars) {
      ols_coefs[[paste0("coef_", var)]][ols_coefs$grade == g] <- coef(mod)[[paste(var)]]
    }
  }
  
  return(ols_coefs)
}

######################################
## STAR Estimates 
######################################

# Feed in df, grades, and aux vars to get experimental estimates 

run_star_estimates <- function(star_data, intermediate_grades, w, aux_vars = NULL) {
  # Remove rows with NA in relevant columns
  relevant_cols <- c(w, "newsch", aux_vars)
  star_data <- star_data %>% 
    filter(if_all(all_of(relevant_cols), ~ !is.na(.)))
  # Prepare auxiliary variables for estimation by removing "miss_lunch_new"
  aux_vars_for_estimation <- setdiff(aux_vars, "miss_lunch_new")
  # Initialize coefficients dataframe with treatment and auxiliary variables
  star_coefs <- data.frame(grade = intermediate_grades,
                           star_estimates = NA)
  for (var in aux_vars) {
    star_coefs[[paste0("coef_", var)]] <- NA
  }
  for (g in intermediate_grades) {
    col_name <- paste0("avsum", g)
    star_data[[col_name]] <- standardize(star_data[[col_name]])
    # Construct formula with auxiliary variables (excluding "miss_lunch_new")
    formula_str <- paste0(col_name, " ~", w)
    if (length(aux_vars_for_estimation) > 0) {
      formula_str <- paste0(formula_str, " + ", paste(aux_vars_for_estimation, collapse = " + "))
    }
    formula_str <- paste0(formula_str, " | newsch")
    # Estimate the model using feols
    model <- feols(as.formula(formula_str), data = star_data)
    # Store coefficients: treatment and auxiliary variables from estimation.
    star_coefs$star_estimates[star_coefs$grade == g] <- coef(model)[w]
    for (var in aux_vars_for_estimation) {
      star_coefs[[paste0("coef_", var)]][star_coefs$grade == g] <- coef(model)[var]
    }
    # For miss_lunch_new, if present in aux_vars, set its coefficient to 0.
    if ("miss_lunch_new" %in% aux_vars) {
      star_coefs[["coef_miss_lunch_new"]][star_coefs$grade == g] <- 0
    }
  }
  return(star_coefs)
}


run_star_new <- function(star_data, intermediate_grades, w, aux_vars = NULL) {
  # Remove rows with NA in relevant columns
  relevant_cols <- c(w, "newsch", aux_vars)
  star_data <- star_data %>% 
    filter(if_all(all_of(relevant_cols), ~ !is.na(.)))
  # Prepare auxiliary variables for estimation by removing "miss_lunch_new"
  aux_vars_for_estimation <- setdiff(aux_vars, "miss_lunch_new")
  # Initialize coefficients dataframe with treatment and auxiliary variables
  star_coefs <- data.frame(grade = intermediate_grades,
                           star_estimates = NA)
  for (var in aux_vars) {
    star_coefs[[paste0("coef_", var)]] <- NA
  }
  for (g in intermediate_grades) {
    col_name <- paste0("mathz", g)
    star_data[[col_name]] <- standardize(star_data[[col_name]])
    # Construct formula with auxiliary variables (excluding "miss_lunch_new")
    formula_str <- paste0(col_name, " ~", w)
    if (length(aux_vars_for_estimation) > 0) {
      formula_str <- paste0(formula_str, " + ", paste(aux_vars_for_estimation, collapse = " + "))
    }
    formula_str <- paste0(formula_str, " | newsch")
    # Estimate the model using feols
    model <- feols(as.formula(formula_str), data = star_data)
    # Store coefficients: treatment and auxiliary variables from estimation.
    star_coefs$star_estimates[star_coefs$grade == g] <- coef(model)[w]
    for (var in aux_vars_for_estimation) {
      star_coefs[[paste0("coef_", var)]][star_coefs$grade == g] <- coef(model)[var]
    }
    # For miss_lunch_new, if present in aux_vars, set its coefficient to 0.
    if ("miss_lunch_new" %in% aux_vars) {
      star_coefs[["coef_miss_lunch_new"]][star_coefs$grade == g] <- 0
    }
  }
  return(star_coefs)
}

######################################
## ESC Estimates  
######################################

# Feed in df, exp coefs, grades, y_s, w, aux_vars 

run_lu_estimates <- function(nyc_data, star_coefs, grades, y_s, w, aux_vars = NULL){
  if (length(aux_vars) > 0) {
    nyc_data <- nyc_data %>%
      drop_na(all_of(aux_vars))
  }
  
  # Extract the coefficients for W --> 3rd grade scores 
  coefs <- as.matrix(star_coefs[1,-1])
  
  # Define df for storing esc estimates 
  esc_coefs <- data.frame(
    grade = grades, 
    esc_estimates = NA_real_
  )
  
  formula_str <- paste0(y_s, "~", w)
  if (length(aux_vars) > 0) {
    formula_str <- paste0(y_s, "~", w, "+", paste(aux_vars, collapse = "+"))
  }
  formula <- as.formula(formula_str)
  
  # Reshape NYC Data to include the desired aux vars 
  X <- model.matrix(
    formula, 
    data = nyc_data
  )
  X <- X[, -1]
  
  # Create predicted values and the esc control term 
  nyc_data[[paste0("y_hat")]] <- as.numeric((X) %*% t(coefs))
  nyc_data[[paste0("resid")]] <- nyc_data[[paste0(y_s)]] - nyc_data[[paste0("y_hat")]] 
  
  
  for (g in grades) {
    formula_str <- paste0("z_score", g, "~", w)
    # Append the selection control term to the regression formula 
    if (length(aux_vars) > 0) {
      formula_str <- paste0("z_score", g, "~", w, "+", paste0(aux_vars, collapse = "+"))
    }
    formula <- as.formula(paste0(formula_str, " + resid"))
    
    mod <- lm(formula, data = nyc_data)
    esc_coefs$esc_estimates[esc_coefs$grade == g] <- coef(mod)[[w]]
    esc_coefs$resid[esc_coefs$grade == g] <- coef(mod)[["resid"]]
  }
  
  return(esc_coefs)
}

run_lu_estimates <- function(nyc_data, star_coefs, grades, y_s, w, aux_vars = NULL){
  if (length(aux_vars) > 0) {
    nyc_data <- nyc_data %>%
      drop_na(all_of(aux_vars))
  }
  
  # Extract the coefficients for W --> 3rd grade scores 
  coefs <- as.matrix(star_coefs[1,-1])
  
  # Define df for storing esc estimates 
  esc_coefs <- data.frame(
    grade = grades, 
    esc_estimates = NA_real_
  )
  
  formula_str <- paste0(y_s, "~", w)
  if (length(aux_vars) > 0) {
    formula_str <- paste0(y_s, "~", w, "+", paste(aux_vars, collapse = "+"))
  }
  formula <- as.formula(formula_str)
  
  # Reshape NYC Data to include the desired aux vars 
  X <- model.matrix(
    formula, 
    data = nyc_data
  )
  X <- X[, -1]
  
  # Create predicted values and the esc control term 
  nyc_data[[paste0("y_hat")]] <- as.numeric((X) %*% t(coefs))
  nyc_data[[paste0("resid")]] <- nyc_data[[paste0(y_s)]] - nyc_data[[paste0("y_hat")]] 
  
  
  for (g in grades) {
    formula_str <- paste0("z_score", g, "~", w)
    # Append the selection control term to the regression formula 
    if (length(aux_vars) > 0) {
      formula_str <- paste0("z_score", g, "~", w, "+", paste0(aux_vars, collapse = "+"))
    }
    formula <- as.formula(paste0(formula_str, " + resid"))
    
    mod <- lm(formula, data = nyc_data)
    esc_coefs$esc_estimates[esc_coefs$grade == g] <- coef(mod)[[w]]
    esc_coefs$resid[esc_coefs$grade == g] <- coef(mod)[["resid"]]
  }
  
  return(esc_coefs)
}

######################################
## Surrogate Estimates
######################################

run_surr_estimates <- function(df, coefs, grades, y_s, w, aux_vars = NULL){
  if (length(aux_vars) > 0) {
    df <- df %>%
      drop_na(all_of(aux_vars))
  }
  
  coefs <- as.matrix(coefs[1,-1])
  
  surr_coefs <- data.frame(
    grade = grades
  )
  
  for (g in grades) {
    # Regress z-score for each grade on y_s and aux_vars 
    formula <- as.formula(paste0(paste0("z_score", g), "~", paste0(c(y_s, aux_vars), collapse = "+")))
    
    model <- lm(formula, data = df)
    
    beta <- coef(model)[[y_s]]
    gamma <- coefs[1]
    
    surr_coefs$surr_estimates[surr_coefs$grade == g] <- beta*gamma
  }
  
  return(surr_coefs)
}


#############################


