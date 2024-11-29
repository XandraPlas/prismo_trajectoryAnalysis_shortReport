##----------------------------------------------------------------------------------------------- ##
#                             TRAJECTORIES DEPRESSION - FUNCTIONS
##----------------------------------------------------------------------------------------------- ##

# function to create dummy variables and scale numeric variables
to_dummy_and_scaled <- function(df_pred) {
  
  # Save outcome variable
  df_outcome <- df_pred %>% select(trajectory)
  
  # Create dummy variables
  df_pred$gender <- ifelse(df_pred$gender == "1", 0, 1) 
  
  df_pred <- df_pred %>%
    mutate(gender = as.integer(gender),
           rank_cat = factor(rank_cat),
           education_cat = factor(education_cat),
           work_function = factor(work_function),
           yr_deployment_cat = factor(yr_deployment_cat),
           Prev_deployment_dummy = as.integer(Prev_deployment_dummy)) 
  
  binary_cols <- df_pred %>%
    select(gender, Prev_deployment_dummy)
  
  dummy_cols <- df_pred %>%
    select(rank_cat, education_cat, work_function, yr_deployment_cat) %>%
    fastDummies::dummy_cols(remove_first_dummy = FALSE, remove_selected_columns = TRUE)  # Convert factors to dummies
  
  df_encoded <- cbind(binary_cols, dummy_cols)
  
  # Scale numerical variables
  # ----------------------------------------------- #
  df_scaled <- df_pred %>%
    select(-c(trajectory, gender, rank_cat, education_cat, work_function, yr_deployment_cat, Prev_deployment_dummy)) %>%
    scale(center = TRUE, scale = TRUE) %>%
    as.matrix()
  
  # Combine categorical with numerical
  # ----------------------------------------------- #
  df_x_prep <- cbind(df_encoded, df_scaled, df_outcome)
  
  return(df_x_prep)
}











# Function to calculate odds ratios
calculate_odds_ratios <- function(model, lambda) {
  # Extract coefficients at the selected lambda
  coeffs <- coef(model, lambda)
  odds_ratios_list <- list()
  
  # Loop through each class to compute odds ratios
  for (class_name in names(coeffs)) {
    # Get coefficients for the current class
    class_coeffs <- coeffs[[class_name]]
    variable_names <- rownames(class_coeffs)
    coefficients <- as.numeric(class_coeffs)
    
    # Compute odds ratios
    odds_ratios <- exp(coefficients)
    
    # Create interpretation
    interpretation <- ifelse(
      odds_ratios > 1, "Increases odds",
      ifelse(odds_ratios < 1, "Decreases odds", "No effect")
    )
    
    # Combine into a data frame
    class_df <- data.frame(
      Variable = variable_names,
      Coefficient = coefficients,
      Odds_Ratio = odds_ratios,
      Interpretation = interpretation,
      row.names = NULL
    )
    
    # Store the data frame in the list
    odds_ratios_list[[class_name]] <- class_df
  }
  
  return(odds_ratios_list)
}
