# Path to data storage 
external_path <- "/Users/juhui/Desktop/research/data/"

# Package Installments 
packages <- c("dplyr", "tidyr", "haven", "dplyr", "ggplot2", "fixest", "sn", "moments")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

check_na <- function(data, var_list) {
  # Pre-allocate the data frame with appropriate number of rows
  na_table <- data.frame(
    Variable = character(length = length(var_list)),
    NAs = numeric(length = length(var_list)),
    Not_NAs = numeric(length = length(var_list)),
    Total_n = numeric(length = length(var_list)),
    stringsAsFactors = FALSE
  )
  
  # Loop through each variable in the var_list
  for (i in seq_along(var_list)) {
    var <- var_list[i]
    na_count <- sum(is.na(data[[var]]))  # Count NA values
    not_na_count <- sum(!is.na(data[[var]]))  # Count non-NA values
    total_n <- nrow(data)  # Total rows in the data
    
    # Assign the results directly to the pre-allocated table
    na_table$Variable[i] <- var
    na_table$NAs[i] <- na_count
    na_table$Not_NAs[i] <- not_na_count
    na_table$Total_n[i] <- total_n
  }
  
  return(na_table)
}

