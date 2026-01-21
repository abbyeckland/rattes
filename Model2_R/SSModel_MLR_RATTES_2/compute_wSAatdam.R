# Import packages

library(dplyr)

# update version of code to save to files being saved:
version <- "010626" 


### IMPORT DATA

# set wd to use in relative paths
setwd(".")

# import SA at dam data for each year, 1699-present
SA_data <- read.csv("SYfiles_Jan2026/MLR_SAatDamInitial122925.csv", header=TRUE, stringsAsFactors = FALSE)
colnames(SA_data) <- sub("^X", "", colnames(SA_data))
colnames(SA_data)[1] <- "SID"  # Rename the first column to 'SID'
print(colnames(SA_data))

# import SS data to get yrp
ss_data <- read.csv("SYfiles_Jan2026/MLR_Input122925.csv", header=TRUE, stringsAsFactors = FALSE)
#ss_data <- ss_data %>% select(-elev) # remove this line once rerun all code, shouldn't conflict anymore

# merge yrp column to SA_data
column_names <- c("SID", "yrp")
# Subset the data frame by column names
yrp <- ss_data[, column_names, drop=FALSE]
SA_data <- merge(yrp, SA_data, by="SID", all=FALSE)


# Function to check for NaNs, NAs, or negative values in a matrix
check_invalid_values_matrix <- function(matrix) {
  any_nan <- any(is.nan(matrix))
  any_na <- any(is.na(matrix))
  any_negative <- any(matrix < 0, na.rm = TRUE)
  
  if (any_nan || any_na || any_negative) {
    cat("Invalid values detected in matrix:\n")
    if (any_nan) cat("- NaN values present.\n")
    if (any_na) cat("- NA values present.\n")
    if (any_negative) cat("- Negative values present.\n")
  } else {
    cat("No NaN, NA, or negative values detected in matrix.\n")
  }
}

# Convert SA_data to matrix, excluding SID and yrp columns
SA_matrix <- as.matrix(SA_data[, -c(1, 2)])
yrp_values <- SA_data$yrp

# Check for invalid values in the matrix
check_invalid_values_matrix(SA_matrix)


### COMPUTE WSAATDAM FOR ALL SITES

# Extract column years
col_years <- as.numeric(colnames(SA_matrix))

# Initialize cumulative mean matrix with same dimensions as SA_matrix
cumulative_mean_matrix <- SA_matrix

# Loop to compute time-weighted wSAatdam
for (i in 1:nrow(SA_matrix)) {
  yrc_value <- yrp_values[i]
  
  # Identify columns that are within or after the 'yrc' value
  valid_cols <- col_years >= yrc_value
  
  # Extract valid values for calculation
  valid_values <- SA_matrix[i, valid_cols]
  
  # Calculate cumulative mean for valid values
  if (length(valid_values) > 0) {
    cumulative_means <- cumsum(valid_values) / seq_along(valid_values)
    
    # Replace values after yrc with cumulative means in the result matrix
    all_values <- SA_matrix[i, ]
    all_values[valid_cols] <- cumulative_means
    
    cumulative_mean_matrix[i, ] <- all_values
  }
}

# Check for invalid values in the result matrix
check_invalid_values_matrix(cumulative_mean_matrix)

# Add back SID and yrp columns for saving to file
cumulative_mean_df <- data.frame(SID = SA_data$SID, yrp = yrp_values, cumulative_mean_matrix)

# Save transformed data to folder
file_path <- paste0("wSAatdam_predSites_notTransformed_", version, ".csv")
write.csv(cumulative_mean_df, file = file_path, row.names = FALSE)

