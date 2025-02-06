# Load libraries
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Define underwriting years
uw_years <- 2003:2023

# Define maximum development years
max_dev_year <- 21

# Function to generate cumulative claims for a single underwriting year
generate_claims <- function(uw_year, max_dev_year) {
  # Starting number of claims for dev_year 1
  base_claims <- sample(10:50, 1)

  # Define a growth rate for development years (10% to 30%)
  growth_rate <- runif(1, 0.1, 0.3)

  # Generate claim numbers for each development year
  claim_numbers <- cumsum(base_claims * (1 + growth_rate) ^ (1:max_dev_year))

  # Create a dataframe for the underwriting year
  data.frame(
    uw_year = uw_year,
    dev_year = 1:max_dev_year,
    claim_number = round(claim_numbers)
  )
}

# Generate the claim development triangle by applying the function to all underwriting years
triangle <- bind_rows(lapply(uw_years, generate_claims, max_dev_year = max_dev_year))

# Sort
triangle <- triangle %>%
  arrange(uw_year, dev_year)

# Remove the bottom right of the triangle

triangle_data <- triangle %>%
  filter(uw_year + dev_year <= 2024)

# output
setwd("data-raw")
usethis::use_data(triangle_data, overwrite = TRUE)

