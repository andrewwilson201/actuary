
library(tidyverse)

# set seed for reproducibility
set.seed(123)

# define the number of accident and development periods
n_acc <- 10
n_dev <- 10

# base parameters for accident-specific decay:
# k_base is the decay parameter for the first accident year,
# delta is the amount by which k increases for each subsequent accident year.
k_base <- 0.5
delta <- 0.1

# generate triangle data
triangle_df <- tibble(
  accident = 1:n_acc,
  base_claim = runif(n_acc, min = 100, max = 500)
) |>
  # each accident year gets its own decay parameter.
  mutate(accident_k = k_base + (accident - 1) * delta)

# generate a common first development factor (for period 1) between 4 and 5.
first_growth <- runif(1, min = 4, max = 5)

# create the triangle
triangle_long <- triangle_df |>
  # Cceate all accident-development combinations.
  crossing(tibble(development = 1:n_dev)) |>
  arrange(accident, development) |>
  group_by(accident) |>
  # for each accident year, compute:
  # - the incremental development factor for each development period, and
  # - the cumulative development factor (as the cumulative product).
  mutate(dev_factor = 1 + (first_growth - 1) * exp(-accident_k * (development - 1)),
         cum_factor = cumprod(dev_factor),
         value = base_claim * cum_factor,
         # only the first (n_dev - accident + 1) development periods are observed.
         value = if_else(development > (n_dev - accident + 1), NA_real_, value)) |>
  ungroup() |>
  # filter out the bottom right
  filter(!is.na(value))

# rename columns and remove some
triangle_data <- triangle_long |>
  transmute(uw_year = accident,
            dev_year = development,
            claim_number = value)

# output
setwd("data-raw")
usethis::use_data(triangle_data, overwrite = TRUE)

