library(dplyr)
library(purrr)
library(tidyr)

# parameters

poisson_mean <- 5
lognorm_param1 <- 16
lognorm_param2 <- 1.5
num_sims <- 20000
territory_perils <- c("US_W", "US_Q", "JAP_W", "JAP_Q")

# create dataset with monte carlo sims from frequency / severity params above
# and randomly assign a territory peril

losses <- tibble(num = rpois(num_sims, poisson_mean)) |> 
  mutate(year = 1:nrow(pick(num))) |> 
  mutate(year = map2(year, num, function(x, y) rep(x, y))) |> 
  unnest(year) |> 
  select(year) |> 
  rowwise() |> 
  mutate(amount = rlnorm(1, lognorm_param1, lognorm_param2)) |> 
  mutate(territory_peril = sample(territory_perils, nrow(pick(year)), replace = TRUE)) |> 
  ungroup()

# output
setwd("data-raw")
usethis::use_data(losses, overwrite = TRUE)
