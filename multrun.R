# devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "1.1.0")
library(PoPS)
library(terra)
library(folderfun)

setff("In", "Q:/Shared drives/Data/Raster/Regional/SLF_1km/")
setff("cals", "Q:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/")

infected_file <- ffIn("slf_2020.tif")
host_file <- ffIn("toh.tif")
total_populations_file <- ffIn("all_plants.tif")
temp <- TRUE
temperature_coefficient_file <- ffIn("new_weather/temp5years.tif")
precip <- FALSE
precipitation_coefficient_file <- ""
model_type <- "SI"
latency_period <- 0
time_step <- "month"
season_month_start <- 5
season_month_end <- 11
start_date <- '2020-01-01'
end_date <- '2024-12-31'
use_lethal_temperature <- TRUE
temperature_file <- ffIn("new_weather/crit_temp5years.tif")
lethal_temperature <- -15
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
mortality_frequency <- 'year'
mortality_frequency_n <- 1
management <- FALSE
treatment_dates <- c('2020-12-24')
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
anthropogenic_dir <- "NONE"

number_of_iterations <- 100
number_of_cores <- 10

pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
random_seed <- NULL
output_frequency <- "year"
output_frequency_n <- 1
movements_file = ""
use_movements <- FALSE
start_exposed <- FALSE
generate_stochasticity <- TRUE
establishment_stochasticity <- TRUE
movement_stochasticity <- TRUE
deterministic <- FALSE
establishment_probability <- 0.5
dispersal_percentage <- 0.99
quarantine_areas_file <- ""
use_quarantine <- FALSE
use_spreadrates <- FALSE
use_overpopulation_movements <- TRUE
overpopulation_percentage <- 0.75
leaving_percentage <- 0.50
leaving_scale_coefficient <- 3
exposed_file <- ""
mask <- NULL
write_outputs = "summary_outputs"
output_folder_path = "C:/Users/cmjone25/Desktop/SLF_5year"

means <- read.csv(ffcals("2019_means.csv"))
means[5:6,] <- 0
parameter_means <- t(means)
parameter_means <- parameter_means[1,]
parameter_cov_matrix <- read.csv(ffcals("2019_cov_matrix.csv"))
parameter_cov_matrix[5:6, ] <- 0
parameter_cov_matrix[ , 5:6] <- 0

data <- pops_multirun(infected_file,
                      host_file,
                      total_populations_file,
                      parameter_means,
                      parameter_cov_matrix,
                      temp,
                      temperature_coefficient_file,
                      precip,
                      precipitation_coefficient_file,
                      model_type,
                      latency_period,
                      time_step,
                      season_month_start,
                      season_month_end,
                      start_date,
                      end_date,
                      use_lethal_temperature,
                      temperature_file,
                      lethal_temperature,
                      lethal_temperature_month,
                      mortality_on,
                      mortality_rate,
                      mortality_time_lag,
                      mortality_frequency,
                      mortality_frequency_n,
                      management,
                      treatment_dates,
                      treatments_file,
                      treatment_method,
                      natural_kernel_type,
                      anthropogenic_kernel_type,
                      natural_dir,
                      anthropogenic_dir,
                      number_of_iterations,
                      number_of_cores,
                      pesticide_duration,
                      pesticide_efficacy,
                      random_seed,
                      output_frequency,
                      output_frequency_n,
                      movements_file,
                      use_movements,
                      start_exposed,
                      generate_stochasticity,
                      establishment_stochasticity,
                      movement_stochasticity,
                      deterministic,
                      establishment_probability,
                      dispersal_percentage,
                      quarantine_areas_file,
                      use_quarantine,
                      use_spreadrates,
                      use_overpopulation_movements,
                      overpopulation_percentage,
                      leaving_percentage,
                      leaving_scale_coefficient,
                      exposed_file,
                      mask,
                      write_outputs,
                      output_folder_path)
