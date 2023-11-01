library(deSolve)
library(rootSolve)

rm(list=ls())

#LV competitive model (dynamics of three member within community)
# Function to calculate population changes over time

calculatePopulationChanges <- function(time, populations, params) {
  with(as.list(c(populations, params)), {
    
    # Function to calculate the impact of the Allee effect
    calculateAlleeEffectImpact <- function(alle_effect_coeff, allee_effect, base) {
      return(alle_effect_coeff * allee_effect + base) # Returns the calculated impact
    }
    
    # Function to calculate the cross protection effect
    crossProtect <- function(death_rate_B_coeff, allee_effect, death_rate_B) {
      result <- death_rate_B_coeff * allee_effect + death_rate_B # Returns the cross protection effect
      return(result)
    }
    
    # Calculating the impact of the Allee effect
    allee_effect_impact <- calculateAlleeEffectImpact(alle_effect_coeff, allee_effect, allee_effect_base)
    # Calculating the cross protection effect
    cross_protection_effect <- crossProtect(death_rate_B_coeff, allee_effect, death_rate_B)
    
    # General function for interaction between populations
    interaction_function <- function(rate_1, rate_2, capacity_1, capacity_2, population_1, population_2, coeff_1_2, effect_impact = 0, coeff_1_effect = 0) {
      return(rate_1 * population_1 * (1 - population_1 / capacity_1 - coeff_1_2 * rate_2 / rate_1 * population_2 / capacity_1 - coeff_1_effect * effect_impact / rate_1 * population_1 / capacity_1))
    }
    
    # Calculate changes in populations
    if(time < antibiotic_start_time || time > (antibiotic_start_time + 24)) {
      # Without antibiotic effect
      change_in_a <- interaction_function(growth_rate_A, growth_rate_B, carrying_capacity_A, carrying_capacity_B, a, b, interaction_coeff_A_B, allee_effect_impact, interaction_coeff_A_C)
      change_in_b <- interaction_function(growth_rate_B, growth_rate_A, carrying_capacity_B, carrying_capacity_A, b, a, interaction_coeff_B_A, allee_effect_impact, interaction_coeff_B_C)
      change_in_c <- interaction_function(allee_effect_impact, growth_rate_B, carrying_capacity_C, carrying_capacity_B, c, b, interaction_coeff_C_B, growth_rate_A, interaction_coeff_C_A)
    } else {
      # With antibiotic effect
      change_in_a <- interaction_function(growth_rate_A, growth_rate_B, carrying_capacity_A, carrying_capacity_B, a, b, interaction_coeff_A_B, allee_effect_impact, interaction_coeff_A_C) - death_rate_A * a
      change_in_b <- interaction_function(growth_rate_B, growth_rate_A, carrying_capacity_B, carrying_capacity_A, b, a, interaction_coeff_B_A, allee_effect_impact, interaction_coeff_B_C) - cross_protection_effect * b
      change_in_c <- interaction_function(allee_effect_impact, growth_rate_B, carrying_capacity_C, carrying_capacity_B, c, b, interaction_coeff_C_B, growth_rate_A, interaction_coeff_C_A) - death_rate_C * c
    }
    
    return(list(c(change_in_a, change_in_b, change_in_c))) # Return changes in populations
  })
}


# Function to initialize growth-related parameters
initializeGrowthParameters <- function() {
  return(list(
    growth_rate_A = 0.5,
    growth_rate_B = 0.75,
    carrying_capacity_A = 1.0,
    carrying_capacity_B = 1.0,
    carrying_capacity_C = 1.0
  ))
}

# Function to initialize death-related parameters
initializeDeathParameters <- function() {
  return(list(
    death_rate_A = -0.6,
    death_rate_B = 5.5,
    death_rate_C = 0.8,
    death_rate_B_coeff = -0.62
  ))
}

# Function to initialize interaction-related parameters
initializeInteractionParameters <- function() {
  return(list(
    interaction_coeff_A_B = 1.74,
    interaction_coeff_A_C = 0.8,
    interaction_coeff_B_A = 1.9,
    interaction_coeff_B_C = 1.2,
    interaction_coeff_C_A = 0.8,
    interaction_coeff_C_B = 1.2
  ))
}

# Function to initialize Allee effect parameters
initializeAlleeParameters <- function() {
  return(list(
    alle_effect_coeff = 0.07,
    allee_effect = 6,
    allee_effect_base = 0.3
  ))
}

# Function to initialize miscellaneous parameters
initializeMiscellaneousParameters <- function() {
  return(list(
    antibiotic_start_time = 72
  ))
}

# Combine all parameters
parameters <- c(
  initializeGrowthParameters(),
  initializeDeathParameters(),
  initializeInteractionParameters(),
  initializeAlleeParameters(),
  initializeMiscellaneousParameters()
)

# Initial populations
initial_populations <- c(a = 1, b = 1, c = 1)

# Duration and time resolution of the simulation
simulation_duration = 168 
time_resolution = 1 / 100

# Running the simulation and storing the results
population_simulation <- ode(
  y = initial_populations,
  times = seq(0, simulation_duration, by = time_resolution),
  func = calculatePopulationChanges,
  parms = parameters
)


