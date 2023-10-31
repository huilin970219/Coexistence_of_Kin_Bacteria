## This script was designed to achieve modeling the temperal dynamics among the three bacteria
## Author: Hui Lin; Donglin Wang
## Date: Nov. 2023

library(deSolve)
library(ggplot2)
library(phaseR)
library(rootSolve)

rm(list=ls())

#LV competitive model (dynamics of three member within community)


calculatePopulationChanges <- function(time, populations, params) {
  with(as.list(c(populations, params)), {
    
    calculateAlleeEffectImpact <- function(alle_effect_coeff, allee_effect) {
      alle_effect_coeff * allee_effect + 0.3
    }
    
    crossProtect <- function(death_rate_B_coeff, allee_effect, death_rate_B) {
      death_rate_B_coeff * allee_effect + death_rate_B
    }
    
    allee_effect_impact <- calculateAlleeEffectImpact(alle_effect_coeff, allee_effect)
    cross_protection_effect <- crossProtect(death_rate_B_coeff, allee_effect, death_rate_B)
    
    if(time < antibiotic_start_time | time > (antibiotic_start_time + 24)) {
      change_in_a <- growth_rate_A * a * (1 - a / carrying_capacity_A - 1.74 * growth_rate_B / growth_rate_A * b / carrying_capacity_A - 0.8 * allee_effect_impact / growth_rate_A * c / carrying_capacity_A)
      change_in_b <- growth_rate_B * b * (1 - b / carrying_capacity_B - 1.9 * growth_rate_A / growth_rate_B * a / carrying_capacity_B - 1.2 * allee_effect_impact / growth_rate_B * c / carrying_capacity_B)
      change_in_c <- allee_effect_impact * c * ((1 - c / carrying_capacity_C) - 0.8 * growth_rate_A / allee_effect_impact * a / carrying_capacity_C - 1.2 * growth_rate_B / allee_effect_impact * b / carrying_capacity_C)
    } else {
      change_in_a <- growth_rate_A * a * (1 - a / carrying_capacity_A - 1.74 * growth_rate_B / growth_rate_A * b / carrying_capacity_A - 0.8 * allee_effect_impact / growth_rate_A * c / carrying_capacity_A) - death_rate_A * a
      change_in_b <- growth_rate_B * b * (1 - b / carrying_capacity_B - 1.9 * growth_rate_A / growth_rate_B * a / carrying_capacity_B - 1.2 * allee_effect_impact / growth_rate_B * c / carrying_capacity_B) - cross_protection_effect * b
      change_in_c <- allee_effect_impact * c * ((1 - c / carrying_capacity_C) - 0.8 * growth_rate_A / allee_effect_impact * a / carrying_capacity_C - 1.2 * growth_rate_B / allee_effect_impact * b / carrying_capacity_C) - death_rate_C * c
    }
    list(c(change_in_a, change_in_b, change_in_c))
  })
}

# Parameters values
simulation_duration = 168  # Simulated hours
time_resolution = 1 / 100 # Time resolution of the output

parameters <- c(
  growth_rate_A = 0.5,
  growth_rate_B = 0.75,
  alle_effect_coeff = 0.07,
  allee_effect = 6,
  antibiotic_start_time = 72,
  death_rate_A = -0.6,
  death_rate_B = 5.5,
  death_rate_C = 0.8,
  death_rate_B_coeff = -0.62,
  carrying_capacity_A = 1.0,
  carrying_capacity_B = 1.0,
  carrying_capacity_C = 1.0
)

initial_populations <- c(a = 0.3, b = 0.4, c = 0.3)

# Run a simulation and plot the time series of species abundances
population_simulation <- ode(
  y = initial_populations,
  times = seq(0, simulation_duration, by = time_resolution),
  func = calculatePopulationChanges,
  parms = parameters
)


