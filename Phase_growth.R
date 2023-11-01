library(deSolve)
library(phaseR)
library(rootSolve)

rm(list=ls())


#LV competitive model (dynamics of three member within community)

calculatePopulationChanges <- function(time, populations, params) {
  with(as.list(c(populations, params)), {
    
    calculateAlleeEffectImpact <- function(alle_effect_coeff, allee_effect) {
      alle_effect_coeff * allee_effect + allee_effect_base
    }
    
    crossProtect <- function(death_rate_B_coeff, allee_effect, death_rate_B) {
      death_rate_B_coeff * allee_effect + death_rate_B
    }
    
    allee_effect_impact <- calculateAlleeEffectImpact(alle_effect_coeff, allee_effect)
    cross_protection_effect <- crossProtect(death_rate_B_coeff, allee_effect, death_rate_B)
    
    if(time < antibiotic_start_time || time > (antibiotic_start_time + 24)) {
      change_in_a <- growth_rate_A * a * (1 - a / carrying_capacity_A - interaction_coeff_A_B * growth_rate_B / growth_rate_A * b / carrying_capacity_A - interaction_coeff_A_C * allee_effect_impact / growth_rate_A * c / carrying_capacity_A)
      change_in_b <- growth_rate_B * b * (1 - b / carrying_capacity_B - interaction_coeff_B_A * growth_rate_A / growth_rate_B * a / carrying_capacity_B - interaction_coeff_B_C * allee_effect_impact / growth_rate_B * c / carrying_capacity_B)
      change_in_c <- allee_effect_impact * c * ((1 - c / carrying_capacity_C) - interaction_coeff_C_A * growth_rate_A / allee_effect_impact * a / carrying_capacity_C - interaction_coeff_C_B * growth_rate_B / allee_effect_impact * b / carrying_capacity_C)
    } else {
      change_in_a <- growth_rate_A * a * (1 - a / carrying_capacity_A - interaction_coeff_A_B * growth_rate_B / growth_rate_A * b / carrying_capacity_A - interaction_coeff_A_C * allee_effect_impact / growth_rate_A * c / carrying_capacity_A) - death_rate_A * a
      change_in_b <- growth_rate_B * b * (1 - b / carrying_capacity_B - interaction_coeff_B_A * growth_rate_A / growth_rate_B * a / carrying_capacity_B - interaction_coeff_B_C * allee_effect_impact / growth_rate_B * c / carrying_capacity_B) - cross_protection_effect * b
      change_in_c <- allee_effect_impact * c * ((1 - c / carrying_capacity_C) - interaction_coeff_C_A * growth_rate_A / allee_effect_impact * a / carrying_capacity_C - interaction_coeff_C_B * growth_rate_B / allee_effect_impact * b / carrying_capacity_C) - death_rate_C * c
    }
    list(c(change_in_a, change_in_b, change_in_c))
  })
}

# Function to get growth-related parameters
getGrowthParams <- function() {
  list(
    growth_rate_A = 0.5,
    growth_rate_B = 0.75,
    carrying_capacity_A = 1.0,
    carrying_capacity_B = 1.0,
    carrying_capacity_C = 1.0
  )
}

# Function to get interaction coefficients
getInteractionCoeffs <- function() {
  list(
    interaction_coeff_A_B = 1.74,
    interaction_coeff_A_C = 0.8,
    interaction_coeff_B_A = 1.9,
    interaction_coeff_B_C = 1.2,
    interaction_coeff_C_A = 0.8,
    interaction_coeff_C_B = 1.2
  )
}

# Function to get allee effect-related parameters
getAlleeEffectParams <- function() {
  list(
    alle_effect_coeff = 0.07,
    allee_effect = 6,
    allee_effect_base = 0.3
  )
}

# Function to get antibiotic-related parameters
getAntibioticParams <- function() {
  list(
    antibiotic_start_time = 72
  )
}

# Function to get death rate-related parameters
getDeathRateParams <- function() {
  list(
    death_rate_A = -0.6,
    death_rate_B = 5.5,
    death_rate_C = 0.8,
    death_rate_B_coeff = -0.62
  )
}

# Combine all parameters into a single list
params <- c(
  getGrowthParams(),
  getInteractionCoeffs(),
  getAlleeEffectParams(),
  getAntibioticParams(),
  getDeathRateParams()
)



# Parameters values
initial_populations <- c(a = 1, b = 1, c = 1)
simulation_duration = 168 
time_resolution = 1 / 100

# Data frame to store the results
phaseX <- data.frame(g_ratio = rep(1.0, 500), d_ratio = rep(1.0, 500), xf = rep(1.0, 500), yf = rep(1.0, 500), zf = rep(1.0, 500))

# Simulations to generate the phase space

Nsims = 500


updateParamA <- function(index, totalSims) {
  return(1 + (index - 1) * 6 / totalSims)
}


updateParamDa3 <- function(index, totalSims) {
  return(-1 + (index - 1) * (-1) / totalSims)
}

# Function to perform a single simulation
performSimulation <- function(initialPopulations, simulationDuration, timeResolution, params) {
  ode(
    y = initialPopulations,
    times = seq(0, simulationDuration, by = timeResolution),
    func = calculatePopulationChanges,
    parms = params
  )
}

# Function to calculate fraction
calculateFraction <- function(numerator, total) {
  return(numerator / total)
}

for (i in 1:Nsims) {
  print(i)  # Print the current simulation number
  params$A <- updateParamA(i, Nsims)
  for (j in 1:Nsims) {
    params$da3 <- updateParamDa3(j, Nsims)
    LV.simulation <- performSimulation(initial_populations, simulation_duration, time_resolution, params)
    finalState <- LV.simulation[nrow(LV.simulation), 2:4]
    phaseX[i * Nsims + j - Nsims,] <- c(params$A, params$da3, finalState)
  }
}

# Calculating fractions for phaseX, phaseY and phaseZ
phaseX$Fraction1 <- apply(phaseX[, c("xf", "yf", "zf")], 1, function(row) calculateFraction(row["xf"], sum(row)))

phaseY$Fraction1 <- apply(phaseY[, c("xf", "yf", "zf")], 1, function(row) calculateFraction(row["yf"], sum(row)))

phaseZ$Fraction1 <- apply(phaseZ[, c("xf", "yf", "zf")], 1, function(row) calculateFraction(row["zf"], sum(row)))



