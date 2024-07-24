# Clear environment
rm(list = ls())

# Load required libraries
library(tidyverse)
library(progress)

# Specify random number generation
set.seed(123)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Source model and results functions
source(file="./00_IBMfunctions.R")

# Read parameter values from calibrated model
traces <- read.csv(file = "./param_traces.csv")
traces <- traces[, c(2:10)]
colnames(traces) <- c("mort", "inf_out", "inf_min", "inf_sub", "min_out",
                      "min_sub", "sub_min", "sub_clin", "clin_sub")

# Set cohort parameter values 
nreps   <- 1000              # Number of repetitions
nsteps  <- 121               # Number of time steps (months)
nyears  <- (nsteps - 1)/12   # Number of time steps (years)
npeople <- 10000             # Number of people in each cohort
treatment <- 0.70            # Treatment rate (Removal of clinical disease)

# Run cohort simulation to endpoints
if(!dir.exists("data")) dir.create(file.path("data"))

if(!dir.exists("./data/IBM")) dir.create(file.path("./data/IBM"))
pb <- progress_bar$new(format = "[:bar] :percent", total = nreps)

for(i in 1:nreps){
   IBM_temp <- full_simulation_inf(nsteps, npeople, 1, traces, "random", treatment)
   save(IBM_temp, file = paste("./data/IBM/IBM_", i, ".Rdata", sep = ""))
   pb$tick()
}
cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Completed: Cohort simulation to endpoints\n"))


# Subset to those who develop any disease
if(!dir.exists("./data/disease")) dir.create(file.path("./data/disease"))
pb <- progress_bar$new(format = "[:bar] :percent", total = nreps)

for(i in 1:nreps){
   disease <- as.data.frame(array(0, dim = c(nsteps, 0)))
   load(paste("./data/IBM/IBM_", i, ".Rdata", sep = ""))
   for (j in 1:npeople){
      if (length(which(IBM_temp[, j] %in% c('m', 's'))) > 0)
         disease <- cbind(disease, IBM_temp[, j])
   }
   save(disease,file=paste("./data/disease/disease_", i, ".Rdata", sep = ""))
   rm("IBM_temp", "disease")
   pb$tick()
}
cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Completed: Subsetting to those who develop any disease\n"))

