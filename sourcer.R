# scripts to source
source("utility_script.R")
source("plotting_utility.R")

# required packages
packs <- c("dplyr", "reshape2", "runjags", "MCMCpack", "mcmcplots",
           'parallel', "coda")

package_load(packs)