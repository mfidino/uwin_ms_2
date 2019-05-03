
# scripts to source
source("utility_script.R")
source("plotting_utility.R")
source("prepare_data.R")

# required packages
packs <- c("dplyr", "reshape2", "runjags", "MCMCpack", "mcmcplots",
           'parallel', "coda", "stringr", "data.table", "DMwR")

package_load(packs)
