
args <- commandArgs(T)
outdir <- args[1]
jid <- as.numeric(args[2])

# read in all the functions
source("~/functions.R")

# Load in the model 
source("~/makemodel.R")

# create parameters
set.seed(jid*123456)
parameters <- simParams(network)

# Run analysis
result <- runSim(inits=inits, timesteps=t,parameters=c(parameters$"value",K.values))

# write output file
save(list(result,parameters),file=paste(outdir,"results_",jid,".Rdata"))