
args <- commandArgs(T)
outdir <- args[1]
jid <- as.numeric(args[2])

library(deSolve)

# load the model
source("/clusterdata/uqgheman/jummy/root-solve-ode-parallel/R/model.R")

# read in the data
load("/clusterdata/uqgheman/jummy/root-solve-ode-parallel/R/NoNegs.Rdata")
this <- results.winners[[jid]]
seed <- this[[1]]
runSim.result <- this[[2]]
params <- this[[3]]

inits <- rep(1,length(runSim.result$y))
t <- seq(1,10000,0.1)
K.values <- c(10,10)
parameters <- c(params$value,K.values)

# Run analysis
results <- ode(inits,t,runModel,parameters)
output <- list(seed,runSim.result,params,results)

# write output file
filename <- paste(outdir,"runODE_t2_",seed,".Rdata",sep="")
save(output, file=filename)
