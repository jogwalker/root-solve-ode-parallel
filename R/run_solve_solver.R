
args <- commandArgs(T)
outdir <- args[1]
jid <- as.numeric(args[2])

library(deSolve)

# load the model
source("/clusterdata/uqgheman/jummy/qualitative-modeling/R/model.R")

# read in the data
load("/clusterdata/uqgheman/jummy/qualitative-modeling/R/NoNegs.Rdata")
this <- results.winners[[jid]]
seed <- this[[1]]
runSim.result <- this[[2]]
params <- this[[3]]

inits <- rep(1,length(runSim.result$y))
t <- seq(1,10000,1)
K.values <- c(10,10)
parameters <- c(params$value,K.values)

# Run analysis
results <- ode(inits,t,runModel,parameters)
output <- list(seed,runSim.result,params,results)

# write output file
filename <- paste("runODE_",seed,".Rdata",sep="")
save(output, file=filename)
