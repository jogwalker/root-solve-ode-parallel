
args <- commandArgs(T)
outdir <- args[1]
jid <- as.numeric(args[2])

# read in all the functions
source("/clusterdata/uqgheman/jummy/qualitative-modeling/R/functions.R")

# Load in the model 
source("/clusterdata/uqgheman/jummy/qualitative-modeling/R/makemodel.R")

# create parameters
set.seed(jid*123456)
parameters <- simParams(network)

# Run analysis
result <- runSim(inits=inits, timesteps=t,parameters=c(parameters$"value",K.values))

# write output file
output <- list(result=result,parameters=parameters)
filename <- ifelse(is.na(result[1]),
                   paste(outdir,"results_",jid,".Rdata",sep=""),
                   paste(outdir,"results_",jid,"success.Rdata",sep="")
)
save(output, file=filename)
