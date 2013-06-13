
args <- commandArgs(T)

outdir <- args[1]
jid <- as.numeric(args[2])

# read in all the functions
source("~/functions.R")

# Load in network data 
load("~/network.RData")
if (class(network)=="data.frame") {node.names <- unique(network[,1])}
else if (class(network)=="matrix") {node.names <- colnames(network)}

# define K limited values
K.limited <- c("grass","trees")
K.values <- c(100,100) # in kg per m^2
K.logical <- rep(FALSE,length(node.names))
K.logical[match(K.limited,node.names)] <- TRUE

# define inits
inits <- rep(1,length(node.names))

# define max timesteps
t <- 100000


# Create parameters
source("~/repo/ibd_ibs/bias_check_simulations/R/parameters.R")
params$path <- hsqPathnames(params, outdir)

# Run analysis
runSimulation(geno, fam, ibs, gcta, params[jid, ])
