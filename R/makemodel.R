# preparations to run model to run each time

# Load in network data
load("/clusterdata/uqgheman/jummy/qualitative-modeling/R/network.rdata")

# get node names
if (class(network)=="data.frame") {node.names <- unique(network[,1])}
if (class(network)=="matrix") {node.names <- colnames(network)}


# define K limited values
K.limited <- c("grass","trees")
K.values <- c(10,10) # in 10 kg per m^2
K.logical <- rep(FALSE,length(node.names))
K.logical[match(K.limited,node.names)] <- TRUE

# this will have already been created
# createModel(LVcreate(node.names,K.logical))
source("/clusterdata/uqgheman/jummy/qualitative-modeling/R/model.R")

# define inits
inits <- rep(1,length(node.names))

# define max timesteps
t <- 100000

