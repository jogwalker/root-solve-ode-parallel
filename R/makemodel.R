# make the model once for running in parallel

# get node names
if (class(network)=="data.frame") {node.names <- unique(network[,1])}
if (class(network)=="matrix") {node.names <- colnames(network)}


# define K limited values
K.limited <- c("grass","trees")
K.values <- c(10,10) # in 10 kg per m^2
K.logical <- rep(FALSE,length(node.names))
K.logical[match(K.limited,node.names)] <- TRUE

createModel(LVcreate(node.names,K.logical))

# define inits
inits <- rep(1,length(node.names))

# define max timesteps
t <- 100000