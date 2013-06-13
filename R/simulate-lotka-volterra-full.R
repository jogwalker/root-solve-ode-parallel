# 3 June 2013
# Josephine Walker
# new attempt after simulate-lotka-volterra.R
# generalize !

#-----------------------------------------------------------------------------

# Load required packages
library(deSolve) # might not actually need this now?
library(reshape)
library(rootSolve)

#-------------

# Step 1: define the network (this is done in simulate_network.R)
# network can either be square matrix or long (melted) data frame (with columns "to", "from","class")

#-------------

load("network.rdata")
if (class(network)=="data.frame") {node.names <- unique(network[,1])}
else if (class(network)=="matrix") {node.names <- colnames(network)}

#-------------

# Step 2: create a function to define the form of a lotka-volterra equation system for use in the ode solver 

#-------------

# This function takes a character vector of node names, and a logical vector of equal length which determines whether that node has logistic growth (K limited) or not
# the function returns a key to the parameters and a list of equation strings
LVcreate <- function(node.names,logistic) {
  nodes <- node.names
  n <- length(nodes)
  logistic <- logistic
  equations <- rep(NA,n)
  index <- expand.grid(1:n,1:n)
  index$list <- as.numeric(rownames(index))
  for (i in 1:n) {  
    temp <- index[which(index$Var2==i),]
    x <- which(temp$Var1==temp$Var2)
    parts <- rep(NA,n)
    for (j in 1:n){
      # diagonal interactions
      if (j == x) {
        if (logistic[j]==TRUE) {parts[j] <- paste("p",temp$list[j],"*",nodes[temp$Var1[j]],"*(1-",nodes[temp$Var1[j]],"/K_",nodes[temp$Var1[j]],")",sep="")
        }
        else {parts[j] <- paste("p",temp$list[j],"*",nodes[temp$Var2[j]],sep="")}
      } 
      # off-diagonal interactions
      else {parts[j] <- paste("p",temp$list[j],"*",nodes[temp$Var1[j]],"*",nodes[temp$Var2[j]],sep="")}
    }
    # put it together
    equations[i] <- paste(parts[],collapse="+")
  }
  
  # create a key of what the parameters refer to
  key <- expand.grid(nodes,nodes)
  names(key) <- c("from","to")
  key$parameter <- paste("p",index$list,sep="")
  key$Kname <- NA
  ind <- key$from == key$to & logistic[match(key$to,nodes)]
  key$Kname[ind] <-  paste("K_", key$to[ind], sep="")
  
  # return key and equations
  result <- list(nodes,key,equations)
  return (result)
}

#-------------

# Step 3: create the model for running with differential equation solver using output from LVcreate 

#-------------

createModel <- function(LVcreate.result) {
  nodes <- LVcreate.result[[1]]
  key <- LVcreate.result[[2]]
  Kname <- key$Kname[which(!is.na(key$Kname))]
  equations <- LVcreate.result[[3]]
  
  filename <- tempfile()
  
  cat("runModel <- function(t,y,p,K=NULL) {\n",file=filename)
  
  # define initial values for each node
  for (i in 1:length(nodes)) {
    cat(c(nodes[i]," <- y[",i,"]\n"),sep="",file=filename,append=TRUE)
  }
  
  # define all parameters
  for (i in 1:nrow(key)) {
    cat(c(key[i,"parameter"]," <- p[",i,"]\n"),sep="",file=filename,append=TRUE)
  }
  
  if(length(Kname) > 0) {
    for (i in 1:length(Kname)) {
      cat(c(Kname[i]," <- K[",i,"]\n"),sep="",file=filename,append=TRUE)
    }
  }  
  
  # define the equations
  cat("ydot <- rep(NA,length(y))\n",file=filename,append=TRUE)
  for (i in 1:length(equations)) {
    cat(c("ydot[",i,"] <- ", equations[i]," \n"),sep="",file=filename,append=TRUE)
  }
  
  cat("list(ydot)\n}"
      ,file=filename,append=TRUE)
  
  source(filename,keep.source=TRUE)
  paste("Function called \'runModel(t,y,p,K)\' created")
}

#-------------

# Step 4: create function to simulate all the parameters for each run based on defined network

#-------------

simParams <- function(n,network) {
  n <- n # number of simulation runs
  .network <- network
  if(class(.network)=="matrix") {
    .network <- melt(.network) 
    names(.network)  <- c("to","from","class")
  }
  
  #create a data frame to hold the parameter estimates
  parameters <- data.frame(matrix(nrow=nrow(.network),ncol=n,dimnames=list(NULL,seq(1,n,1))))
  
  # create a list of indices for each link class based on network
  link.class <- as.character(unique(.network$class))
  n.class <- length(link.class)
  link.index <- data.frame(matrix(ncol=n.class,nrow=nrow(parameters)))
  names(link.index) <- link.class
  for (i in 1:n.class) {
    link.index[,link.class[i]] <- as.logical(match(.network$class,link.class[i]))
  }
  link.index[is.na(link.index)] <- FALSE # clean up NAs
  
  # simulate values based on each class for each iteration
  # here they are defined explicitly for each class, so this has to change if your classes change
  # "ab" "uk" "dd" "co" "py" "na" "mu" "pr" "in"
  for (i in 1:n) {
    # standard beta distribution (positive or negative)
    parameters[which(link.index$ab),i] <- rbeta(sum(link.index$ab),0.5,1)
    parameters[which(link.index$co),i] <- -rbeta(sum(link.index$co),0.5,1)
    parameters[which(link.index$py),i] <- -rbeta(sum(link.index$py),0.5,1)
    parameters[which(link.index$mu),i] <- rbeta(sum(link.index$mu),0.5,1)
    parameters[which(link.index$pr),i] <- rbeta(sum(link.index$pr),0.5,1)
    parameters[which(link.index$"in"),i] <- -rbeta(sum(link.index$"in"),0.5,1)
    #other distributions
    parameters[which(link.index$uk),i] <- unknown.links(sum(link.index$uk))
    parameters[which(link.index$dd),i] <- -runif(sum(link.index$dd),0,1)
    parameters[which(link.index$na),i] <- 0
    
    ## confirm - what do I want to do with dd?.... round to fewer digits?
  }  
  parameters<- (cbind(.network,parameters)) 
  return(parameters)
}

#-------------

# Step 5: Function to determine how unknown links are handled
# input is length of vector of parameters to be simulated

#-------------

unknown.links <- function(n) {
  n <- n
  param.vector <- rep(NA,n)
  p <- runif(1,0,1) # how many will be zero
  sign <- 0.5 # what proportion positive 
  p.vector <- rbinom(n,1,p)
  sign.vector <- rbinom(n,1,sign)
  sign.vector[which(sign.vector==0)] <- -1
  param.vector <- rbeta(n,0.5,1)*p.vector*sign.vector
  return(param.vector)
}

#-------------

# Step 6: decide which nodes will have logistic growth and define K values
# should check whether equilibrium values are sensitive to K
# could also do competition model and see what values lead to some % tree cover

#-------------

K.limited <- c("grass","trees")
K.values <- c(10,10) # in 100 kg per m^2

K.logical <- rep(FALSE,length(node.names))
K.logical[match(K.limited,node.names)] <- TRUE


#-------------

# Step 7: define range of initial values
# what are the units? individuals per km^2?
# sources will be makgadikgadi report, and central statistics office...
#-------------

# simple
inits.simple <- data.frame(node.names)
inits.simple$value <- 1

# or more complicated 
inits <- data.frame(matrix(nrow=length(node.names),ncol=3))
names(inits) <- c("node","low","high")
inits$node <- node.names
inits[which(inits$node=="water"),2:3] <- c(1,1)
inits[which(inits$node=="grass"),2:3] <- c(1,1)
inits[which(inits$node=="trees"),2:3] <- c(1,1)
inits[which(inits$node=="crops"),2:3] <- c(1,1)
inits[which(inits$node=="dom_rum"),2:3] <- c(1,1)
inits[which(inits$node=="dom_eq"),2:3] <- c(1,1)
inits[which(inits$node=="wild_rum"),2:3] <- c(1,1)
inits[which(inits$node=="wild_eq"),2:3] <- c(1,1)
inits[which(inits$node=="par_eq_dom"),2:3] <- c(1,1)
inits[which(inits$node=="par_eq_wild"),2:3] <- c(1,1)
inits[which(inits$node=="par_rum_dom"),2:3] <- c(1,1)
inits[which(inits$node=="par_rum_wild"),2:3] <- c(1,1)
inits[which(inits$node=="par_eq_inf"),2:3] <- c(1,1)
inits[which(inits$node=="par_rum_inf"),2:3] <- c(1,1)
inits[which(inits$node=="carnivores"),2:3] <- c(1,1)
inits[which(inits$node=="elephants"),2:3] <- c(1,1)
inits[which(inits$node=="people"),2:3] <- c(1,1)

#-------------

# Step 8: define range of reasonable equilibrium values (for some)
# currently made up values
#-------------

eq.range <- inits
eq.range[,c("low")] <- 0
eq.range[,c("high")] <- NA
eq.range[which(eq.range$node=="wild_rum"),2:3] <- c(.5,1.5)
eq.range[which(eq.range$node=="wild_eq"),2:3] <- c(.5,1.5)
eq.range[which(eq.range$node=="carnivores"),2:3] <- c(.05,.15)
eq.range[which(eq.range$node=="elephants"),2:3] <- c(.1,.4)
# also, all runsteady output $y should be > 0

#-------------

# Step 9: function to run solver (using function runsteady())

#-------------
# takes vector of inits, length of timesteps, and vector of parameters 
# uses function try() to see if runsteady() finds an equilibrium,
# returns list of runsteady() output if it works, NA otherwise

runSim <- function(inits, timesteps, parameters) {
  run.try <- try(runsteady(y=inits, times = c(0,timesteps),func=runModel,parms=parameters,mf=22),silent=TRUE)
  temp <- NA
  if (class(run.try)[1]!="try-error") {temp <- run.try}
  return(temp)
}

#-------------

# Step 10: run all the loops

#-------------
runSimLoop <- function(n,inits,timesteps,network,keep.params=TRUE){
  # simulate all the parameters up front
  parameters.labeled <- simParams(n,network)
  parameters <- data.frame(parameters.labeled[,-c(1:3)])
  
  # initialize results list and counter
  results <- list()
  successes <- list()
  
  # run the loop
  for (i in 1:n) {
    run.try <- try(
      runsteady(y=inits, times = c(0,timesteps),func=runModel,parms=parameters[,i],mf=22)
      ,silent=TRUE)
    if (class(run.try)[1]!="try-error") {
      results[[i]] <- run.try
      successes <- c(successes,i)
    }
    else {results[[i]] <- NA}
    print(paste(i," out of ",n, " runs completed \n", length(successes)," successes"))
  }
  if(keep.params==FALSE) {
    return(c(successes,results))
  }
  else {return (c(parameters.labeled,successes,results))}
}


#-------------

# Step 11: filter results

#-------------


#-------------


# define reasonable initial values (or range to test) 
# figure out if each simulated model has a stable equilibrium (how?)
# is stable equilibrium sensitive to initial values?
# figure out if stable equilibrium is reasonable
# calculate Jacobian from simulated values
# press perturbation!
# what is the new equilibrium? or if I know the old equilibrium values can I quantitatively press perturb to get new equilibrium?
# ----> also!! write up what the model is. this will help make those arbitrary decisions because they will be recorded.

######################
# another way to make the model... as a matrix






#######################


# test

node.names
logs <- rep(FALSE,length(node.names))
logs

createModel(LVcreate(node.names,logs))

inits <- rep(1,17)
t <- 10000
pars <- simParams(1,network)
print(system.time(
  x <- runSim(inits, t, pars[,4])
))

#-----------------------

# test 2

x <- letters[1:10]
log <- rep(FALSE,length(letters))
network.10 <- expand.grid(x,x)
network.10[,3] <- "uk"
names(network.10) <- c("to","from","class")

createModel(LVcreate(x,log))
inits <- rep(1,length(x))
t <- 10000
# pars <- simParams(1,network.10) # need to actually define the classes of the nodes
pars <- simParams(1,network)[1:100,4]
x <- runSim(inits,t,pars)

#----------------------

# test 3

n <- 2 # iterations
t <- 10000
y0 <- rep(1,length(node.names))
logs <- rep(FALSE,length(node.names))
createModel(LVcreate(node.names,logs))
x <- runSimLoop(n=n,inits=y0,timesteps=t,network=network,keep.params=FALSE)
x[[2]]

# test 4

names.4 <- letters[1:4]
network.4 <- expand.grid(names.4,names.4)
network.4[,3] <- "uk"
names(network.4) <- c("to","from","class")
network.4[which(network.4$to==network.4$from),"class"] <- "dd"
network.4[2,3] <- "pr"
network.4[5,3] <-  "py"
#.....