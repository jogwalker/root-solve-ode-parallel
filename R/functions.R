# functions for parallel run_solve_simulations

# Load packages
library(deSolve) 
library(reshape)
library(rootSolve)

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

### create the Model

createModel <- function(LVcreate.result) {
  nodes <- LVcreate.result[[1]]
  key <- LVcreate.result[[2]]
  Kname <- key$Kname[which(!is.na(key$Kname))]
  equations <- LVcreate.result[[3]]
  
  filename <- "/clusterdata/uqgheman/jummy/qualitative-modeling/R/model.R"
  
  cat("runModel <- function(t,y,p) {\n",file=filename)
  
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
      cat(c(Kname[i]," <- p[",nrow(key)+i,"]\n"),sep="",file=filename,append=TRUE)
    }
  }  
  
  # define the equations
  cat("ydot <- rep(NA,length(y))\n",file=filename,append=TRUE)
  for (i in 1:length(equations)) {
    cat(c("ydot[",i,"] <- ", equations[i]," \n"),sep="",file=filename,append=TRUE)
  }
  
  cat("list(ydot)\n}"
      ,file=filename,append=TRUE)
  
  #source(filename,keep.source=TRUE)
  #paste("Function called \'runModel(t,y,p,K)\' created")
}

### function to simulate the parameters for all

simParamsAll <- function(n,network) {
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

# function to simulate just one set of parameters

simParams <- function(network) {
  .network <- network
  if(class(.network)=="matrix") {
    .network <- melt(.network) 
    names(.network)  <- c("to","from","class")
  }
  
  #create a data frame to hold the parameter estimates
  parameters <- data.frame(vector(length=nrow(.network)))
  names(parameters) <- "value"
  
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
    # standard beta distribution (positive or negative)
    parameters[which(link.index$ab),] <- rbeta(sum(link.index$ab),0.5,1)
    parameters[which(link.index$co),] <- -rbeta(sum(link.index$co),0.5,1)
    parameters[which(link.index$py),] <- -rbeta(sum(link.index$py),0.5,1)
    parameters[which(link.index$mu),] <- rbeta(sum(link.index$mu),0.5,1)
    parameters[which(link.index$pr),] <- rbeta(sum(link.index$pr),0.5,1)
    parameters[which(link.index$"in"),] <- -rbeta(sum(link.index$"in"),0.5,1)
    #other distributions
    parameters[which(link.index$uk),] <- unknown.links(sum(link.index$uk))
    parameters[which(link.index$dd),] <- -runif(sum(link.index$dd),0,1)
    parameters[which(link.index$na),] <- 0

  parameters<- (cbind(.network,parameters)) 
  return(parameters)
}


###Function to determine how unknown links are handled
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

#### function to run a simulation

runSim <- function(inits, timesteps, parameters) {
  run.try <- try(runsteady(y=inits, times = c(0,timesteps),func=runModel,parms=parameters,mf=22),silent=TRUE)
  temp <- NA
  if (class(run.try)[1]!="try-error") {temp <- run.try}
  return(temp)
}