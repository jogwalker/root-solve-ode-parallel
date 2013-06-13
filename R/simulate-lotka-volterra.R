# 30 May 2013
# Josephine Walker
# After meeting with Mike Bode
# 'qualitative modelling' project
# In this script I plan to simulate the equilibrium values of the Makgadikgadi Pans ecosystem
# This will involve setting up the Lotka-Volterra equations with simulated values and running them to equilibrium


# First get it working with just a few interacting species
# ordinary differential equations or difference equations?

#11 June trying to add in 'try' or 'tryCatch' to continue loop despite error
# what I want to do is suppress Warnings, ignore errors (next loop)
# also try rootSolve package - runsteady
#------------------------------------------------------------------------------

library(deSolve)
library(rootSolve)

# define the function to use in deSolve

model <- function(t,y,p) { # y are initial conditions, p are parameters, t is time
  
  # initial values
  A <- y[1]
  B <- y[2]
  C <- y[3]
  
  # parameters
  p1 <- p[1]
  p2 <- p[2]
  p3 <- p[3]
  p4 <- p[4]
  p5 <- p[5]
  p6 <- p[6]
  p7 <- p[7]
  p8 <- p[8]
  p9 <- p[9]
  
  #initialize vector to hold results
  ydot = rep(NA,length(y))
  
  # the equations
  ydot[1] <- p1*A+p2*A*B+p3*A*C
  ydot[2] <- p4*B+p5*B*A+p6*B*C
  ydot[3] <- p7*C+p8*C*A+p9*C*B
  
  # return for next calculation
  list(ydot)
}

#-----------------------------------------------------------------------------

# define parameters

# model run time
t = seq(0,1000,1)

# Initial conditions
Ainit = 1000
Binit = 1000
Cinit = 1000

inits = c(Ainit,Binit,Cinit)

# simulate values for parameters based on classifications

simParams9 <- function() {
  
  parameters <- vector("numeric",9)
  parameters[1] <- round(runif(1,0,1),digits=2)
  parameters[2] <- round(runif(1,-1,0),digits=2)
  parameters[3] <- 0
  parameters[4] <- round(runif(1,-1,0),digits=2)
  parameters[5] <- round(runif(1,0,1),digits=2)
  parameters[6] <- round(runif(1,-1,0),digits=2)
  parameters[7] <- round(runif(1,-1,0),digits=2)
  parameters[8] <- 0
  parameters[9] <- round(runif(1,0,1),digits=2)    
  
  return(parameters)
  
}
#----------------------------------------------------------------------------

# run model

runOnce <- function() {
  
  #odeup = try(ode(inits,t,model,simParams9()))
  print(system.time(
    run <- try(runsteady(y = inits,times= c(0,1000000),func=model,parms=simParams9(),mf=22),silent=TRUE) 
  ))
  x <- NA
  if (class(run)[1]!="try-error") {x <- run}
  # get results
  #output <- data.frame(run)
  #names(output) <- c("time","A","B","C")
  #output$total = apply(output[,-1],1,sum)
  return(x)
}


# plot results

plotRun <- function(x) {
  plot(1,1,xlim=c(0,length(t)-1),ylim=c(0,5000), type="n",xlab="time(days)",ylab="population size")
  lines(x$time,x$A,col="black") # points or lines
  lines(x$time,x$B,col="orange")
  lines(x$time,x$C,col="red")
}

# run as functions
#sink(file='TEMP.TXT')
x <- runOnce()
#sink()

plotRun(x)

#------------------------------------------------------------------------------

# simulate parameters many times and check if any result in non-zero equilibrium
# for some reason the loop doesn't work

runLooped <- function (N,timesteps) { 
  
  Ntimes <- N
  
  # model run time
  #time = seq(0,timesteps,1)
  timesteps <- timesteps
  
  # Initial conditions
  Ainit = 1000
  Binit = 1000
  Cinit = 1000
  
  inits = c(Ainit,Binit,Cinit)
  
  
  # run odeup for each set of parameters
  #all.output <- array(dim=c(length(time),length(inits),Ntimes))
  #all.output <- array(dim=c(length(inits),Ntimes))
  output <- list()
  
  for(i in 1:Ntimes) {
    inits = c(Ainit,Binit,Cinit)
    t = time 
    param = simParams9()
    #temp = try(ode(inits,t,model,param))
    temp = try(runsteady(y = inits,times= c(0,timesteps),func=model,parms=param,mf=22),silent=TRUE)
    if(class(temp)[1] != "try-error") {
      output[[i]] <- temp
      #all.output[,1,i] <- odeup[,2]
      #all.output[,2,i] <- odeup[,3]
      #all.output[,3,i] <- odeup[,4] 
    }
    else {output[[i]] <- NA}
  }
  
  # count non-zero equilibrium values
  #index.nonzero <- vector("logical",Ntimes)
  #if (all(all.output[dim(all.output)[1],2:4,i])>1) {index.nonzero[i] <- TRUE}
  
  return(output)
}

Ntimes <- 5
t <- 10000
xx <- runLooped(N = Ntimes,timesteps = t)
dim(xx)

################################################################################
#------------------------------------------------------------------------------#
################################################################################
# try running as difference equations instead of ode

# 1. generate random parameters and reasonable initial values
# 2. iterate each time step of system of equations
# 3. at each time step, check if any equations are <= 0, if so stop
# 4. repeat until some stable equilibrium is reached
# 5. check if equilibrium values match reasonable range
# 6. keep set of reasonable equilibrium and parameter values, discard others
# 7. calculate Jacobian matrix for each and confirm stability
# 8. simulate management: press perturbations and equilibrium changes

# try with three species 
N <- 10 # number of iterations
n.params = 9
p <- matrix(nrow=N,ncol=n.params)
timesteps <- 5
results <- list()

for (j in 1:N) {
  
  # simulate parameters
  
  p[j,] <- round(runif(n.params,0,1),digits=2)
  p[j,c(2:4,6:7)] <- -p[j,c(2:4,6:7)]
  
  # p[1] growthrate A
  # p[2] coefficient A*B 
  # p[3] coefficient A*C
  # p[4] growthrate B
  # p[5] coefficient B*A 
  # p[6] coefficient B*C
  # p[7] growthrate C
  # p[8] coefficient C*A 
  # p[9] coefficient C*B
  
  
  A1 <- 1000
  B1 <- 1000
  C1 <- 1000
  
  A2 <- 0
  B2 <- 0
  C2 <- 0
  
  t <- 0
  
  # for (i in 1:timesteps) {
  #   A2 <- A1[i] + p[1]*A1[i] + p[2]*A1[i]*B1[i] + p[3]*A1[i]*C1[i]
  #   B2 <- B2[i] + p[4]*B1[i] + p[5]*B1[i]*A1[i] + p[6]*B1[i]*C1[i]
  #   C2 <- C2[i] + p[7]*C1[i] + p[8]*C1[i]*A1[i] + p[9]*C1[i]*B1[i]
  #   cat(c(A2,B2,C2),"\n", i) 
  #   
  #   t[i] <- i
  #   if (i < timesteps) {A1[i+1]=A2; B1[i+1]=B2; C1[i+1]=C2}
  #   temp <- c(A1[i],B1[i],C1[i])
  #   if (any(temp < 0, na.rm = T)) break
  # }
  # results <- c(results,cbind(A1,B1,C1))
  # }
  
  for (i in 1:timesteps) {
    A2 <- A1[i] + p[1]*A1[i] + p[2]*A1[i]*B1[i] 
    B2 <- B2[i] + p[4]*B1[i] + p[5]*B1[i]*A1[i] 
    
    cat(c(A2,B2),"\n", i) 
    
    t[i] <- i
    if (i < timesteps) {A1[i+1]=A2; B1[i+1]=B2}
    temp <- c(A1[i],B1[i])
    if (any(temp < 0, na.rm = T)) break
  }
  results <- c(results,cbind(A1,B1))
}

# l <- list()
# l[[1]] <- 1
# l[[2]] <- 2
# l[[3]] <- "a"
# l[[4]] <- 4
# 
# o <- list()
# 
# for(i in 1:4)
# {
#   print(i)
#   o[[i]] <- try(log(l[[i]]))
#   print(o[[i]])
#   #if(class(o[[i]]) == "try-error")
#   #{
#   #  print("this is a problem")
#   #}
# }
