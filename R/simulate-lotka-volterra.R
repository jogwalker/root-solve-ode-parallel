# 30 May 2013
# Josephine Walker
# After meeting with Mike Bode
# 'qualitative modelling' project
# In this script I plan to simulate the equilibrium values of the Makgadikgadi Pans ecosystem
# This will involve setting up the Lotka-Volterra equations with simulated values and running them to equilibrium


# First get it working with just a few interacting species
# ordinary differential equations or difference equations?
#------------------------------------------------------------------------------

library(deSolve)

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

                       

#----------------------------------------------------------------------------
  
# run model
odeup = dede(inits,t,model,parameters)

# get results
output <- data.frame(odeup)
names(output) <- c("time","A","B","C")
output$total = apply(output[,-1],1,sum)

# plot results

plot(1,1,xlim=c(0,length(t)-1),ylim=c(0,5000), type="n",xlab="time(days)",ylab="population size")
lines(output$time,output$A,col="black") # points or lines
lines(output$time,output$B,col="orange")
lines(output$time,output$C,col="red")

