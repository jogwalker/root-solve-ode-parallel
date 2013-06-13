# 6 May 2013 After conversation with Mike Bode and Chris, try simulating values for the system
# 7 May 2013 - managed to get stable iterations, plot result by proportion, and calculate mean and confidence intervals for simulated response variables. next steps: implement management strategy and see how it changes...compared to what? or see what economic effect is...? compare degree of suppression? how do I calculate the suppression amount from the inverse matrix? (partial derivative?...)
# 8 May 2013
#-------------------------------------------------------------------------------
# Define the system: + for positive interaction, - for negative interaction, 0 for no interaction, NA for unknown interaction

# Nodes: water, grass, trees, crops, dom_rum, dom_eq, wild_rum, wild_eq, par_eq_dom, par_eq_wild, par_rum_dom, par_rum_wild, par_eq_inf, par_rum_inf, carnivores, elephants, people

# potential management actions: add tourists, add fence, add anthelmintics

# Define links:

nodes <- c("water", "grass", "trees", "crops", "dom_rum", "dom_eq", "wild_rum", "wild_eq", "par_eq_dom", "par_eq_wild", "par_rum_dom", "par_rum_wild", "par_eq_inf", "par_rum_inf", "carnivores", "elephants", "people")

n <- length(nodes)

network <- matrix(nrow=n,ncol=n,dimnames=list(nodes,nodes))

# define by type of interaction to reduce uncertainty (effect of column on row)
# ab = abiotic/incidental increase, uk = unknown, mu = mutualistic, pr = predator, py = prey, co = competition, dd = density-dependent (diagonal), na = none, in = incidental reduction
network[,"water"] <- 'ab'
# network["water","water"] <- 'na' #'dd'
network[,"grass"] <- c('uk','dd','co','uk','py','py','py','py','uk','uk','uk','uk','ab','ab','na','uk','uk')
network[,"trees"] <- c('uk','co','dd','uk','uk','uk','uk','uk','na','na','na','na','na','na','uk','py','uk')
network[,"crops"] <- c('uk','uk','uk','dd','uk','uk','na','na','na','na','na','na','na','na','na','py','mu')
network[,"dom_rum"] <- c('uk','pr','na','uk','dd','co','co','co','na','na','py','na','in','in','py','na','mu')
network[,"dom_eq"] <- c('uk','pr','na','uk','co','dd','co','co','py','na','na','na','in','in','py','na','mu')
network[,"wild_rum"] <- c('uk','pr','na','na','co','co','dd','co','na','na','na','py','in','in','py','na','na') 
network[,"wild_eq"] <- c('uk','pr','na','na','co','co','co','dd','na','py','na','na','in','in','py','na','na')
network[,"par_eq_dom"] <- c('uk','na','na','na','na','pr','na','na','dd','na','na','na','mu','na','na','na','na')
network[,"par_eq_wild"] <- c('uk','na','na','na','na','na','na','pr','na','dd','na','na','mu','na','na','na','na')
network[,"par_rum_dom"] <- c('uk','na','na','na','pr','na','na','na','na','na','dd','na','na','mu','na','na','na')
network[,"par_rum_wild"] <- c('uk','na','na','na','na','na','pr','na','na','na','na','dd','na','mu','na','na','na')
network[,"par_eq_inf"] <- c('uk','na','na','na','na','na','na','na','mu','mu','na','na','dd','na','na','na','na')
network[,"par_rum_inf"] <- c('uk','na','na','na','na','na','na','na','na','na','mu','mu','na','dd','na','na','na')
network[,"carnivores"] <- c('uk','na','na','na','pr','pr','pr','pr','na','na','na','na','na','na','dd','uk','na')
network[,"elephants"] <- c('uk','uk','pr','pr','na','na','uk','uk','na','na','na','na','na','na','uk','dd','na')
network[,"people"] <- c('uk','uk','uk','mu','mu','mu','na','na','na','na','na','na','na','na','in','na','dd')


save(network, file='network.rdata')

#-------------------------------------------------------------------------------
# Some thoughts, limitations, and assumptions of the network

# assume no hunting??? (no interaction between people and wild grazers)
# assume parasites have no effect on grass - could have effect through change in eating habits of host?
# do infective parasites have dd effect on themselves?
# what is water affected by?
# what known effects can I use to check the validity of the model?
# doesn't say anything about structure of groupings (eg wild ruminants - what species?)

# Thoughts and assumptions on the simulation method

# can define predator prey interaction - pr = .1py
# Raymond defines unknown interactions as existing or not existing before running simulations
# alternatively, could just draw from distribution between -1 and 1 or whatever.
# for some unknowns - they will be linked, so if ij is positive ji should be negative???
# could define different distributions to draw from
# runif does not choose extreme values (eg 0,1)
# for unknowns, might want to increase chance it is 0 if not simulating that first
# how many decimal points to take?
# -- to get stable iterations, reduce the proportion of unknown that are nonzero (to a point);
# make waterxwater na or dd, reduce range of values for each cell, or set dd = -1

# confirm that the NEGATIVE inverse indicates a press perturbation by reduction - Raymond takes negative inverse and then takes the negative of it......

# some potential validation criteria -
# domestic animals compete more with each other than with wild? or equids compete more with each other than with wild? or do they?

# how to include management strategies? and compare to base model...
# press perturbation on ruminant and equid domestic parasites. do any of the mean responses change significantly? (account for multiple testing...)
# change something in the system, compare new system to old. how to account for uncertainty?


#-------------------------------------------------------------------------------

# Simulate values for the matrix (see Raymond et al 2010)
# Keep all matrices for now, with separate vector of stability

itN <- 100000 #number of iterations

# initialize matrix to hold simulation results
mat.all <- array(dim=c(dim(network),itN))
stable <- vector(length=itN)
dets <- vector(length=itN)
kappas  <- vector(length=itN)

maxval <- 1

for (k in 1:itN) {
  # simulate values based on each network cell definition
  for (i in 1:nrow(network)) {
    for (j in 1:ncol(network)) {
       if (network[i,j] == "ab") {mat.all[i,j,k] <- round(runif(1,0,maxval),digits=2)}
       else if (network[i,j] == "uk") {if (rbinom(1,1,.33) == 1) {mat.all[i,j,k] <- round(runif(1,-maxval,maxval),digits=2)} else mat.all[i,j,k] <- 0} 
       # reducing the probability of giving a value to unknown increases proportion stable
       # else if (network[i,j] == "uk") {mat.all[i,j,k] <- round(runif(1,-1,1),digits=2)} 
       else if (network[i,j] == "mu") {mat.all[i,j,k] <- round(runif(1,0,maxval),digits=2)}
       else if (network[i,j] == "pr") {mat.all[i,j,k] <- round(runif(1,-maxval,0),digits=2)}
       else if (network[i,j] == "py") {mat.all[i,j,k] <- round(runif(1,0,maxval),digits=2)}
       else if (network[i,j] == "co") {mat.all[i,j,k] <- round(runif(1,-maxval,0),digits=2)}
       else if (network[i,j] == "dd") {mat.all[i,j,k] <- -1} #round(runif(1,-maxval,-.25),digits=2)}
       else if (network[i,j] == "na") {mat.all[i,j,k] <- 0} 
       else if (network[i,j] == "in") {mat.all[i,j,k] <- round(runif(1,-maxval,0),digits=2)}
    }
  }
  # check for Lyapunov stability (this is based on Raymond)
  if (all(Re(eigen(mat.all[,,k],only.values=T)$values)<0)) {stable[k] <- TRUE}
  
  # calculate the determinant and kappa
  dets[k] <- det(mat.all[,,k])
  kappas[k] <- kappa(mat.all[,,k])
}

# how many iterations are stable?
sum(stable)

# how many are non-invertible?
x <- length((which(dets==0)))
which(dets==0) %in% stable 

# take subset of stable matrices and work with them 
mat.stable <- mat.all[,,which(stable)]
dets.stable <- dets[which(stable)]
kappas.stable <- kappas[which(stable)]

# what do the stable matrices look like?
summary(dets.stable) # really small determinants! (all negative)
summary(kappas.stable)

# invert the matrices
mat.inv <- array(dim=dim(mat.stable), dimnames=list(nodes,nodes))

# take the negative inverse of the community matrix
for (k in 1:sum(stable)) {
  mat.inv[,,k] <- -solve(mat.stable[,,k])
} # set small values to zero?

#-------------------------------------------------------------------------------
# some preliminary analysis 

# should check if matrices fit some validation criteria...
# meanwhile...

# keep track of proportion of simulations with pos, neg, zero effect, and mean and 95% CI of effect
vals <- c("pos","neg","zero","tot","ppos","pneg","pzero","mean","sd","error","CI-low","CI-high","p")
resp.sign.all <- sign(mat.inv)
alpha = 0.05


# Response to reduction of each variable
resp.neg.all <- array(dim=c(length(nodes),length(nodes),length(vals)),dimnames=list(nodes,nodes,vals))
# Response to addition of each variable 
# (no need to really do these separately because they are just sign opposites of each other)
resp.pos.all <- array(dim=c(length(nodes),length(nodes),length(vals)),dimnames=list(nodes,nodes,vals))

for (i in 1:length(nodes)) {
  for (j in 1:length(nodes)) {
    resp.neg.all[i,j,"pos"] <- length(which(-resp.sign.all[i,j,]==1))
    resp.neg.all[i,j,"neg"] <- length(which(-resp.sign.all[i,j,]==-1))
    resp.neg.all[i,j,"zero"] <- length(which(-resp.sign.all[i,j,]==0))
    resp.neg.all[i,j,"tot"] <- sum(resp.neg.all[i,j,c("pos","neg","zero")])
    resp.neg.all[i,j,"ppos"] <- resp.neg.all[i,j,"pos"]/resp.neg.all[i,j,"tot"]
    resp.neg.all[i,j,"pneg"] <- resp.neg.all[i,j,"neg"]/resp.neg.all[i,j,"tot"]
    resp.neg.all[i,j,"pzero"] <- resp.neg.all[i,j,"zero"]/resp.neg.all[i,j,"tot"]
    resp.neg.all[i,j,"mean"] <- mean(-mat.inv[i,j,])
    resp.neg.all[i,j,"sd"] <- sd(-mat.inv[i,j,])
    resp.neg.all[i,j,"error"] <- qnorm(1-alpha/2)*resp.neg.all[i,j,"sd"]/sqrt(resp.neg.all[i,j,"tot"])
    resp.neg.all[i,j,"CI-low"] <- resp.neg.all[i,j,"mean"]-resp.neg.all[i,j,"error"]
    resp.neg.all[i,j,"CI-high"] <- resp.neg.all[i,j,"mean"]+resp.neg.all[i,j,"error"]
    resp.neg.all[i,j,"p"] <- t.test(-mat.inv[i,j,],alternative="two.sided",mu=0)$p.value
    
    resp.pos.all[i,j,"pos"] <- length(which(resp.sign.all[i,j,]==1))
    resp.pos.all[i,j,"neg"] <- length(which(resp.sign.all[i,j,]==-1))
    resp.pos.all[i,j,"zero"] <- length(which(resp.sign.all[i,j,]==0))
    resp.pos.all[i,j,"tot"] <- sum(resp.pos.all[i,j,c("pos","neg","zero")])
    resp.pos.all[i,j,"ppos"] <- resp.pos.all[i,j,"pos"]/resp.pos.all[i,j,"tot"]
    resp.pos.all[i,j,"pneg"] <- resp.pos.all[i,j,"neg"]/resp.pos.all[i,j,"tot"]
    resp.pos.all[i,j,"pzero"] <- resp.pos.all[i,j,"zero"]/resp.pos.all[i,j,"tot"]
    resp.pos.all[i,j,"mean"] <- mean(mat.inv[i,j,])
    resp.pos.all[i,j,"sd"] <- sd(mat.inv[i,j,])
    resp.pos.all[i,j,"error"] <- qnorm(1-alpha/2)*resp.pos.all[i,j,"sd"]/sqrt(resp.pos.all[i,j,"tot"])
    resp.pos.all[i,j,"CI-low"] <- resp.pos.all[i,j,"mean"]-resp.pos.all[i,j,"error"]
    resp.pos.all[i,j,"CI-high"] <- resp.pos.all[i,j,"mean"]+resp.pos.all[i,j,"error"]
    resp.pos.all[i,j,"p"] <- t.test(mat.inv[i,j,],alternative="greater",mu=0)$p.value
  }
}

#resp.neg.all[,"carnivores",]
#resp.pos.all[,"carnivores",]
#resp.pos.all[,,c("CI-low","mean","CI-high","p")]
length(which(resp.pos.all[,,"CI-low"]>=0))
length(which(resp.pos.all[,,"CI-high"]<=0))
length(which(resp.pos.all[,,"p"]<=0.2))

pt(abs(resp.pos.all[1,1,"mean"])/resp.pos.all[1,1,"sd"],400,low=F) > (0.05/17^2)

resp.pos.all[,,"ppos"]
resp.pos.all[,,"p"]
# plot responses

#barplot(as.matrix(resp.neg.all[,"",c("ppos","pneg")]),cex.names=0.7,xlab="Proportion",horiz=TRUE,las=1)

# plot response to negative perturbation of domestic ruminant parasites 
sub <- resp.neg.all[,"par_rum_dom",c("ppos","pneg","pzero")]
barplot(t(sub),cex.names=0.7,xlab="Proportion",horiz=TRUE,las=1)
title("Press perturbation on par_rum_dom",line=3)

# plot response to negative perturbation of domestic equid parasites 
sub2 <- resp.neg.all[,"par_eq_dom",c("ppos","pneg","pzero")]
barplot(t(sub2),cex.names=0.7,xlab="Proportion",horiz=TRUE,las=1)
title("Press perturbation on par_eq_dom",line=3)

# plot response to negative perturbation of domestic equid and ruminant parasites together
both <- -mat.inv[,c("par_eq_dom","par_rum_dom"),]
both.sum <- apply(both,c(1,3),sum) #sum effects of variables of interest on all nodes (1) for each iteration (3)

both.sign <- sign(both.sum)
both <- matrix(nrow=nrow(both.sign),ncol=3,dimnames=list(nodes, c("pos","neg","zero")))

for (i in 1:nrow(both.sign)) {
  both[i,"pos"] <- length(which(both.sign[i,]==1))
  both[i, "neg"] <- length(which(both.sign[i,]==-1))
  both[i,"zero"] <- length(which(both.sign[i,]==0))
}
both.p <- both[,1:3]/ncol(both.sign)
barplot(t(both.p),cex.names=0.7,xlab="Proportion",horiz=TRUE,las=1)
title("Press perturbation on par_eq_dom and par_rum_dom",line=3)

#-------------------------------
# test what would happen if reduced number of people and domestic ruminants (together)
resp <- -mat.inv[,c("people","dom_rum"),]
resp.sum <- apply(resp,c(1,3),sum) #sum effects of variables of interest on all nodes (1) for each iteration (3)

resp.sign <- sign(resp.sum)
resp <- matrix(nrow=nrow(resp.sign),ncol=3,dimnames=list(nodes, c("pos","neg","zero")))

for (i in 1:nrow(resp.sign)) {
  resp[i,"pos"] <- length(which(resp.sign[i,]==1))
  resp[i, "neg"] <- length(which(resp.sign[i,]==-1))
  resp[i,"zero"] <- length(which(resp.sign[i,]==0))
}

resp <- data.frame(resp)
resp$tot <- resp$pos+resp$neg+resp$zero
resp$ppos  <- resp$pos/resp$tot
resp$pneg <- resp$neg/resp$tot
resp$pzero <- resp$zero/resp$tot

resp.plot <- t(as.matrix(cbind(resp$ppos,resp$pneg,resp$pzero)))
dimnames(resp.plot) <- list(c("pos","neg","zero"),nodes)

# plot the responses
barplot(resp.plot,cex.names=0.7,xlab="Proportion",horiz=TRUE,las=1)
title("People/domestic ruminants suppression",line=3)
#--------------------------------------




