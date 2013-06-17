# load the list of filenames (without)

all.filenames <- read.table("C:/Users/Josephine/Dropbox/Dropbox Documents/Bristol_PhD/R/qualitative modeling/results/successes.txt")

results <- list()

# read in success one at a time
# unpack output, name with run number

for (i in 1:nrow(all.filenames)) {
  file <- paste("C:/Users/Josephine/Dropbox/Dropbox Documents/Bristol_PhD/R/qualitative modeling/results/successes/results_",all.filenames[i,"V1"],"success.RData",sep="")
    load(file)
  
  number[i] <- all.filenames[i,"V1"]
  results[[i]] <- list(number[i], output$result, output$parameters)
}

# check if result is "steady"
# check if any result$y <= 0
summary <- data.frame(matrix(nrow=length(results),ncol=3))
names(summary) <- c("id","steady","all.positive")
summary$id <- number

for (i in 1:nrow(summary)) {
  summary$steady[i] <- attr(results[[i]][2][[1]],"steady")
  ifelse(all(results[[i]][2][[1]]$y > 0),summary$all.positive[i] <- TRUE,summary$all.positive[i] <- FALSE)
}

table(summary$steady)
table(summary$all.positive)

# for any steady ones run ode 
# if the number is >= 20,000 set.seed(number), otherwise, set.seed(number*123456)

results.winners <- subset(results,summary$all.positive==TRUE)
winners.ind <- which(summary$all.positive==TRUE)
winners <- summary$id[winners.ind]

save(results.winners, file="NoNegs.Rdata")

####################################

inits <- rep(1,length(results.winners[[1]][2][[1]]$y))
t <- seq(1,1000,1)

for (i in 1:length(winners.ind)) {
  parameters <- c(results.winners[[i]][[3]]$value,K.values)
  results.winners[[i]][[4]] <- ode(inits,t,runModel,parameters)
}

#a <- runsteady(y=inits, times = c(0,1000),func=runModel,parms=parameters,mf=22)






