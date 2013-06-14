# load the list of filenames (without)

all.filenames <- read.table("C:/Users/Josephine/Dropbox/Dropbox Documents/Bristol_PhD/R/qualitative modeling/results/successnum.txt")

results <- list()

# read in success one at a time
# unpack output, name with run number

for (i in 1:nrow(all.filenames)) {
  file <- paste("C:/Users/Josephine/Dropbox/Dropbox Documents/Bristol_PhD/R/qualitative modeling/results/results_",all.filenames[i,"V1"],"success.RData",sep="")
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

# for any steady ones run ode 
# if the number is >= 20,000 set.seed(number), otherwise, set.seed(number*123456)

results.winners <- subset(results,summary$all.positive==TRUE)
winners.ind <- which(summary$all.positive==TRUE)
winners <- summary$id[winners.ind]

for (i in 1:length(winners)) {
  
}








