# load the list of filenames (without)

all.filenames <- read.table("successnum.txt")

results <- list()

# read in success one at a time

for (i in 1:length(all.filenames)) {
  load(paste("C:/Users/Josephine/Dropbox/Dropbox Documents/Bristol_PhD/R/qualitative modeling/results/results_",all.filenames[i],"success.RData",sep=""))
  
  number <- all.filenames[i]
  results[[i]] <- list(number, output$result, output$parameters)
}


# unpack output, name with run number

# check if result is "steady"
attr(result,"steady")

# check if any result$y <= 0

# keep list of steady, no zero IDs

# for any steady ones run ode 
# if the number is >= 20,000 set.seed(number), otherwise, set.seed(number*123456)