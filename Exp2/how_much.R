# Sets current work directory to the one containing the csv file with the data
dataDir <- "./Task/data/"
dataFile <- "532_cost3_11022019_main_log.csv"

# Loads the csv file into R as a data frame
data <- read.csv(paste(dataDir, dataFile, sep = ""))

# Choices holds all the yes and no decisions everyone made, and rewards holds every offer corresponding to each trial. 
b <- data.frame(choices <- data[4], 
          rewards <- data[3])

# This produces the row numbers (trial numbers) in which people accepted an offer
yes=which(b[[1]] == 1, arr.ind=TRUE)

# Produces a dataframe which shows only the trials in which the person accepted the offer 
yesTrial=b[yes,]

# Sums up the total amount of cents this participant earned 
cents=sum(yesTrial[[2]])

# Converts to dollars for better presentation
dollars=12+(cents/100)
  
# Prints a value the participant earned to the screen, updating depending on the value
sprintf("This participant earned %s dollars", dollars)  







