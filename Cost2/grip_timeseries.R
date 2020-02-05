# explore the relationship between reward/handling and grip strength for a single subject

# script mostly to start thinking of Biopac-friendly functions
# a problem with this analysis is that the data seem cutoff.
# maybe because biopac stopped running or the computer went to sleep
# I'm limiting this test script to just half of the data

# general log info
# psychopy logs each trial either when the reward is achieved or when they quit
# so gripping starts Xs earlier than the Experiment.Time (x = handling time)
# create a vector that's Experiment.Time - RT to get approx when the trial starts

# to do: get every 100th milisecond of smoothed and divide trialStart (i.e. downsample); average pressing per trial

###---------- library and functions ---------
library(tidyverse)

select_signal <- function(physData, start, end, thresh = 2, smoothing = 0) {
  # this function will use the experimental timing from a log to get rid of non-experimental recordings
  # this will be recordings prior to the beginning of the session and when the trials stopped
  # it currently assumes milisecond rate biopac acquisition
  # variables:
  # - physData: the vector with grip strength
  # - start/end: just the time bounds (in ms) containing relevant data 
  # - thresh: sustained threshold for the grip press (can be estimated by plotting the time series)
  
  # the result is a time series vector of grip strength that is aligned with the behavioral log
  
  # where does the first trial start?
  firstTrial <- which(physData > thresh)[1]
  
  # index where to cut pre-experimental data
  preExp <- firstTrial - start 
  postExp <- end + firstTrial
  
  # isolate the signal of interest
  physData <- physData[preExp:postExp]
  
  # smooth the data using a moving average
  # probably unnecessary for this data
  if (smoothing > 0) {
    physData <- stats::filter(physData, sides = 1, filter = rep(1/3, smoothing))
  }
  
  return(physData)
  
}

trial_gripping <- function(physData, trialStart, trialLength, average = TRUE) {
  # this function will generate a data frame that allows to compare the grip strength per choice
  # -
  # - trialLength: in ms, can be handling, or how many samples you want to select (scalar)
  
  # if a single sample set is to be retrieved per trial, create a vector to index below
  if (length(trialLength) == 1) {
    trialLength <- rep(trialLength, length(trialStart))
  }
  
  # empty variables to accrue windows of data
  allSections <- numeric()
  allTrials <- numeric()
  
  # for each trial, get a gripping section (averaged if needed), and an equivalent trial number vector
  for (trial in seq_along(trialStart)) {
    # select this trial's sample
    start <- trialStart[trial]
    end <- start + trialLength[trial]
    section <- physData[start:end]
    
    # average if needed
    if (average) {
      section <- mean(section, na.rm = T)
      trialN <- trial
    } else {
      trialN <- rep(trial, length(section))
    }
    
    # append section
    allSections <- c(allSections, section)
    allTrials <- c(allTrials, trialN)
    
  }
  
  # put it together
  df <- data.frame(Trial = allTrials,
                   GripStrength = allSections)
  
  return(df)
  
}


###---------- analysis ---------
# load data
biopacData <- read.table('378_021318.txt')$V1
behavData <- read.csv('378_phys_13022018_log.csv')

# cut the session in half (due to incomplete n of trials from biopac)
# set a vector with trial-start times
# and round to the nearest ms 
halfData <- behavData %>% 
  filter(Block < 4) %>%
  mutate(trialStart = Experiment.Time - RT,
         trialStart_ms = round(trialStart, digits = 3) * 1000,
         trialEnd_ms = round(Experiment.Time, digits = 3) * 1000)

# isolate relevant grip time series and smooth if desired (filter > 0)
physData <- select_signal(biopacData, halfData$trialStart_ms[1], max(halfData$trialEnd_ms))

# isolate only completed trials and get their grip strength values
trialData <- halfData %>% 
  filter(Choice == 1) %>% 
  select(Handling, Offer, trialStart_ms) %>%
  mutate(Handling_ms = Handling * 1000)

# mean gripping strength per chosen trial
temp <- trial_gripping(physData, trialData$trialStart_ms, trialData$Handling_ms, average = T)

# combine to look at relationships
trialData <- cbind(trialData, temp)

# plot stuff
plot(as.factor(trialData$Handling), trialData$GripStrength, xlab = "Handling Time", ylab = "Grip Strength")






