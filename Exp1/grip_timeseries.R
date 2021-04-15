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
# how to best track potential "fatigue" over time?
# plot average time-wise grip strength across trials per reward, 
# or mean grip for first and second halves of each trial.Also divide by reward

###---------- library and functions ---------
library(tidyverse)
library(data.table)
library(parallel)

select_signal <- function(physData, start, end, thresh = 0, smoothing = 0, id = NA, trialLength) {
  # this function will use the experimental timing from a log to get rid of non-experimental recordings
  # this will be recordings prior to the beginning of the session and when the trials stopped
  # it currently assumes milisecond rate biopac acquisition
  # variables:
  # - physData: the vector with grip strength
  # - start/end: just the time bounds (in ms) containing relevant data 
  # - thresh: sustained threshold for the grip press (can be estimated by plotting the time series)
  # - smoothing: if the data is noisy, apply smoothing
  # - id: if a subjid is provided, it returns a df
  # - trialLength: to ensure that the first grip is from a trial, provide how long the trial should be
  
  # the result is a time series vector of grip strength that is aligned with the behavioral log
  
  # where does the first trial start?
  # adding the first handling helps ensure that only trials are selected (some p's practice gripping)
  if (missing(trialLength)) {
    firstTrial <- which(physData > thresh)[1]
  } else {
    # this is messy
    # basically, count the consecutive number of epochs below and above 0
    # then select the first positive one that matches the trialLength
    # then sum the number of epochs prior to that, and that's where the trial begins
    thrsh <- physData > thresh
    cnt <- rle(thrsh)
    lengths <- tibble(len = cnt$lengths, 
                      vals = cnt$values, 
                      n = seq(length(cnt[[1]]))) %>% 
                      mutate(len2 = floor(len / 1000)) %>% 
                      filter(len2 %in% seq(trialLength, trialLength + 2), vals == TRUE)
    
    firstTrial <- sum(cnt$lengths[seq(lengths$n[1] - 1)]) # sum the n of epochs prior to trial 1
  }
  
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
  
  if (!is.na(id)) {
    physData = data.frame(SubjID = rep(id, length(physData)),
                          signalStrength = physData)
  }
  
  return(physData)
  
}

trial_gripping <- function(physData, trialStart, trialLength, average = TRUE, extras, resample = 0, preData = 0) {
  # this function will generate a data frame that allows to compare the grip strength per choice
  # 
  # - trialLength: in ms, can be handling, or how many samples you want to select (scalar)
  # - average: get the mean grip strength per trial
  # - extras: data frame. extra trial characteristics to add (e.g. reward amount, length, etc) 
  # - resample: if the ts is too long, set adelta to capture (i.e. every 100th row)
  
  # if a single sample set is to be retrieved per trial, create a vector to index below
  if (!missing(trialLength)) {
    trialLength <- rep(trialLength, length(trialStart))
  }
  
  # empty variables to accrue windows of data
  allSections <- numeric()
  allTrials <- numeric()
  allTime <- numeric()
  
  # for each trial, get a gripping section (averaged if needed), and an equivalent trial number vector
  for (trial in seq_along(trialStart)) {
    # select this trial's sample
    start <- trialStart[trial]
    end <- start + trialLength[trial] -1
    section <- physData[(start - preData):end]
    
    # average if requested
    if (average) {
      section <- mean(section, na.rm = T)
      trialN <- trial
      trialTime <- length(section)
    } else {
      trialN <- rep(trial, length(section))
      trialTime <- seq(length(section))
    }
    
    # append section
    allSections <- c(allSections, section)
    allTrials <- c(allTrials, trialN)
    allTime <- c(allTime, trialTime)
    
  }
  
  # put it together
  df <- data.frame(Trial = allTrials,
                   TrialTime = allTime,
                   GripStrength = allSections)
  
  
  # if you have additional metrics to add
  if (! missing(extras)) {
    # get the length of each trial
    counts <- dplyr::count(df, Trial)$n
    
    # repeat each metric to match the length of each trial
    extras <- extras %>% uncount(counts)
    
    # append to grip df
    df <- cbind(extras, df)
  }
  
  # resample if needed
  if (resample > 0) {
    df <- df[seq(1, nrow(df), by = resample), ]
    df <- df %>% 
      group_by(Handling, Offer, Trial) %>%
      mutate(ResampledTime = seq(length(Trial)))
  }
  
  return(df)
  
}


###---------- analysis ---------
### data loading and cleaning up
# recent addition: turning forced travels to acceptances, as they reflect a preference for that trial
# rawChoice will be used to compute the # of mistakes
# The RT is upper-bounded because a glitch in the code made two 14s trials last longer (among all p's) 
setwd("./data")
files <- dir(pattern = '_log.csv')

behavData <- tibble(SubjID = files) %>%
  mutate(contents = map(SubjID, ~ read_csv(., col_types = cols())))  %>%
  mutate(Cost = substring(SubjID, 5, 8),
         Cost = case_when(Cost == "wait" ~ "Wait",
                          Cost == "cogT" ~ "Cognitive",
                          Cost == "phys" ~ "Physical",
                          Cost == "phea" ~ "Easy"),
         SubjID = as.integer(substring(SubjID, 0, 3))) %>%
  unnest() %>%
  mutate(rawChoice = Choice,
         RT = ifelse(RT > 14.1, 14, RT),
         Choice = ifelse(Choice == 2, 1, Choice),
         Half = ifelse(Block < 4, "Half_1", "Half_2"),
         Cost = factor(Cost, levels = c("Physical", "Cognitive", "Wait", "Easy")),
         optimal = case_when(
           (Handling == 10 & Offer < 8) ~ 0,
           (Handling == 14 & Offer < 20) ~ 0,
           TRUE ~ 1
         )) 

# get a simple subject list and the number of subjects
subjList <- unique(behavData$SubjID)
nSubjs <- length(subjList)

# load the biopac data and create a shorter list of subjects with grip data
setwd('../../biopac/cost2/')
files <- dir(pattern = '.txt')
shortList <-sapply(files, substr, 1, 3)

biopacData <- tibble(SubjID = files) %>%
  mutate(contents = map(SubjID, ~ fread(.)))  %>%
  mutate(SubjID = as.integer(substring(SubjID, 0, 3))) %>%
  unnest()

# keep only the behavioral data for participants with biopac data
# cut the session in half (due to incomplete n of trials from biopac)
# set a vector with trial-start times
# and round to the nearest ms 
halfData <- behavData %>% 
  filter(SubjID %in% shortList,
         Block < 4) %>%
  mutate(trialStart = Experiment.Time - RT,
         trialStart_ms = round(trialStart, digits = 3) * 1000,
         trialEnd_ms = round(Experiment.Time, digits = 3) * 1000)

# isolate relevant grip time series and smooth if desired (filter > 0)
#physData <- select_signal(biopacData, halfData$trialStart_ms[1], max(halfData$trialEnd_ms))

# create a smaller dataset to append the start and end of each session to biopacData
temp <- halfData %>% 
  group_by(SubjID) %>% 
  summarise(sessionStart = min(trialStart_ms),
            sessionEnd = max(trialEnd_ms),
            firstHandling = Handling[1])

# select gripping signal from the session
biopacData <- biopacData %>%
  left_join(temp, by = "SubjID") %>%
  filter(SubjID != 419) %>% # some subjects' data are atypical, with low gripping values. Add those here during testing.
  plyr::dlply("SubjID", identity)

biopacData <- mclapply(biopacData, function(data) {
  select_signal(data$V1, 
                unique(data$sessionStart), 
                unique(data$sessionEnd), 
                trialLength = unique(data$firstHandling))
})


# isolate only completed trials and get their grip strength values
trialData <- halfData %>% 
  filter(Choice == 1,
         SubjID != 419) %>% 
  select(SubjID, Handling, Offer, trialStart_ms) %>%
  mutate(Handling_ms = Handling * 1000) %>%
  plyr::dlply("SubjID", identity)

# gripping strength per chosen trial (mean gripping possible)
trialGrip <- mclapply(seq_along(biopacData), function(i) {
  trial_gripping(biopacData[[i]], 
                 trialData[[i]]$trialStart_ms, 
                 trialData[[i]]$Handling_ms, 
                 average = F, 
                 extras = trialData[[i]], 
                 resample = 100,
                 preData = 3000)
})

# get mean grip per offer x handling x subject into 1 df
trialGrip_all <- mclapply(trialGrip, function(data) {
  data %>%
    group_by(SubjID, Handling, Offer, ResampledTime) %>%
    summarise(mGrip = mean(GripStrength)) %>%
    ungroup()
  })
trialGrip_all <- do.call(rbind, trialGrip_all)


### plot checks
# mean gripping per handling x offer across subjects
# trialGrip_all %>%
#   group_by(Handling, Offer, TrialTime) %>%
#   summarise(meanGrip = mean(mGrip, na.rm = T),
#             se = sd(mGrip, na.rm = T) / nSubjs) %>%
#   ggplot(aes(TrialTime, meanGrip, group = as.factor(Offer), color = as.factor(Offer), fill = as.factor(Offer))) +
#     #geom_smooth() +
#     geom_line() +
#     geom_line(aes(y = meanGrip - se), linetype = "dashed") +
#     geom_line(aes(y = meanGrip + se), linetype = "dashed") +
#     #geom_ribbon(ymin = ymin, ymax = ymax, alpha = 0.3) +
#     facet_wrap(vars(Handling)) +
#     theme_classic()


# individuals
# i <- 9
# 
# trialGrip[[i]] %>%
#   mutate(ResampledTime = seq(nrow(.)),
#          Offer = as.factor(Offer)) %>%
#   #filter(Handling == 10) %>%
#   ggplot(aes(TrialTime, GripStrength, group = Offer, color = Offer, fill = Offer)) +
#   geom_point(size = 0.5) +
#   geom_smooth() +
#   #ylim(0, NA) +
#   facet_wrap(vars(Handling)) +
#   theme_classic()

# good examples: 1, 3, 4, 7, 9, 10, 11, 15, 19, 21, 24, 30
# bad: 6, 8, 12, 14, 16, 17, 18, 20, 22, 23, 25, 26, 27, 28, 29, 31

# selected


# good subs
subs <- c(1, 3, 4, 7, 9, 10, 11, 15, 19, 21, 24, 30)



trialGrip_all %>%
  filter(SubjID %in% unique(SubjID)[subs]) %>%
  mutate(Offer = as.factor(Offer)) %>%
  filter(Handling == 14,
         Offer != 4) %>%
  group_by(Handling, Offer, ResampledTime) %>%
  summarise(meanGrip = mean(mGrip, na.rm = T),
            se = sd(mGrip, na.rm = T) / length(unique(SubjID))) %>%
  ggplot(aes(ResampledTime, meanGrip, group = Offer, color = Offer, fill = Offer)) +
    #geom_point(pch = 21, color = "black", alpha = 0.3) +
    geom_line(size = 1.5) +
    geom_vline(xintercept = 30, linetype = "dashed", alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, size = 1.5) +
    # geom_smooth(aes(ResampledTime, mGrip, color = Offer, fill = Offer, group = Offer)) +
    geom_line(aes(y = meanGrip - se), linetype = "dashed", size = 1) +
    geom_line(aes(y = meanGrip + se), linetype = "dashed", size = 1) +
    #facet_wrap(vars(Handling)) +
    labs(x = "Time in Trial (sec)", y = "Mean Gripping Strength (a.u.)") +
    scale_x_continuous(breaks = seq(0, 170, by = 10), labels = seq(-3, 14)) +
    ylim(NA, 2.5) +
    theme(legend.position = c(0.9, 0.25),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size = 22))






