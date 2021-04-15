# This script will browse the attention logs and print out who doesn't meet criteria 
# Per pre-registration, anyone taking over 3 secs to respond is considered inattentive
# and should be removed from the participant pool

# Quietly load tidyverse
suppressMessages(library(tidyverse))

# change directory and find the attention files
setwd('./data/')
files <- dir(pattern = "attention")

# load them, find subjects with 2 or more attention trials longer than 3 secs, and print to terminal
t1 <- tibble(SubjID = files) %>% 
  mutate(content = map(SubjID, ~read_csv(., col_names = F, col_types = cols())),
         SubjID = substring(SubjID, 11, 13)) %>%
  unnest() %>% 
  group_by(SubjID) %>%
  summarize(inattention = sum(floor(X2) > 3)) %>%
  filter(inattention > 1) %>%
  print(SubjID)

# print out the current number of valid participants
finalN <- length(files) - nrow(t1)

sprintf("Current valid n of subjects: %s", finalN)