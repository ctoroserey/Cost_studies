# simulating decisions from Garrett and Daw 2019
# this is just a playground to estimate how models affect their experiment
# this was helpful in understanding their approach
# and really, I did it to motivate myself to do it with my own data

# there are 4 trial types, spanning low to high reward rates
trialType <- list(HDLR = c(rep(0, 8), 20),
                  HDHR = c(rep(0, 8), 80),
                  LDLR = c(rep(0, 2), 20),
                  LDHR = c(rep(0, 2), 80))

# earning rates (pts/time) per trial
trialTimes <- list(HDLR = 8, LDHR = 2, HDHR = 8, LDLR = 2)
trialRewards <- list(HDLR = 20, LDHR = 80, HDHR = 80, LDLR = 20)
  
# blocks had different environmental richness based on rate of trial occurrences (per 7 trials)
richTrials <- c("HDLR", "HDHR", "LDLR", rep("LDHR", 4))
poorTrials <- c("LDHR", "HDHR", "LDLR", rep("HDLR", 4))

# optimal policy for each block
richAccept <- list(HDLR = F, LDHR = T, HDHR = F, LDLR = F)
poorAccept <- list(HDLR = F, LDHR = T, HDHR = T, LDLR = T)

# starting parameters for their Rescorla-Wagner-type model
# the model updates per second
# p = is the expected reward value
# alpha is the adaptation rate
p <- 0
alpha <- 0.005

# this is just to store a session's evolution of the reward estimate
pStore <- numeric()
evol <- data.frame(Reward = NA,
                   Delay = NA,
                   p = NA,
                   trialType = NA,
                   Choice = NA,
                   richOptimal = NA,
                   poorOptimal = NA)

# start simulating a series of 'n' trials in an 'env' environment
n <- 100
env <- poorTrials
for (run in seq(n)) {
  # randomly select the next 7 trials based on environmental richness
  sequence <- sample(env)
  
  # go through trials
  for (trial in sequence) {
    # if the trial's rate is higher than the estimated OC (average reward rate p * delay)
    # accept <- rich/poorAccept[[trial]] for optimality
    OC <- trialTimes[[trial]] * p
    accept <- ifelse(trialRewards[[trial]] > OC, T, F)
    
    # updated p at each second within the trial (one 0 = no update during offer window)
    if (accept) {
      for (r in trialType[[trial]]) {
        p <- p + alpha * (r - p) # could turn this into a function, so multiple models can be tested
        evol <- rbind(evol, c(r, trialTimes[[trial]], p, trial, accept, richAccept[[trial]], poorAccept[[trial]]))
      }
    } else {
      p <- p + alpha * (0 - p) # 0 is for the 1 second in the offer window
      evol <- rbind(evol, c(r, trialTimes[[trial]], p, trial, accept, richAccept[[trial]], poorAccept[[trial]]))
    }
  }
}

# plot the evolution of p
evol %>%
  filter(!is.na(p)) %>%
  mutate(p = as.numeric(p),
         ReceiptHDHR = ifelse((Reward > 0 & Choice == T & trialType == "HDHR"), -1.5, -6),
         ReceiptHDLR = ifelse((Reward > 0 & Choice == T & trialType == "HDLR"), -2.5, -6),
         ReceiptLDLR = ifelse((Reward > 0 & Choice == T & trialType == "LDLR"), -2, -6),
         ReceiptLDHR = ifelse((Reward > 0 & Choice == T & trialType == "LDHR"), -1, -6),
         optimalReceipt = ifelse((Reward > 0 & poorOptimal != Choice), -3, -6), 
         Time = seq(nrow(.))) %>%
  ggplot(aes(Time, p)) +
  geom_line() +
  ylim(-3, NA) +
  labs(title = paste("MVT reward rate evolution, alpha =", alpha)) +
  geom_point(aes(Time, optimalReceipt), pch = 21, color = "black", fill = "grey", size = 1.5) +
  geom_point(aes(Time, ReceiptLDHR, fill = trialType), pch = 21, color = "black", size = 1.5) +
  geom_point(aes(Time, ReceiptHDHR, fill = trialType), pch = 21, color = "black", size = 1.5) +
  geom_point(aes(Time, ReceiptLDLR, fill = trialType), pch = 21, color = "black", size = 1.5) +
  geom_point(aes(Time, ReceiptHDLR, fill = trialType), pch = 21, color = "black", size = 1.5) +
  theme_classic()



