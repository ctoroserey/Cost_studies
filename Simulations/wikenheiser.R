
library(tidyverse)
library(data.table)

# basic calculation
offer <- c(4, 8, 20)
handling <- c(14, 10, 2) + 2 # add 2 to H and T to account for offer and reward disclosure windows
travel <- c(2, 6, 14) + 2
offerRates_basic <- sapply(3:1, function(thresh) {(sum(tail(offer, n = thresh))) / ((handling * thresh) + (travel * 3))}) # good


# new stuff
pmat <- rbind(c(1, 1, 1),
              c(0, 1, 1),
              c(0, 0, 1))
A <- seq(0, 1.5, length.out = 4)


## single example comparing the two at a single acceptance strategy
p <- pmat[1, ]

# basic numerator
bn <- (offer %*% p)

# wikenheiser version
wn <- (offer %*% p) - ((offer %*% (1 - p)) * A)  

# basic denominator
bd <- ((handling * sum(p)) + (travel * 3)) #

# adapted wikenheiser numerator (since our participants don't select among delays, calculate per handling time)
wd <- (1 + (handling * sum(p)))


## loop through every possible combination 
# eventually mapply()
modeledRates <- list()

for (iter in seq_along(A)) {
  
  a <- A[iter]
  
  rstore <- matrix(0, 3, 3)
    
  for (prow in 1:nrow(pmat)) {
    
    # acceptances (binary per reward)
    p <- pmat[prow, ]
    
    # numerator from wikenheiser, denominator from us
    wn <- (offer %*% p) - ((offer %*% (1 - p)) * a)  
    bd <- ((handling * sum(p)) + (travel * 3))
    
    # model
    rstore[, prow] <- wn[1, 1] / bd

  }
  
  # clean up into a dataframe with info
  rstore <- as.data.frame(rstore)
  colnames(rstore) <- offer
  rstore$H <- handling
  modeledRates[[iter]] <- rstore
  
}

offerRates_wik <- do.call(rbind, modeledRates) %>%
  mutate(A = rep(A, each = 3)) %>%
  gather(O, Rate, -c(H:A)) %>%
  mutate(A = as.character(A),
         O = factor(O, levels = list(4, 8, 20)),
         H = factor(H, levels = list(4, 12, 16)))


## plots

# earnings per acceptance threshold
# as A increases, the earnings or OC per sec decrease
offerRates_wik %>%  ggplot(aes(O, Rate, group = H, color = H)) +
    geom_point(size = 3) +
    geom_line(size = 1.5) +
    scale_color_manual(values = c("purple", "grey50", "darkgoldenrod3")) +
    ylim(0, 1.2) +
    labs(y = "Expected reward per sec.", x = "Reward acceptance threshold") +
    facet_wrap(vars(A)) +
    theme(#legend.position = c(0.85, 0.2),
          legend.key = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 22))


# acceptances at each A leve
# fairly similar to our data at 0.5, but not quite
# we basically don't gain anything
offerRates_wik %>% 
  group_by(H, A) %>% 
  mutate(maxRate = max(Rate)) %>%
  ungroup() %>%
  mutate(wikChoice = ifelse(as.numeric(O) > (maxRate * as.numeric(H)), 1, 0)) %>%#,
         #OCChoice = ifelse(as.numeric(O) > (as.numeric(A) * as.numeric(H)), 1, 0)) %>%
  #gather(choiceType, Choice, wikChoice:OCChoice) %>%
  ggplot(aes(interaction(O, H), wikChoice)) + 
    geom_point(size = 3) + 
    geom_line(aes(group = interaction(H, A)), size = 1) +
    labs(x = "Offer.Handling", y = "Accepted") +
    facet_wrap(vars(A), ncol = 2) +
    theme(legend.key = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size = 12))










