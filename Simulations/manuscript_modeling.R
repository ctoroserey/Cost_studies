# OC-specific computation of the negative log likelihood for model optimization
negLogLik <- function(params, choice, handling, reward) {
  gamma <- exp(params[2]) * handling
  model <- exp(params[1]) * (reward - gamma)
  p = 1 / (1 + exp(-model))
  p[p == 1] <- 0.999
  p[p == 0] <- 0.001
  tempChoice <- rep(NA, length(choice))
  tempChoice[choice == 1] <- log(p[choice == 1])
  tempChoice[choice == 0] <- log(1 - p[choice == 0]) # log of probability of choice 1 when choice 0 occurred
  negLL <- -sum(tempChoice)
  return(negLL)
}

model <- expr(exp(params$temp) * (reward - (params$gamma * handling)))

negLogLik <- function(params, choice, handling, reward) {
  gamma <- exp(params[2]) * handling
  #model <- exp(params[1]) * (reward - gamma)
  p = 1 / (1 + exp(-eval(model)))
  p[p == 1] <- 0.999
  p[p == 0] <- 0.001
  tempChoice <- rep(NA, length(choice))
  tempChoice[choice == 1] <- log(p[choice == 1])
  tempChoice[choice == 0] <- log(1 - p[choice == 0]) # log of probability of choice 1 when choice 0 occurred
  negLL <- -sum(tempChoice)
  return(negLL)
}



# optimize the OC model
optimizeOCModel <- function(Data, Algorithm = "NLOPT_LN_NEWUOA", simplify = F, model = negLogLik, params = ) {

  # Data: The participant's log
  # Algorithm: probably let be
  # model: an external function to minimize (in this case OC, separately defined as negloglik)

  # Prep data
  handling <- Data$Handling
  reward <- Data$Offer
  choice <- Data$Choice

  # Prep list of results to be returned
  out <- list()
  out$percentQuit <- mean(choice == 0) * 100
  out$percentAccept <- mean(choice == 1) * 100
  miss <- (choice != 1) & (choice != 0)
  out$percentMiss <- mean(miss)  * 100
  choice <- choice[!miss]
  reward <- as.numeric(reward[!miss])
  handling <- as.numeric(handling[!miss])

  # Establish lower and upper bounds
  LB <- round(log((min(reward)/max(handling)) * 0.99), digits = 4)
  UB <- round(log((max(reward)/min(handling)) * 1.01), digits = 4) # in reality this should be the second largest, since no one would reject the highest val

  # Begin defining parameters
  # If choices are one-sided (i.e. all accepted), assign the upper or lower bound
  if ((sum(choice) == length(choice)) | (sum(choice) == 0)) {
    ifelse(sum(choice) == length(choice), out$Gamma <- exp(LB), out$Gamma <- exp(UB))
    out$temperature <- 1 # it was NA, but in theory a temperature of 1 also indicates noiseless estimates, and allows for easier fit computations
    out$LL <- 0
  } else {
    # Create a feasible region (search space)
    params <- as.matrix(expand.grid(temperature = c(-1, 1), gamma = seq(LB, UB, length = 3)))
    # Create a list to check the minimization of the negative log lik.
    info <- list()
    info$negLL <- Inf
    # Define the options to be used during optimization
    # Consider looking into other optimization algorithms and global minima
    opts <- list("algorithm" = Algorithm,
                 "xtol_rel" = 1.0e-8)
    # Optimize the sOC over all possible combinations of starting points
    for (i in seq(nrow(params))) {
      tempInfo <- nloptr(x0 = params[i, ],
                         eval_f = model,
                         lb = c(log(0.001), LB),
                         ub = c(-log(0.001), UB),
                         opts = opts,
                         choice = choice,
                         handling = handling,
                         reward = reward)
      if (tempInfo$objective < info$negLL) {
        #print("Minimized")
        info$negLL <- tempInfo$objective
        info$params <- tempInfo$solution
      }
    }
    out$Gamma <- exp(info$params[2])
    out$temperature <- exp(info$params[1])
    out$LL <- -info$negLL
  }

  # Summarize the outputs
  out$LL0 <- log(0.5) * length(choice)
  out$Rsquared <- 1 - (out$LL/out$LL0) # pseudo r-squred, quantifying the proportion of deviance reduction vs chance
  out$subjOC <- out$Gamma * handling
  out$probAccept <- 1 / (1 + exp(-(out$Scale*(reward - out$subjOC))))
  out$predicted <- reward > out$subjOC
  out$predicted[out$predicted == TRUE] <- 1
  out$percentPredicted <- mean(out$predicted == choice)

  # adjust the probabilities in case of extreme gammas
  if (out$Gamma <= exp(LB)) {
    out$Gamma <- exp(LB) # temporary condition because Nlopt is not respecting the lower bound
    out$probAccept <- rep(1, length(choice))
  } else if (out$Gamma == exp(UB)) {
    out$probAccept <- rep(0, length(choice))
  }

  # if doing this with dplyr::do(), return a simplified data.frame instead with the important parameters
  if (simplify) {
    out <- data.frame(out[c(seq(8), 12)])
  }

  return(out)

}

summaryOC <- list()
summaryOC$all <- dataWth %>%
  mutate(Block = case_when(
    Block %in% c(1, 2) ~ 1,
    Block %in% c(3, 4) ~ 2,
    Block %in% c(5, 6) ~ 3
  )) %>%
  group_by(SubjID, Cost, Half) %>% #group_by(SubjID, Cost, Block) %>%
  do(optimizeOCModel(., simplify = T)) %>%
  ungroup()

# plot
(summaryOC$plot <- ggplot(summaryOC$all, aes(Cost, Gamma, fill = Cost)) +
  geom_hline(yintercept = 0.7, alpha = 0.9, color = "gray40", size = 1, linetype = "dashed") +
  geom_jitter(pch = 21, size = 3, show.legend = F) +
  ylim(0, 1.5) +
  labs(x = "") +
  scale_fill_manual(values = colsWth) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 16)))

summaryOC$all %>%
  select(SubjID, Cost, Gamma, Half) %>%
  spread(Half, Gamma) %>%
  ggplot(aes(Half_1, Half_2, fill = Cost)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point(pch = 21, color = "black", size = 3) +
    ylim(0, 2) +
    xlim(0, 2) +
    theme_minimal()




# create a list with possible starting values for model parameters
params <- list(temp = c(-1, 1), 
          gamma = seq(0.25, 1.5, length.out = 5), 
          time = seq(0, 1, length.out = 5))

# combine into every possible combination
params_grid <- as.matrix(expand.grid(params))






























