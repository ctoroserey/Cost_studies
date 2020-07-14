
# using NLOPTR is pointless here
# instead I'm just going to feed any sensible combination of parameters to the model and minimize the negLL manually
# simpler, can evaluate models easier, and faster
# for tomorrow: have negloglik return the parameters per iteration
# the idea being that you can enter a model and data, and the function will return the lowest LL and associated parameters
# maybe also R-squares and things in the current optimization function

# simpler form of optimization that allows inputting any model expression into a single function call
# only good for static models
optimize_model <- function(subjData, params, model, simplify = F) {
  # this function finds the combination of parameter values that minimizes the neg log likelihood of a logistic regression
  # used to rely on NLOPTR, but it's too cumbersome for the low-dimensional estimates I'm performing.
  #
  # subjData: a participant's log
  # params: a list of vectors. Each vector is the possible values a given parameter can take. Names in list must match model expression
  # model: using `expr()`, define the model (use <param>[[1]] for free parameters to be estimated. R limitation.). Ex: expr(temp[[1]] * (reward - (gamma[[1]] * handling)))
  
  # extract basic choice information
  handling <- subjData$Handling
  reward <- subjData$Offer
  choice <- subjData$Choice
  rt <- subjData$RT
  cost <- subjData$Cost
  trial <- subjData$TrialN
  rawChoice <- subjData$rawChoice
  expTime <- round(subjData$ExpTime) + 0.1 # to avoid dividing 0 by 0
  blockTime <- subjData$blockTime
    
  # combine parameters into every possible combination
  params <- expand.grid(params)
  
  # Prep list of results to be returned
  out <- list()
  out$percentQuit <- mean(choice == 0) * 100
  out$percentAccept <- mean(choice == 1) * 100 
  
  LLs <- sapply(seq(nrow(params)), function(i) {
    
    # isolate the parameters for this iteration
    # and then store them as variables
    # FIGURE OUT HOW TO NOT STORE THEM AS DATAFRAMES
    pars <- params[i, ]
    lapply(seq_along(pars), function(variable) {assign(colnames(pars)[variable], pars[variable], envir = .GlobalEnv)})
    
    # estimate the probability of acceptance per the model
    p = 1 / (1 + exp(-eval(model)))
    p[p == 1] <- 0.999
    p[p == 0] <- 0.001
    
    # get the likelihood of the observations based on the model
    tempChoice <- rep(NA, length(choice))
    tempChoice[choice == 1] <- log(p[choice == 1])
    tempChoice[choice == 0] <- log(1 - p[choice == 0]) # log of probability of choice 1 when choice 0 occurred
    negLL <- -sum(tempChoice)
  })
  
  # chosen parameters  
  out$LL <- min(LLs)
  chosen_params <- params[which(LLs == out$LL), ]
  lapply(seq_along(chosen_params), function(variable) {assign(colnames(chosen_params)[variable], chosen_params[variable], envir = .GlobalEnv)})
  
  # Summarize the outputs
  out$LL0 <- -(log(0.5) * length(choice))
  out$Rsquared <- 1 - (out$LL / out$LL0) # pseudo r-squared, quantifying the proportion of deviance reduction vs chance
  out$loglikSpace <- LLs # in case you want to examine the concavity of the likelihood space
  out$probAccept <- 1 / (1 + exp(-eval(model)))
  out$Params <- chosen_params
  #out$predicted <- reward > out$subjOC
  #out$predicted[out$predicted == TRUE] <- 1
  #out$percentPredicted <- mean(out$predicted == choice) 
  
  # if doing this with dplyr::do(), return a simplified data.frame instead with the important parameters
  if (simplify) {
    out <- round(data.frame(out[-c(6, 7)]), digits = 2)
    colnames(out) <- c("percentQuit",
                       "percentAccept",
                       "LL",
                       "LL0",
                       "Rsq",
                       colnames(chosen_params))
  }
  
  return(out)
}

# for a trial-wise acceptance version. Add model expression eventually
optimize_model_dyn <- function(subjData, params, simplify = F) {
  # get every combination of parameters
  params <- expand.grid(params)
  
  # relevant behavior elements
  o <- subjData$Offer
  h <- subjData$Handling
  c <- subjData$Choice
  a <- o * c # accepted offers
  time <- subjData$ExpTime
  tau <- time - dplyr::lag(time, default = 0) # how long has it been since the last update?
  
  # Prep list of results to be returned
  out <- list()
  out$percentQuit <- mean(c == 0) * 100
  out$percentAccept <- mean(c == 1) * 100 
  
  # iterate through possible parameters and get the LL
  LLs <- sapply(seq(nrow(params)), function(i) {
    # isolate the parameters for this iteration
    # and then store them as variables
    tempr <- params[i, 1]
    alpha <- params[i, 2]
    
    # update rule (inspired by Constantino and Daw, 2015)
    gamma <- rep(0, nrow(subjData))
    
    # calculate gammas
    for (i in seq(nrow(subjData))) {
      if (i == 1) {
        gamma[i] <- 0.25
      } else {
        
        # delta <- (o[i] / h[i]) - gamma[i]
        # gamma[i + 1] <- gamma[i] + (1 - (1 - alpha) ^ h[1]) * delta
        gamma[i] <- (((1 - alpha) ^ tau[i]) * (a[i - 1] / tau[i])) + (1 - (1 - alpha) ^ tau[i]) * gamma[i - 1] # maybe expr(model)
      }
    }
    
    # estimate the probability of acceptance per the model
    p = 1 / (1 + exp(-(tempr * (o - (gamma * h)))))
    p[p == 1] <- 0.999
    p[p == 0] <- 0.001
    
    # get the likelihood of the observations based on the model
    tempChoice <- rep(NA, length(c))
    tempChoice[c == 1] <- log(p[c == 1])
    tempChoice[c == 0] <- log(1 - p[c == 0]) # log of probability of choice 1 when choice 0 occurred
    negLL <- -sum(tempChoice)
  } 
  )
  
  # chosen parameters  
  out$LL <- min(LLs)
  chosen_params <- params[which(LLs == out$LL), ]
  lapply(seq_along(chosen_params), function(variable) {assign(colnames(chosen_params)[variable], chosen_params[variable], envir = .GlobalEnv)})
  
  # Summarize the outputs
  out$LL0 <- -(log(0.5) * length(c))
  out$Rsquared <- 1 - (out$LL / out$LL0) # pseudo r-squared, quantifying the proportion of deviance reduction vs chance
  out$loglikSpace <- LLs # in case you want to examine the concavity of the likelihood space
  
  # get the optimized gamma to export the probability of acceptance
  # update rule (inspired by Constantino and Daw, 2015)
  gamma <- rep(0, nrow(subjData))
  
  # calculate gammas
  for (i in seq(nrow(subjData))) {
    if (i == 1) {
      gamma[i] <- 0.25
    } else {
      
      # delta <- (o[i] / h[i]) - gamma[i]
      # gamma[i + 1] <- gamma[i] + (1 - (1 - alpha) ^ h[1]) * delta
      gamma[i] <- (((1 - alpha[[1]]) ^ tau[i]) * (a[i - 1] / tau[i])) + (1 - (1 - alpha[[1]]) ^ tau[i]) * gamma[i - 1] # maybe expr(model)
    }
  }
  
  out$rate <- gamma
  out$probAccept <- 1 / (1 + exp(-(tempr[[1]] * (o - (gamma * h)))))
  out$Params <- chosen_params
  #out$predicted <- reward > out$subjOC
  #out$predicted[out$predicted == TRUE] <- 1
  #out$percentPredicted <- mean(out$predicted == choice) 
  
  # if doing this with dplyr::do(), return a simplified data.frame instead with the important parameters
  if (simplify) {
    out <- round(data.frame(out[-c(6, 7, 8)]), digits = 2)
    colnames(out) <- c("percentQuit",
                       "percentAccept",
                       "LL",
                       "LL0",
                       "Rsq",
                       colnames(chosen_params))
  }
  
  return(out)
  
}

# remove break time and start counting from 0
standardize_time <- function(subjData) {
  # remove the break time (variable across subjects) and start counting time from 0 (otherwise it can add physical effort calibration)
  breakTime <- min(subjData$ExpTime[subjData$Block == 4]) - max(subjData$ExpTime[subjData$Block == 3])
  subjData$ExpTime[which(subjData$Block > 3)] <- subjData$ExpTime[which(subjData$Block > 3)] - breakTime + 16
  subjData$ExpTime <- subjData$ExpTime - min(subjData$ExpTime)
  
  return(subjData)
}

# visually compare the values for a given parameter result across costs
param_compare_plot <- function(summary, param = "gamma", color = colsBtw, meanRate = 0.74) {
  plot <- ggplot(summary, aes_string("Cost", param, fill = "Cost")) +
    geom_hline(yintercept = meanRate, alpha = 0.9, size = 1, linetype = "dashed") + # mean of highest earning rates across blocks
    geom_boxplot(show.legend = F) +
    geom_jitter(width = 0.1, alpha = 0.5, show.legend = F, pch = 21, size = 3) +
    labs(x = "") +
    scale_fill_manual(values = color) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 16))
  
  return(plot)
}

## which models to run?
baseOC_nloptr <- F
bOC <- T
baseLogistic <- F # to test whether the brute search converges to a conventional logistic through glm()
fwOC <- F
dOC <- F
twOC <- T

# which experimental dataset?
data <- dataBtw %>% 
  filter(Cost != "Easy") %>%
  group_by(SubjID) %>% 
  do(standardize_time(.)) %>% 
  ungroup()



## ORIGINAL OC NLOPTR
if (baseOC_nloptr) {
  print("Running base OC model on NLOPTR...")
  
  summaryOC <- list()
  summaryOC$all <- data %>%
    group_by(Cost, SubjID) %>%
    do(optimizeOCModel(., simplify = T)) %>%
    ungroup()
}


## ORIGINAL OC
# great correspondence with NLOPTR, just much slower since it surveys the whole parameter space
if (bOC) {
  print("Running base OC model through grid search...")
  
  # model to be fit
  # make sure that you specify the inverse temperature
  # extra parameters as dfs for now, that's why the `[[1]]`
  model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
  
  # create a list with possible starting values for model parameters
  # parameter names must match model ones
  spaceSize <- 30
  params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                 gamma = seq(0.25, 1.5, length.out = spaceSize))
  
  # fit to each subject
  baseOC <- data %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model(., params, model_expr, simplify = T)) %>%
    ungroup()
}


## Fawcett OC (2012)
# nope.. the definition doesn't work for multi-choice foraging
if (fwOC) {
  print("Running Fawcett et al (2012) OC model through grid search...")
  
  # model to be fit
  # make sure that you specify the inverse temperature
  # extra parameters as dfs for now, that's why the `[[1]]`
  model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * (rev(expTime) - handling))))
  
  # create a list with possible starting values for model parameters
  # parameter names must match model ones
  spaceSize <- 30
  params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                 gamma = seq(0, 2, length.out = spaceSize))
  
  # fit to each subject
  fawcettOC <- data %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model(., params, model_expr, simplify = T)) %>%
    ungroup()
}


## BASIC LOGISTIC
# the results mostly match what glm() outputs
if (baseLogistic) {
  print("Running a simple logistic (Handling + Reward) through brute search...")
  # model to be fit
  model_expr <- expr(intercept[[1]] + (betaRwd[[1]] * reward) + (betaHand[[1]] * handling))
  
  # create a list with possible starting values for model parameters
  spaceSize <- 20
  params <- list(intercept = seq(-1, 1, length.out = spaceSize), 
                 betaRwd = seq(-5, 5, length.out = spaceSize),
                 betaHand = seq(-5, 5, length.out = spaceSize))
  
  # fit to each subject
  baseLogistic <- data %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model(., params, model_expr, simplify = T)) %>%
    ungroup()
}


## Dundon, Garrett, et al (2020)
# for a single subject, estimating the global gamma as usual is better than an evolving 
# one using eq 3 on their paper for cost2.
# the estimate using their equation overharvests
# transfer this into a functon to apply to others
# but also try on cost3 subjects
# if anything, it can be used as a mold to try out evolving models

# concluding: MVT-style rate tracking doesn't seem to apply here
# subjects have a clear idea of the rate. What's curious is that for 058 gamma predicts perfect choice
# though the number of trials correctly predicted is indeed higher than the Dundon ver
# this also seems to apply with within-subject peeps.
# fix the break time though. That shouldn't be there.

## Tracking the ongoing rate based on choices
# to track dynamicOC without choice, remove the lag and choice from the model
if (dOC) {
  print("Running ongoing OC (based on choice history) using a grid search...")

  # model to be fit
  # make sure that you specify the inverse temperature
  # extra parameters as dfs for now, that's why the `[[1]]`
  model_expr <- expr(tempr[[1]] * (reward - ((dplyr::lag(cumsum(reward * choice), default = 0) / (expTime ^ alpha[[1]])) * handling)))

  # create a list with possible starting values for model parameters
  # parameter names must match model ones
  spaceSize <- 30
  params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                 alpha = seq(0.25, 2, length.out = spaceSize))
  
  # fit to each subject
  dynamicOC <- data %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model(., params, model_expr, simplify = T)) %>%
    ungroup()
}

# param_compare_plot(dynamicOC, param = "alpha", meanRate = NaN)
# 
# # plot the ongoing results
# id <- "105"
# 
# # plot dynamicOC with choice 
# # to plot dynamicOC without choice, remove the lag and choice from cumulativeRate
# data %>%
#   filter(SubjID == id) %>%
#   left_join(dynamicOC, by = "SubjID") %>%
#   mutate(#cumulativeRate = (dplyr::lag(cumsum(Offer * Choice), default = 0) / ExpTime) ^ eta,
#     #cumulativeRate = ifelse(is.nan(cumulativeRate), 0, cumulativeRate),
#     cumulativeRate = gamma * (Handling / eta),
#     trialRate = Offer, #trialRate = ifelse(trialRate > 3, 3, trialRate)
#     fitChoice = ifelse(trialRate > cumulativeRate, 0.5, -5),
#     newChoice = ifelse(Choice == 0, -5, Choice)) %>% 
#   ggplot(aes(TrialN, cumulativeRate)) +
#   geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
#   geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
#   geom_hline(yintercept = filter(baseOC, SubjID == id)$gamma) + # single gamma estimated for an individual
#   geom_hline(yintercept = 0.74, linetype = "dashed", color = "grey30") + # mean optimal rate across blocks
#   geom_line(aes(color = Handling), size = 0.5) +
#   geom_point(aes(color = Handling), size = 1.2) +
#   geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
#   geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
#   annotate("text", x = 220, y = filter(baseOC, SubjID == id)$gamma + 0.25, label = "Fitted \n Gamma", size = 5) +
#   annotate("text", x = 220, y = 0.55, label = "Optimal", size = 5, color = "grey30") +
#   scale_fill_discrete(name = "Offer") +
#   scale_color_continuous(breaks = c(2, 10, 14), labels = c(2, 10, 14)) +
#   ylim(0, NA) +
#   labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 16))



## trial-wise updating of gamma
# odds: 426, 420
if (twOC) {
  print("Running trial-wise OC model through grid search...")
  
  # create a list with possible starting values for model parameters
  # parameter names must match model ones
  spaceSize <- 30
  params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                 alpha = seq(0, 1, length.out = spaceSize))
  
  # fit to each subject
  trialwiseOC <- data %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model_dyn(., params, simplify = T)) %>%
    ungroup()
}

param_compare_plot(trialwiseOC, "tempr")

#plot



