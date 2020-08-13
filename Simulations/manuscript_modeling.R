
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
  out$pAccept <- 1 / (1 + exp(-eval(model)))
  out$params <- chosen_params
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
optimize_model_dyn <- function(subjData, params, simplify = F, gammaStart = 0) {
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
    for (i in seq(nrow(subjData) - 1)) {
      if (i == 1) {
        gamma[i] <- gammaStart
      } else {
        
        delta <- (a[i - 1] / tau[i]) - gamma[i - 1]
        gamma[i] <- gamma[i - 1] + (1 - (1 - alpha) ^ tau[i]) * delta
        gamma[i] <- (((1 - alpha) ^ tau[i]) * (a[i - 1] / tau[i])) + (1 - (1 - alpha) ^ tau[i]) * gamma[i - 1] # maybe expr(model)
        gamma[i] <- ((1 - (1 - alpha) ^ tau[i]) * (a[i - 1] / tau[i])) + ((1 - alpha) ^ tau[i]) * gamma[i - 1] # switched version that means higher alpha = more learning
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
  for (j in seq(nrow(subjData))) {
    if (j == 1) {
      gamma[j] <- gammaStart
    } else {
      delta <- (a[j - 1] / tau[j]) - gamma[j - 1]
      gamma[j] <- gamma[j - 1] + (1 - (1 - alpha[[1]]) ^ tau[j]) * delta
      gamma[j] <- (((1 - alpha[[1]]) ^ tau[j]) * (a[j - 1] / tau[j])) + (1 - (1 - alpha[[1]]) ^ tau[j]) * gamma[j - 1] # maybe expr(model)
      gamma[j] <- ((1 - (1 - alpha[[1]]) ^ tau[j]) * (a[j - 1] / tau[j])) + ((1 - alpha[[1]]) ^ tau[j]) * gamma[j - 1] # maybe expr(model)
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

# variants of vanilla C&D
optimize_model_dyn_us <- function(subjData, params, simplify = F, gammaStart = 0) {
  # get every combination of parameters
  params <- expand.grid(params)
  
  # relevant behavior elements
  o <- subjData$Offer
  h <- subjData$Handling
  t <- 20 - h
  c <- subjData$Choice
  a <- o * c # accepted offers
  time <- subjData$ExpTime
  
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
    k <- params[i, 3]
    #tau <- (time - dplyr::lag(time, default = 0)) ^ k
    tau <- lag((h ^ k * c) + t) # we dont expect cognitive to affect total time, just the handling time
    
    # update rule (inspired by Constantino and Daw, 2015)
    gamma <- rep(0, nrow(subjData))
    
    # calculate gammas
    for (i in seq(nrow(subjData) - 1)) {
      if (i == 1) {
        gamma[i] <- gammaStart
      } else {
        
        # delta <- (o[i] / h[i]) - gamma[i]
        # gamma[i + 1] <- gamma[i] + (1 - (1 - alpha) ^ h[1]) * delta
        #gamma[i] <- (((1 - alpha) ^ tau[i]) * (a[i - 1] / tau[i])) + (1 - (1 - alpha) ^ tau[i]) * gamma[i - 1]
        gamma[i] <- ((1 - (1 - alpha) ^ tau[i]) * (a[i - 1] / tau[i])) + ((1 - alpha) ^ tau[i]) * gamma[i - 1] # switched version that means higher alpha = more learning
      }
    }
    
    # estimate the probability of acceptance per the model
    p = 1 / (1 + exp(-(tempr * (o - (gamma * h ^ k)))))
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
  #tau <- (time - dplyr::lag(time, default = 0)) ^ k[[1]]
  tau <- lag((h ^ k[[1]] * c) + t)
  
  # calculate gammas
  for (j in seq(nrow(subjData))) {
    if (j == 1) {
      gamma[j] <- gammaStart
    } else {
      #gamma[j] <- (((1 - alpha[[1]]) ^ tau[j]) * (a[j - 1] / tau[j])) + (1 - (1 - alpha[[1]]) ^ tau[j]) * gamma[j - 1]
      gamma[j] <- ((1 - (1 - alpha[[1]]) ^ tau[j]) * (a[j - 1] / tau[j])) + ((1 - alpha[[1]]) ^ tau[j]) * gamma[j - 1]
    } 
  }
  
  out$rate <- gamma
  out$probAccept <- 1 / (1 + exp(-(tempr[[1]] * (o - (gamma * h ^ k[[1]])))))
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

# plot the results from the dynamic model for a single individual
plot_dyn <- function(id = "58", exp = "btw", gammaOne = 0, showChoices = -0.6) {
  
  if (exp == "btw") {
    # choose subject + params
    sub <- filter(dataBtw, SubjID == id)
    spaceSize <- 50
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   alpha = seq(0, 1, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn(sub, params, simplify = F, gammaStart = gammaOne)
    
    # get baseline gamma from a simple fit
    # model to be fit
    # make sure that you specify the inverse temperature
    # extra parameters as dfs for now, that's why the `[[1]]`
    model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
    
    # create a list with possible starting values for model parameters
    # parameter names must match model ones
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    
    # fit to subject
    baseOC <- optimize_model(sub, params, model_expr, simplify = T)
    
    # plot
    ratePlot <- sub %>%
      mutate(earningRate = dplyr::lag(cumsum(Offer * Choice), default = 0) / ExpTime,
             trialRate = Offer / Handling,
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$probAccept,
             g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
        geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
        geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
        geom_hline(yintercept = baseOC$gamma) + # single gamma estimated for an individual
        geom_hline(yintercept = 0.74, linetype = "dashed", color = "grey30") + # mean optimal rate across blocks
        geom_line(aes(color = Handling), size = 0.5) +
        geom_point(aes(color = Handling), size = 1.2) +
        geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        #geom_point(aes(TrialN, pChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        annotate("text", x = max(sub$TrialN) + 8, y = baseOC$gamma + 0.25, label = "Fitted \n Gamma", size = 5) +
        annotate("text", x = max(sub$TrialN) + 8, y = 0.55, label = "Optimal", size = 5, color = "grey30") +
        annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
        annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
        scale_fill_discrete(name = "Offer") +
        scale_color_continuous(breaks = c(2, 10, 14), labels = c(2, 10, 14)) +
        #geom_line(aes(TrialN, earningRate)) +
        ylim(showChoices, NA) +
        labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 16))
    
  } else if (exp == "wth") {
    # choose subject + params
    sub <- filter(dataWth, SubjID == id)
    spaceSize <- 100
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   alpha = seq(0, 1, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn(sub, params, simplify = F, gammaStart = gammaOne)
    
    # plot
    ratePlot <- sub %>%
      mutate(earningRate = dplyr::lag(cumsum(Offer * Choice), default = 0) / ExpTime,
             trialRate = Offer / Handling,
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$probAccept,
             g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
        geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
        geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
        geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey30") + # mean optimal rate across blocks
        geom_line(size = 0.5, color = "grey50") +
        geom_point(aes(color = Cost), size = 1.2) +
        geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        annotate("text", x = max(sub$TrialN) + 8, y = 0.63, label = "Optimal", size = 5, color = "grey30") +
        annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
        annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
        scale_fill_discrete(name = "Offer") +
        scale_color_manual(values = colsWth) +
        #geom_line(aes(TrialN, earningRate)) +
        ylim(showChoices, NA) +
        labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 16))
  }
  
  suppressWarnings(print(ratePlot))
}
plot_dyn_us <- function(id = "58", exp = "btw", gammaOne = 0) {
  
  if (exp == "btw") {
    # choose subject + params
    sub <- filter(dataBtw, SubjID == id)
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   alpha = seq(0, 1, length.out = spaceSize),
                   k = seq(-1, 1, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn_us(sub, params, simplify = F, gammaStart = gammaOne)
    
    # get baseline gamma from a simple fit
    # model to be fit
    # make sure that you specify the inverse temperature
    # extra parameters as dfs for now, that's why the `[[1]]`
    model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
    
    # create a list with possible starting values for model parameters
    # parameter names must match model ones
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    
    # fit to subject
    baseOC <- optimize_model(sub, params, model_expr, simplify = T)
    
    # plot
    ratePlot <- sub %>%
      mutate(trialRate = Offer / Handling,
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$probAccept,
             g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
        geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
        geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
        geom_hline(yintercept = baseOC$gamma) + # single gamma estimated for an individual
        geom_hline(yintercept = 0.74, linetype = "dashed", color = "grey30") + # mean optimal rate across blocks
        geom_line(aes(color = Handling), size = 0.5) +
        geom_point(aes(color = Handling), size = 1.2) +
        geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        #geom_point(aes(TrialN, pChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        annotate("text", x = max(sub$TrialN) + 8, y = baseOC$gamma + 0.25, label = "Fitted \n Gamma", size = 5) +
        annotate("text", x = max(sub$TrialN) + 8, y = 0.55, label = "Optimal", size = 5, color = "grey30") +
        annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
        annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
        scale_fill_discrete(name = "Offer") +
        scale_color_continuous(breaks = c(2, 10, 14), labels = c(2, 10, 14)) +
        ylim(0, NA) +
        labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 16))
    
  } else if (exp == "wth") {
    # choose subject + params
    sub <- filter(dataWth, SubjID == id)
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   alpha = seq(0, 1, length.out = spaceSize),
                   k = seq(-1, 1, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn_us(sub, params, simplify = F, gammaStart = gammaOne)
    
    # plot
    ratePlot <- sub %>%
      mutate(trialRate = Offer / Handling,
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$probAccept,
             g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
        geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
        geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
        geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey30") + # mean optimal rate across blocks
        geom_line(size = 0.5, color = "grey50") +
        geom_point(aes(color = Cost), size = 1.2) +
        geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
        annotate("text", x = max(sub$TrialN) + 8, y = 0.63, label = "Optimal", size = 5, color = "grey30") +
        annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
        annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
        scale_fill_discrete(name = "Offer") +
        scale_color_manual(values = colsWth) +
        ylim(0, NA) +
        labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 16))
  }
  
  suppressWarnings(print(ratePlot))
}

# to test the effect of parameters 
plot_alphas <- function(alphas, k = 1, exp = "btw", gammaStart = 0.5) {
  if (exp == "btw") {
    # choose sample subject + params
    sub <- filter(dataBtw, SubjID == "58")
    
    # relevant behavior elements
    o <- sub$Offer
    h <- sub$Handling
    t <- 20 - h
    
    # store subj data with gammas
    iters <- list()
    
    # iterate over the alphas provided
    for (param in seq_along(alphas)) {
      
      # set up param
      alpha <- alphas[param]
      
      # get the trial-wise gamma to export the probability of acceptance
      c <- rep(0, nrow(sub))
      gamma <- rep(0, nrow(sub))
      
      # calculate gammas
      for (j in seq(nrow(sub))) {
        if (j == 1) {
          gamma[j] <- gammaStart
          c[j] <- ifelse(o[j] / (h[j] ^ k) > gamma[j], 1, 0)
        } else {
          tau <- (h[j - 1] ^ k * c[j - 1]) + t[j - 1]
          a <- o[j - 1] * c[j - 1]
          gamma[j] <- ((1 - (1 - alpha) ^ tau) * (a / tau)) + ((1 - alpha) ^ tau) * gamma[j - 1]
          c[j] <- ifelse(o[j] / (h[j] ^ k) > gamma[j], 1, 0)
        } 
      }
      
      #  store subject data under a specific alpha
      sub$g <- gamma
      sub$alpha <- alpha
      sub$c <- c
      iters[[param]] <- sub
    }
    
    # combine data frames with gammas per alpha
    temp <- do.call(rbind, iters)
    
    # create a list with possible starting values for model parameters
    # parameter names must match model ones
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
    
    # fit to subject
    baseOC <- optimize_model(sub, params, model_expr, simplify = T)
    
    # plot
    ratePlot <- temp %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / (Handling ^ k),
      fitChoice = ifelse(trialRate > g, -0.25, -5),
      newChoice = ifelse(Choice == 1, -0.5, -5),
      trialRate = ifelse(trialRate > 3, 3, trialRate),
      g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g, group = alpha)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      #geom_hline(yintercept = baseOC$gamma) + # single gamma estimated for an individual
      geom_line(aes(TrialN, optimalRate), alpha = 0.8, color = "goldenrod") + # mean optimal rate across blocks
      geom_line(aes(color = alpha), size = 0.5) +
      geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #annotate("text", x = max(sub$TrialN) + 8, y = baseOC$gamma + 0.25, label = "Fitted \n Gamma", size = 5) +
      scale_fill_discrete(name = "Offer") +
      #scale_color_continuous(breaks = c(2, 10, 14), labels = c(2, 10, 14)) +
      ylim(0, NA) +
      labs(x = "Trial Number", y = "Earning rate") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
    
    
  } else if (exp == "wth") {
    # choose sample subject + params
    sub <- filter(dataWth, SubjID == "109")
    
    # relevant behavior elements
    o <- sub$Offer
    h <- sub$Handling
    t <- 20 - h
    
    # store subj data with gammas
    iters <- list()
    
    # iterate over the alphas provided
    for (param in seq_along(alphas)) {
      
      # set up param
      alpha <- alphas[param]
      
      # get the trial-wise gamma to export the probability of acceptance
      c <- rep(0, nrow(sub))
      gamma <- rep(0, nrow(sub))
      
      # calculate gammas
      for (j in seq(nrow(sub))) {
        if (j == 1) {
          gamma[j] <- gammaStart
          c[j] <- ifelse(o[j] / h[j] ^ k> gamma[j], 1, 0)
        } else {
          tau <- (h[j - 1] ^ k * c[j - 1]) + t[j - 1]
          a <- o[j - 1] * c[j - 1]
          gamma[j] <- ((1 - (1 - alpha) ^ tau) * (a / tau)) + ((1 - alpha) ^ tau) * gamma[j - 1]
          c[j] <- ifelse(o[j] / h[j] ^ k > gamma[j], 1, 0)
        } 
      }
      
      #  store subject data under a specific alpha
      sub$g <- gamma
      sub$alpha <- alpha
      sub$c <- c
      iters[[param]] <- sub
    }
    
    # combine data frames with gammas per alpha
    temp <- do.call(rbind, iters)
    
    # plot
    ratePlot <- temp %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / Handling ^ k,
      fitChoice = ifelse(trialRate > g, -0.25, -5),
      newChoice = ifelse(Choice == 1, -0.5, -5),
      trialRate = ifelse(trialRate > 3, 3, trialRate),
      g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g, group = alpha)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      geom_line(aes(TrialN, optimalRate), alpha = 0.8, color = "goldenrod") + # mean optimal rate across blocks
      geom_line(aes(color = alpha), size = 0.5) +
      scale_fill_discrete(name = "Offer") +
      #scale_color_manual(values = colsWth) +
      ylim(0, NA) +
      labs(x = "Trial Number", y = "Earning rate") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
  }
  
  suppressWarnings(print(ratePlot))
  
  acc <- temp %>% 
    group_by(alpha) %>%
    summarise(accuracy = mean(c == Choice))
  
  return(acc)
}

# preferred function, which works for all types of models considered so far
# like, k = 1 and alpha = 0 returns the single parameter model
# generative version: model calculates choices per fitted previous choices, not just based on observed sub choices
optimize_model_dyn_us3 <- function(subjData, params, simplify = F, gammaStart = 0) {
  # get every combination of parameters
  params <- expand.grid(params)
  
  # relevant behavior elements
  o <- subjData$Offer
  h <- subjData$Handling
  t <- 20 - h
  obs_c <- subjData$Choice
  
  # Prep list of results to be returned, and keep track per iteration
  out <- list()
  out$percentQuit <- mean(obs_c == 0) * 100
  out$percentAccept <- mean(obs_c == 1) * 100 
  out$params <- as.numeric() # winning parameters
  out$rate <- as.numeric() # to store the winning gammas
  out$pAccept <- as.numeric() # probability of acceptance for winning parameters
  out$fit_c <- as.numeric() # to store the winning fitted choices
  out$LLs <- rep(0, nrow(params))
  out$LL <- Inf
  
  # iterate through possible parameters and get the LL
  #LLs <- sapply(seq(nrow(params)), function(i) {
  for (param in seq(nrow(params))) {
    # isolate the parameters for this iteration
    # and then store them as variables
    tempr <- params[param, "tempr"]
    alpha <- params[param, "alpha"]
    if ("k" %in% colnames(params)) {
      # use a vector so I can use a common indexing for fixed and to-be-fitted values of k
      k <- rep(params[param, "k"], nrow(subjData)) # if there is a free k parameter, use its fit, otherwise use the values from exp1
    } else {
      k <- subjData$mK
    }
    
    # get the trial-wise gamma to export the probability of acceptance
    c <- rep(0, nrow(subjData))
    gamma <- rep(0, nrow(subjData))
    
    # calculate gammas iterating over trials
    # latex versions: gamma[i]: \gamma_t = (1 - (1 - \alpha) ^ {\tau_t}) \dfrac{R_{t-1} A_{t-1}}{\tau_t} +  (1 - \alpha) ^ {\tau_t} \gamma_{t-1}
    # ugly tau: \tau_{t} = (H_{t-1} ^ {k_{t-1}}  A_{t - 1}) + T_{t - 1}
    # prettier, vectorized tau: \tau = (h ^ k  a) + t
    # gamma for pretty tau: \gamma_t = (1 - (1 - \alpha) ^ {\tau_{t-1}}) \dfrac{r_{t-1} a_{t-1}}{\tau_{t-1}} +  (1 - \alpha) ^ {\tau_{t-1}} \gamma_{t-1}
    # k_t = \dfrac {1} {N}\sum_{cost}^{t - 1} k_{cost}
    # P(A): P(A)_t = \dfrac{1}{1 + exp^{-(\beta[r_t - \gamma_th_t^{k_t}])}}
    for (i in seq(nrow(subjData))) {
      if (i == 1) {
        # the environmental rate that participants start with
        # if we are fitting a single gamma, then add gamma prior with alpha = 0 and k = 1 (check that it reproduces results)
        # otherwise this functions as a prior bias on the environmental rate
        if ("gammaPrior" %in% colnames(params)) {
          gamma[i] <- params[param, "gammaPrior"]
        } else {
          gamma[i] <- gammaStart
        }
        
        # choose if the trial rate > env. rate
        c[i] <- ifelse(o[i] / h[i] > gamma[i], 1, 0)
      } else {
        # non-linear estimate of the elapsed time since the last choice
        tau <- (h[i - 1] ^ k[i - 1] * c[i - 1]) + t[i - 1]
        
        # was the previous offer accepted?
        a <- o[i - 1] * c[i - 1]
        
        # gamma is updated by how much weight is given to the recently experienced reward rate (i.e. left part of eq)
        gamma[i] <- ((1 - (1 - alpha) ^ tau) * (a / tau)) + ((1 - alpha) ^ tau) * gamma[i - 1]
        
        # choose if the prospect's reward rate, non-linearly discounted as above > env. rate
        # in other words, is the local-focus on handling time being affected, or a global environmental rate? (or something in between?)
        c[i] <- ifelse(o[i] / (h[i] ^ k[i]) > gamma[i], 1, 0)
      } 
    }
    
    # estimate the probability of acceptance based on the difference between the offer rate vs global rate
    # this is a rehash of the eq updating c[i] above
    p = 1 / (1 + exp(-(tempr * (o - (gamma * h ^ k)))))
    p[p == 1] <- 0.999
    p[p == 0] <- 0.001
    
    # get the likelihood of the observations based on the model
    tempChoice <- rep(NA, length(obs_c))
    tempChoice[obs_c == 1] <- log(p[obs_c == 1])
    tempChoice[obs_c == 0] <- log(1 - p[obs_c == 0]) # log of probability of choice 1 when choice 0 occurred
    negLL <- -sum(tempChoice)
    
    out$LLs[param] <- negLL
    
    # if these parameters improve the fit, store the rate and choices
    if (negLL < out$LL) {
      diff <- out$LL - negLL
      out$LL <- negLL
      out$rate <- gamma
      out$fit_c <- c
      out$params <- params[param, ]
      out$pAccept <- p
      
      # break if the reduction is too small to be significant
      # based on visual assessment of the likelihood space, which looks convex
      # if (abs(diff) <= 0.0001) {
      #   break
      # }
    }
  } 
  
  # Summarize the outputs
  out$LL0 <- -(log(0.5) * length(obs_c))
  out$Rsquared <- 1 - (out$LL / out$LL0) # pseudo r-squared, quantifying the proportion of deviance reduction vs chance
  
  # if doing this with dplyr::do(), return a simplified data.frame instead with the important parameters
  if (simplify) {
    out <- round(data.frame(out[-c(4, 5, 6, 7)]), digits = 2)
    colnames(out) <- c("percentQuit",
                       "percentAccept",
                       colnames(params),
                       "LL",
                       "LL0",
                       "Rsq")
  }
  
  return(out)
}
optimize_model_dyn_us3 <- function(subjData, params, simplify = F, gammaStart = 0.5) {
  # get every combination of parameters
  params <- expand.grid(params)
  
  # relevant behavior elements
  o <- subjData$Offer
  h <- subjData$Handling
  t <- 20 - h
  obs_c <- subjData$Choice
  
  # Prep list of results to be returned, and keep track per iteration
  out <- list()
  out$percentQuit <- mean(obs_c == 0) * 100
  out$percentAccept <- mean(obs_c == 1) * 100 
  out$params <- as.numeric() # winning parameters
  out$rate <- as.numeric() # to store the winning gammas
  out$pAccept <- as.numeric() # probability of acceptance for winning parameters
  out$fit_c <- as.numeric() # to store the winning fitted choices
  out$LLs <- rep(0, nrow(params))
  out$LL <- Inf
  
  # iterate through possible parameters and get the LL
  #LLs <- sapply(seq(nrow(params)), function(i) {
  for (param in seq(nrow(params))) {
    # isolate the parameters for this iteration
    # and then store them as variables
    # vectors are initialized by a single value, but they get updated
    tempr <- params[param, "tempr"]
    alpha <- params[param, "alpha"]
    ifelse("alpha_k" %in% colnames(params), alpha_k <- params[param, "alpha_k"], alpha_k <- 0.15) # evolving time estimate learning. turn into free parameter
    ifelse("k" %in% colnames(params), k <- rep(params[param, "k"], nrow(subjData)), k <- subjData$mK)
    ifelse("gammaPrior" %in% colnames(params), gamma <- rep(params[param, "gammaPrior"], nrow(subjData)), gamma <- rep(gammaStart, nrow(subjData)))
    
    # vectors to store evolving variables
    a <- rep(0, nrow(subjData))
    tau <- rep(0, nrow(subjData))
    K <- k # experienced k, mainly for updating rule
    
    # calculate gammas iterating over trials
    # latex versions: gamma[i]: \gamma_t = (1 - (1 - \alpha) ^ {\tau_t}) \dfrac{R_{t-1} A_{t-1}}{\tau_t} +  (1 - \alpha) ^ {\tau_t} \gamma_{t-1}
    # ugly tau: \tau_{t} = (H_{t-1} ^ {k_{t-1}}  A_{t - 1}) + T_{t - 1}
    # prettier, vectorized tau: \tau = (h ^ k  a) + t
    # gamma for pretty tau: \gamma_t = (1 - (1 - \alpha) ^ {\tau_{t-1}}) \dfrac{r_{t-1} a_{t-1}}{\tau_{t-1}} +  (1 - \alpha) ^ {\tau_{t-1}} \gamma_{t-1}
    # k_t = \dfrac {1} {N}\sum_{cost}^{t - 1} k_{cost}
    # P(A): P(A)_t = \dfrac{1}{1 + exp^{-(\beta[r_t - \gamma_th_t^{k_t}])}}
    i <- 1
    while (i < nrow(subjData)) {
      # choose if the prospect's reward rate, non-linearly discounted as above > env. rate
      # in other words, is the local-focus on handling time being affected, or a global environmental rate? (or something in between?)
      a[i] <- ifelse(o[i] / (h[i] ^ k[i]) > gamma[i], 1, 0)
      
      # amount earned
      ao <- o[i] * a[i]
      
      # k update versions
      # k[i + 1] <- theta[i] * k[i] * c[i] + (1 - theta[i]) * k[i - 1]
      # K[l + 1] <- theta[l] * ks[l] + (1 - theta[l]) * K[l - 1]
      # K[l + 1] <- (theta[l] ^ (-1 / l)) * ks[l] + (1 - (theta[l] ^ (-1 / l))) * K[l]
      # K[l + 1] <- K[l] + 1/l * (ks[l] - K[l]) # Sutton & Barto, page 37. Add choice to ks[l]?
      left <- ((1 - alpha_k) ^ i) * k[1]
      right <- alpha_k * (1 - alpha_k) ^ (i - seq(i)) * (K[seq(i)] ^ a[seq(i)])
      k[i + 1] <-  left + sum(right)
      
      # non-linear estimate of the elapsed time since the last choice
      tau[i] <- (h[i] ^ k[i] * a[i]) + t[i]
      
      # gamma is updated by how much weight is given to the recently experienced reward rate (i.e. left part of eq)
      gamma[i + 1] <- ((1 - (1 - alpha) ^ tau[i]) * (ao / tau[i])) + ((1 - alpha) ^ tau[i]) * gamma[i]
      
      i <- i + 1
    }
    
    # estimate the probability of acceptance based on the difference between the offer rate vs global rate
    # this is a rehash of the eq updating c[i] above
    p = 1 / (1 + exp(-(tempr * (o - (gamma * h ^ k)))))
    p[p == 1] <- 0.999
    p[p == 0] <- 0.001
    
    # get the likelihood of the observations based on the model
    tempChoice <- rep(NA, length(obs_c))
    tempChoice[obs_c == 1] <- log(p[obs_c == 1])
    tempChoice[obs_c == 0] <- log(1 - p[obs_c == 0]) # log of probability of choice 1 when choice 0 occurred
    negLL <- -sum(tempChoice)
    
    out$LLs[param] <- negLL
    
    # if these parameters improve the fit, store the rate and choices
    if (negLL < out$LL) {
      diff <- out$LL - negLL
      out$LL <- negLL
      out$rate <- gamma
      out$fit_c <- c
      out$params <- params[param, ]
      out$pAccept <- p
      
      #plot(k, type = "b")
      
      # break if the reduction is too small to be significant
      # based on visual assessment of the likelihood space, which looks convex
      # if (abs(diff) <= 0.0001) {
      #   break
      # }
    }
  } 
  
  # Summarize the outputs
  out$LL0 <- -(log(0.5) * length(obs_c))
  out$Rsquared <- 1 - (out$LL / out$LL0) # pseudo r-squared, quantifying the proportion of deviance reduction vs chance
  
  # if doing this with dplyr::do(), return a simplified data.frame instead with the important parameters
  if (simplify) {
    out <- round(data.frame(out[-c(4, 5, 6, 7)]), digits = 2)
    colnames(out) <- c("percentQuit",
                       "percentAccept",
                       colnames(params),
                       "LL",
                       "LL0",
                       "Rsq")
  }
  
  return(out)
}  # these produce diff results for basic model. see why
plot_dyn_us3 <- function(id = "58", exp = "btw", gammaOne = 0.5, showChoices = -0.6) {
  
  if (exp == "btw") {
    # choose subject + params
    sub <- filter(dataBtw, SubjID == id)
    spaceSize <- 30
    params <- list(tempr = seq(0, 2, length.out = spaceSize), 
                   alpha = seq(0, 0.5, length.out = spaceSize),
                   k = seq(0, 2, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn_us3(sub, params, simplify = F, gammaStart = gammaOne)
    print(temp$params)
    
    # get baseline gamma from a simple fit
    # model to be fit
    # make sure that you specify the inverse temperature
    # extra parameters as dfs for now, that's why the `[[1]]`
    model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
    
    # create a list with possible starting values for model parameters
    # parameter names must match model ones
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    
    # fit to subject
    baseOC <- optimize_model(sub, params, model_expr, simplify = T)
    
    k <- as.numeric(temp$params["k"])
    
    # plot
    ratePlot <- sub %>%
      mutate(trialRate = Offer / (Handling ^ k),
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$probAccept,
             g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      #geom_hline(yintercept = baseOC$gamma) + # single gamma estimated for an individual
      #geom_hline(yintercept = 0.74, linetype = "dashed", color = "grey30") + # mean optimal rate across blocks
      geom_line(aes(color = Handling), size = 0.5) +
      geom_point(aes(color = Handling), size = 1.2) +
      geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #geom_point(aes(TrialN, pChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #annotate("text", x = max(sub$TrialN) + 8, y = baseOC$gamma + 0.25, label = "Fitted \n Gamma", size = 5) +
      #annotate("text", x = max(sub$TrialN) + 8, y = 0.55, label = "Optimal", size = 5, color = "grey30") +
      annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
      annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
      scale_fill_discrete(name = "Offer") +
      scale_color_continuous(breaks = c(2, 10, 14), labels = c(2, 10, 14)) +
      ylim(showChoices, NA) +
      labs(x = "Trial Number", y = "Earning rate") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
    
  } else if (exp == "wth") {
    # choose subject + params
    sub <- filter(dataWth, SubjID == id)
    spaceSize <- 30
    params <- list(tempr = seq(0, 2, length.out = spaceSize), 
                   alpha = seq(0, 0.5, length.out = spaceSize),
                   alpha_k = seq(0, 0.5, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn_us3(sub, params, simplify = F, gammaStart = gammaOne)
    print(temp$params)
    
    # plot
    ratePlot <- sub %>%
      mutate(trialRate = Offer / Handling ^ mK,
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$probAccept,
             g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      #geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey30") + # mean optimal rate across blocks
      geom_line(size = 0.5, color = "grey50") +
      geom_point(aes(color = Cost), size = 1.2) +
      geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #annotate("text", x = max(sub$TrialN) + 8, y = 0.63, label = "Optimal", size = 5, color = "grey30") +
      annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
      annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
      scale_fill_discrete(name = "Offer") +
      scale_color_manual(values = colsWth) +
      ylim(showChoices, NA) +
      labs(x = "Trial Number", y = "Earning rate") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
  }
  
  suppressWarnings(print(ratePlot))
}

# reproduce main experimental plot findings using fitted model values
# one for each experiment, since the plots are different
recover_results_btw <- function(fitsList, binary = F) {
  # extract fits
  fits <- do.call(base::c, sapply(fitsList, "[", "pAccept"))
  if (binary) {
    # to get stochastic-less choices
    fits <- ifelse(fits > 0.5, 1, 0)
  }
  data <- dataBtw %>%
    filter(Cost != "Easy") %>%
    mutate(fitChoice = fits) #rbinom(length(fits), 1, fits))
  
  # plot
  data %>%
    group_by(SubjID, Cost, Handling, Offer, optimal) %>%
    summarise(pAccept = mean(fitChoice)) %>%
    group_by(Cost, Handling, Offer, optimal) %>%
    summarise(meanComplete = mean(pAccept),
              SE = sd(pAccept) / sqrt(length(pAccept))) %>%
    ggplot(aes(interaction(Offer, Handling), meanComplete, color = Cost)) + 
    geom_point(size = 3) + 
    geom_errorbar(aes(ymin = meanComplete - SE, ymax = meanComplete + SE), width = 0.2, size = 1) +
    geom_line(aes(group = interaction(Handling, Cost)), size = 1) +
    geom_point(aes(y = round(optimal)), shape = 21, fill = "grey75", color = "grey30", size = 3, show.legend = F) +
    labs(x = "Offer.Handling", y = "Proportion Accepted") +
    scale_color_manual(values = colsBtw) +
    theme(legend.key = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size = 12))
  
}
recover_results_wth <- function(fitsList, binary = F, matrix = T) {
  # use model fits to reproduce observed behavioral plots
  # the reported pre-post results, including mixed effects, use half 1 and half 2 raw, but the results are pretty much conserved either way.
  # here the middle blocks are eliminated to denote prevent the uneven assignment of the middle blocks to bias stuff either way
  
  # extract fits
  fits <- do.call(c, sapply(fitsList, "[", "pAccept"))
  if (binary) {
    # to get stochastic-less choices
    fits <- ifelse(fits > 0.5, 1, 0)
  }
  data <- dataWth %>%
    mutate(fitChoice = fits) #rbinom(length(fits), 1, fits))
  
  # plot prop accepted
  plot <- data %>%
    group_by(SubjID) %>%
    mutate(Block = case_when(
      Block %in% c(1, 2) ~ 1,
      Block %in% c(3, 4) ~ 2,
      Block %in% c(5, 6) ~ 3
    )) %>%
    filter(Block != 2) %>%
    mutate(Block = ifelse(Block == 1, "First Two Blocks", "Last Two Blocks")) %>%
    group_by(Cost, Block, Offer) %>%
    summarise(propAccept = mean(fitChoice),
              SE = sd(fitChoice) / sqrt(nSubjs_wth)) %>%
    ungroup() %>%
    mutate(Offer = ifelse(Offer == 20, 12, Offer)) %>% 
    ggplot(aes(Offer, propAccept, color = Cost)) +
    geom_point(size = 3, show.legend = T) +
    geom_line(size = 1, show.legend = F) +
    scale_color_manual(values = colsWth) +
    scale_fill_manual(values = colsWth) +
    geom_errorbar(aes(ymin = propAccept - SE, ymax = propAccept + SE), width = 0.3, size = 1, show.legend = T) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(4, 8, 12), labels = c(4, 8, 20)) +
    labs(x = "Reward", y = "Proportion Accepted") +
    facet_wrap(vars(Block)) +
    theme(legend.key = element_blank(),
          legend.position = c(0.9, 0.25),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size = 16))
  
  print(plot)
  
  # plot the pre-post coefficient matrix, but only for the binary version (based on GLM def)
  if (binary & matrix) {
    mixLogis_pre <- list()
    mixData <- data %>%
      mutate(Choice = fitChoice,
             Cost = factor(Cost, levels = list("Cognitive", "Physical", "Wait-C", "Wait-P")),
             Block = case_when(
               Block %in% c(1, 2) ~ 1,
               Block %in% c(3, 4) ~ 2,
               Block %in% c(5, 6) ~ 3
             )) %>%
      filter(Block != 2,
             Half == "Half_1") %>%
      group_by(SubjID, Cost, Offer) %>%
      summarize(totalChoices = length(Choice),
                totalAccepted = sum(Choice),
                totalQuits = sum(Choice == 0),
                propAccepted = mean(Choice))
    
    mixLogis_pre$Cognitive <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData)
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Wait-C"))
    mixLogis_pre$`Wait-C` <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData)
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Wait-P"))
    mixLogis_pre$`Wait-P` <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData)
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Physical"))
    mixLogis_pre$Physical <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData)
    
    
    # mixed logistic for the second half
    mixLogis_post <- list()
    mixData <- data %>%
      mutate(Choice = fitChoice,
             Cost = factor(Cost, levels = list("Cognitive", "Physical", "Wait-C", "Wait-P")),
             Block = case_when(
               Block %in% c(1, 2) ~ 1,
               Block %in% c(3, 4) ~ 2,
               Block %in% c(5, 6) ~ 3
             )) %>%
      filter(Block != 2,
             Half == "Half_2") %>%
      group_by(SubjID, Cost, Offer) %>%
      summarize(totalChoices = length(Choice),
                totalAccepted = sum(Choice),
                totalQuits = sum(Choice == 0),
                propAccepted = mean(Choice))
    
    mixLogis_post$Cognitive <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData)
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Wait-C"))
    mixLogis_post$`Wait-C` <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData)
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Wait-P"))
    mixLogis_post$`Wait-P` <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData)
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Physical"))
    mixLogis_post$Physical <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData)
    
    ## now produce a summary matrix
    # get the beta and pvalue matrices
    betasPre <- betaMatrix(mixLogis_pre)
    betasPost <- betaMatrix(mixLogis_post)
    
    # and combine them
    betaMat <- matrix(NA, nrow = nrow(betasPre$Betas), ncol = nrow(betasPre$Betas))
    betaMat[lower.tri(betaMat)] <- betasPre$Betas[lower.tri(betasPre$Betas)]
    betaMat[upper.tri(betaMat)] <- betasPost$Betas[upper.tri(betasPost$Betas)]
    dimnames(betaMat) <- list(names(mixLogis_pre), names(mixLogis_pre))
    
    pvalMat <- matrix(NA, nrow = nrow(betasPre$Pvals), ncol = nrow(betasPre$Pvals))
    pvalMat[lower.tri(pvalMat)] <- betasPre$Pvals[lower.tri(betasPre$Pvals)]
    pvalMat[upper.tri(pvalMat)] <- betasPost$Pvals[upper.tri(betasPost$Pvals)]
    dimnames(pvalMat) <- list(names(mixLogis_pre), names(mixLogis_pre))
    
    # remove uninteresting comparisons
    diag(betaMat) <- 0
    betaMat[rbind(c(1, 3), c(2, 4), c(3, 1), c(4,2))] <- NA
    pvalMat[rbind(c(1, 3), c(2, 4), c(3, 1), c(4,2))] <- NA
    
    # aand plot
    col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))
    corrplot(betaMat,
             is.corr = F,
             p.mat = pvalMat,
             outline = T,
             #insig = "p-value",
             sig.level = 0.05,
             na.label = "square",
             na.label.col = "grey",
             method = "color",
             tl.col = "black",
             tl.cex = 0.8,
             col = col2(200))
  }
  
}

# function for idea: people adapt their Ks with decreasing strength
# meaning that early in the experiment, they keep track of the perceived time per cost, but they slowly settle into previous history towards the end
# similar result to the cumulative mean, but more flexible early in the experiment
# theta vector fixed for now. Think of how to make it free
recalibrate_k <- function(ks, thetaRange = c(1, 0)) {
  # rule
  theta <- seq(thetaRange[1], thetaRange[2], length.out = length(ks))
  
  learnedK <- rep(0, length(ks))
  for (l in seq_along(ks)) {
    if (l == 1) {
      learnedK[l] <- ks[1]
    } else {
      learnedK[l] <- theta[l] * ks[l] + (1 - theta[l]) * learnedK[l - 1]
    }
  }
  
  return(learnedK)
}
recalibrate_k <- function(ks, thetaRange, choice) {
  # rule
  theta <- seq(thetaRange[1], thetaRange[2], length.out = length(ks))
  
  # vector to store the evolving estimates of K 
  K <- ks
  K[1] <- 1
  
  a <- theta[1]
  l <- 1
  while (l < length(ks)) {
    # K[l + 1] <- theta[l] * ks[l] * choice[l] + (1 - theta[l]) * K[l]
    # K[l + 1] <- K[l] + (1/l * (ks[l] - K[l])) * choice[l] # Sutton & Barto, page 37. Add choice to ks[l]?
    left <- ((1 - a) ^ l) * K[1]
    right <- a * (1 - a) ^ (l - seq(l)) * ks[seq(l)] ^ choice[seq(l)]
    K[l + 1] <-  left + sum(right)
    
    l <- l + 1
  }
  
  return(K)
}




## which models to run?
baseOC_nloptr <- F
bOC <- T
baseLogistic <- F # to test whether the brute search converges to a conventional logistic through glm()
fwOC <- F
dOC <- F
twOC <- F
recovery <- T

# which experimental dataset?
data <- dataBtw %>% 
  filter(Cost != "Easy")


## ORIGINAL OC NLOPTR
if (baseOC_nloptr) {
  print("Running base OC model on NLOPTR...")
  
  summaryOC <- data %>%
    group_by(Cost, SubjID) %>%
    do(optimizeOCModel(., simplify = T)) %>%
    ungroup()
}

## ORIGINAL OC
# great correspondence with NLOPTR, just much slower since it surveys the whole parameter space
if (bOC) {
  print("Running base OC model through grid search...")
  
  # # model to be fit
  # # make sure that you specify the inverse temperature
  # # extra parameters as dfs for now, that's why the `[[1]]`
  # model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
  # 
  # # create a list with possible starting values for model parameters
  # # parameter names must match model ones
  # spaceSize <- 30
  # params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
  #                gamma = seq(0.25, 1.5, length.out = spaceSize))
  # 
  # # fit to each subject
  # baseOC <- data %>%
  #   group_by(Cost, SubjID) %>%
  #   do(optimize_model(., params, model_expr, simplify = T)) %>%
  #   ungroup()
  
  # alternative using the big model function
  spaceSize <- 30
  params <- list(tempr = seq(0, 2, length.out = spaceSize), 
                 gammaPrior = seq(0.25, 1.5, length.out = spaceSize),
                 alpha = 0,
                 k = 1,
                 alpha_k = 0)
  
  # fit to each subject
  baseOC <- dataBtw %>%
    filter(Cost != "Easy") %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model_dyn_us3(., params, simplify = T)) %>%
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


## trial-wise updating of gamma
if (twOC) {
  print("Running trial-wise OC model through grid search...")
  
  ## vanilla C&W (slightly adapted for prey selection)
  # create a list with possible starting values for model parameters
  spaceSize <- 30
  params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                 alpha = seq(0, 1, length.out = spaceSize))
  
  # between subject exp
  trialwiseOC_btw <- dataBtw %>%
    filter(Cost != "Easy") %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model_dyn(., params, simplify = T)) %>%
    ungroup() %>%
    distinct(SubjID, .keep_all = T)
  
  
  # within subject exp
  trialwiseOC_wth <- dataWth %>% 
    group_by(SubjID) %>%
    do(optimize_model_dyn(., params, simplify = T)) %>%
    ungroup() %>%
    distinct(SubjID, .keep_all = T)
  
  
  ## modified with non linear time
  # create a list with possible starting values for model parameters
  params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                 alpha = seq(0.001, 1, length.out = spaceSize),
                 k = seq(0, 2, length.out = spaceSize))
  
  # between subject exp
  trialwiseOC_btw_us <- dataBtw %>%
    filter(Cost != "Easy") %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model_dyn_us(., params, simplify = T)) %>%
    ungroup() %>%
    distinct(SubjID, .keep_all = T)
  
  
  # between subject exp
  trialwiseOC_wth_us <- dataWth %>% 
    group_by(SubjID) %>%
    do(optimize_model_dyn_us(., params, simplify = T)) %>%
    ungroup() %>%
    distinct(SubjID, .keep_all = T)
  
}


### test
# taking the cumulative average of the subjective time perception makes subjective overall time seem slower.
# a free parameter could be added to index the degree to which participants cared about this (replaces k as a free parameter, since k is adopted from exp 1)
# weird that the higher the weight on cum_means, the lower the OC ie less selective. Think about it.
# fit btw exp
# alpha will often default to the next lowest to 0
# so I tested a number of iterations to ensure that the results persist
# at 30 and 50 it's the same. As alpha is reduced to ~0, k increases for cog.
# but ~0.02 alpha seems sensible.
spaceSize <- 30
params <- list(tempr = seq(0, 2, length.out = spaceSize), 
               alpha = seq(0, 0.5, length.out = spaceSize),
               k = seq(0, 2, length.out = spaceSize),
               alpha_k = seq(0, 0.5, length.out = spaceSize))

trialwiseOC_btw_us_new <- dataBtw %>%
  filter(Cost != "Easy") %>%
  group_by(Cost, SubjID) %>%
  do(optimize_model_dyn_us3(., params, simplify = T)) %>%
  ungroup() %>%
  distinct(SubjID, .keep_all = T)

param_compare_plot(trialwiseOC_btw_us_new, "k", meanRate = 1)

#btw ks to apply to wth
ks <- trialwiseOC_btw_us_new %>% 
  group_by(Cost) %>% 
  summarise(mK = median(k)) %>%
  rename(simpleCost = Cost) # to merge without replacing

# reset data: 
tryCatch(dataWth <- dataWth %>% select(-mK), error = function(e) {print("Oops, no need to remove mK")})
if (! "mK" %in% colnames(dataWth)) {
  dataWth <- dataWth %>%
    mutate(simpleCost = ifelse(Cost %in% c("Wait-C", "Wait-P"), "Wait", as.character(Cost))) %>% 
    left_join(ks, by = "simpleCost") #%>%
    # group_by(SubjID) %>%
    # mutate(laggedK = dplyr::lag(mK, default = 1),
    #        mK = mK * Choice,
    #        mK = ifelse(mK == 0, laggedK, mK),
    #        mK = recalibrate_k(mK, thetaRange = c(0.9, 0))) %>%
    # ungroup()
}

# # check that the k evolution looks ok
# dataWth %>%
#   filter(SubjID %in% c("109", "461")) %>%
#   ggplot(aes(TrialN, mK, color = SubjID, group = SubjID)) +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   geom_line(show.legend = F) +
#   theme_minimal()

# fit wth
params <- list(tempr = seq(0, 2, length.out = spaceSize), 
               alpha = seq(0, 0.5, length.out = spaceSize),
               alpha_k = seq(0, 0.5, length.out = spaceSize))

trialwiseOC_wth_us_new <- dataWth %>%
  group_by(SubjID) %>%
  do(optimize_model_dyn_us3(., params, simplify = T)) %>%
  ungroup() %>%
  distinct(SubjID, .keep_all = T)



### Result recovery
if (recovery) {
  ### Between subjects
  # attempt to recover the observed results from 16.2.2, reproducing the prop-completed plot
  # so far both single gamma and alpha + k perform well, though the latter is preferred due to explanatory power and generalization to exp 2
  
  ## recover results using basic model
  #model_expr <- expr(tempr[[1]] * (reward - (gammaPrior[[1]] * handling)))
  
  # create a list with possible starting values for model parameters
  # # parameter names must match model ones
  # spaceSize <- 30
  # params <- list(tempr = seq(0, 2, length.out = spaceSize), 
  #                gammaPrior = seq(0.25, 1.5, length.out = spaceSize),
  #                alpha = 0,
  #                k = 1)
  # 
  # # fit to each btw subject
  # temp_bOC_fits <- dataBtw %>%
  #   filter(Cost != "Easy") %>%
  #   plyr::dlply("SubjID", identity) %>%
  #   lapply(., optimize_model_dyn_us3, params, simplify = F)
  # 
  # recover_results_btw(temp_bOC_fits, binary = T)
  # 
  # ## alpha only
  # # create a list with possible starting values for model parameters
  # spaceSize <- 30
  # params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
  #                alpha = seq(0, 0.5, length.out = spaceSize),
  #                k = 1)
  # 
  # # between subject exp
  # temp_dOC_fits <- dataBtw %>%
  #   filter(Cost != "Easy") %>%
  #   plyr::dlply("SubjID", identity) %>%
  #   lapply(., optimize_model_dyn_us3, params, simplify = F)
  
  # fit btw exp final model (alpha + nonlinear k)
  params <- list(tempr = seq(0, 2, length.out = spaceSize),
                 alpha = seq(0, 0.5, length.out = spaceSize),
                 k = seq(0, 2, length.out = spaceSize),
                 alpha_k = seq(0, 0.5, length.out = spaceSize))
  
  temp_tw_fits_btw <- dataBtw %>%
    filter(Cost != "Easy") %>%
    plyr::dlply("SubjID", identity) %>%
    lapply(., optimize_model_dyn_us3, params, simplify = F, gammaStart = 0)
  
  
  recover_results_btw(temp_tw_fits_btw, binary = F)
  
  ### Within subjects (reproduce 16.3.3)
  # partial replication! without temperature I get the first block fine
  # with temperature I get the last one perfect (i.e. based on probabilities). First block still matches kind of.
  # the reported pre-post results, including mixed effects, use half 1 and half 2 raw, but the results are pretty much conserved either way.
  # here the middle blocks are eliminated to denote prevent the uneven assignment of the middle blocks to bias stuff either way
  
  # fit wth
  spaceSize <- 30
  params <- list(tempr = seq(0, 2, length.out = spaceSize),
                 alpha = seq(0, 0.5, length.out = spaceSize),
                 alpha_k = seq(0, 0.5, length.out = spaceSize))
  
  temp_tw_fits_wth <- dataWth %>%
    plyr::dlply("SubjID", identity) %>%
    lapply(., optimize_model_dyn_us3, params, simplify = F)
  
  recover_results_wth(temp_tw_fits_wth, binary = T)
}







