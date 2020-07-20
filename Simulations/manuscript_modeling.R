
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
plot_dyn <- function(id = "58", exp = "btw", gammaOne = 0) {
  
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
        ylim(-0.6, NA) +
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
        ylim(-0.6, NA) +
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
      ylim(-0.6, NA) +
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
      ylim(-0.6, NA) +
      labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
  }
  
  suppressWarnings(print(ratePlot))
}

# plot gamma for stipulated alpha and k parameters
plot_alphas <- function(alpha, k = 1, exp = "btw", gammaStart = 0) {
  # parameter setting
  
  
  if (exp == "btw") {
    # choose sample subject + params
    sub <- filter(dataBtw, SubjID == "58")
    
    # relevant behavior elements
    o <- sub$Offer
    h <- sub$Handling
    t <- 20 - h
    c <- sub$Choice
    a <- o * c # accepted offers
    time <- sub$ExpTime
    tau <- lag((h ^ k * c) + t)
    
    # get the trial-wise gamma to export the probability of acceptance
    gamma <- rep(0, nrow(sub))
    
    # calculate gammas
    for (j in seq(nrow(sub))) {
      if (j == 1) {
        gamma[j] <- gammaStart
      } else {
        gamma[j] <- ((1 - (1 - alpha) ^ tau[j]) * (a[j - 1] / tau[j])) + ((1 - alpha) ^ tau[j]) * gamma[j - 1]
      } 
    }
    
    # create a list with possible starting values for model parameters
    # parameter names must match model ones
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    
    # fit to subject
    baseOC <- optimize_model(sub, params, model_expr, simplify = T)
    
    # plot
    ratePlot <- sub %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / Handling,
      g = gamma,
      fitChoice = ifelse(trialRate > g, -0.25, -5),
      newChoice = ifelse(Choice == 1, -0.5, -5),
      trialRate = ifelse(trialRate > 3, 3, trialRate),
      g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      geom_hline(yintercept = baseOC$gamma) + # single gamma estimated for an individual
      geom_line(aes(TrialN, optimalRate), alpha = 0.8, color = "goldenrod") + # mean optimal rate across blocks
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
      ylim(-0.6, NA) +
      labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
    
  } else if (exp == "wth") {
    # choose subject + params
    sub <- filter(dataWth, SubjID == "109")
    
    # relevant behavior elements
    o <- sub$Offer
    h <- sub$Handling
    t <- 20 - h
    c <- sub$Choice
    a <- o * c # accepted offers
    time <- sub$ExpTime
    tau <- lag((h ^ k * c) + t)
    
    # get the trial-wise gamma to export the probability of acceptance
    gamma <- rep(0, nrow(sub))
    
    # calculate gammas
    for (j in seq(nrow(sub))) {
      if (j == 1) {
        gamma[j] <- gammaStart
      } else {
        gamma[j] <- ((1 - (1 - alpha) ^ tau[j]) * (a[j - 1] / tau[j])) + ((1 - alpha) ^ tau[j]) * gamma[j - 1]
      } 
    }
    
    # plot
    ratePlot <- sub %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / Handling,
      g = gamma,
      fitChoice = ifelse(trialRate > g, -0.25, -5),
      newChoice = ifelse(Choice == 1, -0.5, -5),
      trialRate = ifelse(trialRate > 3, 3, trialRate),
      g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      geom_line(aes(TrialN, optimalRate), alpha = 0.8, color = "goldenrod") + # mean optimal rate across blocks
      geom_line(size = 0.5, color = "grey50") +
      geom_point(aes(color = Cost), size = 1.2) +
      geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      annotate("text", x = max(sub$TrialN) + 8, y = 0.63, label = "Optimal", size = 5, color = "grey30") +
      annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
      annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
      scale_fill_discrete(name = "Offer") +
      scale_color_manual(values = colsWth) +
      ylim(-0.6, NA) +
      labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
  }
  
  suppressWarnings(print(ratePlot))
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
                 alpha = seq(0, 1, length.out = spaceSize),
                 k = seq(-2, 2, length.out = spaceSize))
  
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

#plot
param_compare_plot(trialwiseOC_btw_us, "k", meanRate = 1) 


### test
# taking the cumulative average of the subjective time perception makes subjective overall time seem slower.
# a free parameter could be added to index the degree to which participants cared about this (replaces k as a free parameter, since k is adopted from exp 1)
# weird that the higher the weight on cum_means, the lower the OC ie less selective. Think about it.
optimize_model_dyn_us2 <- function(subjData, params, simplify = F, gammaStart = 0) {
  # get every combination of parameters
  params <- expand.grid(params)
  
  # relevant behavior elements
  o <- subjData$Offer
  h <- subjData$Handling
  t <- 20 - h
  c <- subjData$Choice
  a <- o * c # accepted offers
  k <- subjData$mK
  time <- subjData$ExpTime
  
  # Prep list of results to be returned
  out <- list()
  out$percentQuit <- mean(c == 0) * 100
  out$percentAccept <- mean(c == 1) * 100 
  
  # iterate through possible parameters and get the LL
  LLs <- sapply(seq(nrow(params)), function(i) {
    # isolate the parameters for this iteration
    # and then store them as variables
    tempr <- params[i, "tempr"]
    alpha <- params[i, "alpha"]
    if ("k" %in% colnames(params)) {
      k <- params[i, "k"] # if there is a free k parameter, use its fit, otherwise use the values from exp1
    }
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
  if ("k" %in% colnames(params)) {
    k <- k[[1]] 
  } 
  tau <- lag((h ^ k * c) + t)
  
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
  out$probAccept <- 1 / (1 + exp(-(tempr[[1]] * (o - (gamma * h ^ k)))))
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

plot_dyn_us2 <- function(id = "58", exp = "btw", gammaOne = 0) {
  
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
      ylim(-0.6, NA) +
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
                   alpha = seq(0, 1, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn_us2(sub, params, simplify = F, gammaStart = gammaOne)
    
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
      ylim(-0.6, NA) +
      labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
  }
  
  suppressWarnings(print(ratePlot))
}

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
  out$LL <- 100000000000
  
  # iterate through possible parameters and get the LL
  #LLs <- sapply(seq(nrow(params)), function(i) {
  for (i in seq(nrow(params))) {
    # isolate the parameters for this iteration
    # and then store them as variables
    tempr <- params[i, "tempr"]
    alpha <- params[i, "alpha"]
    if ("k" %in% colnames(params)) {
      # use a vector so I can use a common indexing for fixed and to-be-fitted values of k
      k <- rep(params[i, "k"], nrow(subjData)) # if there is a free k parameter, use its fit, otherwise use the values from exp1
    } else {
      k <- subjData$mK
    }
    
    # get the trial-wise gamma to export the probability of acceptance
    c <- rep(0, nrow(sub))
    gamma <- rep(0, nrow(sub))
    
    # calculate gammas iterating over trials
    for (i in seq(nrow(sub))) {
      if (i == 1) {
        gamma[i] <- gammaStart
        c[i] <- ifelse(o[i] / h[i] > gamma[i], 1, 0)
      } else {
        tau <- (h[i - 1] ^ k[i - 1] * c[i - 1]) + t[i - 1]
        a <- o[i - 1] * c[i - 1]
        gamma[i] <- ((1 - (1 - alpha) ^ tau) * (a / tau)) + ((1 - alpha) ^ tau) * gamma[i - 1]
        c[i] <- ifelse(o[i] / (h[i] ^ k[i]) > gamma[i], 1, 0)
      } 
    }
    
    # estimate the probability of acceptance per the model
    p = 1 / (1 + exp(-(tempr * (o - (gamma * h ^ k)))))
    p[p == 1] <- 0.999
    p[p == 0] <- 0.001
    
    # get the likelihood of the observations based on the model
    tempChoice <- rep(NA, length(obs_c))
    tempChoice[obs_c == 1] <- log(p[obs_c == 1])
    tempChoice[obs_c == 0] <- log(1 - p[obs_c == 0]) # log of probability of choice 1 when choice 0 occurred
    negLL <- -sum(tempChoice)
    
    # if these parameters improve the fit, store the rate and choices
    if (negLL < out$LL) {
      #write(negLL, stdout())
      out$LL <- negLL
      out$rate <- gamma
      out$fit_c <- c
      out$params <- params[i, ]
      out$pAccept <- p
    }
  } 
  
  # Summarize the outputs
  out$LL0 <- -(log(0.5) * length(obs_c))
  out$Rsquared <- 1 - (out$LL / out$LL0) # pseudo r-squared, quantifying the proportion of deviance reduction vs chance
  
  # if doing this with dplyr::do(), return a simplified data.frame instead with the important parameters
  if (simplify) {
    out <- round(data.frame(out[-c(4, 5, 6)]), digits = 2)
    colnames(out) <- c("percentQuit",
                       "percentAccept",
                       colnames(params),
                       "LL",
                       "LL0",
                       "Rsq")
  }
  
  return(out)
}


plot_dyn_us3 <- function(id = "58", exp = "btw", gammaOne = 0) {
  
  if (exp == "btw") {
    # choose subject + params
    sub <- filter(dataBtw, SubjID == id)
    spaceSize <- 30
    params <- list(tempr = seq(0, 2, length.out = spaceSize), 
                   alpha = seq(0, 1, length.out = spaceSize),
                   k = seq(-1, 1, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn_us3(sub, params, simplify = F, gammaStart = gammaOne)
    
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
      ylim(-0.6, NA) +
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
                   alpha = seq(0, 1, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_dyn_us3(sub, params, simplify = F, gammaStart = gammaOne)
    
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
      ylim(-0.6, NA) +
      labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 16))
  }
  
  suppressWarnings(print(ratePlot))
}


#base choices off evolving gamma, not past observed participant choices
#this might be an issue
plot_alphas <- function(alpha, k = 1, exp = "btw", gammaStart = 0) {
  # parameter setting
  
  
  if (exp == "btw") {
    # choose sample subject + params
    sub <- filter(dataBtw, SubjID == "58")
    
    # relevant behavior elements
    o <- sub$Offer
    h <- sub$Handling
    t <- 20 - h
    c <- sub$Choice
    a <- o * c # accepted offers
    time <- sub$ExpTime
    tau <- lag((h ^ k * c) + t)
    
    # get the trial-wise gamma to export the probability of acceptance
    gamma <- rep(0, nrow(sub))
    
    # calculate gammas
    for (j in seq(nrow(sub))) {
      if (j == 1) {
        gamma[j] <- gammaStart
      } else {
        gamma[j] <- ((1 - (1 - alpha) ^ tau[j]) * (a[j - 1] / tau[j])) + ((1 - alpha) ^ tau[j]) * gamma[j - 1]
      } 
    }
    
    # create a list with possible starting values for model parameters
    # parameter names must match model ones
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    
    # fit to subject
    baseOC <- optimize_model(sub, params, model_expr, simplify = T)
    
    # plot
    ratePlot <- sub %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / Handling,
      g = gamma,
      fitChoice = ifelse(trialRate > g, -0.25, -5),
      newChoice = ifelse(Choice == 1, -0.5, -5),
      trialRate = ifelse(trialRate > 3, 3, trialRate),
      g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      geom_hline(yintercept = baseOC$gamma) + # single gamma estimated for an individual
      geom_line(aes(TrialN, optimalRate), alpha = 0.8, color = "goldenrod") + # mean optimal rate across blocks
      geom_line(aes(color = Handling), size = 0.5) +
      #geom_point(aes(color = Handling), size = 1.2) +
      #geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #geom_point(aes(TrialN, pChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      annotate("text", x = max(sub$TrialN) + 8, y = baseOC$gamma + 0.25, label = "Fitted \n Gamma", size = 5) +
      #annotate("text", x = max(sub$TrialN) + 8, y = 0.55, label = "Optimal", size = 5, color = "grey30") +
      #annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
      #annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
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
    sub <- filter(dataWth, SubjID == "109")
    
    # relevant behavior elements
    o <- sub$Offer
    h <- sub$Handling
    t <- 20 - h
    c <- sub$Choice
    a <- o * c # accepted offers
    time <- sub$ExpTime
    tau <- lag((h ^ k * c) + t)
    
    # get the trial-wise gamma to export the probability of acceptance
    gamma <- rep(0, nrow(sub))
    
    # calculate gammas
    for (j in seq(nrow(sub))) {
      if (j == 1) {
        gamma[j] <- gammaStart
      } else {
        gamma[j] <- ((1 - (1 - alpha) ^ tau[j]) * (a[j - 1] / tau[j])) + ((1 - alpha) ^ tau[j]) * gamma[j - 1]
      } 
    }
    
    # plot
    ratePlot <- sub %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / Handling,
      g = gamma,
      fitChoice = ifelse(trialRate > g, -0.25, -5),
      newChoice = ifelse(Choice == 1, -0.5, -5),
      trialRate = ifelse(trialRate > 3, 3, trialRate),
      g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      geom_line(aes(TrialN, optimalRate), alpha = 0.8, color = "goldenrod") + # mean optimal rate across blocks
      geom_line(size = 0.5, color = "grey50") +
      #geom_point(aes(color = Cost), size = 1.2) +
      #geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #annotate("text", x = max(sub$TrialN) + 8, y = 0.63, label = "Optimal", size = 5, color = "grey30") +
      #annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
      #annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
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
  
  fitChoice <- ifelse(o / h > gamma, 1, 0)
  return(sum(fitChoice == sub$Choice) / nrow(sub))
}

plot_alphas2 <- function(alpha, k = 1, exp = "btw", gammaStart = 0) {
  if (exp == "btw") {
    # choose sample subject + params
    sub <- filter(dataBtw, SubjID == "58")
    
    # relevant behavior elements
    o <- sub$Offer
    h <- sub$Handling
    t <- 20 - h
    
    # get the trial-wise gamma to export the probability of acceptance
    c <- rep(0, nrow(sub))
    gamma <- rep(0, nrow(sub))
    
    # calculate gammas
    for (j in seq(nrow(sub))) {
      if (j == 1) {
        gamma[j] <- gammaStart
        c[j] <- ifelse(o[j] / h[j] > gamma[j], 1, 0)
      } else {
        tau <- (h[j - 1] ^ k * c[j - 1]) + t[j - 1]
        a <- o[j - 1] * c[j - 1]
        gamma[j] <- ((1 - (1 - alpha) ^ tau) * (a / tau)) + ((1 - alpha) ^ tau) * gamma[j - 1]
        c[j] <- ifelse(o[j] / h[j] > gamma[j], 1, 0)
      } 
    }
    
    # create a list with possible starting values for model parameters
    # parameter names must match model ones
    spaceSize <- 30
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    
    # fit to subject
    baseOC <- optimize_model(sub, params, model_expr, simplify = T)
    
    # plot
    ratePlot <- sub %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / Handling,
      g = gamma,
      fitChoice = ifelse(trialRate > g, -0.25, -5),
      newChoice = ifelse(Choice == 1, -0.5, -5),
      trialRate = ifelse(trialRate > 3, 3, trialRate),
      g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      geom_hline(yintercept = baseOC$gamma) + # single gamma estimated for an individual
      geom_line(aes(TrialN, optimalRate), alpha = 0.8, color = "goldenrod") + # mean optimal rate across blocks
      geom_line(aes(color = Handling), size = 0.5) +
      #geom_point(aes(color = Handling), size = 1.2) +
      #geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #geom_point(aes(TrialN, pChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      annotate("text", x = max(sub$TrialN) + 8, y = baseOC$gamma + 0.25, label = "Fitted \n Gamma", size = 5) +
      #annotate("text", x = max(sub$TrialN) + 8, y = 0.55, label = "Optimal", size = 5, color = "grey30") +
      #annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
      #annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
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
    # choose sample subject + params
    sub <- filter(dataWth, SubjID == "109")
    
    # relevant behavior elements
    o <- sub$Offer
    h <- sub$Handling
    t <- 20 - h
    
    # get the trial-wise gamma to export the probability of acceptance
    c <- rep(0, nrow(sub))
    gamma <- rep(0, nrow(sub))
    
    # calculate gammas
    for (j in seq(nrow(sub))) {
      if (j == 1) {
        gamma[j] <- gammaStart
        c[j] <- ifelse(o[j] / h[j] > gamma[j], 1, 0)
      } else {
        tau <- (h[j - 1] ^ k * c[j - 1]) + t[j - 1]
        a <- o[j - 1] * c[j - 1]
        gamma[j] <- ((1 - (1 - alpha) ^ tau) * (a / tau)) + ((1 - alpha) ^ tau) * gamma[j - 1]
        c[j] <- ifelse(o[j] / h[j] > gamma[j], 1, 0)
      } 
    }
    
    # plot
    ratePlot <- sub %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / Handling,
      g = gamma,
      fitChoice = ifelse(trialRate > g, -0.25, -5),
      newChoice = ifelse(Choice == 1, -0.5, -5),
      trialRate = ifelse(trialRate > 3, 3, trialRate),
      g = ifelse(g > 3, 3.2, g)) %>%
      ggplot(aes(TrialN, g)) +
      geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
      geom_point(aes(TrialN, trialRate, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black", size = 1) +
      geom_line(aes(TrialN, optimalRate), alpha = 0.8, color = "goldenrod") + # mean optimal rate across blocks
      geom_line(size = 0.5, color = "grey50") +
      #geom_point(aes(color = Cost), size = 1.2) +
      #geom_point(aes(TrialN, fitChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #geom_point(aes(TrialN, newChoice, fill = factor(Offer, levels = c(4, 8, 20))), pch = 21, color = "black") +
      #annotate("text", x = max(sub$TrialN) + 8, y = 0.63, label = "Optimal", size = 5, color = "grey30") +
      #annotate("text", x = max(sub$TrialN), y = -0.3, label = "Predicted choices", size = 3) +
      #annotate("text", x = max(sub$TrialN), y = -0.55, label = "Observed choices", size = 3) +
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
  
  return(sum(c == sub$Choice) / nrow(sub))
}

# ma <- function(x, n = 5){
#   as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))
# }

#btw ks
ks <- trialwiseOC_btw_us %>% 
  group_by(Cost) %>% 
  summarise(mK = median(k))

# reset data: dataWth <- dataWth %>% select(-mK)
# also, consider only updating time when the offer was accepted: so cum_mean(mK * Choice) (doesn't matter, since h ^ k is multiplied by choice anyways)
dataWth <- dataWth %>% select(-mK)
if (! "mK" %in% colnames(dataWth)) {
  dataWth <- dataWth %>%
    mutate(Cost = ifelse(Cost %in% c("Wait-C", "Wait-P"), "Wait", as.character(Cost))) %>%
    left_join(ks, by = "Cost") %>%
    group_by(SubjID) %>%
    mutate(laggedK = dplyr::lag(mK, default = 0),
           mK = mK * Choice,
           mK = ifelse(mK == 0, laggedK, mK),
           mK = cumsum(mK) / TrialN,
           mK = mK * 0.2) %>% # maybe this can be the weight given to this average, a free param
    ungroup()
}

spaceSize <- 30
params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
               alpha = seq(0, 1, length.out = spaceSize))

# trialwiseOC_wth_us2 <- dataWth %>%
#   group_by(SubjID) %>%
#   do(optimize_model_dyn_us2(., params, simplify = T)) %>%
#   ungroup() %>%
#   distinct(SubjID, .keep_all = T)

trialwiseOC_wth_us2
id <- "109"

plot_dyn_us2(id, exp = "wth", gammaOne = 0.2)
plot_dyn_us(id, exp = "wth", gammaOne = 0.2)





