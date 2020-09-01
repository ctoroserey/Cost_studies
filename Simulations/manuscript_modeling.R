#!/usr/bin/env Rscript

# using NLOPTR is pointless here
# instead I'm just going to feed any sensible combination of parameters to the model and minimize the negLL manually
# simpler, can evaluate models easier, and faster
# for tomorrow: have negloglik return the parameters per iteration
# the idea being that you can enter a model and data, and the function will return the lowest LL and associated parameters
# maybe also R-squares and things in the current optimization function

#$ -pe omp 28

## libraries
library(tidyverse)
library(ggfortify)
library(knitr)
library(pander) 
library(nloptr)
#library(pwr)
library(lme4)
#library(lmerTest)
library(gridExtra)
library(reshape2)
library(corrplot)
#library(survival)
library(ggfortify)
#library(patchwork)
library(data.table)
library(parallel)

# remove break time and start counting from 0
standardize_time <- function(subjData) {
  # remove the break time (variable across subjects) and start counting time from 0 (otherwise it can add physical effort calibration)
  breakTime <- min(subjData$ExpTime[subjData$Block == 4]) - max(subjData$ExpTime[subjData$Block == 3])
  subjData$ExpTime[which(subjData$Block > 3)] <- subjData$ExpTime[which(subjData$Block > 3)] - breakTime + 16
  subjData$ExpTime <- subjData$ExpTime - min(subjData$ExpTime)
  
  return(subjData)
}

# summary matrices for refrence-changing models
betaMatrix <- function(model, rearrange = NA) {
  # get a similarity matrix of the resulting coefficient pairings for the cost conditions
  # first, do a full_join based on column names on the list of coefficient vectors from each dummy code relevel
  # then match the names of columns and rows so NAs are in the diagonal
  
  # get the names of the reference group per model iteration
  refnames <- names(model)
  
  # coefficient matrix
  temp <- lapply(model, function(data) {coefficients(data)$SubjID[1, 2:4]})
  mixCoeffs <- bind_rows(temp) 
  preln <- ifelse("Cost" %in% substr(names(mixCoeffs), 1, 4), 4, 5) # count how many characters precede the name of each cost (diff across studies)
  dimnames(mixCoeffs) <- list(refnames, substr(names(mixCoeffs), preln + 1, 20))
  mixCoeffs <- as.matrix(mixCoeffs[, match(rownames(mixCoeffs), colnames(mixCoeffs))])
  mixCoeffs[is.na(mixCoeffs)] <- 0
  
  # now the pvals
  temp <- lapply(model, function(data) {as.list(summary(data)$coefficients[2:4, 4])})
  mixPvals <- as.matrix(bind_rows(temp)) 
  dimnames(mixPvals) <- list(refnames, substr(colnames(mixPvals), preln + 1, 20))
  mixPvals <- as.matrix(mixPvals[, match(rownames(mixPvals), colnames(mixPvals))])
  mixPvals[is.na(mixPvals)] <- 1
  
  # if you would like to re-arrange the coefficient order, supply a vector with the desired sequence
  if (length(rearrange) > 1) {
    mixCoeffs <- mixCoeffs[rearrange, rearrange]
    dimnames(mixCoeffs) <- list(rearrange, rearrange)
    mixPvals <- mixPvals[rearrange, rearrange]
    dimnames(mixPvals) <- list(rearrange, rearrange)
  }
  
  # combine matrices into list to return
  out <- list(Betas = round(mixCoeffs, digits = 2),
              Pvals = round(mixPvals, digits = 5))
  
  return(out)
  
}

# simpler form of optimization that allows inputting any model expression into a single function call
# only good for static models
optimize_model_static <- function(subjData, params, model, simplify = F) {
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
  effort <- subjData$Effort
    
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

# variants of vanilla C&D
optimize_model_cd <- function(subjData, params, simplify = F, gammaStart = 0) {
  # get every combination of parameters
  params <- expand.grid(params)
  
  # relevant behavior elements
  o <- subjData$Offer
  h <- subjData$Handling + 2
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
    tempr <- params[i, "tempr"]
    alpha <- params[i, "alpha"]
    #tau <- (time - dplyr::lag(time, default = 0)) ^ s
    tau <- lag((h * c) + t) # we dont expect cognitive to affect total time, just the handling time
    
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
    p = 1 / (1 + exp(-(tempr * (o - (gamma * h ^ s)))))
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
  #tau <- (time - dplyr::lag(time, default = 0)) ^ s[[1]]
  tau <- lag((h * c) + t)
  
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
  out$probAccept <- 1 / (1 + exp(-(tempr[[1]] * (o - (gamma * h ^ s[[1]])))))
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

# to test the effect of parameters 
plot_alphas <- function(alphas, s = 1, exp = "btw", gammaStart = 0.5) {
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
          c[j] <- ifelse(o[j] / (h[j] ^ s) > gamma[j], 1, 0)
        } else {
          tau <- (h[j - 1] ^ s * c[j - 1]) + t[j - 1]
          a <- o[j - 1] * c[j - 1]
          gamma[j] <- ((1 - (1 - alpha) ^ tau) * (a / tau)) + ((1 - alpha) ^ tau) * gamma[j - 1]
          c[j] <- ifelse(o[j] / (h[j] ^ s) > gamma[j], 1, 0)
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
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
    
    # fit to subject
    baseOC <- optimize_model_static(sub, params, model_expr, simplify = T)
    
    # plot
    ratePlot <- temp %>%
      mutate(optimalRate = case_when(
        Handling == 2 ~ 0.53,
        Handling == 10 ~ 0.56,
        Handling == 14 ~ 0.62
      ),
      trialRate = Offer / (Handling ^ s),
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
          c[j] <- ifelse(o[j] / h[j] ^ s> gamma[j], 1, 0)
        } else {
          tau <- (h[j - 1] ^ s * c[j - 1]) + t[j - 1]
          a <- o[j - 1] * c[j - 1]
          gamma[j] <- ((1 - (1 - alpha) ^ tau) * (a / tau)) + ((1 - alpha) ^ tau) * gamma[j - 1]
          c[j] <- ifelse(o[j] / h[j] ^ s > gamma[j], 1, 0)
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
      trialRate = Offer / Handling ^ s,
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
# like, s = 1 and alpha = 0 returns the single parameter model
# generative version: model calculates choices per fitted previous choices, not just based on observed sub choices
optimize_model_adaptive <- function(subjData, params, simplify = F, gammaStart = 0.4) {
  #write(paste("Working on subject", unique(SubjData$SubjID)), stdout())
  
  # get every combination of parameters
  params <- expand.grid(params)
  
  # relevant behavior elements
  o <- subjData$Offer
  h <- subjData$Handling + 2 # if they accept the handling, they experience the reward window
  t <- 20 - h # and the travel includes the 2s offer window from the next trial
  obs_c <- subjData$Choice
  effort <- ifelse(subjData$Cost == "Wait", 0, 1)
  trial <- subjData$TrialN
  #fatigue <- 0.01 * subjData$TrialN * effort
  
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
  out$evolvingS <- as.numeric()
  
  # iterate through possible parameters and get the LL
  #LLs <- sapply(seq(nrow(params)), function(i) {
  for (param in seq(nrow(params))) {
    # isolate the parameters for this iteration
    # and then store them as variables
    # vectors are initialized by a single value, but they get updated
    tempr <- params[param, "tempr"]
    alpha <- params[param, "alpha"]
    ifelse("alpha_s" %in% colnames(params), alpha_s <- params[param, "alpha_s"], alpha_s <- 0) # evolving time estimate learning. turn into free parameter
    ifelse("s" %in% colnames(params), s <- rep(params[param, "s"], nrow(subjData)), s <- subjData$mS)
    ifelse("gammaPrior" %in% colnames(params), gamma <- rep(params[param, "gammaPrior"], nrow(subjData)), gamma <- rep(gammaStart, nrow(subjData)))
    
    # vectors to store evolving variables
    a <- rep(0, nrow(subjData))
    tau <- rep(0, nrow(subjData))
    S <- s # experienced s, mainly for updating rule
    #if (! "s" %in% colnames(params)) {s[1] <- 1}

    # calculate gammas iterating over trials
    # latex versions: gamma[i]: \gamma_{t+1} = (1 - (1 - \alpha) ^ {\tau_t}) \dfrac{R_{t} A_{t}}{\tau_t} +  (1 - \alpha) ^ {\tau_t} \gamma_{t}
    # \tau_{t} = H_{t} ^ {s_{t}}  A_{t} + T_{t}
    # s_{t + 1} = (1 - \alpha)^t S_1 + \sum^{t}_{i = 1} \alpha(1 - \alpha)^{t - i} (S_t A_t + s_{t - 1}(1 - A_t))
    # s_{t + 1} = s_t + A_t \left (\frac{1}{t}[S_t - s_t] \right)
    # w_t = \left (\frac{t}{max(t))} \right ) ^\alpha 
    # mS_t = (s_t w_t) + (S_t (1 - w_t))
    # P(A): P(A)_t = \dfrac{1}{1 + exp^{-(\beta[R_t - \gamma_tH_t^{s_t}])}}
    # P(A): P(A)_t = \dfrac{1}{1 + exp^{-(\beta[R_t - \gamma_tH_t^{(s_t w_t) + (S_t (1 - w_t))}])}}
    # P(A)_t = \dfrac{1}{1 + exp^{- \left (\beta \left [R_t - \gamma_tH_t^{(s_t w_t) + (S_t (1 - w_t))} \right] \right)}}
    i <- 1
    while (i < nrow(subjData)) {
      # choose if the prospect's reward rate, non-linearly discounted as above > env. rate
      # in other words, is the local-focus on handling time being affected, or a global environmental rate? (or something in between?)
      a[i] <- ifelse(o[i] / (h[i] ^ s[i])  > gamma[i], 1, 0)
      
      # amount earned
      ao <- o[i] * a[i]
      
      # s update versions
      # s[i + 1] <- alpha_s[i] * S[i] * c[i] + (1 - alpha_s[i]) * s[i]
      # s[i + 1] <- (alpha_s / i) * S[i] * c[i] + (1 - (alpha_s / i)) * s[i]
      # s[i + 1] <- s[i] + 1/i * (S[i] ^ a[i] - s[i]) # Sutton & Barto, page 37. Added choice to S[i], such that no-experienced bias estimates back to 1 (i.e. nominal time)
      s[i + 1] <- s[i] + (1/i * (S[i] - s[i])) * a[i]
      #s[l + 1] <- s[l] + (alpha_s/l * (S[l] - s[l])) * a[l]
      #left <- ((1 - alpha_s) ^ i) * s[1]
      #right <- alpha_s * (1 - alpha_s) ^ (i - seq(i)) * (S[seq(i)] ^ a[seq(i)])
      #right <- alpha_s * (1 - alpha_s) ^ (i - seq(i)) * ((S[seq(i)] * a[seq(i)]) + (dplyr::lag(s[seq(i)], default = 0)) * (1 - a[seq(i)])) # if trial isn't experienced, retain the i-1 s
      #s[i + 1] <-  left + sum(right)
      
      # non-linear estimate of the elapsed time since the last choice
      tau[i] <- (h[i] ^ s[i] * a[i]) + t[i]
      
      # gamma is updated by how much weight is given to the recently experienced reward rate (i.e. left part of eq)
      gamma[i + 1] <- ((1 - (1 - alpha) ^ tau[i]) * (ao / tau[i])) + ((1 - alpha) ^ tau[i]) * gamma[i]
      
      i <- i + 1
    }
    
    # estimate the probability of acceptance based on the difference between the offer rate vs global rate
    # this is a rehash of the eq updating c[i] above
    w <- (trial / max(trial)) ^ alpha_s
    ms <- ((s * w) + (S * (1 - w)))
    #ms <- (s * alpha_s) + (S * (1 - alpha_s))

    p = 1 / (1 + exp(-(tempr * (o - (gamma * h ^ ms)))))
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
      out$fit_c <- a
      out$params <- params[param, ]
      out$pAccept <- p
      out$evolvingS <- s
      out$mtrialS <- ms
      
      #plot(s, type = "b")
      
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
plot_adaptive_sub <- function(id = "58", exp = "btw", gammaOne = 0.5, showChoices = -0.6) {
  
  if (exp == "btw") {
    # choose subject + params
    sub <- filter(dataBtw, SubjID == id)
    spaceSize <- 30
    params <- list(tempr = seq(0, 2, length.out = spaceSize), 
                   alpha = seq(0, 0.5, length.out = spaceSize),
                   s = seq(0, 2, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_adaptive(sub, params, simplify = F, gammaStart = gammaOne)
    print(temp$params)
    
    # get baseline gamma from a simple fit
    # model to be fit
    # make sure that you specify the inverse temperature
    # extra parameters as dfs for now, that's why the `[[1]]`
    model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
    
    # create a list with possible starting values for model parameters
    # parameter names must match model ones
    params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                   gamma = seq(0.25, 1.5, length.out = spaceSize))
    
    # fit to subject
    baseOC <- optimize_model(sub, params, model_expr, simplify = T)
    
    s <- as.numeric(temp$params["s"])
    
    # plot
    ratePlot <- sub %>%
      mutate(trialRate = Offer / (Handling ^ s),
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$pAccept,
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
                   alpha_s = seq(0, 0.2, length.out = spaceSize))
    
    # run model
    temp <- optimize_model_adaptive(sub, params, simplify = F, gammaStart = gammaOne)
    print(temp$params)
    
    # plot
    ratePlot <- sub %>%
      mutate(trialRate = Offer / (Handling ^ mS),
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$pAccept,
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
plot_adaptive_fitsub <- function(id = "58", exp = "btw", showChoices = -0.6) {
  
  if (exp == "btw") {
    # choose subject + params
    sub <- filter(dataBtw, SubjID == id)
    
    # get fits
    temp <- adaptiveOC_btw[[id]]
    s <- as.numeric(temp$params["s"])
    s
    # plot
    ratePlot <- sub %>%
      mutate(Handling = Handling + 2,
             trialRate = Offer / (Handling ^ s),
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$pAccept,
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
    
    # get fits
    temp <- adaptiveOC_wth[[id]]
    
    # plot
    ratePlot <- sub %>%
      mutate(trialRate = Offer / (Handling ^ mS),
             g = temp$rate,
             fitChoice = ifelse(trialRate > g, -0.25, -5),
             newChoice = ifelse(Choice == 1, -0.5, -5),
             trialRate = ifelse(trialRate > 3, 3, trialRate),
             pChoice = temp$pAccept,
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

# simple function meant to create different choice sequences as a function of fixed parameters
# almost good enough to use in NLOPTR as the optimization function
generate_data_wth <- function(subjData, tempr = 0.5, alpha = 0.01, alpha_s = 0.2, gammaStart = 0.4) {
  # simple function meant to create different choice sequences as a function of fixed parameters
  # almost good enough to use in NLOPTR as the optimization function
  
  # relevant behavior elements
  o <- subjData$Offer
  h <- subjData$Handling + 2 # if they accept the handling, they experience the reward window
  t <- 20 - h # and the travel includes the 2s offer window from the next trial
  
  # iterate through possible parameters and get the LL
  # and then store them as variables
  # vectors are initialized by a single value, but they get updated
  s <- subjData$mS
  gamma <- rep(gammaStart, nrow(subjData))
  
  # vectors to store evolving variables
  a <- rep(0, nrow(subjData))
  tau <- rep(0, nrow(subjData))
  S <- s # experienced s, mainly for updating rule
  #if (! "s" %in% colnames(params)) {s[1] <- 1}
  
  # calculate gammas iterating over trials
  # latex versions: gamma[i]: \gamma_{t+1} = (1 - (1 - \alpha) ^ {\tau_t}) \dfrac{R_{t} A_{t}}{\tau_t} +  (1 - \alpha) ^ {\tau_t} \gamma_{t}
  # \tau_{t} = H_{t} ^ {s_{t}}  A_{t} + T_{t}
  # s_{t + 1} = (1 - \alpha)^t S_1 + \sum^{t}_{i = 1} \alpha(1 - \alpha)^{t - i} (S_t A_t + s_{t - 1}(1 - A_t))
  # P(A): P(A)_t = \dfrac{1}{1 + exp^{-(\beta[R_t - \gamma_tH_t^{s_t}])}}
  i <- 1
  while (i < nrow(subjData)) {
    # choose if the prospect's reward rate, non-linearly discounted as above > env. rate
    # in other words, is the local-focus on handling time being affected, or a global environmental rate? (or something in between?)
    a[i] <- ifelse(o[i] / (h[i] ^ s[i])  > gamma[i], 1, 0)
    
    # amount earned
    ao <- o[i] * a[i]
    
    # s update versions
    # s[i + 1] <- alpha_s[i] * S[i] * c[i] + (1 - alpha_s[i]) * s[i]
    # s[i + 1] <- (alpha_s / i) * S[i] * c[i] + (1 - (alpha_s / i)) * s[i]
    # s[i + 1] <- s[i] + 1/i * (S[i] ^ a[i] - s[i]) # Sutton & Barto, page 37. Added choice to S[i], such that no-experienced bias estimates back to 1 (i.e. nominal time)
    #s[l + 1] <- s[l] + (1/l * (S[l] - s[l])) * a[l]
    #s[l + 1] <- s[l] + (alpha_s/l * (S[l] - s[l])) * a[l]
    left <- ((1 - alpha_s) ^ i) * s[1]
    right <- alpha_s * (1 - alpha_s) ^ (i - seq(i)) * (S[seq(i)] ^ a[seq(i)])
    right <- alpha_s * (1 - alpha_s) ^ (i - seq(i)) * ((S[seq(i)] * a[seq(i)]) + (dplyr::lag(s[seq(i)], default = 0)) * (1 - a[seq(i)])) # if trial isn't experienced, retain the i-1 s
    s[i + 1] <-  left + sum(right)
    
    # non-linear estimate of the elapsed time since the last choice
    tau[i] <- (h[i] ^ s[i] * a[i]) + t[i]
    
    # gamma is updated by how much weight is given to the recently experienced reward rate (i.e. left part of eq)
    gamma[i + 1] <- ((1 - (1 - alpha) ^ tau[i]) * (ao / tau[i])) + ((1 - alpha) ^ tau[i]) * gamma[i]
    
    i <- i + 1
  }
  
  # estimate the probability of acceptance based on the difference between the offer rate vs global rate
  # this is a rehash of the eq updating c[i] above
  p = 1 / (1 + exp(-(tempr * (o - (gamma * h ^ s)))))
  p[p == 1] <- 0.999
  p[p == 0] <- 0.001
  
  # store results and return list
  out <- list()
  out$pAccept <- p
  out$rats <- gamma
  out$fit_c <- a
  out$evolvingS <- s
  
  return(out)
  
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
recover_results_wth <- function(fitsList, binary = F, matrix = T, order = F) {
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
    {if (order) group_by(., Cost, BlockOrder, Block, Offer) else group_by(., Cost, Block, Offer)} %>%
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
      {if (order) facet_wrap(vars(BlockOrder, Block)) else facet_wrap(vars(Block))} +
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
    
    mixLogis_post$Cognitive <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Wait-C"))
    mixLogis_post$`Wait-C` <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Wait-P"))
    mixLogis_post$`Wait-P` <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    mixData <- within(mixData, Cost <- relevel(Cost, ref = "Physical"))
    mixLogis_post$Physical <-  glmer(cbind(totalAccepted, totalQuits) ~ Cost + Offer + (1 | SubjID), family = "binomial", data = mixData, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    
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

# function for idea: people adapt their Ss with decreasing strength
# meaning that early in the experiment, they keep track of the perceived time per cost, but they slowly settle into previous history towards the end
# similar result to the cumulative mean, but more flexible early in the experiment
# theta vector fixed for now. Think of how to make it free
recalibrate_s <- function(S, alpha_s, choice) {
  # vector to store the evolving estimates of s
  s <- S
  #s[1] <- 1
  a <- choice
  
  l <- 1
  while (l < length(S)) {
    # s update verslons
    #s[l + 1] <- alpha_s[l] * S[l] * a[l] + (1 - alpha_s[l]) * s[l]
    #s[l + 1] <- (alpha_s / l) * S[l] * a[l] + (1 - (alpha_s / l)) * s[l]
    #s[l + 1] <- s[l] + 1/l * (S[l] ^ a[l] - s[l]) # Sutton & Barto, page 37. Added cholce to S[l], such that no-experlenced blas estlmates back to 1 (l.e. nomlnal tlme)
    #s[l + 1] <- s[l] + (1/l * (S[l] - s[l])) * a[l]
    #s[l + 1] <- s[l] + (alpha_s/l * (S[l] - s[l])) * a[l]
    left <- ((1 - alpha_s) ^ l) * s[1] ^ choice[1]
    right <- alpha_s * (1 - alpha_s) ^ (l - seq(l)) * (S[seq(l)] ^ a[seq(l)])
    right <- alpha_s * (1 - alpha_s) ^ (l - seq(l)) * ((S[seq(l)] * a[seq(l)]) + (dplyr::lag(s[seq(l)], default = 0)) * (1 - a[seq(l)])) # if trial isn't experienced, retain the i-1 s
    s[l + 1] <-  left + sum(right)
    
    l <- l + 1
  }
  
  return(s)
}

# get a summary from a group's fits from the dynamical fits
simplify_results <- function(fitLists, exp = "btw") {
  # get the relevant parameters
  temp <- sapply(fitLists, function(sub) {sub[c("percentQuit", "percentAccept", "params", "LL", "LL0", "Rsquared")]})
  
  # expand them
  df <- unnest(as.tibble(t(temp)))
  
  # append subject id and cost (if relevant)
  if (exp == "btw") {
    subList <- dataBtw %>% 
      filter(Cost != "Easy") %>% 
      distinct(SubjID, .keep_all = T) %>% 
      select(SubjID, Cost)
  } else if (exp == "wth") {
    subList <- dataWth %>% 
      distinct(SubjID, .keep_all = T) %>% 
      select(SubjID)
  }
  
  # final summary
  df <- cbind(subList, df)
  
  return(df)
}

plot_fittedS <- function(fitsList) {
  # extract fits
  fits <- do.call(c, sapply(fitsList, "[", "evolvingS"))
  data <- dataWth %>%
    mutate(evolvingS = fits) %>% #rbinom(length(fits), 1, fits))
    group_by(SubjID) %>%
    mutate(firstCost = Cost[1]) %>%
    ungroup()
  
  
  data %>%
    ggplot(aes(TrialN, evolvingS)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_line(aes(group = SubjID, color = BlockOrder), show.legend = T, alpha = 0.2) +
    geom_smooth(aes(color = BlockOrder), method = "loess") +
    theme_minimal()
}


## aesthetic options
lbls <- c("Wait","Cognitive","Physical","Easy") # between subj
colsBtw = c("#78AB05","#D9541A","deepskyblue4", "darkgoldenrod2") # plot colors (wait, effort)
colsWth <- c("#D9541A", "#78AB05", "dodgerblue4", "deepskyblue3")#"grey30", "grey70") # plot colors (wait, effort)
lthick = 2 # line thickness for plots

######## LOAD DATA
setwd("../Cost2/data")
files <- dir(pattern = '_log.csv')

# load data
dataBtw <- tibble(SubjID = files) %>%
  mutate(contents = map(SubjID, ~ suppressWarnings(read_csv(., col_types = cols()))))  %>%
  mutate(Cost = substring(SubjID, 5, 8),
         Cost = case_when(Cost == "wait" ~ "Wait",
                          Cost == "cogT" ~ "Cognitive",
                          Cost == "phys" ~ "Physical",
                          Cost == "phea" ~ "Easy"),
         SubjID = as.integer(substring(SubjID, 0, 3))) %>%
  unnest() %>%
  rename(TrialN = X1,
         ExpTime = Experiment.Time) %>%
  mutate(rawChoice = Choice, 
         RT = ifelse(RT > 14.1, 14, RT),
         Choice = ifelse(Choice == 2, 1, Choice), # forced travels (2) become acceptances (1)
         Completed = ifelse(rawChoice == 2, 0, rawChoice),
         Half = ifelse(Block < 4, "Half_1", "Half_2"),
         Cost = factor(Cost, levels = c("Physical", "Cognitive", "Wait", "Easy")),
         optimal = case_when(
           (Handling == 10 & Offer < 8) ~ 0,
           (Handling == 14 & Offer < 20) ~ 0,
           TRUE ~ 1
         )
  ) %>%
  group_by(SubjID) %>% 
  do(standardize_time(.)) %>% 
  group_by(SubjID, Block) %>%
  mutate(blockTime = ExpTime - min(ExpTime),
         blockElapsed = blockTime - dplyr::lag(blockTime, default = 0)) %>% # how much time elapsed between trials, counting per block
  ungroup()

# get a simple subject list and the number of subjects
subjList_btw <- unique(dataBtw$SubjID)
nSubjs_btw <- length(subjList_btw)


# First looks at the new data
# The RT is upper-bounded because a glitch in the code made one 10s last 14s
setwd('../../Cost3/data/')
files <- dir(pattern = 'main_log.csv')

# load the data and remove extreme subjects
dataWth <- data_frame(SubjID = files) %>% 
  mutate(contents = map(SubjID, ~ suppressWarnings(read_csv(., col_types = cols()))))  %>%
  mutate(SubjID = substring(SubjID, 0, 3)) %>%
  unnest() %>%
  mutate(RT = ifelse(RT > 10.1, 10, RT),
         Half = ifelse(Block < 4, "Half_1", "Half_2"),
         Btype = BlockType) %>%
  unite(Cost, Cost, BlockType) %>%
  rename(TrialN = X1) %>%
  mutate(rawChoice = Choice,
         Choice = ifelse(Choice == 2, 1, Choice),
         Completed = ifelse(rawChoice == 2, 0, rawChoice),
         Cost = case_when(Cost == "WAIT_0" ~ "Wait-C",
                          Cost == "COGNITIVE_0" ~ "Cognitive",
                          Cost == "GRIP_1" ~ "Physical",
                          Cost == "WAIT_1" ~ "Wait-P"),
         Cost = as.factor(Cost)) %>%
  group_by(SubjID) %>%
  mutate(BlockOrder = ifelse(Btype[1] == 0, "Cognitive1st", "Physical1st")) %>%
  group_by(SubjID) %>% 
  do(standardize_time(.)) %>% 
  group_by(SubjID, Block) %>%
  mutate(blockTime = ExpTime - min(ExpTime)) %>%
  ungroup()

# load the cognitive task performance logs
files <- dir(pattern = 'coglog')
colname <- c("Handling", "Offer", "Outcome", "RT", "Trial_Time", "ExpTime", "Trial_outcome","Type", "Setup")

dataWth_coglogs <- data_frame(SubjID = files) %>% 
  mutate(contents = map(SubjID, ~ suppressWarnings(read_csv(., col_names = colname, col_types = cols()))),
         SubjID = substring(SubjID, 8, 10)) %>%
  unnest() %>%
  mutate(Trial_time = ifelse(RT > 10.1, 10, RT),
         Trial_Time = round(Trial_Time)) %>%
  group_by(SubjID) %>%
  mutate(Half = ifelse(ExpTime < (max(ExpTime) / 2), "Half_1", "Half_2")) %>%
  ungroup()


# Get just the subject list and number of subjects
subjList_wth <- unique(dataWth$SubjID)
nSubjs_wth <- length(subjList_wth)

### MODELS
setwd('../..')

# what makes this run unique?
qualifier <- "publishable_allModels"
write(paste("Run description:", qualifier), stdout())

# how big should the parameter space be?
spaceSize <- 50
write(paste("n of possibilities per parameter:", spaceSize), stdout())

write(paste("n of cores =", detectCores()), stdout())

## which models to run?
bOC <- T
baseLogistic <- T # to test whether the brute search converges to a conventional logistic through glm()
fwOC <- F
dOC <- T
constantinoC <- T


### ADAPTIVE MODELING
# taking the cumulative average of the subjective time perception makes subjective overall time seem slower.
# a free parameter could be added to index the degree to which participants cared about this (replaces k as a free parameter, since s is adopted from exp 1)
# weird that the higher the weight on cum_means, the lower the OC ie less selective. Think about it.

# FIT BETWEEN SUBJECTS
# alpha will often default to the next lowest to 0
# so I tested a number of iterations to ensure that the results persist
# at 30 and 50 it's the same. As alpha is reduced to ~0, s increases for cog.
# but ~0.02 alpha seems sensible.
params <- list(tempr = seq(0, 2, length.out = spaceSize), 
               alpha = seq(0.001, 0.2, length.out = spaceSize),
               s = seq(0.5, 1.5, length.out = spaceSize),
               alpha_s = 1) # in the between subjects version all S_costs are the same, so there is no update. Just enforcing that here to save on computation

write("Fitting between-subject data", stdout())

# fit the model to each individual
system.time(adaptiveOC_btw <- dataBtw %>%
  filter(Cost != "Easy") %>%
  plyr::dlply("SubjID", identity) %>%
  mclapply(., optimize_model_adaptive, params, simplify = F, mc.cores = detectCores()))

# summarise
adaptiveOC_btw_summary <- simplify_results(adaptiveOC_btw)

# plot comparisons
param_compare_plot(adaptiveOC_btw_summary, "s", meanRate = 1)

# standard stats
summary(aov(s ~ Cost, adaptiveOC_btw_summary))
t.test(filter(adaptiveOC_btw_summary, Cost == "Cognitive")$s, mu = 1)

# plot result recovery
recover_results_btw(adaptiveOC_btw, binary = T)

#btw ss to apply to wth
ss <- adaptiveOC_btw_summary %>% 
  group_by(Cost) %>% 
  summarise(mS = median(s)) %>%
  rename(simpleCost = Cost) # to merge without replacing

# reset data: 
tryCatch(dataWth <- dataWth %>% select(-mS), error = function(e) {print("Oops, no need to remove mS")})
if (! "mS" %in% colnames(dataWth)) {
  suppressWarnings(dataWth <- dataWth %>%
    mutate(simpleCost = ifelse(Cost %in% c("Wait-C", "Wait-P"), "Wait", as.character(Cost))) %>% 
    left_join(ss, by = "simpleCost")) 
}

# FIT WITHIN SUBJECTS

write("Fitting within-subject data", stdout())

params <- list(tempr = seq(0, 2, length.out = spaceSize), 
               alpha = seq(0.001, 0.2, length.out = spaceSize),
               alpha_s = seq(0.001, 1, length.out = spaceSize))

# fit per individual
system.time(adaptiveOC_wth <- dataWth %>%
  plyr::dlply("SubjID", identity) %>%
  mclapply(., optimize_model_adaptive, params, simplify = F, mc.cores = detectCores()))

# summarise
adaptiveOC_wth_summary <- simplify_results(adaptiveOC_wth, exp = "wth")

# plot comparisons
alphaCors <- cor.test(adaptiveOC_wth_summary$alpha, adaptiveOC_wth_summary$alpha_s)
plot(adaptiveOC_wth_summary$alpha, adaptiveOC_wth_summary$alpha_s, 
     xlab = "Global alpha (OC)",
     ylab = "Focal alpha (S)",
     main = paste("Correlation = ", round(alphaCors$estimate[[1]], digits = 2), "p = ", round(alphaCors$p.value[[1]], digits = 4)),
     xaxt = "n")
axis(side = 1, at = unique(adaptiveOC_wth_summary$alpha), labels = FALSE)

# plot result recovery
recover_results_wth(adaptiveOC_wth, binary = T, order = T)



# MORE BASIC MODELS
## ORIGINAL OC computed using the dynamic model, showing that the adaptive version can be reduced to a single gamma fit
# great correspondence with NLOPTR, just much slower since it surveys the whole parameter space
if (bOC) {
  write("Running basic model", stdout())
  
  # model to be fit
  # make sure that you specify the inverse temperature
  # extra parameters as dfs for now, that's why the `[[1]]`
  model_expr <- expr(tempr[[1]] * (reward - (gamma[[1]] * handling)))
  
  # create a list with possible starting values for model parameters
  # parameter names must match model ones
  params <- list(tempr = seq(0, 2, length.out = spaceSize),
                 gamma = seq(0.25, 1.5, length.out = spaceSize))
  
  # fit per sub
  system.time(baseOC_btw <- dataBtw %>%
                filter(Cost != "Easy") %>%
                group_by(Cost, SubjID) %>%
                plyr::dlply("SubjID", identity) %>%
                mclapply(., optimize_model_static, params, model_expr, simplify = F, mc.cores = detectCores()))
  
  # # alternative using the big model function
  # spaceSize <- 30
  # params <- list(tempr = seq(0, 2, length.out = spaceSize), 
  #                gammaPrior = seq(0.25, 1.5, length.out = spaceSize),
  #                alpha = 0,
  #                s = 1,
  #                alpha_s = 0)
  # 
  # # fit to each subject
  # system.time(baseOC <- dataBtw %>%
  #               filter(Cost != "Easy") %>%
  #               group_by(Cost, SubjID) %>%
  #               plyr::dlply("SubjID", identity) %>%
  #               mclapply(., optimize_model_adaptive, params, simplify = F, mc.cores = detectCores()))
  
  baseOC_btw_summary <- simplify_results(baseOC_btw)
  #param_compare_plot(baseOC_summary, "gamma", meanRate = 0.7)

}

## basic logistic regression
if (baseLogistic) {
  
  write("Running logistic regressions", stdout())
  
  baseLogistic_btw <- dataBtw %>%
    filter(Cost != "Easy") %>%
    plyr::dlply("SubjID", identity) %>%
    lapply(function(data) {glm(Choice ~ Handling + Offer, data = data, family = "binomial")}) 
  
  
  baseLogistic_wth <- dataWth %>%
    plyr::dlply("SubjID", identity) %>%
    lapply(function(data) {glm(Choice ~ Cost + Offer, data = data, family = "binomial")}) 
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
  params <- list(tempr = seq(-1, 1, length.out = spaceSize), 
                 gamma = seq(0, 2, length.out = spaceSize))
  
  # fit to each subject
  fawcettOC <- data %>%
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
  
  write("Running nominal Dundon, Garrett, et al. version", stdout())
  
  # model to be fit
  # make sure that you specify the inverse temperature
  # extra parameters as dfs for now, that's why the `[[1]]`
  model_expr <- expr(tempr[[1]] * (reward - ((dplyr::lag(cumsum(reward * choice), default = 0) / (expTime ^ alpha[[1]])) * handling)))
  
  # create a list with possible starting values for model parameters
  # parameter names must match model ones
  params <- list(tempr = seq(0, 2, length.out = spaceSize), 
                 alpha = seq(0.001, 0.2, length.out = spaceSize))
  
  # fit to each subject
  dynamicOC_btw <- dataBtw %>%
    group_by(Cost, SubjID) %>%
    do(optimize_model_static(., params, model_expr, simplify = T)) %>%
    ungroup()
  
  dynamicOC_wth <- dataWth %>%
    group_by(SubjID) %>%
    do(optimize_model_static(., params, model_expr, simplify = T)) %>%
    ungroup()
}


## trial-wise updating of gamma
if (constantinoC) {
  
  write("Fitting the simple Constantino & Daw (2015) adaptive model")
  
  # Between subjects
  params <- list(tempr = seq(0, 2, length.out = spaceSize), 
                 alpha = seq(0.001, 0.2, length.out = spaceSize),
                 s = 1,
                 alpha_s = 1) # in the between subjects version all S_costs are the same, so there is no update. Just enforcing that here to save on computation
  
  # fit the model to each individual
  constantinoOC_btw <- dataBtw %>%
    filter(Cost != "Easy") %>%
    plyr::dlply("SubjID", identity) %>%
    mclapply(., optimize_model_adaptive, params, simplify = F, mc.cores = detectCores())
  
  # summarise
  constantinoOC_btw_summary <- simplify_results(adaptiveOC_btw)
  
  
  
  # Within subjects
  params <- list(tempr = seq(0, 2, length.out = spaceSize), 
                 alpha = seq(0.001, 0.2, length.out = spaceSize),
                 s = 1,
                 alpha_s = 1)
  
  # fit per individual
  constantinoOC_wth <- dataWth %>%
    plyr::dlply("SubjID", identity) %>%
    mclapply(., optimize_model_adaptive, params, simplify = F, mc.cores = detectCores())
  
  # summarise
  constantinoOC_wth_summary <- simplify_results(adaptiveOC_wth, exp = "wth")
  
}

### testing grounds

# plot_fittedmS <- function(fitsList) {
#   # extract fits
#   fits <- do.call(c, sapply(fitsList, "[", "mtrialS"))
#   data <- dataWth %>%
#     mutate(evolvingS = fits) %>% #rbinom(length(fits), 1, fits))
#     group_by(SubjID) %>%
#     mutate(firstCost = Cost[1]) %>%
#     ungroup()
#   
#   
#   data %>%
#     ggplot(aes(TrialN, evolvingS)) +
#     geom_hline(yintercept = 1, linetype = "dashed") +
#     geom_line(aes(group = SubjID, color = BlockOrder), show.legend = T, alpha = 0.2) +
#     geom_smooth(aes(color = BlockOrder), method = "loess") +
#     theme_minimal()
# }
# 
# recalibrate_s <- function(S, alpha_s, choice) {
#   # vector to store the evolving estimates of s
#   s <- S
#   #s[1] <- 1
#   a <- choice
#   
#   l <- 1
#   while (l < length(S)) {
#     # s update verslons
#     #s[l + 1] <- alpha_s[l] * S[l] * a[l] + (1 - alpha_s[l]) * s[l]
#     #s[l + 1] <- (alpha_s / l) * S[l] * a[l] + (1 - (alpha_s / l)) * s[l]
#     s[l + 1] <- s[l] + 1/l * (S[l] ^ a[l] - s[l]) # Sutton & Barto, page 37. Added cholce to S[l], such that no-experlenced blas estlmates back to 1 (l.e. nomlnal tlme)
#     s[l + 1] <- s[l] + (1/l * (S[l] - s[l])) * a[l]
#     s[l + 1] <- s[l] + (alpha_s/l * (S[l] - s[l])) * a[l]
#     left <- ((1 - alpha_s) ^ l) * (s[1] ^ choice[1])
#     right <- alpha_s * (1 - alpha_s) ^ (l - seq(l)) * (S[seq(l)] ^ a[seq(l)])
#     right <- alpha_s * (1 - alpha_s) ^ (l - seq(l)) * ((S[seq(l)] * a[seq(l)]) + (dplyr::lag(s[seq(l)], default = 0)) * (1 - a[seq(l)])) # if trial isn't experienced, retain the i-1 s
#     #s[l + 1] <-  left + sum(right)
#     
#     l <- l + 1
#   }
#   
#   return(s)
# }
# 
# generate_data_wth <- function(subjData, tempr = 0.5, alpha = 0.01, alpha_s = 0.2, gammaStart = 0.4) {
#   # simple function meant to create different choice sequences as a function of fixed parameters
#   # almost good enough to use in NLOPTR as the optimization function
#   
#   # relevant behavior elements
#   o <- subjData$Offer
#   h <- subjData$Handling + 2 # if they accept the handling, they experience the reward window
#   t <- 20 - h # and the travel includes the 2s offer window from the next trial
#   trial <- subjData$TrialN
#   
#   # iterate through possible parameters and get the LL
#   # and then store them as variables
#   # vectors are initialized by a single value, but they get updated
#   s <- subjData$mS
#   gamma <- rep(gammaStart, nrow(subjData))
#   
#   # vectors to store evolving variables
#   a <- rep(0, nrow(subjData))
#   tau <- rep(0, nrow(subjData))
#   S <- s # experienced s, mainly for updating rule
#   #s[1] <- 1
#   
#   # calculate gammas iterating over trials
#   # latex versions: gamma[i]: \gamma_{t+1} = (1 - (1 - \alpha) ^ {\tau_t}) \dfrac{R_{t} A_{t}}{\tau_t} +  (1 - \alpha) ^ {\tau_t} \gamma_{t}
#   # \tau_{t} = H_{t} ^ {s_{t}}  A_{t} + T_{t}
#   # s_{t + 1} = (1 - \alpha)^t S_1 + \sum^{t}_{i = 1} \alpha(1 - \alpha)^{t - i} (S_t A_t + s_{t - 1}(1 - A_t))
#   # P(A): P(A)_t = \dfrac{1}{1 + exp^{-(\beta[R_t - \gamma_tH_t^{s_t}])}}
#   i <- 1
#   while (i < nrow(subjData)) {
#     # choose if the prospect's reward rate, non-linearly discounted as above > env. rate
#     # in other words, is the local-focus on handling time being affected, or a global environmental rate? (or something in between?)
#     a[i] <- ifelse(o[i] / (h[i] ^ s[i])  > gamma[i], 1, 0)
#     
#     # amount earned
#     ao <- o[i] * a[i]
#     
#     # s update versions
#     #s[i + 1] <- alpha_s[i] * S[i] * a[i] + (1 - alpha_s[i]) * s[i]
#     #s[i + 1] <- (alpha_s / i) * S[i] * a[i] + (1 - (alpha_s / i)) * s[i]
#     #s[i + 1] <- s[i] + 1/i * (S[i] ^ a[i] - s[i]) # Sutton & Barto, page 37. Added choice to S[i], such that no-experienced bias estimates back to 1 (i.e. nominal time)
#     s[i + 1] <- s[i] + (1/i * (S[i] - s[i])) * a[i]
#     #s[i + 1] <- s[i] + (alpha_s/i * (S[i] - s[i])) * a[i]
#     #s[i + 1] <- s[i] + (alpha_s * (S[i] - s[i])) * a[i]
#     left <- ((1 - alpha_s) ^ i) * s[1] ^ a[1]
#     right <- alpha_s * (1 - alpha_s) ^ (i - seq(i)) * (S[seq(i)] ^ a[seq(i)])
#     right <- alpha_s * (1 - alpha_s) ^ (i - seq(i)) * ((S[seq(i)] * a[seq(i)]) + (dplyr::lag(s[seq(i)], default = 1)) * (1 - a[seq(i)])) # if trial isn't experienced, retain the i-1 s
#     #s[i + 1] <-  left + sum(right)
#     
#     # non-linear estimate of the elapsed time since the last choice
#     tau[i] <- (h[i] ^ s[i] * a[i]) + t[i]
#     
#     # gamma is updated by how much weight is given to the recently experienced reward rate (i.e. left part of eq)
#     gamma[i + 1] <- ((1 - (1 - alpha) ^ tau[i]) * (ao / tau[i])) + ((1 - alpha) ^ tau[i]) * gamma[i]
#     
#     i <- i + 1
#   }
#   
#   # estimate the probability of acceptance based on the difference between the offer rate vs global rate
#   # this is a rehash of the eq updating c[i] above
#   w <- (trial / max(trial)) ^ alpha_s
#   ms <- (S + s) / 2
#   ms <- ((s * w) + (S * (1 - w)))
#   
#   p = 1 / (1 + exp(-(tempr * (o - (gamma * h ^ ms)))))
#   p[p == 1] <- 0.999
#   p[p == 0] <- 0.001
#   
#   # store results and return list
#   out <- list()
#   out$pAccept <- p
#   out$rats <- gamma
#   out$fit_c <- a
#   out$evolvingS <- s
#   out$mtrialS <- ms
#   
#   return(out)
#   
# } 
# 
# 
# 
# 
# tempfits <- dataWth %>%
#   plyr::dlply("SubjID", identity) %>%
#   lapply(., generate_data_wth, alpha_s = 2, alpha = 0.001, tempr = 0.78)
# 
# plot_fittedS(tempfits)
# plot_fittedmS(tempfits) #+ ylim(0.8, 1.2) 
# recover_results_wth(tempfits, binary = F, matrix = F, order = T)








save.image(paste("/restricted/projectnb/cd-lab/Claudio/Cost_studies/data_", Sys.Date(), "_", qualifier, ".RData", sep = ""))

































