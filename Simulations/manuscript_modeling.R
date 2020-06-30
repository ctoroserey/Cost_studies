
# using NLOPTR is pointless here
# instead I'm just going to feed any sensible combination of parameters to the model and minimize the negLL manually
# simpler, can evaluate models easier, and faster
# for tomorrow: have negloglik return the parameters per iteration
# the idea being that you can enter a model and data, and the function will return the lowest LL and associated parameters
# maybe also R-squares and things in the current optimization function

# simpler form of optimization that allows inputting any model expression into a single function call
optimizeModel <- function(subjData, params, model, simplify = F) {
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
  expTime <- round(subjData$ExpTime)
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

# a dynamic estimation of gamma, with resulting plots for MVT-style tracking versus single-parammeter fits
dynamic_gamma <- function(sub, estimatedGamma = 0, alpha = 0.95, meanRate = 0) {
  # this function will treat prey foraging in an MVT fashion, using the equation from Dundon et al 2020
  # sub: a single subject's dataset
  # estimatedGamma: their single gamma, estimated with the original model
  # alpha: for the current model, the discounting of the effect of recent time in the updating of gamma (1 = basic Dundon)
  # meanRate: where to initialize gamma (either 0 or the actual group mean of 0.59 work)
  
  # remove the break time (variable across subjects) and start counting time from 0 (otherwise it can add physical effort calibration)
  breakTime <- min(sub$ExpTime[sub$Block == 4]) - max(sub$ExpTime[sub$Block == 3])
  sub$ExpTime[which(sub$Block > 3)] <- sub$ExpTime[which(sub$Block > 3)] - breakTime
  sub$ExpTime <- sub$ExpTime - min(sub$ExpTime)
  
  # calculate gammas as they evolve per trial
  g <- rep(0, nrow(sub))
  s <- sub$ExpTime
  r <- sub$Offer
  c <- sub$Choice
  
  for (trial in seq(nrow(sub))) {
    if (trial == 1) {
      g[trial] <- meanRate
    } else {
      # base dundon
      #g[trial] <- ((g[trial - 1] * s[trial - 1]) + (r[trial - 1] * c[trial - 1])) / s[trial]
      
      # adapted
      g[trial] <- ((g[trial - 1] * s[trial - 1]^alpha) + (r[trial - 1] * c[trial - 1])) / (s[trial]^alpha)
    }
  }
  
  # get the single gamma (base) and coarsely compare the choice fits
  sub$rate <- g
  summary <- sub %>%
    mutate(value = Offer - (rate * Handling),
           rateChoice = ifelse(Offer > (rate * Handling), 1, 0),
           gammaChoice = ifelse(Offer > (estimatedGamma * Handling), 1, 0),
           rateAcc = Choice == rateChoice,
           gammaAcc = Choice == gammaChoice)
  
  
  # plot the evolving gammas
  ratePlot <- summary %>%
    mutate(trialRate = ifelse(Offer / Handling > 2, 1.5, Offer / Handling),
           offerAccept = case_when(
             (Offer == 4 & rateChoice == 1) ~ -0.05,
             (Offer == 8 & rateChoice == 1) ~ -.075,
             (Offer == 20 & rateChoice == 1) ~ -0.1,
             TRUE ~ -6
           ),
           ratebasedChoice = ifelse(rateChoice == 1, -0.125, -7),
           actualChoice = ifelse(Choice == 1, -0.15, -8),
           optimalChoice = ifelse(optimal == 1, -0.175, -9)) %>%
    ggplot(aes(TrialN, rate)) +
    geom_line(aes(TrialN, trialRate), linetype = "dashed", size = 0.2) +
    geom_point(aes(TrialN, trialRate), size = 0.5) +
    geom_hline(yintercept = estimatedGamma) +
    geom_line(aes(color = Handling), size = 0.5) +
    geom_point(aes(color = Handling), size = 1.2) +
    geom_point(aes(TrialN, offerAccept, fill = factor(Offer, levels = c(4, 8, 20))), color = "black", pch = 21) +
    geom_point(aes(TrialN, ratebasedChoice), pch = 21, color = "black", fill = "grey20") +
    geom_point(aes(TrialN, actualChoice), pch = 21, color = "black", fill = "grey50") +
    geom_point(aes(TrialN, optimalChoice), pch = 21, color = "black", fill = "grey80") +
    scale_fill_discrete(name = "Offer") +
    scale_color_continuous(breaks = c(2, 10, 14), labels = c(2, 10, 14)) +
    ylim(-0.2, NA) +
    labs(x = "Trial Number", y = "Ongoing Opportunity Cost (gamma)") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 16))
  
  
  # compare the choices estimated with a single gamma and the rate-based gamma
  choicePatterns <- summary %>%
    group_by(Handling, Offer) %>%
    summarise(pAccept = mean(Choice),
              prateChoice = mean(rateChoice),
              pgammaChoice = mean(gammaChoice))
  
  # what proportion of choices matched participant behavior?
  accuracy <- round(c(rate = mean(summary$rateAcc), single = mean(summary$gammaAcc)), digits = 2)
  
  # get the negative log likelihood of the observations based on the model
  # this is very coarse. Since the model is normative, I did the likelihoods based on either probability 1 or 0 for choices
  # this is not like in the single parameter model, where the sigmoidal allows for more continuous probabilities
  negLL <- summary %>% 
    mutate(p = ifelse(rateChoice == 1, 0.999, 0.001), 
           tempChoice = ifelse(Choice == 1, log(p), log(1 - p))) %>%
    summarise(negLL = -sum(tempChoice))
  
  # combine outputs and return 
  results <- list(summary = summary, 
                  overallChoices = choicePatterns, 
                  accuracy = accuracy, 
                  negLL = negLL$negLL,
                  plot = ratePlot)
  
  
  return(results)
}



## which models to run?
baseOC_nloptr <- F
bOC <- T
baseLogistic <- F # to test whether the brute search converges to a conventional logistic through glm()
fwOC <- F

# which experimental data?
data <- dataBtw


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
    group_by(SubjID) %>%
    do(optimizeModel(., params, model_expr, simplify = T)) %>%
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
    group_by(SubjID) %>%
    do(optimizeModel(., params, model_expr, simplify = T)) %>%
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
    group_by(SubjID, Cost) %>%
    do(optimizeModel(., params, model_expr, simplify = T)) %>%
    ungroup()
}



##  ADD OTHERS

# Dundon, Garrett, et al (2020)
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

# one possibility for this is that by tracking an evolving rate, sequential quits would reduce the OC, prompting lower reward acceptances
# instead, our participants don't seem to be influenced by that, having a good and stable sense of the environmental rate


# get the base rate of the environment from the final group average
# starting from 0 converges to a very similar result
meanRate <- dataBtw %>% 
  filter(Choice == 1) %>% 
  group_by(SubjID, Cost, Block) %>%
  summarise(mEarn = sum(Offer) / max(blockTime)) %>%
  ungroup() %>%
  summarise(m = mean(mEarn))

# which subject to run?
subjs <- dataBtw %>% plyr::dlply("SubjID", identity)

# fit
#r <- dynamic_gamma(subjs[[1]], 0.98)
allResults <- lapply(as.character(subjList_btw), function(sub) {dynamic_gamma(subjs[[sub]], 
                                                                              filter(baseOC, SubjID == sub)$gamma, 
                                                                              meanRate = meanRate$m, 
                                                                              alpha = 1.01)})

allResults[[10]]$plot

# compare the accuracy and negLL between fits 
# for a fixed alpha of 0.95, it looks like a single gamma works best
# see 275 for a poor MVT fit from a cognitive participant
# accuracy
fits <- sapply(allResults, "[[", "accuracy")
plot(fits["rate",], fits["single",])
abline(a = 0, b = 1)

# negLL (note that this type of model comparison might not be well suited)
plot(sapply(allResults, "[[", "negLL"), baseOC$LL)
abline(a = 0, b = 1)

# let's fit this alpha parameter
as <- seq(0.5, 1.5, length.out = 20)
temp <- lapply(as, function(alpha) {dynamic_gamma(subjs[["58"]], 
                                                  filter(baseOC, SubjID == "58")$gamma, 
                                                  meanRate = meanRate$m, 
                                                  alpha = alpha)})







