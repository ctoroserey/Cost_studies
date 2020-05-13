
# using NLOPTR is pointless here
# instead I'm just going to feed any sensible combination of parameters to the model and minimize the negLL manually
# simpler, can evaluate models easier, and faster
# for tomorrow: have negloglik return the parameters per iteration
# the idea being that you can enter a model and data, and the function will return the lowest LL and associated parameters
# maybe also R-squares and things in the current optimization function



  
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
  out$Rsquared <- 1 - (out$LL / out$LL0) # pseudo r-squred, quantifying the proportion of deviance reduction vs chance
  out$probAccept <- 1 / (1 + exp(-eval(model)))
  out$Params <- chosen_params
  #out$predicted <- reward > out$subjOC
  #out$predicted[out$predicted == TRUE] <- 1
  #out$percentPredicted <- mean(out$predicted == choice) 
  
  # if doing this with dplyr::do(), return a simplified data.frame instead with the important parameters
  if (simplify) {
    out <- round(data.frame(out[-6]), digits = 2)
    colnames(out) <- c("percentQuit",
                       "percentAccept",
                       "LL",
                       "LL0",
                       "Rsq",
                       colnames(chosen_params))
  }
  
  return(out)
}



# model to be fit
# make sure that you specify the inverse temperature
# extra parameters as dfs for now, that's why the `[[1]]`
model_expr <- expr(temp[[1]] * (reward - (gamma[[1]] * handling)))

# create a list with possible starting values for model parameters
spaceSize <- 100
params <- list(temperature = seq(-1, 1, length.out = spaceSize), 
               Gamma = seq(0.25, 1.5, length.out = spaceSize))

# fit to each subject
tempOC <- dataBtw %>%
  group_by(SubjID, Cost) %>%
  do(optimizeModel(., params, model_expr, simplify = T)) %>%
  ungroup()


# plot
ggplot(tempOC, aes(Cost, Gamma, fill = Cost)) +
  geom_hline(yintercept = 0.7, alpha = 0.9, color = "gray40", size = 1, linetype = "dashed") +
  geom_jitter(pch = 21, size = 3, show.legend = F) +
  ylim(0, 1.5) +
  labs(x = "") +
  scale_fill_manual(values = colsWth) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 16))

# summaryOC$all %>%
#   select(SubjID, Cost, Gamma, Half) %>%
#   spread(Half, Gamma) %>%
#   ggplot(aes(Half_1, Half_2, fill = Cost)) +
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#     geom_point(pch = 21, color = "black", size = 3) +
#     ylim(0, 2) +
#     xlim(0, 2) +
#     theme_minimal()



































