
# using NLOPTR is pointless here
# instead I'm just going to feed any sensible combination of parameters to the model and minimize the negLL manually
# simpler, can evaluate models easier, and faster
# for tomorrow: have negloglik return the parameters per iteration
# the idea being that you can enter a model and data, and the function will return the lowest LL and associated parameters
# maybe also R-squares and things in the current optimization function



  
negLogLik <- function(subjData, params, model) {
  # this function finds the combination of parameter values that minimizes the neg log likelihood of a logistic regression
  # used to rely on NLOPTR, but it's too cumbersome for the low-dimensional estimates I'm performing.
  #
  # subjData: a participant's log
  # params: a list of vectors. Each vector is the possible values a given parameter can take. Names in list must match model expression
  # model: using `expr()`, define the model (use <param>[[1]] for free parameters to be estimated. R limitation.)
  
  
  # extract basic choice information
  handling <- subjData$Handling
  reward <- subjData$Offer
  choice <- subjData$Choice
  rt <- subjData$RT
  cost <- subjData$Cost
  trial <- subjData$TrialN
  
  # combine parameters into every possible combination
  params <- expand.grid(params)
  
  LLs <- sapply(seq(nrow(params)), function(i) {
    
    # isolate the parameters for this iteration
    # and then store them as variables
    # FIGURE OUT HOW TO NOT STORE THEM AS DATAFRAMES
    pars <- params_grid[i, ]
    lapply(seq_along(pars), function(variable) {assign(colnames(pars)[variable], pars[variable], envir = .GlobalEnv)})
    
    # estimate the probability of acceptance per the model
    p = 1 / (1 + exp(-eval(model_expr)))
    p[p == 1] <- 0.999
    p[p == 0] <- 0.001
    
    # get the likelihood of the observations based on the model
    tempChoice <- rep(NA, length(choice))
    tempChoice[choice == 1] <- log(p[choice == 1])
    tempChoice[choice == 0] <- log(1 - p[choice == 0]) # log of probability of choice 1 when choice 0 occurred
    negLL <- -sum(tempChoice)
  })
    
    
  return(LLs)
}

# model to be fit
# make sure that you specify the inverse temperature
# extra parameters as dfs for now, that's why the `[[1]]`
model_expr <- expr(temp[[1]] * (reward - (gamma[[1]] * handling)))

# create a list with possible starting values for model parameters
params <- list(temp = seq(-1, 1, length.out = 100), 
               gamma = seq(0.25, 1.5, length.out = 100))

# combine into every possible combination
params_grid <- expand.grid(params)

subj <- 190
subjData <- dataBtw %>%
  filter(SubjID == subj)

test <- negLogLik(subjData, params_grid, model_expr)




# plot the likelihood over the parameter space (surprisingly convex)
plot(t)

min(t)
params_grid[which(t == min(t)), ]
summaryOC$all %>% filter(SubjID == subj)


tempOC <- dataBtw %>%
  group_by(SubjID) %>%
  summarise(min_negLL = min(sapply(seq(nrow(params_grid)), function(i) negLogLik(params_grid[i, ], model_expr, Choice, Handling, Offer))))




# summaryOC <- list()
# summaryOC$all <- dataWth %>%
#   mutate(Block = case_when(
#     Block %in% c(1, 2) ~ 1,
#     Block %in% c(3, 4) ~ 2,
#     Block %in% c(5, 6) ~ 3
#   )) %>%
#   group_by(SubjID, Cost, Half) %>% #group_by(SubjID, Cost, Block) %>%
#   do(optimizeOCModel(., simplify = T)) %>%
#   ungroup()
# 
# # plot
# (summaryOC$plot <- ggplot(summaryOC$all, aes(Cost, Gamma, fill = Cost)) +
#   geom_hline(yintercept = 0.7, alpha = 0.9, color = "gray40", size = 1, linetype = "dashed") +
#   geom_jitter(pch = 21, size = 3, show.legend = F) +
#   ylim(0, 1.5) +
#   labs(x = "") +
#   scale_fill_manual(values = colsWth) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 16)))
# 
# summaryOC$all %>%
#   select(SubjID, Cost, Gamma, Half) %>%
#   spread(Half, Gamma) %>%
#   ggplot(aes(Half_1, Half_2, fill = Cost)) +
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#     geom_point(pch = 21, color = "black", size = 3) +
#     ylim(0, 2) +
#     xlim(0, 2) +
#     theme_minimal()



































