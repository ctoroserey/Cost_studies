
# using NLOPTR is pointless here
# instead I'm just going to feed any sensible combination of parameters to the model and minimize the negLL manually
# simpler, can evaluate models easier, and faster
# for tomorrow: have negloglik return the parameters per iteration
# the idea being that you can enter a model and data, and the function will return the lowest LL and associated parameters
# maybe also R-squares and things in the current optimization function

negLogLik <- function(subjData, params, model) {
  
  handling <- subjData$Handling
  reward <- subjData$Offer
  choice <- subjData$Choice
  
  LLs <- apply(params, 2, function(i) {
    p = 1 / (1 + exp(-eval(model)))
    p[p == 1] <- 0.999
    p[p == 0] <- 0.001
    tempChoice <- rep(NA, length(choice))
    tempChoice[choice == 1] <- log(p[choice == 1])
    tempChoice[choice == 0] <- log(1 - p[choice == 0]) # log of probability of choice 1 when choice 0 occurred
    negLL <- -sum(tempChoice)
  })
    
    
  return(negLL)
}

# model to be fit
model_expr <- expr(params["temp"] * (reward - (params["gamma"] * handling)))

# create a list with possible starting values for model parameters
params <- list(temp = seq(-1, 1, length.out = 100), 
               gamma = seq(0.25, 1.5, length.out = 100))

# combine into every possible combination
params_grid <- as.matrix(expand.grid(params))

subj <- 190
Data <- dataBtw %>%
  filter(SubjID == subj)


t <- sapply(seq(nrow(params_grid)), function(i) negLogLik(params_grid[i, ], model_expr, choice, handling, reward))

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



































