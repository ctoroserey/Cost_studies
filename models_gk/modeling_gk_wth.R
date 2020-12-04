#########################
### Setup
#########################

# set wd
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

library(data.table)
library(foreach)
library(rextras)
library(zeallot)
library(modelfitr)
library(GaussExMax)

load("models_gk/results/btw.Rdata")
load("models_gk/results/wth.Rdata")
source('models_gk/load_data.R')
source('models_gk/likelihoods.R')
Rcpp::sourceCpp('models_gk/likelihoods.cpp')


#########################
### Load Data
#########################

c(dat_btw, dat_wth) %<-% load_data_full()
dat_wth[, Subject := sub]

#########################
### plot within subjects behavior
#########################

mean_sub_wth = dat_wth[Block %in% c(1,2,5,6), .(Choice = mean(Choice)), .(sub, cond, half, first, BlockType, Handling, Offer)]

plot_wth(mean_sub_wth)

# glm1 = lme4::glmer(data=dat_wth, Choice ~ 1 + Offer*cond*half*first + (1 | sub), family=binomial)
# summary(glm1)
#
# phia::testInteractions(glm1, pairwise="cond")
# phia::testInteractions(glm1, pairwise="half")
# phia::testInteractions(glm1, pairwise="first")
# phia::testInteractions(glm1, across="cond", pairwise="half")
# phia::testInteractions(glm1, across="cond", pairwise="first")

#########################
### Define Models
#########################

models = list()

models[["Base"]] = list(
  par_names = "beta",
  mus = c(beta = 1),
  sigma = 10,
  options = list(package="optim",
                 method="Brent",
                 lower=-100,
                 upper=100)
)

models[["OC"]] = list(
  par_names = c("gamma", "beta"),
  mus = c(gamma=0.55, beta=1),
  sigma = c(10, 10),
  options = list(method="nlminb")
)

models[["Adaptive OC"]] = list(
  par_names = c("gamma", "alpha", "beta"),
  mus = c(gamma=0.55, alpha=.01, beta=1),
  sigma = c(10, 10, 10),
  options = list(method="nlminb")
)

models[["Biased OC"]] = list(
  par_names = c("cost_1", "cost_2", "cost_3", "beta"),
  mus = c(cost_1=1, cost_2=1, cost_3=1, beta=1),
  sigma = c(10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Time Sensitivity"]] = list(
  par_names = c("s_1", "s_2", "s_3", "beta"),
  mus = c(s_1=1, s_2=1, s_3=1, beta=1),
  sigma = c(10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Adaptive OC + Biased OC"]] = list(
  par_names = c("gamma", "alpha", "cost_1", "cost_2", "cost_3", "beta"),
  mus = c(gamma=0.55, alpha=.01, cost_1=1, cost_2=1, cost_3=1, beta=1),
  sigma = c(10, 10, 10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Adaptive OC + Time Sensitivity"]] = list(
  par_names = c("gamma", "alpha", "s_1", "s_2", "s_3", "beta"),
  mus = c(gamma=0.55, alpha=.01, s_1=1, s_2=1, s_3=1, beta=1),
  sigma = c(10, 10, 10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Adaptive Biased OC"]] = list(
  par_names = c("cost_1", "cost_2", "cost_3", "C", "w_decay", "beta"),
  mus = c(cost_1=1, cost_2=1, cost_3=1, C=1, w_decay=1, beta=1),
  sigma = c(10, 10, 10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Adaptive Time Sensitivity"]] = list(
  par_names = c("s_1", "s_2", "s_3", "S", "w_decay", "beta"),
  mus = c(s_1=1, s_2=1, s_3=1, S=1, w_decay=1, beta=1),
  sigma = c(10, 10, 10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Adaptive OC + Adaptive Biased OC"]] = list(
  par_names = c("gamma", "alpha", "cost_1", "cost_2", "cost_3", "C", "w_decay", "beta"),
  mus = c(gamma=0.55, alpha=0.01, cost_1=1, cost_2=1, cost_3=1, C=1, w_decay=1, beta=1),
  sigma = c(10, 10, 10, 10, 10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Adaptive OC + Adaptive Time Sensitivity"]] = list(
  par_names = c("gamma", "alpha", "s_1", "s_2", "s_3", "S", "w_decay", "beta"),
  mus = c(gamma=0.55, alpha=0.01, s_1=1, s_2=1, s_3=1, S=1, w_decay=1, beta=1),
  sigma = c(10, 10, 10, 10, 10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Between Time Sensitivity"]] = list(
  par_names = c("w_decay", "beta"),
  mus = c(w_decay=1, beta=1),
  sigma = c(10, 10),
  options = list(method="nlminb"),
  fixed_pars = exp(c(s_1 = fits_btw[["Time Sensitivity"]][["wait"]]$mus[1],
                     s_2 = fits_btw[["Time Sensitivity"]][["cogTask"]]$mus[1],
                     s_3 = fits_btw[["Time Sensitivity"]][["phys"]]$mus[1],
                     S = log(1)))
)

models[["Adaptive Between Biased OC"]] = list(
  par_names = c("C", "w_decay", "beta"),
  mus = c(C=1, w_decay=1, beta=1),
  sigma = c(10, 10, 10),
  options = list(method="nlminb"),
  fixed_pars = c(cost_1 = fits_btw[["Biased OC"]][["wait"]]$mus[1],
                 cost_2 = fits_btw[["Biased OC"]][["cogTask"]]$mus[1],
                 cost_3 = fits_btw[["Biased OC"]][["phys"]]$mus[1])
)

models[["Adaptive Between Time Sensitivity"]] = list(
  par_names = c("S", "w_decay", "beta"),
  mus = c(S=1, w_decay=1, beta=1),
  sigma = c(10, 10, 10),
  options = list(method="nlminb"),
  fixed_pars = exp(c(s_1 = fits_btw[["Time Sensitivity"]][["wait"]]$mus[1],
                 s_2 = fits_btw[["Time Sensitivity"]][["cogTask"]]$mus[1],
                 s_3 = fits_btw[["Time Sensitivity"]][["phys"]]$mus[1]))
)

models[["Adaptive OC + Adaptive Between Biased OC"]] = list(
  par_names = c("gamma", "alpha", "C", "w_decay", "beta"),
  mus = c(gamma=0.55, alpha=0.01, C=1, w_decay=1, beta=1),
  sigma = c(10, 10, 10, 10, 10),
  options = list(method="nlminb"),
  fixed_pars = c(cost_1 = fits_btw[["Adaptive OC + Biased OC"]][["wait"]]$mus[3],
                 cost_2 = fits_btw[["Adaptive OC + Biased OC"]][["cogTask"]]$mus[3],
                 cost_3 = fits_btw[["Adaptive OC + Biased OC"]][["phys"]]$mus[3])
)

models[["Adaptive OC + Adaptive Between Time Sensitivity"]] = list(
  par_names = c("gamma", "alpha", "S", "w_decay", "beta"),
  mus = c(gamma=0.55, alpha=0.01, S=1, w_decay=1, beta=1),
  sigma = c(10, 10, 10, 10, 10),
  options = list(method="nlminb"),
  fixed_pars = exp(c(s_1 = fits_btw[["Adaptive OC + Time Sensitivity"]][["wait"]]$mus[3],
                 s_2 = fits_btw[["Adaptive OC + Time Sensitivity"]][["cogTask"]]$mus[3],
                 s_3 = fits_btw[["Adaptive OC + Time Sensitivity"]][["phys"]]$mus[3]))
)

#########################
### Fit Models
#########################

models_to_fit = models[names(models)[!(names(models) %in% names(fits_wth))]]

fits_wth_new = foreach(i=1:length(models_to_fit)) %do% {
  
  cat("model =", names(models_to_fit)[i], "\n\n")
  
  m = models_to_fit[[i]]
  if (!is.null(m$fixed_pars)) fixed_pars = m$fixed_pars else fixed_pars = c()
  
  this_em <- do.call(em, c(list(dat=dat_wth,
                                mus=log(m$mus),
                                sigma=m$sigma,
                                lik=cost_lik,
                                bic=TRUE,
                                ndata=dat_wth[, .N],
                                obj_args=list(
                                  par_names=m$par_names,
                                  log_par=T,
                                  fixed_pars=fixed_pars
                                ),
                                perturb_start=0.1,
                                emtol=.001,
                                parallel=T,
                                verbose=T),
                           m$options)
  )
  
  this_em_dt = setDT(get_ind_df(this_em))
  dat_predict = get_all_predictions(copy(dat_wth), this_em_dt, m$par_names, within=T, log_par=T, return_mean=T)
  plot_wth(dat_predict, title=names(models)[i])
  plot_params(this_em_dt, m$par_names, exp_par=T, within=T, title=names(models)[i])
  this_em
}

names(fits_wth_new) = names(models_to_fit)

fits_wth = c(fits_wth, fits_wth_new)
fits_wth = fits_wth[names(models)]

#########################
### Compare models
#########################

fits_bic = melt(as.data.table(lapply(fits_wth, function(x) x$bic)), variable.name="model", value.name="bic")

ggplot(fits_bic, aes(x=model, y=bic, group=0)) +
  geom_point() +
  geom_line() +
  theme_abw(x_angle=45)


#########################
### Plot Model Predictions
#########################

model_to_plot = "Adaptive OC + Adaptive Between Time Sensitivity"

params_dt = setDT(get_ind_df(fits_wth[[model_to_plot]]))
if (!is.null(models[[model_to_plot]]$fixed_pars)) fixed_pars = models[[model_to_plot]]$fixed_pars else fixed_pars = c()
predictions = get_all_predictions(dat_wth, params_dt, par_names=models[[model_to_plot]]$par_names, fixed_pars=fixed_pars, within=T, return_mean=T, log_par=T)
plot_wth(predictions, title=model_to_plot)
plot_params(params_dt, par_names=models[[model_to_plot]]$par_names, exp_par=T, within=T)

#########################
### Save workspace
#########################

save.image("models_gk/results/wth.Rdata")
