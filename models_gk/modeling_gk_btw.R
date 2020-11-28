#########################
### Setup
#########################

# set wd
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

load("models_gk/results/btw.Rdata")

library(data.table)
library(foreach)
library(rextras)
library(zeallot)
library(modelfitr)
library(GaussExMax)

source('models_gk/load_data.R')
source('models_gk/likelihoods.R')
Rcpp::sourceCpp('models_gk/likelihoods.cpp')


#########################
### Load Data
#########################

c(dat_btw, dat_wth) %<-% load_data_full()
dat_btw[, Subject := sub]

#########################
### plot within subjects behavior
#########################

mean_sub_btw = dat_btw[, .(Choice = mean(Choice)), .(sub, cond, Handling, Offer)]
plot_btw(mean_sub_btw)

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

models[["Biased OC"]] = list(
  par_names = c("cost", "beta"),
  mus = c(cost=1, beta=1),
  sigma = c(10, 10),
  options = list(method="nlminb")
)

models[["Time Sensitivity"]] = list(
  par_names = c("s", "beta"),
  mus = c(s=1,  beta=1),
  sigma = c(10, 10),
  options = list(method="nlminb")
)

models[["Adaptive OC + Biased OC"]] = list(
  par_names = c("gamma", "alpha", "cost", "beta"),
  mus = c(gamma=0.55, alpha=.01, cost=1, beta=1),
  sigma = c(10, 10, 10, 10),
  options = list(method="nlminb")
)

models[["Adaptive OC + Time Sensitivity"]] = list(
  par_names = c("gamma", "alpha", "s", "beta"),
  mus = c(gamma=0.55, alpha=.01, s=1, beta=1),
  sigma = c(10, 10, 10, 10),
  options = list(method="nlminb")
)

#########################
### Fit Models
#########################

fits_btw = foreach(i=1:length(models)) %do% {
  
  cat("model =", names(models)[i], "\n\n")
  
  m = models[[i]]
  
  this_em <- do.call(em, c(list(dat=dat_btw,
                                mus=log(m$mus),
                                sigma=m$sigma,
                                lik=cost_lik,
                                bic=TRUE,
                                ndata=dat_btw[, .N],
                                obj_args=list(
                                  par_names=m$par_names,
                                  log_par=T
                                ),
                                perturb_start=0.1,
                                emtol=.001,
                                parallel=T,
                                verbose=T),
                           m$options)
                     )
  
  this_em
}

names(fits_btw) = names(models)


#########################
### Compare models
#########################

fits_bic = melt(as.data.table(lapply(fits_btw, function(x) x$bic)), variable.name="model", value.name="bic")

ggplot(fits_bic, aes(x=model, y=bic, group=0)) +
  geom_point() +
  geom_line() +
  theme_abw(x_angle=45)


#########################
### Plot Model Predictions
#########################

model_to_plot = "Adaptive OC + Biased OC"

params_dt = setDT(get_ind_df(fits_btw[[model_to_plot]]))
predictions = get_all_predictions(dat_btw, params_dt, par_names=models[[model_to_plot]]$par_names, return_mean=T, log_par=T)
plot_btw(predictions)

#########################
### Save workspace
#########################

save.image("models_gk/results/btw.Rdata")
