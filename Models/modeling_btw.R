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

models[["OC"]] = list(
  par_names = c("gamma", "beta"),
  mus = c(gamma=1, beta = 1),
  sigma = c(10, 10),
  options = list(method="nlminb")
)

models[["Adaptive OC"]] = list(
  par_names = c("gamma", "alpha", "beta"),
  mus = c(gamma=1, alpha=.01, beta = 1),
  sigma = c(10, 10, 10),
  options = list(method="nlminb")
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

conds = dat_btw[, unique(cond)]

fits_btw = foreach(i=1:length(models)) %do% {
  
  fits_conds = foreach(c=1:length(conds)) %do% {
  
    cat("cond =", conds[c], "// model =", names(models)[i], "\n\n")
    
    m = models[[i]]
    
    this_em <- do.call(em, c(list(dat=dat_btw[cond==conds[c]],
                                  mus=log(m$mus),
                                  sigma=m$sigma,
                                  lik=cost_lik,
                                  bic=TRUE,
                                  ndata=dat_btw[cond==conds[c], .N],
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
  
  names(fits_conds) = conds
  fits_conds
  
}

names(fits_btw) = names(models)


#########################
### Compare models
#########################

fits_bic = melt(as.data.table(lapply(fits_btw, function(x) sum(sapply(x, function(x) x$bic)))), variable.name="model", value.name="bic")

ggplot(fits_bic, aes(x=model, y=bic, group=0)) +
  geom_point() +
  geom_line() +
  theme_abw(x_angle=45)


#########################
### Plot Model Predictions
#########################

model_to_plot = "Time Sensitivity"

params_list = lapply(fits_btw[[model_to_plot]], function(x) get_ind_df(x))
for(i in 1:length(params_list)) params_list[[i]]$cond = names(params_list)[i]
params_dt = rbindlist(lapply(params_list, setDT))
predictions = get_all_predictions(dat_btw, params_dt, par_names=models[[model_to_plot]]$par_names, return_mean=T, log_par=T)
plot_btw(predictions)
plot_params(params_dt, par_names = models[[model_to_plot]]$par_names, exp_par=T)

#########################
### Save workspace
#########################

save.image("models_gk/results/btw.Rdata")
