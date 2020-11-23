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

source('models_gk/load_data.R')
source('models_gk/likelihoods.R')

#########################
### Load Data
#########################

c(dat_btw, dat_wth) %<-% load_data_full()
fit_cols = c("sub", "cond", "Offer", "Choice", "HandleTime", "TravelTime", "condCode")
dat_btw[, Subject := sub]

#########################
### plot between subjects behavior
#########################

mean_sub_btw = dat_btw[, .(Choice = mean(Choice)), .(sub, cond, Handling, Offer)]
plot_btw(mean_sub_btw)

#########################
### Fit super base model -- optimal values + decision noise
#########################

par_names_base = c("beta")
mus_base = c("beta" = 0.1)
sigma_base = c(1)

system.time(fits_btw_base <- dat_btw[, em(.SD,
                                          log(mus_base),
                                          sigma_base,
                                          cost_lik,
                                          lower = -100,
                                          upper = 100,
                                          package = "optim",
                                          method = "Brent",
                                          bic=TRUE,
                                          ndata=.SD[, .N],
                                          obj_args = list(
                                            par_names=par_names_base,
                                            log_par=T
                                          ),
                                          verbose = F,
                                          return_ind_df=T),
                                     .(cond)])

fits_btw_base_predict = get_all_predictions(dat_btw, fits_btw_base, par_names_base, log_par=T, return_mean=T)
plot_btw(fits_btw_base_predict, title="Base Model (MVT + Noise)")
plot_params(fits_btw_base, par_names_base, exp_par=T)


#########################
### Fit base opportunity cost model
#########################

par_names_oc = c("gamma", "beta")
mus_oc = c("gamma" = 1, "beta" = .1)
sigma_oc = c(1, 1)

system.time(fits_btw_oc <- dat_btw[, em(.SD,
                                        log(mus_oc),
                                        sigma_oc,
                                        cost_lik,
                                        bic=TRUE,
                                        ndata=.SD[, .N],
                                        obj_args = list(
                                          par_names=par_names_oc,
                                          log_par=T
                                        ),
                                        perturb_start=0.1,
                                        verbose=T,
                                        parallel=T,
                                        return_ind_df=T),
                                   .(cond)])

fits_btw_oc_predict = get_all_predictions(dat_btw, fits_btw_oc, par_names_oc, log_par=T, return_mean=T)
plot_btw(fits_btw_oc_predict, title="Opportunity Cost")
plot_params(fits_btw_oc, par_names_oc, exp_par=T)


#########################
### Fit time sensitivity model
#########################

par_names_ts = c("s", "beta")
mus_ts = c("s" = 1, "beta" = .1)
sigma_ts = c(1, 1)

system.time(fits_btw_ts <- dat_btw[, em(.SD,
                                        log(mus_ts),
                                        sigma_ts,
                                        cost_lik,
                                        bic=TRUE,
                                        ndata=.SD[, .N],
                                        obj_args = list(
                                          par_names=par_names_ts,
                                          log_par=T
                                        ),
                                        perturb_start=0.1,
                                        verbose=T,
                                        parallel=T,
                                        return_ind_df=T),
                                   .(cond)])

fits_btw_ts_predict = get_all_predictions(dat_btw, fits_btw_ts, par_names_ts, log_par=T, return_mean=T)
plot_btw(fits_btw_ts_predict, title="Time Sensitivity")
plot_params(fits_btw_ts, par_names_ts, exp_par=T)


#########################
### Compare models
#########################

fits_btw_base[, model := "base"]
fits_btw_oc[, model:= "oc"]
fits_btw_ts[, model := "ts"]
fits_all = rbind(fits_btw_base, fits_btw_oc, fits_btw_ts, fill=T)
fits_bic = fits_all[, .(bic = bic[1]), .(model, cond)]

ggplot(fits_bic, aes(x=model, y=bic, group=0)) +
  facet_wrap(~cond) +
  geom_point() +
  geom_line() +
  theme_abw()

ggplot(fits_bic, aes(x=model, y=bic, group=0)) +
  stat_summary(fun=sum, geom="point") +
  stat_summary(fun=sum, geom="line") +
  theme_abw()

#########################
### Save workspace
#########################

save.image("models_gk/results/btw.Rdata")
