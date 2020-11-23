#########################
### Setup
#########################

# set wd
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

load("models_gk/results/wth.Rdata")

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
### Fit super base model -- optimal values + decision noise
#########################

par_names_base = c("beta")
mus_base = c("beta" = 0.1)
sigma_base = c(1)

system.time(fits_wth_base_em <- em(dat_wth,
                                   log(mus_base),
                                   sigma_base,
                                   cost_lik,
                                   lower = -100,
                                   upper = 100,
                                   package = "optim",
                                   method = "Brent",
                                   bic=TRUE,
                                   ndata=dat_wth[, .N],
                                   obj_args = list(
                                     par_names=par_names_base,
                                     log_par=T
                                   ),
                                   parallel=F,
                                   verbose = F,
                                   return_ind_df=T))

fits_wth_base = setDT(fits_wth_base_em)
fits_wth_base_predict = get_all_predictions(dat_wth, fits_wth_base, par_names_base, within=T, log_par=T, return_mean=T)
plot_wth(fits_wth_base_predict, title="Base Model (MVT + Noise)")
plot_params(fits_wth_base, par_names_base, exp_par=T, within=T)


#########################
### Fit opportunity cost model
#########################

par_names_oc = c("cost_1", "cost_2", "cost_3", "beta")
mus_oc = c("cost_1" = 1, "cost_2" = 1, "cost_3" = 1, "beta" = 0.1)
sigma_oc = c(1, 1, 1, 1)

system.time(fits_wth_oc_em <- em(dat_wth,
                                 log(mus_oc),
                                 sigma_oc,
                                 cost_lik,
                                 bic=TRUE,
                                 ndata=dat_wth[, .N],
                                 obj_args = list(
                                   par_names=par_names_oc,
                                   log_par=T
                                 ),
                                 perturb_start=0.1,
                                 verbose=T,
                                 parallel=T))

fits_wth_oc = setDT(get_ind_df(fits_wth_oc_em))
fits_wth_oc_predict = get_all_predictions(dat_wth, fits_wth_oc, par_names_oc, within=T, log_par=T, return_mean=T)
plot_wth(fits_wth_oc_predict, title="Biased OC")
plot_params(fits_wth_oc, par_names_oc, exp_par=T, within=T)


#########################
### Fit opportunity cost 2
#########################

par_names_oc2 = c("cost_1", "cost_2", "cost_3", "C", "beta", "w_decay")
mus_oc2 = c("cost_1" = 1, "cost_2" = 1, "cost_3" = 1, "C" = 1, "beta" = 0.1, "w_decay" = 1)
sigma_oc2 = c(1, 1, 1, 1, 1, 1)

system.time(fits_wth_oc2_em <- em(dat_wth,
                                  log(mus_oc2),
                                  sigma_oc2,
                                  cost_lik,
                                  bic=TRUE,
                                  ndata=dat_wth[, .N],
                                  obj_args = list(
                                    par_names=par_names_oc2,
                                    log_par=T
                                  ),
                                  perturb_start=0.1,
                                  parallel=T,
                                  verbose=T))

fits_wth_oc2 = setDT(get_ind_df(fits_wth_oc2_em))
fits_wth_oc2_predict = get_all_predictions(dat_wth, fits_wth_oc2, par_names_oc2, within=T, log_par=T, return_mean=T)
plot_wth(fits_wth_oc2_predict, title="Adaptive Biased OC")
plot_params(fits_wth_oc2, par_names_oc2, exp_par=T, within=T)

#########################
### Fit time sensitivity model
#########################

par_names_ts = c("s_1", "s_2", "s_3", "beta")
mus_ts = c("s_1" = 1, "s_2" = 1, "s_3" = 1, "beta" = 0.1)
sigma_ts = c(1, 1, 1, 1)

system.time(fits_wth_ts_em <- em(dat_wth,
                                 log(mus_ts),
                                 sigma_ts,
                                 cost_lik,
                                 bic=TRUE,
                                 ndata=dat_wth[, .N],
                                 obj_args = list(
                                   par_names=par_names_ts,
                                   log_par=T
                                 ),
                                 perturb_start=0.1,
                                 parallel=T,
                                 verbose=T))

fits_wth_ts = setDT(get_ind_df(fits_wth_ts_em))
fits_wth_ts_predict = get_all_predictions(dat_wth, fits_wth_ts, par_names_ts, within=T, log_par=T, return_mean=T)
plot_wth(fits_wth_ts_predict, title="TS")
plot_params(fits_wth_ts, par_names_ts, exp_par=T, within=T)


#########################
### Fit adaptive time sensitivity model
#########################

par_names_ts2 = c("s_1", "s_2", "s_3", "S", "beta", "w_decay")
mus_ts2 = c("s_1" = 1, "s_2" = 1, "s_3" = 1, "S" = 1, "beta" = 0.1, "w_decay" = 1)
sigma_ts2 = c(1, 1, 1, 1, 1, 1)

system.time(fits_wth_ts2_em <- em(dat_wth,
                                  log(mus_ts2),
                                  sigma_ts2,
                                  cost_lik,
                                  bic=TRUE,
                                  ndata=dat_wth[, .N],
                                  obj_args = list(
                                    par_names=par_names_ts2,
                                    log_par=T
                                  ),
                                  perturb_start=0.1,
                                  parallel=T,
                                  verbose=T))

fits_wth_ts2 = setDT(get_ind_df(fits_wth_ts2_em))
fits_wth_ts2_predict = get_all_predictions(dat_wth, fits_wth_ts2, par_names_ts2, within=T, log_par=T, return_mean=T)
plot_wth(fits_wth_ts2_predict, title="Adaptive TS")
plot_params(fits_wth_ts2, par_names_ts2, exp_par=T, within=T)

#########################
### Fit adaptive opportunity cost/adaptive time sensitivity model
#########################

par_names_ts3 = c("s_1", "s_2", "s_3", "S", "beta", "w_decay", "gamma", "alpha")
mus_ts3 = c("s_1" = 1, "s_2" = 1, "s_3" = 1, "S" = 1, "beta" = 0.1, "w_decay" = 1, "gamma" = 0.55, "alpha" = 0.2)
sigma_ts3 = c(1, 1, 1, 1, 1, 1, 1, 1)

system.time(fits_wth_ts3_em <- em(dat_wth,
                                  log(mus_ts3),
                                  sigma_ts3,
                                  cost_lik,
                                  bic=TRUE,
                                  ndata=dat_wth[, .N],
                                  obj_args = list(
                                    par_names=par_names_ts3,
                                    log_par=T
                                  ),
                                  perturb_start=0.1,
                                  parallel=T,
                                  verbose=T))

fits_wth_ts3 = setDT(get_ind_df(fits_wth_ts3_em))
fits_wth_ts3_predict = get_all_predictions(dat_wth, fits_wth_ts3, par_names_ts3, within=T, log_par=T, return_mean=T)
plot_wth(fits_wth_ts3_predict, title="Adaptive OC + Adaptive TS")
plot_params(fits_wth_ts3, par_names_ts3, exp_par=T, within=T)


#########################
### Compare models
#########################

fits_wth_base[, model := "base"]
fits_wth_oc[, model := "oc"]
fits_wth_oc2[, model := "oc2"]
fits_wth_ts[, model := "ts"]
fits_wth_ts2[, model := "ts2"]
fits_all = rbind(fits_wth_base, fits_wth_oc, fits_wth_oc2, fits_wth_ts, fits_wth_ts2, fill=T)
fits_bic = fits_all[, .(bic = bic[1]), .(model)]
fits_bic[, model := factor(model, levels=c("base", "oc", "ts", "oc2", "ts2"))]

ggplot(fits_bic, aes(x=model, y=bic, group=0)) +
  geom_point() +
  geom_line() +
  theme_abw()

#########################
### Save workspace
#########################

save.image("models_gk/results/wth2.Rdata")
