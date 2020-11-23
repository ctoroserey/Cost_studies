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

c(dat_btw, dat_wth) %<-% load_data()

#########################
### Fit opportunity cost model
#########################

# par_values = c(gamma = 1, beta = .1)
# par_names = c("gamma", "beta")
# 
# system.time(fits_btw <- dat_btw[,
#                                 modelfitr::fit_model(
#                                   cost_lik,
#                                   log(par_values),
#                                   hessian=T,
#                                   package="optimx",
#                                   method="nmkb",
#                                   obj_args=list(
#                                     dat=.SD,
#                                     par_names=par_names,
#                                     log_par=T
#                                   ),
#                                   aic=T,
#                                   bic=T,
#                                   n_obs=.N,
#                                   return_df=T),
#                                 .(sub, cond),
#                                 .SDcols=names(dat_btw)
# ])
# 
# for (i in par_names) fits_btw[[i]] = exp(fits_btw[[i]])
# 
# fits_btw_melt = melt(fits_btw, measure.vars=c("gamma", "beta"), variable.name="par", value.name="val")
# ggplot(fits_btw_melt, aes(x=cond, y=val)) +
#   facet_wrap(~par, scales="free") +
#   geom_violin() +
#   theme_abw()
# 
# 
# ### get model predictions for each subject and condition
# 
# for (subj in unique(dat_btw[, sub])) {
#   for (cond in unique(dat_btw[, cond])) {
#     dat_btw[sub==subj & cond ==cond, predict_accept := cost_lik(fits_btw[sub==subj & cond==cond, c(gamma, beta)],
#                                                                 dat=.SD,
#                                                                 return_pacc=T)]
#     
#   }
# }
# 
# 
# mean_dat_btw = dat_btw[, .(Choice = mean(Choice),
#                            predict_accept = mean(predict_accept)),
#                        .(sub, cond, Handling, Offer)]
# 
# ggplot(mean_dat_btw, aes(x=Offer, y=Choice, color=cond, fill=cond)) +
#   facet_wrap(~Handling) +
#   stat_summary(fun=mean, geom="point") +
#   stat_summary(fun.data=mean_se, geom="errorbar", width=2) +
#   stat_summary(aes(y=predict_accept), fun=mean, geom="line") +
#   stat_summary(aes(y=predict_accept), fun.data=mean_se, geom="ribbon", color=NA, alpha=.2) +
#   scale_y_continuous("p(accept)", breaks=seq(0, 1, .2)) +
#   coord_cartesian(ylim=c(0,1)) +
#   theme_abw()
# 

#########################
### Fit opportunity cost model using EM
#########################

dat_btw[, Subject := sub]
subs = dat_btw[, unique(Subject)]

fit_btw_oc = dat_btw[, GaussExMax::em(.SD,
                                      .SD[, unique(Subject)],
                                      log(c(gamma=1, beta=0.1)),
                                      c(5, 5),
                                      cost_lik,
                                      emtol=.001,
                                      par_names=c("gamma", "beta"),
                                      log_par=T,
                                      return_ind_df=T), .(cond)]
fit_btw_oc[, model := "oc"]


fit_btw_ts = dat_btw[, GaussExMax::em(.SD,
                                      .SD[, unique(Subject)],
                                      log(c(gamma=1, beta=0.1, s=1)),
                                      c(5, 5, 5),
                                      cost_lik,
                                      emtol=.001,
                                      par_names=c("gamma", "beta", "s"),
                                      log_par=T,
                                      return_ind_df=T), .(cond)]
fit_btw_ts[, model := "ts"]

fit_btw = rbind(fit_btw, fit_btw_ts, fill=TRUE)

fit_btw_melt = melt(fit_btw,
                    measure.vars=c("gamma", "beta", "s"),
                    variable.name="par",
                    value.name="val")
fit_btw_melt[, val := exp(val)]

ggplot(fit_btw_melt, aes(x=cond, y=val)) +
  facet_grid(par ~ model, scales="free") +
  geom_violin() +
  scale_x_discrete(NULL, breaks=c("cogTask", "pheasy", "phys", "wait"), labels=c("Cog", "PhEasy", "Phys", "Wait")) +
  ylab("value") +
  theme_abw()

### get model predictions for each subject and condition

all_par_names = c("gamma", "beta", "s")
for(mod in unique(fit_btw_melt[, model])) {
  colname = paste0("predict_", mod)
  for (subj in unique(dat_btw[, sub])) {
    for (cond in unique(dat_btw[, cond])) {
      sub_ind = which(dat_btw[, sub==subj & cond==cond])
      sub_pars = as.numeric.dt(fit_btw[Subject==subj & cond==cond & model==mod, ..all_par_names])
      sub_pars = sub_pars[!is.na(sub_pars)]
      set(dat_btw, i=sub_ind, j=colname, value=cost_lik(sub_pars, dat=dat_btw[sub_ind], log_par=T, return_pacc=T))
    }
  }
}

mean_dat_btw = dat_btw[, .(Choice = mean(Choice),
                           predict_oc = mean(predict_oc),
                           predict_ts = mean(predict_ts)),
                       .(sub, cond, Handling, Offer)]

ggplot(mean_dat_btw, aes(x=Offer, y=Choice, color=cond, fill=cond)) +
  facet_wrap(~Handling) +
  stat_summary(fun=mean, geom="point") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=2) +
  stat_summary(aes(y=predict_oc), fun=mean, geom="line") +
  stat_summary(aes(y=predict_oc), fun.data=mean_se, geom="ribbon", color=NA, alpha=.2) +
  scale_y_continuous("p(accept)", breaks=seq(0, 1, .2)) +
  coord_cartesian(ylim=c(0,1)) +
  ggtitle("base oc model") +
  theme_abw() + theme(plot.title = element_text(hjust=0.5))

ggplot(mean_dat_btw, aes(x=Offer, y=Choice, color=cond, fill=cond)) +
  facet_wrap(~Handling) +
  stat_summary(fun=mean, geom="point") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=2) +
  stat_summary(aes(y=predict_ts), fun=mean, geom="line") +
  stat_summary(aes(y=predict_ts), fun.data=mean_se, geom="ribbon", color=NA, alpha=.2) +
  scale_y_continuous("p(accept)", breaks=seq(0, 1, .2)) +
  coord_cartesian(ylim=c(0,1)) +
  ggtitle("time sensitivity model") +
  theme_abw() + theme(plot.title = element_text(hjust=0.5))
