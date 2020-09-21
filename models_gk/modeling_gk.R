#########################
### Setup
#########################

# set wd
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

library(data.table)
library(foreach)

reload = function(){
  library(modelfitr)
  detach("package:modelfitr", unload=T)
  library(modelfitr)
}
reload()

source('models_gk/likelihoods.R')

#########################
### Load Data
#########################

dir_btw = 'Cost2/data'
f_btw = list.files(dir_btw, pattern=".csv", full.names=T)
dat_btw = foreach(f = f_btw, .combine=rbind) %do% {
  info = strsplit(basename(f), "_")[[1]]
  data.table(sub=info[1], cond=info[2], date=info[3], fread(f))
}
dat_btw[, Travel := 16-Handling]

# dir_wth = 'Cost3/data'
# f_wth = list.files(dir_wth, pattern="main_log*.csv", full.names=T)
# dat_wth = foreach(f = f_wth, .combine=rbind) %do% {
#   info = strsplit(basename(f), "_")[[1]]
#   data.table(sub=info[1], cond=info[2], date=info[3], fread(f))
# }

#########################
### Fit time sensitivity model
#########################

par_values = c(gamma = 1, beta = .1)

reload()
system.time(fits_btw <- dat_btw[,
                                modelfitr::fit_model(
                                  get_lik,
                                  log(par_values),
                                  hessian=T,
                                  method="nmkb",
                                  obj_args=list(
                                    lfun=time_sensitivity_lik,
                                    R_i=Offer,
                                    H_i=Handling,
                                    Choice=Choice,
                                    log_par=T
                                  ),
                                  aic=T,
                                  bic=T,
                                  n_obs=.N,
                                  return_df=T),
                                .(sub, cond)
])
fits_btw[, gamma := exp(gamma)]
fits_btw[, beta := exp(beta)]

fits_btw_melt = melt(fits_btw, measure.vars=c("gamma", "beta"), variable.name="par", value.name="val")
ggplot(fits_btw_melt, aes(x=cond, y=val)) +
  facet_wrap(~par, scales="free") +
  geom_violin()

dat_btw[, predict_choice := get_lik(fits_btw[sub==sub[1] & cond==cond[1], c(gamma, beta)],
                                    return_liks=T,
                                    time_sensitivity_lik,
                                    R_i=Offer,
                                    H_i=Handling,
                                    Choice=Choice),
        .(sub, cond)]

dat_btw[, predict_accept := ifelse(Choice==1, predict_choice, 1-predict_choice)]
mean_dat_btw = dat_btw[, .(Choice = mean(Choice),
                           predict_accept = mean(predict_accept)),
                       .(sub, cond, Handling, Offer)]

ggplot(mean_dat_btw, aes(x=Offer, y=Choice, color=cond, linetype=)) +
  facet_wrap(~Handling) +
  stat_summary(fun=mean, geom="point") +
  stat_summary(fun.data=mean_se, geom="errorbar") +
  stat_summary(aes(y=predict_accept), fun=mean, geom="line")
