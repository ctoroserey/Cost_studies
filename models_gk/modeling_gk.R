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

par_values = c(gamma = 1, beta = .1)
par_names = c("gamma", "beta")

system.time(fits_btw <- dat_btw[,
                                modelfitr::fit_model(
                                  cost_lik,
                                  log(par_values),
                                  hessian=T,
                                  package="optimx",
                                  method="nmkb",
                                  obj_args=list(
                                    dat=.SD,
                                    par_names=par_names,
                                    log_par=T
                                  ),
                                  aic=T,
                                  bic=T,
                                  n_obs=.N,
                                  return_df=T,
                                  verbose = T),
                                .(sub, cond),
                                .SDcols=names(dat_btw)
])

### get model predictions for each subject and condition

get_all_predictions(dat_btw, fits_btw, par_names, log_par=T)

mean_dat_btw = dat_btw[, .(Choice = mean(Choice),
                           predict_accept = mean(predict_accept)),
                       .(sub, cond, Handling, Offer)]

ggplot(mean_dat_btw, aes(x=Offer, y=Choice, color=cond, fill=cond)) +
  facet_wrap(~Handling) +
  stat_summary(fun=mean, geom="point") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=2) +
  stat_summary(aes(y=predict_accept), fun=mean, geom="line") +
  stat_summary(aes(y=predict_accept), fun.data=mean_se, geom="ribbon", color=NA, alpha=.2) +
  scale_y_continuous("p(accept)", breaks=seq(0, 1, .2)) +
  coord_cartesian(ylim=c(0,1)) +
  ggtitle("Base opportunity cost") +
  theme_abw() +
  theme(plot.title = element_text(hjust=0.5))


