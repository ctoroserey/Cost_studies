### define likelihood functions

# Rcpp::sourceCpp(file.path(dirname(sys.frame(1)$ofile), "likelihoods.cpp"))

cost_p <- function(dat,
                   cost = 0,
                   s = 1,
                   beta = 1,
                   alpha = 0,
                   C = NA,
                   S = NA,
                   w_decay = 0,
                   gamma = NA,
                   return_pacc = F) {
  
  n_trials = dat[, .N]
  R_i = dat[, Offer]
  H_i = dat[, HandleTime]
  C_i = dat[, Choice]
  T_i = dat[, TravelTime]
  if(is.na(gamma)) {
    gamma = mean(dat[, RR])
  }
  
  if (alpha > 0) {
    for (i in 1:((length(R_i)) - 1)) {
      tau <- H_i[i]^s[i] * C_i[i] + T_i[i]
      gamma[i + 1] <- (1 - (1 - alpha)^tau) * ((R_i[i] * C_i[i]) / tau) + ((1 - alpha)^tau * gamma[i])
    }
  }
  
  w_t = ((n_trials - 1:n_trials) / n_trials) ^ w_decay
  
  if(!is.na(C)) {
    if (length(cost) == 1) cost = rep(cost, n_trials)
    cost = cost * w_t + C * (1 - w_t)
  }
  
  if(!is.na(S)) {
    s = s * w_t + S * (1 - w_t)
  }
  
  p_accept <- 1 / (1 + exp(-1/beta * (R_i - (gamma + cost) * H_i^s)))
  
  if (return_pacc)
    return(p_accept)
  else
    return(ifelse(C_i, p_accept, 1 - p_accept))
}

cost_p_old <- function(dat,
                       gamma = 0,
                       beta = 1,
                       s = 1,
                       alpha = 0,
                       big_s = 0,
                       w_decay = 0,
                       fit_gamma = F,
                       return_pacc = F) {
  
  R_i = dat[, Offer]
  H_i = dat[, HandleTime]
  C_i = dat[, Choice]
  T_i = dat[, TravelTime]
  w_t = 1 ^ w_decay
  if (!fit_gamma){
    gamma = dat[1, RR] + gamma
  }
  
  if (!missing(alpha)) {
    
    if (length(s) == 1) {
      s = rep(s, length(R_i))
    }
    
    if (alpha > 0) {
      for (i in 1:((length(R_i)) - 1)) {
        tau <- H_i[i]^s[i] * C_i[i] + T_i[i]
        gamma[i + 1] <- (1 - (1 - alpha)^tau) * ((R_i[i] * C_i[i]) / tau) + ((1 - alpha)^tau * gamma[i])
        w_t[i+1] = ((i+1) / length(R_i)) ^ w_decay
      }
    }
  }
  
  p_accept <- 1 / (1 + exp(-1/beta * (R_i - gamma * H_i^(s*(w_t) + big_s*(1-w_t)))))
  
  if (return_pacc)
    return(p_accept)
  else
    return(ifelse(C_i, p_accept, 1 - p_accept))
}


log_transform_pars = function(par) {
  par[!grepl("cost", names(par)) & names(par) != "C"] <- exp(par[!grepl("cost", names(par)) & names(par) != "B"])
  par
}

cost_lik <- function(par, dat, par_names = NULL, fixed_pars=c(), log_par = F, min_p = 1e-5, return_pacc = F, debug = F, ...) {
  
  if (debug){
    for(i in 1:length(par_names)) {
      cat(par_names[i], " = ", round(par[i], 3), sep = "")
      cat(", ", sep="")
    }
  }
  
  if (!is.null(par_names)) names(par) <- par_names
  if (log_par) par = log_transform_pars(par)
  par = c(par, fixed_pars)
  
  pass_through = list(...)
  
  # get gamma by condition
  
  cost_par_index = grep("^cost_", names(par))
  if (length(cost_par_index) > 1) {
    cost_ind = as.numeric(sapply(names(par)[cost_par_index], function(x) substring(x, nchar(x), nchar(x))))
    cost_cond = numeric(max(cost_ind))
    for(i in 1:length(cost_par_index)) {
      cost_cond[cost_ind[i]] = par[cost_par_index[i]]
    }
    cost_vec = dat[, cost_cond[condCode]]
    pass_through$cost = cost_vec
  }
  
  # get s by condition
  
  s_par_index = grep("^s_", names(par))
  if (length(s_par_index) > 1) {
    s_ind = as.numeric(sapply(names(par)[s_par_index], function(x) substring(x, nchar(x), nchar(x))))
    s_cond = numeric(max(s_ind))
    for(i in 1:length(s_par_index)) {
      s_cond[s_ind[i]] = par[s_par_index[i]]
    }
    s_vec = dat[, s_cond[condCode]]
    pass_through$s = s_vec
  }
  
  par_to_remove = c(cost_par_index, s_par_index)
  if (length(par_to_remove) > 0) par = par[-par_to_remove]
  
  liks <- do.call(cost_p, c(list(dat=dat), as.list(par), pass_through,return_pacc=return_pacc))
  
  if (return_pacc){
    return(liks)
  }else{
    nll = -sum(log(ifelse(liks < min_p | is.na(liks), min_p, liks)))
    if (debug){
      cat(" // nll = ", round(nll, 3), "\n", sep="")
    }
    return(nll)
  }
}

get_all_predictions = function(all_data, fits, par_names, within=F, return_mean=FALSE, ...) {
  
  for (subj in unique(all_data[, Subject])) {
    
    if(within) {
      
      all_data[Subject==subj, predict_accept := cost_lik(as.numeric.dt(fits[Subject==subj, par_names, with=F]),
                                                         dat=.SD,
                                                         return_pacc=T,
                                                         ...)]
    } else {
      
      for (cond in unique(dat_btw[, cond])) {
        
        all_data[Subject==subj & cond ==cond, predict_accept := cost_lik(as.numeric.dt(fits[Subject==subj & cond==cond, par_names, with=F]),
                                                                         dat=.SD,
                                                                         return_pacc=T,
                                                                         ...)]
        
        
      }
    }
  }
  
  if (return_mean) {
    if (within) {
      all_data[Block %in% c(1,2,5,6), .(Choice = mean(Choice),
                                        predict_accept = mean(predict_accept)),
               .(sub, cond, half, first, Offer, BlockType)]
    } else {
      all_data[, .(Choice = mean(Choice),
                   predict_accept = mean(predict_accept)),
               .(Subject, cond, Handling, Offer)]
    }
  }
  
}

get_exp_pars = function(fits, par_names) {
  fits_exp = copy(fits)
  for(p in par_names){
    if (!grepl("gamma", p) & p != "big_g"){
      fits_exp[[p]] = exp(fits[[p]])
    }
  }
  
  fits_exp
}

plot_btw = function(plot_data, title="") {
  cols = c("#78AB05","#D9541A","deepskyblue4", "darkgoldenrod2")
  
  pl = ggplot(plot_data, aes(x=Offer, y=Choice, color=cond, fill=cond)) +
    facet_wrap(~Handling) +
    stat_summary(fun=mean, geom="point") +
    stat_summary(fun.data=mean_se, geom="errorbar", width=2)
  
  if ("predict_accept" %in% names(plot_data)) {
    pl = pl + 
      stat_summary(aes(y=predict_accept), fun=mean, geom="line") +
      stat_summary(aes(y=predict_accept), fun.data=mean_se, geom="ribbon", color=NA, alpha=.2)
  } else {
    pl = pl + stat_summary(fun=mean, geom="line")
  }
  
  pl = pl + 
    scale_y_continuous("p(accept)", breaks=seq(0, 1, .2)) +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    coord_cartesian(ylim=c(0,1)) +
    ggtitle(title) +
    theme_abw() +
    theme(plot.title = element_text(hjust=0.5))
  
  pl
}

plot_wth = function(plot_data, title="") {
  cols = c("#D9541A", "#78AB05", "dodgerblue4", "deepskyblue3")
  
  pdat = copy(plot_data)
  pdat[, cond_block := factor(interaction(cond, BlockType), levels=c("cogTask.0", "phys.1", "wait.0", "phys.1"), labels=c("Cognitive", "Physical", "Wait-C", "Wait-P"))]

  pl = ggplot(pdat, aes(x=Offer, y=Choice, color=cond_block, fill=cond_block)) +
    facet_grid(first~half) +
    stat_summary(fun=mean, geom="point") +
    stat_summary(fun.data=mean_se, geom="errorbar", width=2)
  
  if ("predict_accept" %in% names(plot_data)) {
    pl = pl + stat_summary(aes(y=predict_accept), fun=mean, geom="line") +
      stat_summary(aes(y=predict_accept), fun.data=mean_se, geom="ribbon", color=NA, alpha=.2)
  } else {
    pl = pl + stat_summary(fun=mean, geom="line")
  }
  
  pl = pl + scale_y_continuous("p(accept)", breaks=seq(0, 1, .2)) +
    coord_cartesian(ylim=c(0,1)) +
    scale_color_manual(NULL, values=cols) +#, labels=c("wait-cog", "cog-cog", "wait-phys", "phys-phys")) +
    scale_fill_manual(NULL, values=cols) +#, labels=c("wait-cog", "cog-cog", "wait-phys", "phys-phys")) +
    ggtitle(title) +
    theme_abw() +
    theme(plot.title = element_text(hjust=0.5))
  
  pl
}


plot_params = function(fits, par_names, exp_par=F, within=F) {
  
  if(exp_par) {
    fits = get_exp_pars(fits, par_names)
  }
  
  mfits = melt(fits, measure.vars=par_names, variable.name="par", value.name="val")
  
  if (within) {
    ggplot(mfits, aes(x=par, y=val)) +
      geom_violin() +
      xlab(NULL) +
      theme_abw()
  } else {
    ggplot(mfits, aes(x=cond, y=val)) +
      facet_wrap(~par, scales="free") +
      geom_violin() +
      theme_abw()
  }
  
}