### define likelihood functions

time_sensitivity_lik = function(R_i, H_i, Choice, gamma, beta) {
  p_accept = 1 / (1 + exp(-beta * (R_i - gamma * H_i)))
  p_choice = ifelse(as.logical(Choice), p_accept, 1-p_accept)
  p_choice
}

time_sensitivity_rate_lik = function(R_i, H_i, T_i, Choice, gamma, beta, Block=NULL) {
  accept_rate = R_i / (gamma * H_i + T_i)
  
  if(is.null(Block)){
    global_rate = sum(R_i * Choice) / sum(gamma * H_i * Choice + T_i)
  } else {
    blocks = unique(Block)
    global_rate = numeric(length(R_i))
    for (i in unique(block)) {
      global_rate[Block==i] = sum(R_i[Block==i] * Choice[Block==i]) / sum(gamma * H_i[Block==i] * Choice[Block==i] + T_i[Block==i])
    }
  }
  
  p_accept = 1 / (1 + exp(-beta * (accept_rate - global_rate)))
  p_choice = ifelse(as.logical(Choice), p_accept, 1-p_accept)
  p_choice
}

get_lik = function(par, lfun, par_names=NULL, log_par=F, min_p=1e-5, return_liks=F, ...) {
  if(log_par) par=exp(par)
  if(!is.null(par_names)) names(par) = par_names
  liks = do.call(lfun, c(as.list(par), list(...)))
  liks[liks < min_p] = min_p
  if (return_liks)
    return(liks)
  else
    return(-sum(log(liks)))
}
