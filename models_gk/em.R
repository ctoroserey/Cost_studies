#' em
#'
#' Run expectation maximization algorithm to fit parameters to each subject of an experiment
#'
#' @param dat data.frame; data
#' @param subs string vector; subject ids
#' @param mus vector; starting mean parameters
#' @param sigma vector; starting sd in parameter
#' @param likfun function; the objective function
#'
#' @export
em = function(dat, subs, mus, sigma, likfun, emtol = 1e-3, maxiter = 100, lower = -Inf, upper = Inf, parallel = F, print = T, ...){
  
  nparam = length(mus)
  
  ### set up design matrix and sigma matrix ###
  dm = designmatrix(nparam, length(subs))
  sigma = Matrix::Diagonal(length(sigma), sigma)
  
  ### first iteration ###
  oldparams = c(mus, diag(sigma))
  fit = estep(dat, subs, mus, sigma, likfun, lower = lower, upper = upper, parallel = parallel, ...)
  par_list = mstep(fit$x, dm, fit$h, sigma)
  mus = par_list$mus
  sigma = par_list$sigma
  newparams = c(mus, diag(sigma))
  
  
  ### repeat until convergence ###
  iter = 1
  while(max(abs((newparams-oldparams)/oldparams)) > emtol & iter < maxiter){
    iter = iter + 1
    oldparams = newparams
    
    fit = estep(dat, subs, mus, sigma, likfun, startx = fit$x, lower = lower, upper = upper, parallel = parallel, ...)
    par_list = mstep(fit$x, dm, fit$h, sigma)
    mus = par_list$mus
    sigma = par_list$sigma
    newparams = c(mus, diag(sigma))
    
    if(print){
      cat("\n")
      cat("iter:", iter, "\n")
      cat("mus:", round(mus, 5), "\n")
      cat("sigma:", round(diag(sigma), 5), "\n")
      cat("change:", round(abs((newparams-oldparams)/oldparams), 5), "\n")
      cat("max:", round(max(abs((newparams-oldparams)/oldparams)), 5), "\n")
    }
  }
  
  list(mus = mus, sigma = sigma, x = fit$x, l = fit$l, h = fit$h)
}

estep = function(dat, subs, mus, sigma, likfun, startx = NULL, lower = -Inf, upper = Inf, parallel = F, ...){
  nsub = length(subs)
  nparam = length(mus)
  if(is.null(startx)) startx = matrix(rep(mus, nsub), ncol = nsub)
  
  h = array(0, c(nparam, nparam, nsub))
  l = numeric(nsub)
  x = matrix(0, nrow = nparam, ncol = nsub)
  
  if(!parallel){
    for(i in 1:nsub){
      cat(i,"..",sep="")
      sub_fit = opt_sub(gaussian_prior, startx[,i], lower = lower, upper = upper, mus = mus, sigma = sigma, dat = dat[dat$Subject == subs[i],], likfun = likfun, ...)
      x[,i] = sub_fit$par
      l[i] = sub_fit$lik
      h[,,i] = sub_fit$hess
    }
  }else{
    dat_list = split(dat, dat$Subject)
    startx_list = simplify2array(apply(startx, 2, list))
    opt_sub_parallel = function(d_sub, stx, ...){
      opt_sub(gaussian_prior, stx, lower = lower, upper = upper, dat = d_sub, ...)
    }
    fit_list = parallel::mcmapply(opt_sub_parallel, d_sub = dat_list, stx = startx_list,
                                  MoreArgs = c(list(mus = mus), sigma = sigma, likfun = likfun, list(...)),
                                  mc.preschedule = F, mc.cores = parallel::detectCores() - 1)
    # fit_list = mapply(opt_sub_parallel, d_sub = dat_list, stx = startx_list,
    #                   MoreArgs = c(list(mus = mus), sigma = sigma, likfun = overharvest_objective, list(...)))
    
    if(length(fit_list) == 0) browser()
    x = simplify2array(fit_list[1,])
    l = as.numeric(fit_list[2,])
    h = simplify2array(fit_list[3,])
  }
  
  list(x = x, l = l, h = h)
}

mstep = function(x, dm, h, sigma, full = F){
  nsub = ncol(x)
  isigma = solve(sigma)
  
  ### find betas ###
  m1 = apply(dm, 3, function(x) t(x) %*% isigma %*% x)
  if(class(isigma) != "matrix"){
    mean_m1 = numeric(nrow(m1[[1]]))
    for(i in 1:nsub){
      mean_m1 = mean_m1 + diag(m1[[i]])
    }
    mean_m1 = mean_m1 / nsub
    mean_m1 = Matrix::Diagonal(length(mean_m1), mean_m1)
  }else{
    m1 = array(m1, c(nrow(sigma), nrow(sigma), nsub))
    mean_m1 = apply(m1, 1:2, mean)
  }
  
  mean_m2 = numeric(nrow(x))
  for(i in 1:nsub){
    mean_m2 = mean_m2 + t(dm[,,i]) %*% isigma %*% x[,i]
  }
  mean_m2 = mean_m2 / nsub
  
  mus = as.numeric(solve(mean_m1) %*% mean_m2)
  
  ### find sigma ###
  sigma = (x-mus) %*% t(x-mus) / nsub + apply(h, 1:2, mean, na.rm = T)
  
  if(!full){
    sigma = Matrix::Diagonal(nrow(sigma), diag(sigma))
  }
  
  list(mus = mus, sigma = sigma)
}

opt_sub = function(obj, start, lower = -Inf, upper = Inf, method = "Nelder-Mead", optControl = list(), ...){
  
  it = 1
  fit = optimx::optimx(start, obj, lower = lower, upper = upper, method = method, hessian = F, control = c(optControl, kkt = F), ...)
  
  while(fit$convcode != 0 & it < 5){
    it = it + 1
    start = rnorm(length(start), start, abs(start))
    new_fit = optimx::optimx(start, obj, lower = lower, upper = upper, method = method, hessian = F, control = c(optControl, kkt = F), ...)
    if(new_fit$value[1] < fit$value[1]) fit = new_fit
  }
  
  par = as.numeric(fit[1:length(start)])
  val = fit$value[1]
  hess = solve(pracma::hessian(obj, par, ...))
  if(det(hess) < 0) hess = solve(numDeriv::hessian(obj, par, ...))
  #if(det(hess) < 0) browser()
  if(det(hess) < 0){
    cat("\n negative hessian! replacing with NA \n")
    hess = matrix(NA, length(start), length(start))
  }
  
  return(list(par = par, lik = val, hess = hess))
}

gaussian_prior = function(param_values, mus, sigma, dat, likfun, ...){
  d = length(param_values)
  lp = -d/2 * log(2*pi) - 1/2 * log(det(sigma)) - 1/2 * t(param_values - mus) %*% solve(sigma) %*% (param_values - mus)
  nll = likfun(param_values, dat, ...)
  map = nll - lp[1]
  if(is.infinite(map) | is.na(map)) map = 1e10
  if(is.nan(map)) browser()
  map
}

designmatrix = function(npar, nsub) array(diag(npar), c(npar, npar, nsub))

ibic = function(x, l, h, mus, sigma, ndata){
  nHypPar = length(c(mus, diag(sigma)))
  lml(x, l, h) + nHypPar/2 * log(ndata)
}

lml = function(x, l, h){
  # laplace approximation to the log marginal likelihood
  nparam = nrow(x)
  nsub = ncol(x)
  
  log_det_hess = apply(h, 3, function(x) log(det(x)))
  -nparam/2 * log(2*pi) * nsub + sum(l) - sum(log_det_hess) / 2
}

fit_all = function(dat, subs, obj, start, lower = -Inf, upper = Inf, parallel = F, ...){
  nparam = length(start)
  nsub = length(subs)
  h = array(0, c(nparam, nparam, nsub))
  l = numeric(nsub)
  x = matrix(0, nrow = nparam, ncol = nsub)
  
  if(!parallel){
    for(i in 1:nsub){
      cat(i,"..",sep="")
      sub_fit = opt_sub(obj, start, lower = lower, upper = upper, dat = dat[dat$Subject == subs[i],], ...)
      x[,i] = sub_fit$par
      l[i] = sub_fit$lik
      h[,,i] = sub_fit$hess
    }
  }else{
    dat_list = split(dat, dat$Subject)
    opt_sub_parallel = function(d_sub, stx, ...){
      opt_sub(gaussian_prior, stx, lower = lower, upper = upper, dat = d_sub, ...)
    }
    fit_list = parallel::mclapply(dat_list, opt_sub_parallel,
                                  MoreArgs = c(stx = start, list(...)),
                                  mc.cores = parallel::detectCores())
    
    x = simplify2array(fit_list[1,])
    l = as.numeric(fit_list[2,])
    h = simplify2array(fit_list[3,])
  }
  
  list(x = x, l = l, h = h)
  
}
