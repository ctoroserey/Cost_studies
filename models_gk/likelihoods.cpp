#include <Rcpp.h>
using namespace Rcpp;


NumericVector vecpow(const NumericVector base, const NumericVector exp) {
  NumericVector out(base.size());
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow);
  return out;
}

// [[Rcpp::export]]
NumericVector cost_p_cpp(DataFrame dat,
                         NumericVector cost=NumericVector::create(0),
                         NumericVector s=NumericVector::create(1),
                         double beta=1,
                         double alpha=0,
                         double C = NA_REAL,
                         double S = NA_REAL,
                         double w_decay = 0,
                         double gamma = NA_REAL,
                         bool return_pacc=false) {
  
  NumericVector R_i = dat["Offer"],
                H_i = dat["HandleTime"],
                C_i = dat["Choice"],
                T_i = dat["TravelTime"];
  if (NumericVector::is_na(gamma)) {
    NumericVector rr = dat["RR"];
    gamma = mean(rr);
  }
  
  int n_trials = R_i.size();
  
  IntegerVector trial_vec = seq_len(n_trials);
  NumericVector w_t = pow((n_trials - as<NumericVector>(trial_vec)) / n_trials, w_decay);

  if (cost.length() == 1) {
    cost = rep(cost, n_trials);
  }
  if(!NumericVector::is_na(C)) {
    cost = cost * w_t + C * (1 - w_t);
  }

  if(s.length() == 1) {
    s = rep(s, n_trials);
  }
  if(!NumericVector::is_na(S)) {
    s = s * w_t + S * (1 - w_t);
  }
  
  NumericVector gvec(n_trials);
  if (alpha > 0) {
    gvec[0] = gamma;
    for (int i=0; i<n_trials-1; i++) {
      double tau = pow(H_i[i],s[i]) * C_i[i] + T_i[i];
      gvec[i+1] = (1 - pow(1 - alpha, tau)) * ((R_i[i] * C_i[i]) / tau) + pow(1 - alpha, tau) * gvec[i];
    }
  } else {
    gvec.fill(gamma);
  }

  NumericVector p_accept = 1 / (1 + exp(-1/beta * (R_i - (gvec + cost) * vecpow(H_i, s))));
    
  if (return_pacc)  {
    return p_accept;
  } else {
    NumericVector p_choice = ifelse(abs(C_i-1) < 1e-5, p_accept, 1-p_accept);
    return p_choice;
  }
  
}

// [[Rcpp::export]]
NumericVector cost_p_cpp_old(DataFrame dat,
                         double gamma=0,
                         double beta=1,
                         NumericVector s=NumericVector::create(1),
                         double alpha=0,
                         NumericVector big_s = 0,
                         double w_decay = 0,
                         bool fit_gamma=false,
                         bool return_pacc=false) {
  
  NumericVector R_i = dat["Offer"],
                H_i = dat["HandleTime"],
                C_i = dat["Choice"],
                T_i = dat["TravelTime"];
  
  if (!fit_gamma) {
    NumericVector rr = dat["RR"];
    gamma = rr[0] + gamma;
  }
  
  int n_trials = R_i.size();
  
  NumericVector gvec(n_trials),
    wvec = rep(pow(1, w_decay), n_trials);
  
  if (s.length() == 1)
    s = rep(s, n_trials);
  
  if (big_s.length() == 1)
    big_s = rep(big_s, n_trials);
  
  if (alpha > 0) {
    gvec[0] = gamma;
    for (int i=0; i<n_trials-1; i++) {
      double tau = pow(H_i[i],s[i]) * C_i[i] + T_i[i];
      gvec[i+1] = (1 - pow(1 - alpha, tau)) * ((R_i[i] * C_i[i]) / tau) + (pow(1 - alpha, tau) * gvec[i]);
      wvec[i+1] = pow((i+2) / n_trials, w_decay);
    }
  } else {
    gvec.fill(gamma);
  }
  
  NumericVector p_accept = 1 / (1 + exp(-1/beta * (R_i - gvec * vecpow(H_i, (s*(wvec) + big_s*(1-wvec))))));
  if (return_pacc)  {
    return p_accept;
  } else {
    NumericVector p_choice = ifelse(abs(C_i-1) < 1e-5, p_accept, 1-p_accept);
    return p_choice;
  }
  
}



