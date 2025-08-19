setwd('/ihome/jbwang/chl471/bartpp')
library(gtools)
library(tidyverse)
library(SoftBart)
library(progress)
library(moderate.mediation)
source('lib/bart_mediate2.R')
source('lib/get_clever_cov.R')
source('lib/get_ps.R')

appdat = readRDS("gpath.RDS")
flip_p = mean(appdat$t)

name = c("white", "male", "educ", "age_at_visit")
result = modmed(data = appdat, 
                treatment = "t", 
                mediator = "m", 
                outcome = "y", 
                moderators.cont = name, 
                m.model = list(intercept = name, treatment = name), 
                y.model = list(intercept = name, treatment = NULL, 
                               mediator = NULL, tm = NULL),
                m.scale = "continuous", 
                y.scale = "continuous", 
                comp.mod.cont.values = NULL, 
                ref.mod.cont.values = NULL, 
                method = "mc",
                nmc = 1000,
                nboot = 1000,
                conf.level = 0.95)
result$m.model
result$y.model

formula_ps = t ~ white + male + educ + age_at_visit
formula_m = m ~ white + male + educ + age_at_visit + t
formula_y = y ~ white + male + educ + age_at_visit + m + t

n_burn_in = 5000; n_after_burn_in = 2500;

### Clever covariates and ps score
clever_cov = get_clever_cov(appdat, appdat, formula_m, 'm', 'y', 't')
m0_hat = clever_cov$m0_hat
m1_hat = clever_cov$m1_hat
pi_hat = get_ps(appdat, appdat, formula_ps)

### Fit
vcbcmf = bart_mediate1(
  data_train = appdat, 
  data_test = appdat, 
  formula_m = formula_m, 
  formula_y = formula_y, 
  pi_hat_train = pi_hat, 
  pi_hat_test = pi_hat,
  m0_hat_train = m0_hat, 
  m0_hat_test = m0_hat, 
  m1_hat_train = m1_hat, 
  m1_hat_test = m1_hat,
  mediator_name = 'm', 
  outcome_name = 'y', 
  treat_name = 't',
  n_iter = n_burn_in + n_after_burn_in, 
  burnin = n_after_burn_in, 
  num_tree = 20, gamma = 0.5,
  chain = 1
)

gen = function(n = 1000, case = 1){
  
  sigma.m = 1; sigma.y = 1; 
  e.m = rnorm(n, 0, sigma.m)
  e.y = rnorm(n, 0, sigma.y)
  
  ### Case 1
  if(case == 1){
    
    c = c(0.44, 0.06, 0.54, 0.42, 0.28, 0.41, 0.63, 0.78, 0.78, 0.85, 0.72, 0.29)
    x1 = rnorm(n); x2 = rnorm(n); 
    ps = plogis(x1 + x2)
    t = rbinom(n, 1, ps)
    
    b0.m = c[1] + c[2] * x1
    b1.m = c[4] + c[3] * x2 
    
    b0.y = c[5] + c[6] * x1
    b1.y = c[8] + c[7] * x2
    b2.y = c[9] + c[10] * x1
    b3.y = c[11] + c[12] * x2 
    
    m = b0.m + b1.m * t + e.m
    y = b0.y + b1.y * t + b2.y * m + b3.y * t * m + e.y
    
    m0 = b0.m
    m1 = b0.m + b1.m
    y1m1 = b0.y + b1.y + b2.y * m1 + b3.y * m1
    y1m0 = b0.y + b1.y + b2.y * m0 + b3.y * m0
    y0m0 = b0.y + b2.y * m0
    
    true_TE = y1m1 - y0m0
    true_NIE = y1m1 - y1m0
    true_NDE = y1m0 - y0m0
    
    data = cbind.data.frame(x1 = x1, x2 = x2, t = t, m = m, y = y, true_NIE = true_NIE, true_NDE = true_NDE)
    
  } else if(case == 2){
    
    c = c(0.44, 0.06, 0.54, 0.42, 0.28, 0.41, 0.63, 0.78, 0.78, 0.85)
    x1 = rnorm(n); x2 = rnorm(n); 
    ps = plogis(x1 + x2)
    t = rbinom(n, 1, ps)
    
    b0.m = c[1] + c[2] * x1^2
    b1.m = c[3] + c[4] * x2 
    
    b0.y = c[5] + c[6] * abs(x1)
    b1.y = c[7] + c[8] * x2
    b2.y = c[9] + c[10] * sin(x1)
    
    m = b0.m + b1.m * t + e.m
    y = b0.y + b1.y * t + b2.y * m + e.y
    
    m0 = b0.m
    m1 = b0.m + b1.m
    y1m1 = b0.y + b1.y + b2.y * m1
    y1m0 = b0.y + b1.y + b2.y * m0
    y0m0 = b0.y + b2.y * m0
    
    true_TE = y1m1 - y0m0
    true_NIE = y1m1 - y1m0
    true_NDE = y1m0 - y0m0
    
    data = cbind.data.frame(x1 = x1, x2 = x2, t = t, m = m, y = y, true_NIE = true_NIE, true_NDE = true_NDE)
    
  } else if(case == 3){
    
    x1 = runif(n); x2 = runif(n); x3 = runif(n); x4 = runif(n); x5 = runif(n);
    ps = plogis(x1 + x2 * x3 - x4 * x5)
    t = rbinom(n, 1, ps)
    
    b0.m = x1^2 + 2 * x2 * x5
    b1.m = 2
    
    b0.y = x5
    b1.y = 2 * x4
    b2.y = 4 * (x3 - 0.5)^2
    
    m = b0.m + b1.m * t + e.m
    y = b0.y + b1.y * t + b2.y * m + e.y
    
    m0 = b0.m
    m1 = b0.m + b1.m
    y1m1 = b0.y + b1.y + b2.y * m1
    y1m0 = b0.y + b1.y + b2.y * m0
    y0m0 = b0.y + b2.y * m0
    
    true_TE = y1m1 - y0m0
    true_NIE = y1m1 - y1m0
    true_NDE = y1m0 - y0m0
    
    data = cbind.data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5,
                            t = t, m = m, y = y, true_NIE = true_NIE, true_NDE = true_NDE)
    
  } else if(case == 4){
    
    x1 = runif(n); x2 = runif(n); x3 = runif(n); x4 = runif(n); x5 = runif(n);
    ps = plogis(x2 * x3 - x4 * x5)
    t = rbinom(n, 1, ps)
    
    b0.m = 1 / (x1 + 1) + x2 
    b1.m = 2 * abs(x3 - 0.5)
    
    b0.y = x5
    b1.y = 2 * x4
    b2.y = 4 * (x5 - 0.5)^2
    b3.y = 2 * sin(pi * x1 * x2)
    
    m = b0.m + b1.m * t + e.m
    y = b0.y + b1.y * t + b2.y * m^2 + b3.y * t * m + e.y
    
    m0 = b0.m 
    m1 = b0.m + b1.m 
    y1m1 = b0.y + b1.y + b2.y * m1^2 + b3.y * m1
    y1m0 = b0.y + b1.y + b2.y * m0^2 + b3.y * m0
    y0m0 = b0.y + b2.y * m0^2
    
    true_TE = y1m1 - y0m0
    true_NIE = y1m1 - y1m0
    true_NDE = y1m0 - y0m0
    
    data = cbind.data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5,
                            t = t, m = m, y = y, true_NIE = true_NIE, true_NDE = true_NDE)
    
  } else if(case == 5){ ### Use MM, original m and y, generate some for the true effects

    
    m_c = result$m.model$Estimate
    y_c = result$y.model$Estimate
    
    subdf = appdat[sample(nrow(appdat), n, replace = TRUE), , drop = FALSE]
    white = subdf$white
    male = subdf$male
    educ = subdf$educ
    age_at_visit = subdf$age_at_visit
    origin_t = subdf$t
    origin_m = subdf$m
    origin_y = subdf$y
    t = ifelse(runif(n) < flip_p, 1 - origin_t, origin_t)
    
    m0_base = m_c[1] + white * m_c[2] + male * m_c[3] + educ * m_c[4] + age_at_visit * m_c[5]
    m1_base = m_c[6] + white * m_c[7] + male * m_c[8] + educ * m_c[9] + age_at_visit * m_c[10]
    
    m0 = ifelse(origin_t == 0, origin_m, origin_m - m1_base)
    m1 = ifelse(origin_t == 1, origin_m, origin_m + m1_base)
    m = ifelse(t == 1, m1, m0)
    
    y0 = ifelse(origin_t == 0, origin_y, origin_y - y_c[6] - m1_base * y_c[7] - origin_m * y_c[8])
    y1 = ifelse(origin_t == 1, origin_y, origin_y + y_c[6] + m1_base * y_c[7] + (origin_m + m1_base) * y_c[8])
    y = ifelse(t == 1, y1, y0)
    
    true_NIE = m1_base * (y_c[7] + y_c[8])
    true_NDE = y_c[6] + m0 * y_c[8]
    
    data = cbind.data.frame(x1 = white, x2 = male, x3 = educ, x4 = age_at_visit,
                            t = t, m = m, y = y, true_NIE = true_NIE, true_NDE = true_NDE)
    
  } else if (case == 6){ ### Use VCBCMF, original m and y, generate some for the true effects
    
    subdf = appdat[sample(nrow(appdat), n, replace = TRUE), , drop = FALSE]
    white = subdf$white
    male = subdf$male
    educ = subdf$educ
    age_at_visit = subdf$age_at_visit
    origin_t = subdf$t
    origin_m = subdf$m
    origin_y = subdf$y
    t = ifelse(runif(n) < flip_p, 1 - origin_t, origin_t)
    newdata = data.frame(white = white, male = male, educ = educ, age_at_visit = age_at_visit, t = t)
    
    X_m = quantile_normalize_bart(
      preprocess_df(
        newdata
      )[[1]]
    )
    
    pi_hat = get_ps(newdata, newdata, formula_ps)
    X_m = cbind(X_m, pi_hat)
    
    mu_m = vcbcmf$forest_mu_m$do_predict(X_m)
    tau = vcbcmf$forest_tau$do_predict(X_m)
    
    mu_m_unscale = mean(appdat$m) + (mu_m * sd(appdat$m))
    tau_unscale = tau * sd(appdat$m)
    
    m0 = ifelse(origin_t == 0, origin_m, origin_m - tau_unscale)
    m1 = ifelse(origin_t == 1, origin_m, origin_m + tau_unscale)
    m = ifelse(t == 1, m1, m0)
    
    newdata$m = m
    newdata$m0_hat = mu_m_unscale
    newdata$m1_hat = mu_m_unscale + tau_unscale
    
    X_y = quantile_normalize_bart(
      preprocess_df(
        newdata
      )[[1]]
    )
    X_y = cbind(X_y, pi_hat)
    
    mu_y = vcbcmf$forest_mu_y$do_predict(X_y)
    zeta = vcbcmf$forest_zeta$do_predict(X_y)
    d = vcbcmf$forest_d$do_predict(X_y)
    
    mu_y_unscale = mean(appdat$y) + (mu_y * sd(appdat$y)) - mean(appdat$m) / sd(appdat$m) * d * sd(appdat$y)
    zeta_unscale = zeta * sd(appdat$y)
    d_unscale = d * sd(appdat$y) / sd(appdat$m)
    
    
    y0 = ifelse(origin_t == 0, origin_y, origin_y - zeta_unscale - tau_unscale * d_unscale)
    y1 = ifelse(origin_t == 1, origin_y, origin_y + zeta_unscale + tau_unscale * d_unscale)
    y = ifelse(t == 1, y1, y0)
    
    true_NIE = d_unscale * tau_unscale
    true_NDE = zeta_unscale
    
    data = cbind.data.frame(x1 = newdata$white, x2 = newdata$male, x3 = newdata$educ, x4 = newdata$age_at_visit,
                            t = newdata$t, m = newdata$m, y = y, true_NIE = true_NIE, true_NDE = true_NDE)
  }
  
  return(data)
}

dat1 = NULL; dat2 = NULL; dat3 = NULL; dat4 = NULL; dat5 = NULL; dat6 = NULL; dat7 = NULL; dat8 = NULL;
for(i in 1:50){
  dat1 = rbind(dat1, gen(case = 1, n = 10000))
  dat2 = rbind(dat2, gen(case = 2, n = 10000))
  dat3 = rbind(dat3, gen(case = 3, n = 10000))
  dat4 = rbind(dat4, gen(case = 4, n = 10000))
  dat5 = rbind(dat5, gen(case = 5, n = 10000))
  dat6 = rbind(dat6, gen(case = 6, n = 10000))
}

saveRDS(dat1, file = "dat1.rds")
saveRDS(dat2, file = "dat2.rds")
saveRDS(dat3, file = "dat3.rds")
saveRDS(dat4, file = "dat4.rds")
saveRDS(dat5, file = "dat5.rds")
saveRDS(dat6, file = "dat6.rds")
