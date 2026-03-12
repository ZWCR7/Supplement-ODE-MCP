library(odeintr)
library(numDeriv)
library(roptim)
library(foreach)
library(doParallel)

#General Structural Change with Normal Errors
#############################################################################################################################################################################################
rm(list = ls())
pars = c("a", "b", "c")

FN.sys = '
    dxdt[0] = c*(x[0] - x[0]*x[0]*x[0]/3 + x[1]);
    dxdt[1] = -(x[0] - a + b*x[1])/c;
    '

compile_sys("ODE", FN.sys, pars, TRUE)

source('~/PDE-MCP/utils/residualFunc.R')
source('~/PDE-MCP/utils/seededInterval.R')
source('~/PDE-MCP/utils/StatGenerator.R')
source('~/PDE-MCP/utils/SNOT-Algorithm.R')
source('~/PDE-MCP/utils/mopsFunc.R')

n = 2000
K = 10

set.seed(1)

Time = 0.1*(1:n)

true_cp = rep(0, K)
for(j in 1:K)
{
  true_cp[j] = j*floor(n/(K + 1)) + runif(1, - floor(n^(1/4)),  floor(n^(1/4)))
}

indexp = c(1, true_cp, n)
m = 40

set.seed(3)
betav = rnorm(K + 1, 0.6, 0.2); initj = c(-1, 1)
Curve_True = matrix(0, n, 2)

for (j in 1:(K + 1))
{
  betaj = betav[j]
  startj = ceiling(indexp[j]); endj = floor(indexp[j+1])
  ODE_set_params(a = 0.2, b = betaj, c = 3)

  if (j == 1)
  {
    init_startj = 0
  }
  else
  {
    init_startj = Time[startj - 1]
  }

  Curve_True[(startj:endj), ] = as.matrix(ODE_at(init = initj, times = c(init_startj, Time[startj:endj]), start = init_startj)[-1, -1])
  initj = Curve_True[endj, ] + rnorm(2, 0, 0.2)
}

signaltude = rep(0, K)
for (j in 2:(K + 1))
{
  startj0 = ceiling(indexp[j - 1]); endj0 = floor(indexp[j]); nj0 = endj0 - startj0 + 1

  startj = ceiling(indexp[j]); endj = floor(indexp[j+1]); nj = endj - startj + 1
  betaj2 = betav[j-1]

  init_startj = Time[startj - 1]; initj2 = Curve_True[startj - 1, ]
  Curve1 = Curve_True[startj:endj, ]

  ODE_set_params(a = 0.2, b = betaj2, c = 3)
  Curve2 = as.matrix(ODE_at(init = initj2, times = c(init_startj, Time[startj:endj]), start = init_startj)[-1, -1])

  signaltude[j-1] = sqrt(sum((Curve1 - Curve2)^2)/(nj0 + nj)/2)
}

signal_strength = min(signaltude)

jpeg('test.jpg', width = 2000, height = 1000)
plot(Time, Curve_True[, 1], type = 'l')
dev.off()

gammav = c(0.5, 0.7, 0.9, 1.1, 1.3)

for (set_s in 1:length(gammav))
{
  gamma = gammav[set_s]
  #load(paste0('Results_FN/Mops-FN-General_n =', n, '_noise sigma_N', set_s, '.RData'))
  
  Simu = function(s)
  {
    set.seed(s)

    pars = c("a", "b", "c")

    FN.sys = '
    dxdt[0] = c*(x[0] - x[0]*x[0]*x[0]/3 + x[1]);
    dxdt[1] = -(x[0] - a + b*x[1])/c;
    '

    compile_sys("ODE", FN.sys, pars, TRUE)

    epsilon = matrix(rnorm(2*n, 0, sd = gamma*signal_strength), n, 2)
    Y = Curve_True + epsilon
    
    #Y = CP_RES[[s]]$Y; Time = CP_RES[[s]]$Time

    indexO = seq(1, n, by = 2)
    YO = Y[indexO, ]; TO = Time[indexO]
    nO = length(indexO)

    intervals = seeded.intervals(nO, m)
    bi_stat = Stat_Generator(YO, TO, intervals, method = 'L-BFGS-B', beta0 = 0.6, betal = 0, betau = 1.1)
    #bi_stat = CP_RES[[s]]$detect_res_ER$bi_stat
    detect_res_ER = SNot(bi_stat, thres_beta = NULL, thres_GLR = NULL, q.max = 40, method = 'upper.bound', Statistic = 'ER', gap = m/2)

    cpts_thresh_ER = sort(detect_res_ER$cpts_list[[1]])
    cpts_ER = cpts_thresh_ER*2 - 1

    if (min(cpts_ER) <= 3) cpts_ER = cpts_ER[-1]
    if (max(cpts_ER) >= n - 2) cpts_ER = cpts_ER[-length(cpts_ER)]
    Mops_ER_Res = PMops(Y, Time, Omega = diag(rep(1, 3)), alpha = 0.2, method = 'L-BFGS-B', candidate_cpt = cpts_ER, beta0 = 0.6, betal = 0, betau = 1.1)
    #####################################################################################################################################################################

    detect_res_GLR = SNot(bi_stat, thres_beta = NULL, thres_GLR = NULL, q.max = 40, method = 'upper.bound', Statistic = 'GLR', gap = m/2)

    cpts_thresh_GLR = sort(detect_res_GLR$cpts_list[[1]])
    cpts_GLR = cpts_thresh_GLR*2 - 1

    if (min(cpts_GLR) <= 3) cpts_GLR = cpts_GLR[-1]
    if (max(cpts_GLR) >= n - 2) cpts_GLR = cpts_GLR[-length(cpts_GLR)]
    Mops_GLR_Res = PMops(Y, Time, Omega = diag(rep(1, 3)), alpha = 0.2, method = 'L-BFGS-B', candidate_cpt = cpts_GLR, beta0 = 0.6, betal = 0, betau = 1.1)

    return(list(Mops_ER_Res = Mops_ER_Res, Mops_GLR_Res = Mops_GLR_Res, detect_res_ER = detect_res_ER, detect_res_GLR = detect_res_GLR, Y = Y, Time = Time))
  }

  nsimu = 100

  aa = Sys.time()
  cl = makeCluster(100)
  registerDoParallel(cl)

  CP_RES = foreach(s = 1:nsimu, .packages = c("odeintr", "numDeriv", "roptim")) %dopar% Simu(s)

  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()

  print(bb - aa)

  save(list = c('CP_RES', 'true_cp'), file = paste0('Results_FN/Mops-FN-General_n =', n, '_K = 20_noise sigma_N', set_s, '.RData'))
}

############################################################################################################################################################################################################################

#General Structural Change with t Errors
#############################################################################################################################################################################################
rm(list = ls())
pars = c("a", "b", "c")

FN.sys = '
    dxdt[0] = c*(x[0] - x[0]*x[0]*x[0]/3 + x[1]);
    dxdt[1] = -(x[0] - a + b*x[1])/c;
    '

compile_sys("ODE", FN.sys, pars, TRUE)

source('~/PDE-MCP/utils/residualFunc.R')
source('~/PDE-MCP/utils/seededInterval.R')
source('~/PDE-MCP/utils/StatGenerator.R')
source('~/PDE-MCP/utils/SNOT-Algorithm.R')
source('~/PDE-MCP/utils/mopsFunc.R')

n = 2000
K = 10

set.seed(1)

Time = 0.1*(1:n)

true_cp = rep(0, K)
for(j in 1:K)
{
  true_cp[j] = j*floor(n/(K + 1)) + runif(1, - floor(n^(1/4)),  floor(n^(1/4)))
}

indexp = c(1, true_cp, n)
m = 40

set.seed(3)
betav = rnorm(K + 1, 0.6, 0.2); initj = c(-1, 1)
Curve_True = matrix(0, n, 2)

for (j in 1:(K + 1))
{
  betaj = betav[j]
  startj = ceiling(indexp[j]); endj = floor(indexp[j+1])
  ODE_set_params(a = 0.2, b = betaj, c = 3)

  if (j == 1)
  {
    init_startj = 0
  }
  else
  {
    init_startj = Time[startj - 1]
  }

  Curve_True[(startj:endj), ] = as.matrix(ODE_at(init = initj, times = c(init_startj, Time[startj:endj]), start = init_startj)[-1, -1])
  initj = Curve_True[endj, ] + rnorm(2, 0, 0.2)
}

signaltude = rep(0, K)
for (j in 2:(K + 1))
{
  startj0 = ceiling(indexp[j - 1]); endj0 = floor(indexp[j]); nj0 = endj0 - startj0 + 1

  startj = ceiling(indexp[j]); endj = floor(indexp[j+1]); nj = endj - startj + 1
  betaj2 = betav[j-1]

  init_startj = Time[startj - 1]; initj2 = Curve_True[startj - 1, ]
  Curve1 = Curve_True[startj:endj, ]

  ODE_set_params(a = 0.2, b = betaj2, c = 3)
  Curve2 = as.matrix(ODE_at(init = initj2, times = c(init_startj, Time[startj:endj]), start = init_startj)[-1, -1])

  signaltude[j-1] = sqrt(sum((Curve1 - Curve2)^2)/(nj0 + nj)/2)
}

signal_strength = min(signaltude)

jpeg('test.jpg', width = 2000, height = 1000)
plot(Time, Curve_True[, 1], type = 'l')
dev.off()

gammav = c(0.5, 0.7, 0.9, 1.1, 1.3)

for (set_s in 1:length(gammav))
{
  gamma = gammav[set_s]

  Simu = function(s)
  {
    set.seed(s)

    pars = c("a", "b", "c")

    FN.sys = '
    dxdt[0] = c*(x[0] - x[0]*x[0]*x[0]/3 + x[1]);
    dxdt[1] = -(x[0] - a + b*x[1])/c;
    '

    compile_sys("ODE", FN.sys, pars, TRUE)

    epsilon = matrix(rt(2*n, 3)*gamma*signal_strength/sqrt(3), n, 2)
    Y = Curve_True + epsilon

    indexO = seq(1, n, by = 2)
    YO = Y[indexO, ]; TO = Time[indexO]
    nO = length(indexO)

    intervals = seeded.intervals(nO, m)
    bi_stat = Stat_Generator(YO, TO, intervals, method = 'L-BFGS-B', beta0 = 0.6, betal = 0, betau = 1.1)

    detect_res_ER = SNot(bi_stat, thres_beta = NULL, thres_GLR = NULL, q.max = 40, method = 'upper.bound', Stat = 'ER', gap = m/2)

    cpts_thresh_ER = sort(detect_res_ER$cpts_list[[1]])
    cpts_ER = cpts_thresh_ER*2 - 1

    if (min(cpts_ER) <= 3) cpts_ER = cpts_ER[-1]
    if (max(cpts_ER) >= n - 2) cpts_ER = cpts_ER[-length(cpts_ER)]
    Mops_ER_Res = PMops(Y, Time, Omega = diag(rep(1, 3)), alpha = 0.2, method = 'L-BFGS-B', candidate_cpt = cpts_ER, beta0 = 0.6, betal = 0, betau = 1.1)
    #####################################################################################################################################################################

    detect_res_GLR = SNot(bi_stat, thres_beta = NULL, thres_GLR = NULL, q.max = 40, method = 'upper.bound', Stat = 'GLR', gap = m/2)

    cpts_thresh_GLR = sort(detect_res_GLR$cpts_list[[1]])
    cpts_GLR = cpts_thresh_GLR*2 - 1

    if (min(cpts_GLR) <= 3) cpts_GLR = cpts_GLR[-1]
    if (max(cpts_GLR) >= n - 2) cpts_GLR = cpts_GLR[-length(cpts_GLR)]
    Mops_GLR_Res = PMops(Y, Time, Omega = diag(rep(1, 3)), alpha = 0.2, method = 'L-BFGS-B', candidate_cpt = cpts_GLR, beta0 = 0.6, betal = 0, betau = 1.1)

    return(list(Mops_ER_Res = Mops_ER_Res, Mops_GLR_Res = Mops_GLR_Res, detect_res_ER = detect_res_ER, detect_res_GLR = detect_res_GLR, Y = Y, Time = Time))
  }

  nsimu = 100

  aa = Sys.time()
  cl = makeCluster(100)
  registerDoParallel(cl)

  CP_RES = foreach(s = 1:nsimu, .packages = c("odeintr", "numDeriv", "roptim")) %dopar% Simu(s)

  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()

  print(bb - aa)

  save(list = c('CP_RES', 'true_cp'), file = paste0('Results_FN/Mops-FN-General_n =', n, '_noise sigma_T', set_s, '.RData'))
}

############################################################################################################################################################################################################################







