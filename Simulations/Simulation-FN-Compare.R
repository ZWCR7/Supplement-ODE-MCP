library(odeintr)
library(numDeriv)
library(roptim)
library(foreach)
library(doParallel)

library(not)
library(nsp)
library(strucchangeRcpp)

pars = c("a", "b", "c")

FN.sys = '
    dxdt[0] = c*(x[0] - x[0]*x[0]*x[0]/3 + x[1]);
    dxdt[1] = -(x[0] - a + b*x[1])/c;
    '

compile_sys("ODE", FN.sys, pars, TRUE)

residual_func = function(Y, Time, start, init, beta)
{
  n = length(Time); Time1 = unique(c(start, Time))
  ODE_set_params(a = 0.2, b = beta, c = 3)
  solve_y = ODE_at(init = init, times = Time1, start = start)
  
  if (length(Time1) == (n + 1)) pred_y = solve_y[-1, 2]
  if (length(Time1) == n) pred_y = solve_y[, 2]
  
  residual = sum((Y - pred_y)^2)
  
  return(residual)
}

seeded.intervals = function(n, m = 10, decay = sqrt(2), unique.int = T)
{
  n	= as.integer(n)
  depth	= log(n, base = decay)
  depth	= ceiling(depth)
  M	= sum(2^(1:depth)-1)
  
  boundary_mtx = matrix(NA, ncol = 2)
  colnames(boundary_mtx) = c("st", "end")
  boundary_mtx[1, ] = c(1, n)
  
  depth	<- log(n, base = decay)
  depth	<- ceiling(depth)
  
  
  for(i in 2:depth){
    int_length	<- n * (1/decay)^(i-1)
    
    n_int		<- ceiling(round(n/int_length, 14))*2-1		# sometimes very slight numerical inaccuracies
    
    boundary_mtx	<- rbind(boundary_mtx,
                          cbind(floor(seq(1, n-int_length, length.out = (n_int))), 
                                ceiling(seq(int_length, n, length.out = (n_int)))))
  }
  
  if(unique.int) boundary_mtx = unique(boundary_mtx);
  
  
  length.interval = boundary_mtx[, 2] - boundary_mtx[, 1]
  index_include = which(length.interval >= m)
  
  boundary_mtx = boundary_mtx[index_include, ]
  return(boundary_mtx)
}

RCUSUM = function(Y, Time, intervals, method = 'L-BFGS-B', beta0 = NULL, betal = NULL, betau = NULL)
{
  num = nrow(intervals); p = length(beta0); p1 = length(beta0) + 2
  
  print(paste0('Number of intervals:', num))
  
  cusum_stat = rep(0, num)
  cusum_pos = rep(0, num)
  beta_diff = rep(0, num)
  
  for (inval in 1:num)
  {
    aa = Sys.time()
    
    si = intervals[inval, 1]; ei = intervals[inval, 2]
    
    Xb = Inf; Zb = -Inf
    for (bi in (si + p1):(ei - p1 - 1))
    {
      Yb1 = Y[si:(bi-1)]; Yb2 = Y[((bi+1):ei)]
      Tb1 = Time[si:(bi-1)]; Tb2 = Time[(bi+1):ei]
      
      start1 = Time[si]; start2 = Time[bi]
      
      residual1_func = function(beta) 
      {
        q = length(beta)
        return(residual_func(Yb1, Tb1, start1, beta[(p+1):q], beta[1:p]))
      }
      
      residual2_func = function(beta) 
      {
        q = length(beta)
        return(residual_func(Yb2, Tb2, start2, beta[(p+1):q], beta[1:p]))
      }
      
      init1 = c(Y[si], -Y[si]); init2 = c(Y[bi], -Y[bi])
      initl = c(-3, -3); initu = c(3, 3)
      lower = c(betal, initl); upper = c(betau, initu)
      #lower2 = c(betal, initl); upper2 = c(betau, initu)
      
      optim_res1 = optim(c(beta0, init1), residual1_func, NULL, method = method, lower = lower, upper = upper)
      optim_res2 = optim(c(beta0, init2), residual2_func, NULL, method = method, lower = lower, upper = upper)
      
      betahat1 = optim_res1$par[1:p]; parahat2 = optim_res2$par
      start_value_hat = optim_res1$par[(p+1):p1]
      
      ODE_set_params(a = 0.2, b = betahat1, c = 3)
      end_value_hat = as.vector(as.matrix(ODE_at(init = start_value_hat, times = c(start1, start2), start = start1))[-1, -1])
      
      parahat1 = c(betahat1, end_value_hat)
      residualhat1 = optim_res1$value; residualhat2 = optim_res2$value
      
      mi = ei - si
      
      scale0 = sqrt((bi - si)*(ei - bi)/mi)
      
      Xtemp = abs(residualhat1 + residualhat2)
      betadiff_temp = (parahat1 - parahat2)*c(1, 1, 0)
      Ztemp = scale0*norm(betadiff_temp, '2')
      
      if (Xtemp < Xb)
      {
        Xb = Xtemp
        cusum_pos[inval] = bi
        cusum_stat[inval] = Xb
        #beta_diff[inval] = scale0*norm(parahat1 - parahat2, '2')
      }
      
      if (Ztemp > Zb)
      {
        Zb = Ztemp
        beta_diff[inval] = Zb
      }
    }
    
    bb = Sys.time()
    
    print(paste0('interval', inval, '-------finished; Time = ', bb - aa))
  }
  
  return(list(candidate_cpt = cusum_pos, candidate_stat = cusum_stat, beta_diff = beta_diff))
}

Not_ERDD = function(Y, Time, m, thres_beta = NULL, thres_residual = NULL, q.max = NULL, method = 'L-BFGS-B', beta0 = NULL, betal = NULL, betau = NULL, gap = 10, decay = sqrt(2))
{
  n = length(Time)
  
  intervals = seeded.intervals(n, m, decay)
  
  bi_stat = RCUSUM(Y, Time, intervals, method, beta0, betal, betau)
  s = 1; e = n
  
  candidate_cpt = bi_stat$candidate_cpt
  candidate_stat = bi_stat$candidate_stat
  beta_diff = bi_stat$beta_diff
  
  if (is.null(thres_beta))
  {
    regular_points = seq(1, n, length.out = q.max + 2)
    regular_points = regular_points[-c(1, q.max + 2)]
    
    infornum = 0
    for (j in 1:nrow(intervals))
    {
      lower = intervals[j, 1]; upper = intervals[j, 2]
      index_info = ((regular_points >= lower) & (regular_points <= upper))
      
      if (sum(index_info != 0))
      {
        infornum = infornum + 1
      }
    }
    
    thres_beta = sort(beta_diff, decreasing = T)[infornum]
  }
  
  cpts_list = list()
  cpts_stat_list = list()
  beta_diff_list = list()
  num_cpts = rep(0, length(thres_beta))
  
  for (i in 1:length(thres_beta))
  {
    index_diff = which(beta_diff >= thres_beta[i])
    intervals_temp = intervals[index_diff, ]
    candidate_cpt_temp = candidate_cpt[index_diff]
    beta_diff_temp = beta_diff[index_diff]
    candidate_stat_temp = candidate_stat[index_diff]
    
    cpts = c(); cpts_stat = c(); cpts_diff = c()
    
    index_stop = 9999; not_num = 1
    while(index_stop > 2)
    {
      print(paste0('Not Step:', not_num))
      
      interval_lengths = intervals_temp[, 2] - intervals_temp[, 1]
      
      min_length = min(interval_lengths)
      length_candidate_index = which(interval_lengths <= min_length + 2)
      
      maxdiff = max(beta_diff_temp[length_candidate_index])
      index_cpt = which(beta_diff_temp == maxdiff)
      stat = candidate_stat_temp[index_cpt]
      diff = beta_diff_temp[index_cpt]
      
      cpt = unique(candidate_cpt_temp[index_cpt])
      
      Belong_Func = function(x) return((x[1] >= cpt + gap) + (x[2] <= cpt - gap))
      Belong_index = apply(intervals_temp, 1, Belong_Func)
      
      index_belong = which(Belong_index == 1)
      intervals_temp = as.matrix(intervals_temp[index_belong, ])
      
      candidate_cpt_temp = candidate_cpt_temp[index_belong]
      candidate_stat_temp = candidate_stat_temp[index_belong]
      beta_diff_temp = beta_diff_temp[index_belong]
      
      cpts = c(cpts, cpt)
      cpts_stat = c(cpts_stat, stat)
      cpts_diff = c(cpts_diff, diff)
      
      index_stop = length(intervals_temp)
      
      not_num = not_num + 1
    }
    
    num_cpts[i] = length(cpts)
    cpts_list = c(cpts_list, list(cpts))
    cpts_stat_list = c(cpts_stat_list, list(cpts_stat))
    beta_diff_list = c(beta_diff_list, list(cpts_diff))
  }
  
  return(list(cpts_list = cpts_list, cpts_stat_list = cpts_stat_list, beta_diff_list = beta_diff_list, bi_stat = bi_stat, num_cpts = num_cpts))
}

R_Mops = function(Y, Time, Omega, alpha, method = 'L-BFGS-B', candidate_cpt, beta0 = NULL, betal = NULL, betau = NULL)
{
  n = length(Time); p = length(beta0); p1 = length(beta0) + 2
  mid_cpt1 = c(0, candidate_cpt); mid_cpt2 = c(candidate_cpt, n + 1)
  
  W = rep(0, length(candidate_cpt))
  for (j in 1:length(candidate_cpt))
  {
    print(paste0('Candidate change point', j))
    
    s0 = mid_cpt1[j] + 1; e0 = mid_cpt2[j] - 1
    index0_O = seq(2*ceiling((s0 - 1)/2) + 1, 2*floor((e0 - 1)/2) + 1, by = 2)
    index0_E = seq(2*ceiling(s0/2), 2*floor(e0/2), by = 2)
    
    s1 = mid_cpt1[j + 1] + 1; e1 = mid_cpt2[j + 1] - 1
    index1_O = seq(2*ceiling((s1 - 1)/2) + 1, 2*floor((e1 - 1)/2) + 1, by = 2)
    index1_E = seq(2*ceiling(s1/2), 2*floor(e1/2), by = 2)
    
    inner = mid_cpt2[j]
    
    YO0 = Y[index0_O]; TO0 = Time[index0_O]
    YE0 = Y[index0_E]; TE0 = Time[index0_E]
    
    YO1 = Y[index1_O]; TO1 = Time[index1_O]
    YE1 = Y[index1_E]; TE1 = Time[index1_E]
    
    n0 = e0 - s0 + 1; n1 = e1 - s1 + 1
    startO0 = TO0[1]; startO1 = Time[inner]
    startE0 = TE0[1]; startE1 = Time[inner]
    
    initO0 = c(YO0[1], -YO0[1]); initO1 = c(Y[inner], -Y[inner])
    initE0 = c(YE0[1], -YE0[1]); initE1 = c(Y[inner], -Y[inner])
    
    residualO0_func = function(beta) 
    {
      q = length(beta)
      return(residual_func(YO0, TO0, startO0, beta[(p+1):q], beta[1:p]))
    }
    
    residualE0_func = function(beta) 
    {
      q = length(beta)
      return(residual_func(YE0, TE0, startE0, beta[(p+1):q], beta[1:p]))
    }
    
    residualO1_func = function(beta) 
    {
      q = length(beta)
      return(residual_func(YO1, TO1, startO1, beta[(p+1):q], beta[1:p]))
    }
    
    residualE1_func = function(beta) 
    {
      q = length(beta)
      return(residual_func(YE1, TE1, startE1, beta[(p+1):q], beta[1:p]))
    }
    
    initl = c(-3, -3); initu = c(3, 3)
    lower = c(betal, initl); upper = c(betau, initu)
    
    optim_resO0 = optim(c(beta0, initO0), residualO0_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
    optim_resE0 = optim(c(beta0, initE0), residualE0_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
    
    optim_resO1 = optim(c(beta0, initO1), residualO1_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
    optim_resE1 = optim(c(beta0, initE1), residualE1_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
    
    betahatO0 = optim_resO0$par[1:p]; parahatO1 = optim_resO1$par
    start_value_hatO0 =  optim_resO0$par[(p+1):p1]
    
    ODE_set_params(a = 0.2, b = betahatO0, c = 3)
    end_value_hatO0 = as.vector(as.matrix(ODE_at(init = start_value_hatO0, times = c(startO0, startO1), start = startO0))[-1, -1])
    parahatO0 = c(betahatO0, end_value_hatO0)
    
    betahatE0 = optim_resE0$par[1:p]; parahatE1 = optim_resE1$par
    start_value_hatE0 =  optim_resE0$par[(p+1):p1]
    
    ODE_set_params(a = 0.2, b = betahatE0, c = 3)
    end_value_hatE0 = as.vector(as.matrix(ODE_at(init = start_value_hatE0, times = c(startE0, startE1), start = startE0))[-1, -1])
    parahatE0 = c(betahatE0, end_value_hatE0)
    
    print(parahatO0); print(parahatO1); print(parahatE0); print(parahatE1)
    
    SObar_diff = parahatO0 - parahatO1
    SEbar_diff = parahatE0 - parahatE1
    
    W[j] = n0*n1/(n0 + n1)*t(SObar_diff) %*% Omega %*% SEbar_diff
  }
  
  emp_level = rep(0, length(W))
  
  for (j in 1:length(W))
  {
    numj = 1 + sum(W <= -abs(W[j]))
    denj = max(1, sum(W >= abs(W[j])))
    
    emp_level[j] = numj/denj
  }
  
  indx_thres = which(emp_level <= alpha)
  if (length(indx_thres) == 0)
  {
    thres = Inf
  }
  else
  {
    thres = min(abs(W[indx_thres]))
  }
  
  Final_cpt = candidate_cpt[which(W >= thres)]
  
  return(Final_cpt)
}


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

pars = c("a", "b", "c")

FN.sys = '
    dxdt[0] = c*(x[0] - x[0]*x[0]*x[0]/3 + x[1]);
    dxdt[1] = -(x[0] - a + b*x[1])/c;
    '

compile_sys("ODE", FN.sys, pars, TRUE)

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
  Curve1 = Curve_True[startj:endj, 1]
  
  ODE_set_params(a = 0.2, b = betaj2, c = 3)
  Curve2 = as.matrix(ODE_at(init = initj2, times = c(init_startj, Time[startj:endj]), start = init_startj)[-1, -1])[, 1]
  
  signaltude[j-1] = sqrt(sum((Curve1 - Curve2)^2)/(nj0 + nj))
}

signal_strength = min(signaltude)

gammav = c(0.5, 0.7, 0.9, 1.1, 1.3)

for (set_s in 1:length(gammav))
{
  gamma = gammav[set_s]
  
  Simu = function(s)
  {
    pars = c("a", "b", "c")
    
    FN.sys = '
    dxdt[0] = c*(x[0] - x[0]*x[0]*x[0]/3 + x[1]);
    dxdt[1] = -(x[0] - a + b*x[1])/c;
    '
    
    compile_sys("ODE", FN.sys, pars, TRUE)
    
    set.seed(s)
    
    epsilon = matrix(rnorm(2*n, 0, sd = gamma*signal_strength), n, 2)
    Y = Curve_True + epsilon
    
    indexO = seq(1, n, by = 2)
    YO = Y[indexO, ]; TO = Time[indexO]
    Y1 = Y[, 1]; YO1 = YO[, 1]
    
    detect_res = Not_ERDD(YO1, TO, m, NULL, NULL, q.max = 20, method = 'L-BFGS-B', beta0 = 0.6, betal = 0, betau = 1.1, gap = m/2)
    
    # index_thres = max(which(detect_res$num_cpts <= 20))
    # cpts_thresh = sort(detect_res$cpts_list[[index_thres]])
    
    cpts_thresh = sort(detect_res$cpts_list[[1]])
    cpts = cpts_thresh*2 - 1
    
    if (min(cpts) <= 3) cpts = cpts[-1]
    if (max(cpts) >= n - 2) cpts = cpts[-length(cpts)]
    
    Mops_Res = R_Mops(Y1, Time, Omega = diag(c(1, 1, 0)), alpha = 0.2, method = 'L-BFGS-B', candidate_cpt = cpts, beta0 = 0.6, betal = 0, betau = 1.1)
    ##################################################################################################################################################################################
    
    not_results1 = not(Y1, method = 'not', contrast = 'pcwsConstMean')
    not_summary1 = features(not_results1, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res1 = not_summary1$cpt
    
    not_results2 = not(Y1, method = 'not', contrast = 'pcwsConstMeanHT')
    not_summary2 = features(not_results2, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res2 = not_summary2$cpt
    
    not_results3 = not(Y1, method = 'not', contrast = 'pcwsLinContMean')
    not_summary3 = features(not_results3, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res3 = not_summary3$cpt
    
    not_results4 = not(Y1, method = 'not', contrast = 'pcwsLinMean')
    not_summary4 = features(not_results4, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res4 = not_summary4$cpt
    
    not_results5 = not(Y1, method = 'not', contrast = 'pcwsQuadMean')
    not_summary5 = features(not_results5, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res5 = not_summary5$cpt
    
    not_results6 = not(Y1, method = 'not', contrast = 'pcwsConstMeanVar')
    not_summary6 = features(not_results6, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res6 = not_summary6$cpt
    
    not_results_O = not(YO1, method = 'not', contrast = 'pcwsQuadMean')
    loc_O = min(which(not_results_O$solution.path$n.cpt >= 20))
    cpts_not = sort(not_results_O$solution.path$cpt[[loc_O]])
    
    cpts_not = cpts_not*2 - 1
    
    Not_Mops_Res = R_Mops(Y1, Time, Omega = diag(c(1, 1, 0)), alpha = 0.2, method = 'L-BFGS-B', candidate_cpt = cpts_not, beta0 = 0.6, betal = 0, betau = 1.1)
    ##################################################################################################################################################################################
    
    Nsp_Res_1 = nsp_poly(Y1, deg = 1, alpha = 0.2); Nsp_Res_2 = nsp_poly(Y1, deg = 2, alpha = 0.2); Nsp_Res_3 = nsp_poly(Y1, deg = 3, alpha = 0.2)
    #Nsp_Self_1 = nsp_poly_selfnorm(Y1, deg = 1, alpha = 0.2); Nsp_Self_2 = nsp_poly_selfnorm(Y1, deg = 2, alpha = 0.2); Nsp_Self_3 = nsp_poly_selfnorm(Y1, deg = 3, alpha = 0.2)
    ##################################################################################################################################################################################
    
    #save(list = c('Mops_Res'), file = paste0('test/FN-Compare_n =', n, '_noise sigma_', set_s, '_seed_', s, '.RData'))
    
    return(list(Mops_Res = Mops_Res, Not_SIC_Res1 = Not_SIC_Res1, Not_SIC_Res2 = Not_SIC_Res2, Not_SIC_Res3 = Not_SIC_Res3, Not_SIC_Res4 = Not_SIC_Res4, Not_SIC_Res5 = Not_SIC_Res5, Not_SIC_Res6 = Not_SIC_Res6, Not_Mops_Res = Not_Mops_Res, 
                Nsp_Res_1 = Nsp_Res_1, Nsp_Res_2 = Nsp_Res_2, Nsp_Res_3 = Nsp_Res_3, detect_res = detect_res, Y = Y, Time = Time))
  }
  
  nsimu = 100
  
  aa = Sys.time()
  cl = makeCluster(100)
  registerDoParallel(cl)
  
  CP_RES = foreach(s = 1:nsimu, .packages = c("odeintr", "numDeriv", "roptim", "not", "nsp", "strucchangeRcpp")) %dopar% Simu(s)
  
  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()
  
  print(bb - aa)
  
  save(list = c('CP_RES', 'true_cp'), file = paste0('Results_FN/Mops-FN-Compare_n = ', n, '_K = ', K, '_noise sigma_N', set_s, '.RData'))
  
}

















