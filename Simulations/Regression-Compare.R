##Simulation for Normal Noises

library(odeintr)
library(numDeriv)
library(roptim)
library(deSolve)
library(foreach)
library(doParallel)

library(not)
library(nsp)
library(strucchangeRcpp)

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

RCUSUM = function(Y, Latent, intervals)
{
  num = nrow(intervals); p1 = 2
  
  print(paste0('Number of intervals:', num))
  
  cusum_stat = rep(0, num)
  cusum_pos = rep(0, num)
  beta_diff = rep(0, num)
  
  for (inval in 1:num)
  {
    aa = Sys.time()
    
    si = intervals[inval, 1]; ei = intervals[inval, 2]
    
    Xb = Inf;  Zb = -Inf
    for (bi in (si + p1 + 1):(ei - p1 - 1))
    {
      Yb1 = Y[si:(bi-1)]; Yb2 = Y[((bi+1):ei)]
      Tb1_reg = Latent[si:(bi-1)]; Tb2_reg = Latent[(bi+1):ei]
    
      optim_res1 = lm(Yb1 ~ Tb1_reg); optim_res2 = lm(Yb2 ~ Tb2_reg)
      
      parahat1 = optim_res1$coefficients[2]; parahat2 = optim_res2$coefficients[2]
      residualhat1 = sum((optim_res1$residuals)^2); residualhat2 = sum((optim_res2$residuals)^2)
      
      mi = ei - si
      
      scale0 = sqrt((bi - si)*(ei - bi)/mi)
      scale1 = scale2 = 1
      
      residual1 = scale1*residualhat1
      residual2 = scale2*residualhat2
      
      # Z1 = cbind(rep(1:length(Tb1_reg)), Tb1_reg); sigma_scale1 = solve(t(Z1) %*% Z1)[2, 2]*sqrt(bi - si)
      # Z2 = cbind(rep(1:length(Tb2_reg)), Tb2_reg); sigma_scale2 = solve(t(Z2) %*% Z2)[2, 2]*sqrt(ei - bi)
      
      Xtemp = abs(residual1 + residual2)
      betadiff_temp = (parahat1 - parahat2)
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

RSBS = function(Y, Latent, m, thres_beta = NULL, thres_residual = NULL, q.max = NULL, gap = 10, decay = sqrt(2))
{
  n = length(Y)
  
  intervals = seeded.intervals(n, m, decay)
  #intervals = random.intervals(n, 5000, m)
  
  bi_stat = RCUSUM(Y, Latent, intervals)
  s = 1; e = n
  
  candidate_cpt = bi_stat$candidate_cpt
  candidate_stat = bi_stat$candidate_stat
  beta_diff = bi_stat$beta_diff
  
  if (is.null(thres_beta))
  {
    #thres_const = quantile(beta_diff, seq(0.5, 0.99, 0.01))
    #thres_const = quantile(beta_diff, seq(100/length(beta_diff), 0.99, 1/length(beta_diff)))
    # num_infor = 2*q.max/(decay - 1) + q.max*log(n, decay) - (q.max + 1)*log(q.max + 1, decay) - q.max*log(m, decay)
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
  
  if (is.null(thres_residual))
  {
    thres_residual = max(candidate_stat)
  }
  
  index_diff = which(beta_diff > thres_beta)
  intervals = intervals[index_diff, ]
  candidate_cpt = candidate_cpt[index_diff]
  beta_diff = beta_diff[index_diff]
  candidate_stat = candidate_stat[index_diff]
  
  # if (is.null(thres_residual))
  # {
  #   thres_residual = sort(candidate_stat)
  # }
  # 
  # thres_residual = thres_residual[-1]
  
  cpts_list = list()
  cpts_stat_list = list()
  beta_diff_list = list()
  num_cpts = rep(0, length(thres_residual))
  
  for (i in 1:length(thres_residual))
  {
    index_residual = which(candidate_stat <= thres_residual[i])
    intervals_temp = intervals[index_residual, ]
    candidate_cpt_temp = candidate_cpt[index_residual]
    beta_diff_temp = beta_diff[index_residual]
    candidate_stat_temp = candidate_stat[index_residual]
    
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

R_Mops = function(Y, Latent, Omega, alpha, candidate_cpt)
{
  n = length(Y); p1 = 2
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
    
    YO0 = Y[index0_O]; TO0_reg = Latent[index0_O]
    YE0 = Y[index0_E]; TE0_reg = Latent[index0_E]
    
    YO1 = Y[index1_O]; TO1_reg = Latent[index1_O]
    YE1 = Y[index1_E]; TE1_reg = Latent[index1_E]
    
    n0 = e0 - s0 + 1; n1 = e1 - s1 + 1
     
    optim_resO0 = lm(YO0 ~ TO0_reg); optim_resE0 = lm(YE0 ~ TE0_reg)
    optim_resO1 = lm(YO1 ~ TO1_reg); optim_resE1 = lm(YE1 ~ TE1_reg)
    
    # betahatO0 = optim_resO0$coefficients[2]; parahatO1 = optim_resO1$coefficients
    # betahatE0 = optim_resE0$coefficients[2]; parahatE1 = optim_resE1$coefficients
    # 
    # inithatO0 = optim_resO0$coefficients[1]; inithatE0 = optim_resE0$coefficients[1]
    # midhatO0 = inithatO0 + betahatO0*(startO1 - startO0)
    # midhatE0 = inithatE0 + betahatE0*(startE1 - startE0)
    # 
    # parahatO0 = c(midhatO0, betahatO0); parahatE0 = c(midhatE0, betahatE0)
    
    parahatO0 = optim_resO0$coefficients[2]; parahatO1 = optim_resO1$coefficients[2]
    parahatE0 = optim_resE0$coefficients[2]; parahatE1 = optim_resE1$coefficients[2]
    
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
  
  print(emp_level)
  
  return(Final_cpt)
}

set.seed(4)

n = 2000
K = 10

Latent = rnorm(n, 1, 1)

true_cp = rep(0, K)
for(j in 1:K)
{
  true_cp[j] = j*floor(n/(K + 1)) + runif(1, - floor(n^(1/4)),  floor(n^(1/4)))
}

indexp = c(1, true_cp, n)

set.seed(1)
betav = rnorm(K + 1, 4, 1)
Curve_True = rep(0, n)

for (j in 1:(K + 1))
{
  startj = ceiling(indexp[j]); endj = floor(indexp[j+1])
  betaj = betav[j]
  Curve_True[(startj:endj)] = betaj*Latent[(startj:endj)] 
} 

# jpeg('test.jpg', width = 2000, height = 1000)
# plot(1:2000, Curve_True, type = 'l')
# dev.off()

signaltude = rep(0, K)
for (j in 2:(K + 1))
{
  startj0 = ceiling(indexp[j - 1]); endj0 = floor(indexp[j]); nj0 = endj0 - startj0 + 1
  
  startj = ceiling(indexp[j]); endj = floor(indexp[j+1]); nj = endj - startj + 1
  betaj2 = betav[j-1]
  
  Curve1 = Curve_True[startj:endj]
  Curve2 = betaj2*Latent[startj:endj]
  
  signaltude[j-1] = sqrt(sum((Curve1 - Curve2)^2)/(nj0 + nj))
}

signal_strength = min(signaltude)

gammav = c(2.5, 2.7, 2.9)

wn003 = sim_max_holder(1000, 500, 0.03)
lambda = as.numeric(quantile(wn003, 0.8))

for (set_s in 1:length(gammav))
{
  gamma = gammav[set_s]
  
  Simu = function(s)
  {
    set.seed(s+5) #The nsp_tvreg() function breaks in seed 4 under Gaussian noise
    
    #epsilon = rt(n, 3)*gamma*signal_strength/sqrt(3)
    epsilon = rnorm(n, 0, sd = gamma*signal_strength)
    Y = Curve_True + epsilon
    
    indexO = seq(1, n, by = 2)
    YO = Y[indexO]; LO = Latent[indexO]
    
    #thres1 = 0; thres2 = 0
    
    thres1 = NULL; thres2 = NULL; q.max = 20; m = 40
    
    detect_res = RSBS(YO, LO, m = m, thres_beta =NULL, thres_residual = NULL, q.max = q.max, gap = m/2)
    cpts_thresh = sort(detect_res$cpts_list[[1]][1:20])
    cpts = cpts_thresh*2 - 1
    
    if (min(cpts) <= 3) cpts = cpts[-1]
    if (max(cpts) >= n - 2) cpts = cpts[-length(cpts)]
    
    Mops_Res = R_Mops(Y, Latent, Omega = 1, alpha = 0.2, candidate_cpt = cpts)
    ##################################################################################################################################################################################
    not_results1 = not(Y, method = 'not', contrast = 'pcwsConstMean')
    not_summary1 = features(not_results1, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res1 = not_summary1$cpt
    
    not_results2 = not(Y, method = 'not', contrast = 'pcwsConstMeanHT')
    not_summary2 = features(not_results2, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res2 = not_summary2$cpt
    
    not_results3 = not(Y, method = 'not', contrast = 'pcwsLinContMean')
    not_summary3 = features(not_results3, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res3 = not_summary3$cpt
    
    not_results4 = not(Y, method = 'not', contrast = 'pcwsLinMean')
    not_summary4 = features(not_results4, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res4 = not_summary4$cpt
    
    not_results5 = not(Y, method = 'not', contrast = 'pcwsQuadMean')
    not_summary5 = features(not_results5, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res5 = not_summary5$cpt
    
    not_results6 = not(Y, method = 'not', contrast = 'pcwsConstMeanVar')
    not_summary6 = features(not_results6, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res6 = not_summary6$cpt
    
    ##################################################################################################################################################################################
    
    Nsp = nsp_tvreg(Y, as.matrix(Latent))
    Nsp_self = nsp_selfnorm(Y, as.matrix(Latent), M = 1000, lambda = lambda)
    ##################################################################################################################################################################################
    
    return(list(Mops_Res = Mops_Res, Not_SIC_Res1 = Not_SIC_Res1, Not_SIC_Res2 = Not_SIC_Res2, Not_SIC_Res3 = Not_SIC_Res3, 
                Not_SIC_Res4 = Not_SIC_Res4, Not_SIC_Res5 = Not_SIC_Res5, Not_SIC_Res6 = Not_SIC_Res6, Nsp = Nsp, Nsp_self = Nsp_self, detect_res = detect_res, Y = Y, Latent = Latent))
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
  
  save(list = c('CP_RES', 'true_cp'), file = paste0('Results_Lin/Regression-Compare_n =', n, '_noise sigma_N', set_s, '.RData'))
  
}

##################################################################################################################################################################################################
##Simulation for T Noises

rm(list = ls())

library(odeintr)
library(numDeriv)
library(roptim)
library(deSolve)
library(foreach)
library(doParallel)

library(not)
library(nsp)
library(strucchangeRcpp)

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

RCUSUM = function(Y, Latent, intervals)
{
  num = nrow(intervals); p1 = 2
  
  print(paste0('Number of intervals:', num))
  
  cusum_stat = rep(0, num)
  cusum_pos = rep(0, num)
  beta_diff = rep(0, num)
  
  for (inval in 1:num)
  {
    aa = Sys.time()
    
    si = intervals[inval, 1]; ei = intervals[inval, 2]
    
    Xb = Inf;  Zb = -Inf
    for (bi in (si + p1 + 1):(ei - p1 - 1))
    {
      Yb1 = Y[si:(bi-1)]; Yb2 = Y[((bi+1):ei)]
      Tb1_reg = Latent[si:(bi-1)]; Tb2_reg = Latent[(bi+1):ei]
      
      optim_res1 = lm(Yb1 ~ Tb1_reg); optim_res2 = lm(Yb2 ~ Tb2_reg)
      
      parahat1 = optim_res1$coefficients[2]; parahat2 = optim_res2$coefficients[2]
      residualhat1 = sum((optim_res1$residuals)^2); residualhat2 = sum((optim_res2$residuals)^2)
      
      mi = ei - si
      
      scale0 = sqrt((bi - si)*(ei - bi)/mi)
      scale1 = scale2 = 1
      
      residual1 = scale1*residualhat1
      residual2 = scale2*residualhat2
      
      # Z1 = cbind(rep(1:length(Tb1_reg)), Tb1_reg); sigma_scale1 = solve(t(Z1) %*% Z1)[2, 2]*sqrt(bi - si)
      # Z2 = cbind(rep(1:length(Tb2_reg)), Tb2_reg); sigma_scale2 = solve(t(Z2) %*% Z2)[2, 2]*sqrt(ei - bi)
      
      Xtemp = abs(residual1 + residual2)
      betadiff_temp = (parahat1 - parahat2)
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

RSBS = function(Y, Latent, m, thres_beta = NULL, thres_residual = NULL, q.max = NULL, gap = 10, decay = sqrt(2))
{
  n = length(Y)
  
  intervals = seeded.intervals(n, m, decay)
  #intervals = random.intervals(n, 5000, m)
  
  bi_stat = RCUSUM(Y, Latent, intervals)
  s = 1; e = n
  
  candidate_cpt = bi_stat$candidate_cpt
  candidate_stat = bi_stat$candidate_stat
  beta_diff = bi_stat$beta_diff
  
  if (is.null(thres_beta))
  {
    #thres_const = quantile(beta_diff, seq(0.5, 0.99, 0.01))
    #thres_const = quantile(beta_diff, seq(100/length(beta_diff), 0.99, 1/length(beta_diff)))
    # num_infor = 2*q.max/(decay - 1) + q.max*log(n, decay) - (q.max + 1)*log(q.max + 1, decay) - q.max*log(m, decay)
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
  
  if (is.null(thres_residual))
  {
    thres_residual = max(candidate_stat)
  }
  
  index_diff = which(beta_diff > thres_beta)
  intervals = intervals[index_diff, ]
  candidate_cpt = candidate_cpt[index_diff]
  beta_diff = beta_diff[index_diff]
  candidate_stat = candidate_stat[index_diff]
  
  # if (is.null(thres_residual))
  # {
  #   thres_residual = sort(candidate_stat)
  # }
  # 
  # thres_residual = thres_residual[-1]
  
  cpts_list = list()
  cpts_stat_list = list()
  beta_diff_list = list()
  num_cpts = rep(0, length(thres_residual))
  
  for (i in 1:length(thres_residual))
  {
    index_residual = which(candidate_stat <= thres_residual[i])
    intervals_temp = intervals[index_residual, ]
    candidate_cpt_temp = candidate_cpt[index_residual]
    beta_diff_temp = beta_diff[index_residual]
    candidate_stat_temp = candidate_stat[index_residual]
    
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

R_Mops = function(Y, Latent, Omega, alpha, candidate_cpt)
{
  n = length(Y); p1 = 2
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
    
    YO0 = Y[index0_O]; TO0_reg = Latent[index0_O]
    YE0 = Y[index0_E]; TE0_reg = Latent[index0_E]
    
    YO1 = Y[index1_O]; TO1_reg = Latent[index1_O]
    YE1 = Y[index1_E]; TE1_reg = Latent[index1_E]
    
    n0 = e0 - s0 + 1; n1 = e1 - s1 + 1
    
    optim_resO0 = lm(YO0 ~ TO0_reg); optim_resE0 = lm(YE0 ~ TE0_reg)
    optim_resO1 = lm(YO1 ~ TO1_reg); optim_resE1 = lm(YE1 ~ TE1_reg)
    
    # betahatO0 = optim_resO0$coefficients[2]; parahatO1 = optim_resO1$coefficients
    # betahatE0 = optim_resE0$coefficients[2]; parahatE1 = optim_resE1$coefficients
    # 
    # inithatO0 = optim_resO0$coefficients[1]; inithatE0 = optim_resE0$coefficients[1]
    # midhatO0 = inithatO0 + betahatO0*(startO1 - startO0)
    # midhatE0 = inithatE0 + betahatE0*(startE1 - startE0)
    # 
    # parahatO0 = c(midhatO0, betahatO0); parahatE0 = c(midhatE0, betahatE0)
    
    parahatO0 = optim_resO0$coefficients[2]; parahatO1 = optim_resO1$coefficients[2]
    parahatE0 = optim_resE0$coefficients[2]; parahatE1 = optim_resE1$coefficients[2]
    
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
  
  print(emp_level)
  
  return(Final_cpt)
}

set.seed(4)

n = 2000
K = 10

Latent = rnorm(n, 1, 1)

true_cp = rep(0, K)
for(j in 1:K)
{
  true_cp[j] = j*floor(n/(K + 1)) + runif(1, - floor(n^(1/4)),  floor(n^(1/4)))
}

indexp = c(1, true_cp, n)

set.seed(1)
betav = rnorm(K + 1, 4, 1)
Curve_True = rep(0, n)

for (j in 1:(K + 1))
{
  startj = ceiling(indexp[j]); endj = floor(indexp[j+1])
  betaj = betav[j]
  Curve_True[(startj:endj)] = betaj*Latent[(startj:endj)] 
} 

# jpeg('test.jpg', width = 2000, height = 1000)
# plot(1:2000, Curve_True, type = 'l')
# dev.off()

signaltude = rep(0, K)
for (j in 2:(K + 1))
{
  startj0 = ceiling(indexp[j - 1]); endj0 = floor(indexp[j]); nj0 = endj0 - startj0 + 1
  
  startj = ceiling(indexp[j]); endj = floor(indexp[j+1]); nj = endj - startj + 1
  betaj2 = betav[j-1]
  
  Curve1 = Curve_True[startj:endj]
  Curve2 = betaj2*Latent[startj:endj]
  
  signaltude[j-1] = sqrt(sum((Curve1 - Curve2)^2)/(nj0 + nj))
}

signal_strength = min(signaltude)

gammav = c(2.5, 2.7, 2.9)

wn003 = sim_max_holder(1000, 500, 0.03)
lambda = as.numeric(quantile(wn003, 0.8))

for (set_s in 1:length(gammav))
{
  gamma = gammav[set_s]
  
  Simu = function(s)
  {
    set.seed(s+5) #The nsp_tvreg() function breaks in seed 4 under Gaussian noise
    
    epsilon = rt(n, 3)*gamma*signal_strength/sqrt(3)
    #epsilon = rnorm(n, 0, sd = gamma*signal_strength)
    Y = Curve_True + epsilon
    
    indexO = seq(1, n, by = 2)
    YO = Y[indexO]; LO = Latent[indexO]
    
    #thres1 = 0; thres2 = 0
    
    thres1 = NULL; thres2 = NULL; q.max = 20; m = 40
    
    detect_res = RSBS(YO, LO, m = m, thres_beta =NULL, thres_residual = NULL, q.max = q.max, gap = m/2)
    cpts_thresh = sort(detect_res$cpts_list[[1]][1:20])
    cpts = cpts_thresh*2 - 1
    
    if (min(cpts) <= 3) cpts = cpts[-1]
    if (max(cpts) >= n - 2) cpts = cpts[-length(cpts)]
    
    Mops_Res = R_Mops(Y, Latent, Omega = 1, alpha = 0.2, candidate_cpt = cpts)
    ##################################################################################################################################################################################
    not_results1 = not(Y, method = 'not', contrast = 'pcwsConstMean')
    not_summary1 = features(not_results1, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res1 = not_summary1$cpt
    
    not_results2 = not(Y, method = 'not', contrast = 'pcwsConstMeanHT')
    not_summary2 = features(not_results2, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res2 = not_summary2$cpt
    
    not_results3 = not(Y, method = 'not', contrast = 'pcwsLinContMean')
    not_summary3 = features(not_results3, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res3 = not_summary3$cpt
    
    not_results4 = not(Y, method = 'not', contrast = 'pcwsLinMean')
    not_summary4 = features(not_results4, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res4 = not_summary4$cpt
    
    not_results5 = not(Y, method = 'not', contrast = 'pcwsQuadMean')
    not_summary5 = features(not_results5, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res5 = not_summary5$cpt
    
    not_results6 = not(Y, method = 'not', contrast = 'pcwsConstMeanVar')
    not_summary6 = features(not_results6, method = 'ic', q.max = 20, penalty = 'sic')
    Not_SIC_Res6 = not_summary6$cpt
    
    ##################################################################################################################################################################################
    
    Nsp = nsp_tvreg(Y, as.matrix(Latent))
    Nsp_self = nsp_selfnorm(Y, as.matrix(Latent), M = 1000, lambda = lambda)
    ##################################################################################################################################################################################
    
    return(list(Mops_Res = Mops_Res, Not_SIC_Res1 = Not_SIC_Res1, Not_SIC_Res2 = Not_SIC_Res2, Not_SIC_Res3 = Not_SIC_Res3, 
                Not_SIC_Res4 = Not_SIC_Res4, Not_SIC_Res5 = Not_SIC_Res5, Not_SIC_Res6 = Not_SIC_Res6, Nsp = Nsp, Nsp_self = Nsp_self, detect_res = detect_res, Y = Y, Latent = Latent))
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
  
  save(list = c('CP_RES', 'true_cp'), file = paste0('Results_Lin/Regression-Compare_n =', n, '_noise sigma_T', set_s, '.RData'))
  
}




