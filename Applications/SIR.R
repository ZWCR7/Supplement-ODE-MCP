library(odeintr)
library(numDeriv)
library(roptim)
library(doParallel)
library(foreach)

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

load('covid19Italy.RData')

exclude = 1:27

ItalyIndex = ItalyIndex[-exclude]
ItalyInfect = ItalyInfect[-exclude]
ItalySuscep = ItalySuscep[-exclude]
ItalyTime = ItalyTime[-exclude]
ItalyRecovery = Population - ItalyInfect -ItalySuscep

ItalyIndex = ItalyIndex - min(ItalyIndex)
Time = ItalyIndex
Y = cbind(ItalySuscep, ItalyInfect)/Population

indexodd = seq(1, length(Time), by = 2)
YO = Y[indexodd, ]; TO = Time[indexodd]

n = length(TO)
m = 16
intervals = seeded.intervals(n, m)

beta0 = c(0.3, 0.3)
betal = c(1, 1)/Population
betau = c(1, 1)

residual_func = function(Y, Time, start, init, beta)
{
  n = length(Time); Time1 = unique(c(start, Time))
  SIR_set_params(alpha = beta[1], gamma = beta[2])
  solve_y = SIR_at(init = init, times = Time1, start = start)
  
  if (length(Time1) == (n + 1)) pred_y = solve_y[-1, -1]
  if (length(Time1) == n) pred_y = solve_y[, -1]
  
  #scale1 = mean(Y[, 1]); scale2 = mean(Y[, 2])
  
  Y1 = Y[, 1]; Y1_est = pred_y[, 1]
  Y2 = Y[, 2]; Y2_est = pred_y[, 2]
  
  #prop2 = (Y2 - Y2_est)/Y2; prop1 = (Y1 - Y1_est)/Y1
  #residual = sum(prop2^2) + sum(prop1^2)
  
  prop2 = (Y2 - Y2_est)*Population
  residual = sum(prop2^2)
  #residual = sum((Y2 - Y2_est)^2)
  
  return(residual)
}

pars = c('alpha', 'gamma')

SIR.sys = '
    dxdt[0] = -alpha*x[1]*x[0];
    dxdt[1] = alpha*x[1]*x[0] - gamma*x[1];
    '

compile_sys("SIR", SIR.sys, pars, TRUE)

WildSeg_Par = function(inval)
{
  pars = c('alpha', 'gamma')
  
  SIR.sys = '
    dxdt[0] = -alpha*x[1]*x[0];
    dxdt[1] = alpha*x[1]*x[0] - gamma*x[1];
    '
  
  compile_sys("SIR", SIR.sys, pars, TRUE)
  
  p = length(beta0)
  p1 = length(c(beta0, YO[1, ]))
  cusum_stat = 0
  cusum_pos = 0
  beta_diff = 0
  
  si = intervals[inval, 1]; ei = intervals[inval, 2]
  
  Xb = Inf; Zb = -Inf
  for (bi in (si + p1):(ei - p1 - 1))
  {
    Yb1 = YO[si:(bi-1), ]; Yb2 = YO[((bi+1):ei), ]
    Tb1 = TO[si:(bi-1)]; Tb2 = TO[(bi+1):ei]
    
    start1 = TO[si]; start2 = TO[bi]
    
    
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
    
    init1 = YO[si, ]; init2 = YO[bi, ]
    initl = c(1, 1)/Population; initu = c(1, 1)
    lower = c(betal, initl); upper = c(betau, initu)
    #lower2 = c(betal, initl); upper2 = c(betau, initu)
    
    optim_res1 = optim(c(beta0, init1), residual1_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
    optim_res2 = optim(c(beta0, init2), residual2_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
    
    betahat1 = optim_res1$par[1:p]; parahat2 = optim_res2$par
    start_value_hat = optim_res1$par[(p+1):p1]
    
    SIR_set_params(alpha = betahat1[1], gamma = betahat1[2])
    end_value_hat = as.vector(as.matrix(SIR_at(init = start_value_hat, times = c(start1, start2), start = start1))[-1, -1])
    
    parahat1 = c(betahat1, end_value_hat)
    residualhat1 = optim_res1$value; residualhat2 = optim_res2$value
    
    
    mi = ei - si
    
    scale0 = sqrt((bi - si)*(ei - bi)/mi)
    #scale1 = sqrt((ei - bi - p1)/(mi - 2*p1)/(bi - si - p1))
    #scale2 = sqrt((bi - si - p1)/(mi - 2*p1)/(ei - bi - p1))
    
    scale1 = scale2 = 1
    
    residual1 = scale1*residualhat1
    residual2 = scale2*residualhat2
    
    Xtemp = abs(residual1 + residual2)
    betadiff_temp = (parahat1 - parahat2)*c(1, 1, 0, 1)
    Ztemp = scale0*norm(betadiff_temp, '2')
    
    if (Xtemp < Xb)
    {
      Xb = Xtemp
      cusum_pos = bi
      cusum_stat = Xb
      #beta_diff[inval] = scale0*norm(parahat1 - parahat2, '2')
    }
    
    if (Ztemp > Zb)
    {
      Zb = Ztemp
      beta_diff = Zb
    }
    print(paste0('Candidate_Point', bi, '---finished'))
  }
  
  cusum_stat = Xb
  
  #save(list = c('cusum_pos', 'cusum_stat', 'beta_diff'), file = paste0('Covid19/Italy_inval_', inval, '.RData'))
  return(list(cusum_pos = cusum_pos, cusum_stat = cusum_stat, beta_diff = beta_diff))
}

aa = Sys.time()
cl = makeCluster(50)
registerDoParallel(cl)

bi_stat = foreach(inval = 1:nrow(intervals), .packages = c("odeintr", "numDeriv")) %dopar% WildSeg_Par(inval)

stopImplicitCluster()
stopCluster(cl)


beta_diff = NULL
candidate_cpt = NULL
candidate_stat = NULL
for (inval in 1:nrow(intervals))
{
  beta_diff = c(beta_diff, bi_stat[[inval]]$beta_diff)
  candidate_cpt = c(candidate_cpt, bi_stat[[inval]]$cusum_pos)
  candidate_stat = c(candidate_stat, bi_stat[[inval]]$cusum_stat)
}

s = 1; e = n; q.max = 15; gap = m/2
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
#thres_residual = sort(candidate_stat)[infornum]
thres_residual = max(candidate_stat)

index_diff = which(beta_diff > thres_beta)
intervals = intervals[index_diff, ]
candidate_cpt = candidate_cpt[index_diff]
beta_diff = beta_diff[index_diff]
candidate_stat = candidate_stat[index_diff]

#thres_residual = thres_residual[-1]

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

pars = c('alpha', 'gamma')

SIR.sys = '
    dxdt[0] = -alpha*x[1]*x[0];
    dxdt[1] = alpha*x[1]*x[0] - gamma*x[1];
    '

compile_sys("SIR", SIR.sys, pars, TRUE)

R_Mops = function(Y, Time, Omega, alpha, method = 'Candidate_Set', candidate_cpt, Candidate_Set = NULL, beta0 = NULL, betal = NULL, betau = NULL)
{
  n = length(Time); p = length(beta0); p1 = length(c(beta0, Y[1, ]))
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
    
    YO0 = Y[index0_O, ]; TO0 = Time[index0_O]
    YE0 = Y[index0_E, ]; TE0 = Time[index0_E]
    
    YO1 = Y[index1_O, ]; TO1 = Time[index1_O]
    YE1 = Y[index1_E, ]; TE1 = Time[index1_E]
    
    n0 = e0 - s0 + 1; n1 = e1 - s1 + 1
    startO0 = TO0[1]; startO1 = Time[inner]
    startE0 = TE0[1]; startE1 = Time[inner]
    
    initO0 = YO0[1, ]; initO1 = Y[inner, ]
    initE0 = YE0[1, ]; initE1 = Y[inner, ]
    
    if (method == 'L-BFGS-B')
    {
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
      
      initl = c(1, 1)/Population; initu = c(1, 1)
      lower = c(betal, initl); upper = c(betau, initu)
      
      optim_resO0 = optim(c(beta0, initO0), residualO0_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
      optim_resE0 = optim(c(beta0, initE0), residualE0_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
      
      optim_resO1 = optim(c(beta0, initO1), residualO1_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
      optim_resE1 = optim(c(beta0, initE1), residualE1_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
      
      betahatO0 = optim_resO0$par[1:p]; parahatO1 = optim_resO1$par
      start_value_hatO0 =  optim_resO0$par[(p+1):p1]
      
      SIR_set_params(alpha = betahatO0[1], gamma = betahatO0[2])
      end_value_hatO0 = as.vector(as.matrix(SIR_at(init = start_value_hatO0, times = c(startO0, startO1), start = startO0))[-1, -1])
      parahatO0 = c(betahatO0, end_value_hatO0)
      
      betahatE0 = optim_resE0$par[1:p]; parahatE1 = optim_resE1$par
      start_value_hatE0 =  optim_resE0$par[(p+1):p1]
      
      SIR_set_params(alpha = betahatE0[1], gamma = betahatE0[2])
      end_value_hatE0 = as.vector(as.matrix(SIR_at(init = start_value_hatE0, times = c(startE0, startE1), start = startE0))[-1, -1])
      parahatE0 = c(betahatE0, end_value_hatE0)
    }
    
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
  
  print(W); print(emp_level)
  
  return(Final_cpt)
}

# index_thres = max(which(num_cpts <= q.max))
# cpts_thresh = sort(cpts_list[[index_thres]])

cpts_thresh = sort(cpts_list[[1]])

cpts_select = cpts_thresh*2 - 1
if (min(cpts_select) <= 10) cpts_select = cpts_select[-1]
if (max(cpts_select) >= length(Time) - 10) cpts_select = cpts_select[-length(cpts_select)]
Omega = diag(c(1, 1, 0, 1))

Mops_Res = R_Mops(Y, Time, Omega = Omega, alpha = 0.2, method = 'L-BFGS-B', candidate_cpt = cpts_select, Candidate_Set = NULL, beta0 = beta0, betal = betal, betau = betau)
Final_cpts = Time[Mops_Res]

cpt_estimates = c(0, Mops_Res, length(Time))
partnum = length(cpt_estimates) - 1

Y_est = NULL
for (i in 1:partnum)
{
  sti = cpt_estimates[i] + 1; edi = cpt_estimates[i + 1]
  
  Yi = Y[sti:edi, ]; Ti = Time[sti:edi]; p = length(beta0)
  
  residuali_func = function(beta) 
  {
    q = length(beta)
    return(residual_func(Yi, Ti, Ti[1], beta[(p+1):q], beta[1:p]))
  }
  
  initl = c(1, 1)/Population; initu = c(1, 1); init0 = c(Yi[1, 1] - 0.05, Yi[1, 2])
  lower = c(betal, initl); upper = c(betau, initu)
  #lower2 = c(betal, initl); upper2 = c(betau, initu)
  
  optim_resi = optim(c(beta0, init0), residuali_func, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
  betai = optim_resi$par
  
  print(betai)
  
  ni = length(Ti)
  SIR_set_params(alpha = betai[1], gamma = betai[2])
  solve_yi = SIR_at(init = betai[3:4], times = Ti, start = Ti[1])
  pred_yi = solve_yi[, -1]
  
  Y_est = rbind(Y_est, pred_yi)
}


global_optimize = function(beta)
{
  q = length(beta)
  return(residual_func(Y, Time, Time[1], beta[(p+1):q], beta[1:p]))
}

init0 = c(Y[1, 1] - 0.1, Y[1, 2])
optim_resG = optim(c(beta0, init0), global_optimize, NULL, method = 'L-BFGS-B', lower = lower, upper = upper)
betaG = optim_resG$par
SIR_set_params(alpha = betaG[1], gamma = betaG[2])
solve_yG = SIR_at(init = betaG[3:4], times = Time, start = Time[1])
Y_G = solve_yG[, -1]

jpeg('SIR_estimates.jpg', width = 1200, height = 900)
plot(Time, Y[, 2], type = 'l', lwd = 3, main = 'Estimations of the Infection Curve', ylab = '', xlab = 'Time', cex.main = 3, ylim = c(0, max(Y_est[, 2], Y[, 2], Y_G[, 2])), 
     cex.lab = 3, cex.axis = 2.5)
lines(Time, Y_est[, 2], col = 'blue', lwd = 3)
lines(Time, Y_G[, 2], col = 'Green', lwd = 3)
points(Time[Mops_Res], Y[Mops_Res, 2], col = 'red', pch = 20, cex = 4)
dev.off()

jpeg(filename = 'covid19InfectedItaly.jpg', width = 1200, height = 900)
plot(Time, Y[, 2], lwd = 3, ylim = c(0, max(Y_est[, 2], Y[, 2], Y_G[, 2])), type = 'l',
     main = 'Dates of Structural Changes', ylab = '', xlab = 'Time', cex.main = 3, cex.lab = 3, cex.axis = 2.5)
text(x = Time[Final_cpts], y = Y[Final_cpts, 2] + 20000/Population, labels = ItalyTime[Final_cpts], cex = 2)
points(Time[Final_cpts], Y[Final_cpts, 2], col = 'red', cex = 4, pch = 20)
dev.off()


