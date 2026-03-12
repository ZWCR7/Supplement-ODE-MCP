Stat_Generator = function(Y, Time, intervals, method = 'L-BFGS-B', beta0 = NULL, betal = NULL, betau = NULL)
{
  num = nrow(intervals); p = length(beta0); p1 = length(c(beta0, Y[1,]))
  
  print(paste0('Number of intervals:', num))
  
  cusum_stat = rep(0, num)
  cusum_pos = rep(0, num)
  beta_diff = rep(0, num)
  
  cusum_statG = rep(0, num)
  cusum_posG = rep(0, num)
  
  for (inval in 1:num)
  {
    aa = Sys.time()
    
    si = intervals[inval, 1]; ei = intervals[inval, 2]
    
    Yinval = Y[si:ei, ]; Tinval = Time[si:ei]
    startG = Time[si]
    
    residualG_func = function(beta) 
    {
      q = length(beta)
      return(residual_func(Yinval, Tinval, startG, beta[(p+1):q], beta[1:p]))
    }
    
    initG = Y[si, ]
    initlG = c(-3, -3); inituG = c(3, 3)
    lowerG = c(betal, initlG); upperG = c(betau, inituG)
    
    optim_resG = optim(c(beta0, initG), residualG_func, NULL, method = method, lower = lowerG, upper = upperG)
    residualG = optim_resG$value
    
    Gb = -Inf; Xb = Inf; Zb = -Inf
    for (bi in (si + p1):(ei - p1 - 1))
    {
      Yb1 = Y[si:(bi-1), ]; Yb2 = Y[(bi+1):ei, ]
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
      
      init1 = Y[si, ]; init2 = Y[bi, ]
      initl = c(-3, -3); initu = c(3, 3)
      lower = c(betal, initl); upper = c(betau, initu)
      #lower2 = c(betal, initl); upper2 = c(betau, initu)
      
      optim_res1 = optim(c(beta0, init1), residual1_func, NULL, method = method, lower = lower, upper = upper)
      optim_res2 = optim(c(beta0, init2), residual2_func, NULL, method = method, lower = lower, upper = upper)
      
      betahat1 = optim_res1$par[1:p]; parahat2 = optim_res2$par
      start_value_hat =  optim_res1$par[(p+1):p1]
      
      ODE_set_params(a = 0.2, b = betahat1, c = 3)
      end_value_hat = as.vector(as.matrix(ODE_at(init = start_value_hat, times = c(start1, start2), start = start1))[-1, -1])
      
      parahat1 = c(betahat1, end_value_hat)
      residualhat1 = optim_res1$value; residualhat2 = optim_res2$value
      
      Gtemp = - residualhat1 - residualhat2 + residualG
      if (Gtemp > Gb)
      {
        Gb = Gtemp
        cusum_posG[inval] = bi
        cusum_statG[inval] = Gb
      }
      
      mi = ei - si
      
      scale0 = sqrt((bi - si)*(ei - bi)/mi)
      # scale1 = scale2 = 1
      # 
      # residual1 = scale1*residualhat1
      # residual2 = scale2*residualhat2
      
      Xtemp = abs(residualhat1 + residualhat2)
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
  
  return(list(Y = Y, Time = Time, candidate_cpt = cusum_pos, candidate_stat = cusum_stat, beta_diff = beta_diff, candidate_cptG = cusum_posG, candidate_statG = cusum_statG, intervals = intervals))
}