PMops = function(Y, Time, Omega, alpha, method = 'L-BFGS-B', candidate_cpt, beta0 = NULL, betal = NULL, betau = NULL)
{
  n = length(Time); p = length(beta0); p1 = length(c(beta0, Y[1,]))
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