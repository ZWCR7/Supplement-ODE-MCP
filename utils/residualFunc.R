residual_func = function(Y, Time, start, init, beta)
{
  n = length(Time); Time1 = unique(c(start, Time))
  ODE_set_params(a = 0.2, b = beta, c = 3)
  solve_y = ODE_at(init = init, times = Time1, start = start)
  
  if (length(Time1) == (n + 1)) pred_y = solve_y[-1, -1]
  if (length(Time1) == n) pred_y = solve_y[, -1]
  
  residual = sum((Y - pred_y)^2)
  
  return(residual)
}