set_s = 5; n = 2000

FDR = matrix(0, 2, set_s); TPR = matrix(0, 2, set_s)
CP_NUM = matrix(0, 2, set_s); Dist = matrix(0, 2, set_s)

for (set in 1:set_s)
{
  load(paste0('Results_FN/Mops-FN-Structural_n = ', n, '_K = 10_noise sigma_N', set, '.RData'))
  #load(paste0('Results_FN/Mops-FN-Structural_n = ', n, '_K = 10_noise sigma_T', set, '.RData'))
  #load(paste0('Results_FN/Mops-FN-General_n = ', n, '_K = 10_noise sigma_N', set, '.RData'))
  #load(paste0('Results_FN/Mops-FN-General_n = ', n, '_K = 10_noise sigma_T', set, '.RData'))
  
  nsimu = length(CP_RES)
  #th = sigma*sqrt(8*log(n))
  
  infornum = matrix(0, 2, nsimu); uninfornum = matrix(0, 2, nsimu)
  distance = matrix(0, 2, nsimu); cpnum = matrix(0, 2, nsimu)
  
  for (i in 1:nsimu)
  {
    for (k in 1:2)
    {
      est_cp = CP_RES[[i]][[k]]
      
      if (length(est_cp) > 0)
      {
        cpnum[k, i] = length(est_cp)
        
        work_cp1 = c(0, est_cp); work_cp2 = c(est_cp, n)
        
        for (j in 1:length(est_cp))
        {
          lower = ceiling((work_cp1[j] + work_cp2[j])/2)
          upper = ceiling((work_cp1[j+1] + work_cp2[j+1])/2)
          
          index_info = ((true_cp >= lower) & (true_cp < upper))
          
          if (sum(index_info != 0))
          {
            infornum[k, i] = infornum[k, i] + 1
          }
          else
          {
            uninfornum[k, i] = uninfornum[k, i] + 1
          }
        }
        
        
        for (j in 1:length(true_cp))
        {
          distance[k, i] = distance[k, i] + min(abs(est_cp - true_cp[j]))
        }
        
        distance[k, i] = distance[k, i]/length(true_cp)
      }
    }
    
  }
  
  fdr_s = uninfornum/cpnum; tpr_s = infornum/length(true_cp)
  
  FDR[, set] = apply(fdr_s, 1, sum, na.rm = T)/nsimu; TPR[, set] = apply(tpr_s, 1, sum, na.rm = T)/nsimu
  CP_NUM[, set] = apply(cpnum, 1, sum, na.rm = T)/nsimu; Dist[, set] = apply(distance, 1, mean, na.rm = T)
}

gamma = c(0.5, 0.7, 0.9, 1.1, 1.3)
jpeg(filename = 'FN-Systematic-N_FDR.jpg', width = 400, height = 300) #Change the name file and main name for different settings
plot(gamma, FDR[1, ], type = 'l', lwd = 3, col = 1, xlab = 'gamma', ylab = '', main = 'FDR-Normal', cex.lab = 2, cex.main = 3, ylim = c(0, 0.2))
lines(gamma, FDR[2, ], type = 'l', lwd = 3, col = 2)
dev.off()

gamma = c(0.5, 0.7, 0.9, 1.1, 1.3)
jpeg(filename = 'FN-Systematic-N_TPR.jpg', width = 400, height = 300) #Change the name file and main name for different settings
plot(gamma, TPR[1, ], type = 'l', lwd = 3, col = 1, xlab = 'gamma', ylab = '', main = 'TPR-Normal', cex.lab = 2, cex.main = 3, ylim = c(0.8, 1))
lines(gamma, TPR[2, ], type = 'l', lwd = 3, col = 2)
legend('bottomleft', legend = c('Proposed', 'GLR'), cex = 1.5, lwd = 3, col = c(1, 2))
dev.off()

###########################################################################################################

rm(list = ls())
set_s = 5; n = 2000

FDR = matrix(0, 11, set_s); TPR = matrix(0, 11, set_s)
CP_NUM = matrix(0, 11, set_s); Dist = matrix(0, 11, set_s)

for (set in 1:set_s)
{
  load(paste0('Results_FN/Mops-FN-Compare_n = ', n, '_K = ', 10, '_noise sigma_N', set, '.RData'))
  
  nsimu = length(CP_RES)
  #th = sigma*sqrt(8*log(n))
  
  infornum = matrix(0, 11, nsimu); uninfornum = matrix(0, 11, nsimu)
  distance = matrix(0, 11, nsimu); cpnum = matrix(0, 11, nsimu)
  
  for (i in 1:nsimu)
  {
    for (k in 1:11)
    {
      est_cp = CP_RES[[i]][[k]]
      
      if (k >= 9)
      {
        est_cp = CP_RES[[i]][[k]]$intervals$midpoints
      }
      
      if (sum(is.na(est_cp)) > 0)
      {
        est_cp = NULL
      }
      
      if (length(est_cp) > 0)
      {
        cpnum[k, i] = length(est_cp)
        
        work_cp1 = c(0, est_cp); work_cp2 = c(est_cp, n)
        
        for (j in 1:length(est_cp))
        {
          lower = ceiling((work_cp1[j] + work_cp2[j])/2)
          upper = ceiling((work_cp1[j+1] + work_cp2[j+1])/2)
          
          index_info = ((true_cp >= lower) & (true_cp < upper))
          
          if (sum(index_info != 0))
          {
            infornum[k, i] = infornum[k, i] + 1
          }
          else
          {
            uninfornum[k, i] = uninfornum[k, i] + 1
          }
        }
        
        
        for (j in 1:length(true_cp))
        {
          distance[k, i] = distance[k, i] + min(abs(est_cp - true_cp[j]))
        }
        
        distance[k, i] = distance[k, i]/length(true_cp)
      }
    }
    
  }
  
  fdr_s = uninfornum/cpnum; tpr_s = infornum/length(true_cp)
  
  FDR[, set] = apply(fdr_s, 1, sum, na.rm = T)/nsimu; TPR[, set] = apply(tpr_s, 1, sum, na.rm = T)/nsimu
  CP_NUM[, set] = apply(cpnum, 1, sum, na.rm = T)/nsimu; Dist[, set] = apply(distance, 1, mean, na.rm = T)
}

#############################################################################################################################

set_s = 3; n = 2000

FDR = matrix(0, 5, set_s); TPR = matrix(0, 5, set_s)
CP_NUM = matrix(0, 5, set_s); Dist = matrix(0, 5, set_s)

for (set in 1:set_s)
{
  load(paste0('Results_Lin/Linear-Compare_n =', n, '_noise sigma_N', set, '.RData'))
  #load(paste0('Results_Lin/Linear-Compare_n =', n, '_noise sigma_T', set, '.RData'))
  
  nsimu = length(CP_RES)
  #th = sigma*sqrt(8*log(n))
  
  infornum = matrix(0, 5, nsimu); uninfornum = matrix(0, 5, nsimu)
  distance = matrix(NA, 5, nsimu); cpnum = matrix(0, 5, nsimu)
  
  for (i in 1:nsimu)
  {
    for (k in 1:5)
    {
      est_cp = CP_RES[[i]][[k]]
      
      if (k >= 4)
      {
        est_cp = CP_RES[[i]][[k]]$intervals$midpoints
      }
      
      
      if (length(est_cp) > 0)
      {
        cpnum[k, i] = length(est_cp)
        
        work_cp1 = c(0, est_cp); work_cp2 = c(est_cp, n)
        
        for (j in 1:length(est_cp))
        {
          lower = ceiling((work_cp1[j] + work_cp2[j])/2)
          upper = ceiling((work_cp1[j+1] + work_cp2[j+1])/2)
          
          # lower = max(est_cp[j] - 100, work_cp1[j])
          # upper = min(est_cp[j] + 100, work_cp2[j + 1])
          
          index_info = ((true_cp >= lower) & (true_cp <= upper))
          
          if (sum(index_info != 0))
          {
            infornum[k, i] = infornum[k, i] + 1
          }
          else
          {
            uninfornum[k, i] = uninfornum[k, i] + 1
          }
        }
        
        distance[k, i] = 0
        for (j in 1:length(true_cp))
        {
          distance[k, i] = distance[k, i] + min(abs(est_cp - true_cp[j]))
        }
        
        distance[k, i] = distance[k, i]/length(true_cp)
      }
    }
    
  }
  
  fdr_s = uninfornum/cpnum; tpr_s = infornum/length(true_cp); pa_s = (tpr_s == 1)
  
  FDR[, set] = apply(fdr_s, 1, sum, na.rm = T)/nsimu; TPR[, set] = apply(tpr_s, 1, sum, na.rm = T)/nsimu
  CP_NUM[, set] = apply(cpnum, 1, sum, na.rm = T)/nsimu; Dist[, set] = apply(distance, 1, mean, na.rm = T)
  Pa = apply(pa_s, 1, sum, na.rm = T)/nsimu
}

##############################################################################################################################################

set_s = 3; n = 2000

FDR = matrix(0, 9, set_s); TPR = matrix(0, 9, set_s)
CP_NUM = matrix(0, 9, set_s); Dist = matrix(0, 9, set_s)

sdFDR = matrix(0, 9, set_s); sdTPR = matrix(0, 9, set_s)
sdCP_NUM = matrix(0, 9, set_s)

for (set in 1:set_s)
{
  load(paste0('Results_Lin/Regression-Compare_n =', n, '_noise sigma_N', set, '.RData'))
  #load(paste0('Results_Lin/Regression-Compare_n =', n, '_noise sigma_T', set, '.RData'))
  
  nsimu = length(CP_RES)
  #th = sigma*sqrt(8*log(n))
  
  infornum = matrix(0, 9, nsimu); uninfornum = matrix(0, 9, nsimu)
  distance = matrix(NA, 9, nsimu); cpnum = matrix(0, 9, nsimu)
  
  for (i in 1:nsimu)
  {
    for (k in 1:9)
    {
      est_cp = CP_RES[[i]][[k]]
      
      if (k >= 8)
      {
        est_cp = CP_RES[[i]][[k]]$intervals$midpoints
      }
      
      if (sum(is.na(est_cp)) > 0)
      {
        est_cp = NULL
      }
      
      if (length(est_cp) > 0)
      {
        cpnum[k, i] = length(est_cp)
        
        work_cp1 = c(0, est_cp); work_cp2 = c(est_cp, n)
        
        for (j in 1:length(est_cp))
        {
          lower = ceiling((work_cp1[j] + work_cp2[j])/2)
          upper = ceiling((work_cp1[j+1] + work_cp2[j+1])/2)
          
          # lower = max(est_cp[j] - 100, work_cp1[j])
          # upper = min(est_cp[j] + 100, work_cp2[j + 1])
          
          index_info = ((true_cp >= lower) & (true_cp <= upper))
          
          if (sum(index_info != 0))
          {
            infornum[k, i] = infornum[k, i] + 1
          }
          else
          {
            uninfornum[k, i] = uninfornum[k, i] + 1
          }
        }
        
        distance[k, i] = 0
        for (j in 1:length(true_cp))
        {
          distance[k, i] = distance[k, i] + min(abs(est_cp - true_cp[j]))
        }
        
        distance[k, i] = distance[k, i]/length(true_cp)
      }
    }
    
  }
  
  fdr_s = uninfornum/cpnum; tpr_s = infornum/length(true_cp)
  
  FDR[, set] = apply(fdr_s, 1, sum, na.rm = T)/nsimu; TPR[, set] = apply(tpr_s, 1, sum, na.rm = T)/nsimu
  CP_NUM[, set] = apply(cpnum, 1, sum, na.rm = T)/nsimu; Dist[, set] = apply(distance, 1, mean, na.rm = T)
  
  sdFDR[, set] = apply(fdr_s, 1, sd, na.rm = T); sdTPR[, set] = apply(tpr_s, 1, sd, na.rm = T)
  sdCP_NUM[, set] = apply(cpnum, 1, sd, na.rm = T)
}
