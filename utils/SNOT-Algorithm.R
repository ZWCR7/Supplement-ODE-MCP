SNot = function(bi_stat, thres_beta = NULL, thres_GLR = NULL, q.max = NULL, method = "upper.bound", Statistic = 'ER', gap = 10)
{
  Y = bi_stat$Y; Time = bi_stat$Time
  
  n = length(Time)
  s = 1; e = n
  
  intervals = bi_stat$intervals
  
  if (Statistic == 'ER')
  {
    candidate_cpt = bi_stat$candidate_cpt
    candidate_stat = bi_stat$candidate_stat
    beta_diff = bi_stat$beta_diff
    
    if (method == 'solution.path')
    {
      thres_beta = 0
      
      cpts_list = list()
      cpts_stat_list = list()
      beta_diff_list = list()
      num_cpts = c()
      
      index_diff = which(beta_diff > thres_beta)
      index_beta = length(index_diff)
      
      intervals_temp = intervals[index_diff, ]
      candidate_cpt_temp = candidate_cpt[index_diff]
      beta_diff_temp = beta_diff[index_diff]
      candidate_stat_temp = candidate_stat[index_diff]
      
      while(index_beta > 1)
      {
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
        
        num_cpts = c(num_cpts, length(cpts))
        cpts_list = c(cpts_list, list(cpts))
        cpts_stat_list = c(cpts_stat_list, list(cpts_stat))
        beta_diff_list = c(beta_diff_list, list(cpts_diff))
        
        thres_beta = min(cpts_diff)
        
        index_diff = which(beta_diff > thres_beta)
        index_beta = length(index_diff)
        
        intervals_temp = intervals[index_diff, ]
        candidate_cpt_temp = candidate_cpt[index_diff]
        beta_diff_temp = beta_diff[index_diff]
        candidate_stat_temp = candidate_stat[index_diff]
      }
    }
    else
    {
      if (method == 'upper.bound')
      {
        if (is.null(q.max)) {return('Error: the upper.bound method needs to specify the parameter q.max')}
        
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
      
      if (method == 'fixed')
      {
        if (is.null(thres_beta)) {return('Error: the fixed method needs to specify the parameter thres_beta')}
        
        thres_beta = thres_beta
      }
      
      cpts_list = list()
      cpts_stat_list = list()
      beta_diff_list = list()
      num_cpts = rep(0, length(thres_beta))
      
      for (i in 1:length(thres_beta))
      {
        index_diff = which(beta_diff > thres_beta[i])
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
    }
    
    return(list(cpts_list = cpts_list, cpts_stat_list = cpts_stat_list, beta_diff_list = beta_diff_list, bi_stat = bi_stat, num_cpts = num_cpts))
  }
  
  if (Statistic == 'GLR')
  {
    candidate_cpt = bi_stat$candidate_cptG
    candidate_stat = bi_stat$candidate_statG
    
    if (method == 'solution.path')
    {
      cpts_list = list()
      cpts_stat_list = list()
      num_cpts = c()
      
      thres_GLR = 0
      index_GLR = which(candidate_stat > thres_GLR)
      index_beta = length(index_GLR)
      
      intervals_temp = intervals[index_GLR, ]
      candidate_cpt_temp = candidate_cpt[index_GLR]
      candidate_stat_temp = candidate_stat[index_GLR]
      
      while(index_beta > 1)
      {
        cpts = c(); cpts_stat = c(); cpts_diff = c()
        
        index_stop = 9999; not_num = 1
        while(index_stop > 2)
        {
          print(paste0('Not Step:', not_num))
          
          interval_lengths = intervals_temp[, 2] - intervals_temp[, 1]
          
          min_length = min(interval_lengths)
          length_candidate_index = which(interval_lengths <= min_length + 2)
          
          maxdiff = max(candidate_stat_temp[length_candidate_index])
          index_cpt = which(candidate_stat_temp == maxdiff)
          stat = candidate_stat_temp[index_cpt]
          
          cpt = unique(candidate_cpt_temp[index_cpt])
          
          Belong_Func = function(x) return((x[1] >= cpt + gap) + (x[2] <= cpt - gap))
          Belong_index = apply(intervals_temp, 1, Belong_Func)
          
          index_belong = which(Belong_index == 1)
          intervals_temp = as.matrix(intervals_temp[index_belong, ])
          
          candidate_cpt_temp = candidate_cpt_temp[index_belong]
          candidate_stat_temp = candidate_stat_temp[index_belong]
          
          cpts = c(cpts, cpt)
          cpts_stat = c(cpts_stat, stat)
          cpts_diff = c(cpts_diff, diff)
          
          index_stop = length(intervals_temp)
          
          not_num = not_num + 1
        }
        
        num_cpts = c(num_cpts, length(cpts))
        cpts_list = c(cpts_list, list(cpts))
        cpts_stat_list = c(cpts_stat_list, list(cpts_stat))
        
        thres_GLR = min(cpts_stat)
        index_GLR = which(candidate_stat > thres_GLR)
        index_beta = length(index_GLR)
        
        intervals_temp = intervals[index_GLR, ]
        candidate_cpt_temp = candidate_cpt[index_GLR]
        candidate_stat_temp = candidate_stat[index_GLR]
      }
    }
    else
    {
      if (method == 'upper.bound')
      {
        if (is.null(q.max)) {return('Error: the upper.bound method needs to specify the parameter q.max')}
        
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
        
        thres_GLR = sort(candidate_stat, decreasing = T)[infornum]
      }
      
      if (method == 'fixed')
      {
        if (is.null(thres_GLR)) {return('Error: the fixed method needs to specify the parameter thres_GLR')}
        thres_GLR = thres_GLR
      }
      
      cpts_list = list()
      cpts_stat_list = list()
      num_cpts = rep(0, length(thres_GLR))
      
      for (i in 1:length(thres_GLR))
      {
        index_GLR = which(candidate_stat >= thres_GLR[i])
        intervals_temp = intervals[index_GLR, ]
        candidate_cpt_temp = candidate_cpt[index_GLR]
        candidate_stat_temp = candidate_stat[index_GLR]
        
        cpts = c(); cpts_stat = c(); cpts_diff = c()
        
        index_stop = 9999; not_num = 1
        while(index_stop > 2)
        {
          print(paste0('Not Step:', not_num))
          
          interval_lengths = intervals_temp[, 2] - intervals_temp[, 1]
          
          min_length = min(interval_lengths)
          length_candidate_index = which(interval_lengths <= min_length + 2)
          
          maxdiff = max(candidate_stat_temp[length_candidate_index])
          index_cpt = which(candidate_stat_temp == maxdiff)
          stat = candidate_stat_temp[index_cpt]
          
          cpt = unique(candidate_cpt_temp[index_cpt])
          
          Belong_Func = function(x) return((x[1] >= cpt + gap) + (x[2] <= cpt - gap))
          Belong_index = apply(intervals_temp, 1, Belong_Func)
          
          index_belong = which(Belong_index == 1)
          intervals_temp = as.matrix(intervals_temp[index_belong, ])
          
          candidate_cpt_temp = candidate_cpt_temp[index_belong]
          candidate_stat_temp = candidate_stat_temp[index_belong]
          
          cpts = c(cpts, cpt)
          cpts_stat = c(cpts_stat, stat)
          cpts_diff = c(cpts_diff, diff)
          
          index_stop = length(intervals_temp)
          
          not_num = not_num + 1
        }
        
        num_cpts[i] = length(cpts)
        cpts_list = c(cpts_list, list(cpts))
        cpts_stat_list = c(cpts_stat_list, list(cpts_stat))
      }
    }
    
    return(list(cpts_list = cpts_list, cpts_stat_list = cpts_stat_list, thres_GLR = thres_GLR, bi_stat = bi_stat, num_cpts = num_cpts))
  }
}