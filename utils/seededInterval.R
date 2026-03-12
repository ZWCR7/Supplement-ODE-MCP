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