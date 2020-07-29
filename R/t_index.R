get_random_population <- function(aoi, values, n_pop){
  out <- raster::sampleRandom(aoi, n_pop, sp= T)
  out <- raster::extract(x=values, y=out, sp=T)
  return(out)
}


kdistr <- function(x, eval.points, ...){
  pdf <- pdfCluster::kepdf(x, eval.points=eval.points, ...)
  pdf@estimate <-  pdf@estimate/sum(pdf@estimate)
  my_pdf <- data.frame(x=pdf@eval.points, y=pdf@estimate)
  return(my_pdf)
}


pkdistr <- function(x_target, x, eval.points = NULL, ...){
  if(is.null(eval.points)) eval.points <- seq(min(x)-0.05, max(x)+0.05, length.out = 100)
  
  x_target <- abs(x_target)
  
  eval.points <-  sort(c(eval.points, -x_target, x_target))
  pdf <- pdfCluster::kepdf(x, eval.points=eval.points, ...)
  pdf@estimate <-  pdf@estimate/sum(pdf@estimate)
  
  p_lhs <- sum(pdf@estimate[1:which(pdf@eval.points == -x_target)])
  p_rhs <- sum(pdf@estimate[which(pdf@eval.points == x_target):length(eval.points)])
  my_p <- sum(p_lhs, p_rhs)#min(p_lhs, p_rhs)
  
  return(my_p)
}


get_critical_IB_values <- function(random_IBs, min_IB, max_IB, length_out = 100, c_values = c(0.32, 0.05)){
  
  my_sequence <- seq(min_IB,
                     max_IB,
                     length.out = length_out)
  
  
  out_list <- list()
  for(c_value in c_values){
    neg_bound <- data.frame(x = my_sequence,
                            y = unlist(lapply(my_sequence,
                                              function(x, xrdm, eval.points)
                                                pkdistr(x_target=x, x=xrdm, eval.points=eval.points),
                                              xrdm=random_IBs, 
                                              eval.points = my_sequence))) %>% 
      dplyr::filter(x<0) %>% 
      mutate(tdfiff = abs(c_value-y)) %>% 
      dplyr::filter(tdfiff == min(tdfiff)) 
    
    pos_bound <- data.frame(x = my_sequence,
                            y = unlist(lapply(my_sequence,
                                              function(x, xrdm, eval.points)
                                                pkdistr(x_target=x, x=xrdm, eval.points=eval.points),
                                              xrdm=random_IBs, 
                                              eval.points = my_sequence))) %>% 
      dplyr::filter(x>0) %>% 
      mutate(tdfiff = abs(c_value-y)) %>% 
      dplyr::filter(tdfiff == min(tdfiff)) 
    out_list[[as.character(c_value)]] = data.frame(c_value = c_value, IB_lower = neg_bound$x, IB_upper = pos_bound$x)
  }
  return(dplyr::bind_rows(out_list))
}


get_range <- function(random_IBs, unknown_IBs){
  my_max <- ifelse(max(c(random_IBs, unknown_IBs)) >= 0.95,
                   1,
                   max(c(random_IBs, unknown_IBs))+0.05)
  my_min <- ifelse(max(c(random_IBs, unknown_IBs)) <= -0.95,
                   -1,
                   min(c(random_IBs, unknown_IBs))-0.05)
  my_range <- list(min=my_min, max=my_max)
  return(my_range)
}