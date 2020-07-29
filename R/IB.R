
#' Compute the Ib index of a sample given a population.
#'
#' Computes the spatial balance of a sample with the Ib index proposed by Tille et al. (2018). As the formula can generate large distance matrices, this function uses the \code{"bigstatsr"} and the \code{"bigdist"} packages.
#'
#' @param population is a matrix representing the population being considered
#' @param samplepoints is a matrix representing sample points from the population
#' @param tmp_dir (optional) indicates a temporary folder where to save the bigmatrices created in the process. If none is provided, the default temporary folder is used.
#' @param tmp_file (optional) indicates a temporary folder where to save the bigmatrices created in the process. If none is provided, a default temporary folder is created.
#' 
#' @details 
#' This index is developped by Tillé et al. (2018) and measure the spreading of a sample drawn from a population. It uses a corrected version of the traditional Moran's I index. Each row of the matrix \bf W should represents a stratum. Each stratum is defined by a particular unit and its neighbouring units. See wpik. The spatial balance measure is equal to
#' I_B =\frac{( \bf s- \bar{s}_w)^\top W ( s- \bar{s}_w)}{\bf √{( s- \bar{s}_w)^\top D ( s- \bar{s}_w) ( s- \bar{s}_w)^\top B ( s- \bar{s}_w)}},
#' where \bf D is the diagonal matrix containing the w_i,
#'\bf \bar{s}_w = 1 \frac{ s^\top W 1}{ 1^\top W 1}
#' and
#' \bf B = W^\top D^{-1} W - \frac{ W^\top 1 1^\top W}{1^\top W 1}.
#'
#' This measure  is more favourable than the original Moran's I index since it is bounded and allows one to discriminate, in absolute terms, between the absence and presence of spatial balance in a sample. The lower bound (-1) indicates maximum spatial balance and the upper bound (+1) indicates a poor spatial balance.
#' 
#'@references Tillé, Y., Dickson, M. M., Espa, G., & Giuliani, D. (2018). Measuring the spatial balance of a sample: A new measure based on Moran’s I index. Spatial Statistics, 23, 182-192.
#'
#' @return IB
#'
#' @export

bigIB <- function(population, samplepoints, tmp_dir, tmp_file, n_cpus=1){
  
  library(bigstatsr)
  library(bigdist)
  library(prodlim)
  
  if(!is.data.frame(population)) population <- as.data.frame(population)
  if(!is.data.frame(samplepoints)) samplepoints <- as.data.frame(samplepoints)
  if(missing(tmp_dir)){
    tmp_dir <- tempdir()
  } else {
    if(!dir.exists(tmp_dir)) stop("tmp_dir does not exist")
  }
  if(missing(tmp_file)){
    tmp_file <- tempfile()
  } else {
    tmp_file <- file.path(tmp_dir, tmp_file)
    if(file.exists(tmp_file)) file.remove(tmp_file)
  }
  
  #Define function to get ww matrix from a fbm
  get_ww_from_bigM <- function(fbm, ind, pik, x){
    rr <- rank(fbm[ind,],ties.method="min")
    ww <- as.integer(rr<=ceiling(1/pik[ind]))
    dec<- as.integer(rr==max(rr*ww))
    if(sum(dec)>0)  ww[dec==1]=ww[dec==1]*(1/pik[ind]-sum(ww-dec))/sum(dec)
    x[ind, ] <- ww
    x[ind, ind]  <- 0
    ww
  }
  
  tmp_fileW <- tempfile()
  
  
  # compute distance matrix using bigdist
  distm <- bigdist::bigdist(mat = as.matrix(population), file = tmp_file)
  
  #Set up vectors and matrices for wfromdpik and IB calculation
  U <- rep(1,nrow(population))
  s <- rep(0, nrow(population)) 
  s[prodlim::row.match(x=samplepoints, table=population)] = 1
  
  pik <- rep(1/2, bigdist::bigdist_size(distm))
  N <- length(pik)
  w <- bigstatsr::FBM(bigdist::bigdist_size(distm), bigdist::bigdist_size(distm), 
                      backingfile = tmp_fileW)
  
  #Compute wfromdpik matrix
  output_nok <- bigstatsr::big_apply(X=distm$fbm, 
                                     a.FUN = get_ww_from_bigM,
                                     ind = bigstatsr::rows_along(distm$fbm),
                                     ncores=n_cpus,
                                     pik = pik,
                                     x = w,
                                     block.size=1)
  #Compute IB
  w <- as.matrix(w[, ])
  wi <- rowSums(w)
  Yb <- sum(wi*s)/sum(w)
  z  <- s-Yb
  B  <- crossprod(w,diag(1/wi))%*%w-(crossprod(w,U))%*%t(U)%*%w/sum(w)
  IB <- c(crossprod(z,w)%*%z  / sqrt(crossprod(z,diag(wi))%*%z * crossprod(z,B)%*%z))
  
  #Remove tmp files
  lapply(list.files(tmp_file), function(x) file.remove(x))
  lapply(list.files(tmp_fileW), function(x) file.remove(x))
  
  #Return IB
  return(IB)
}	