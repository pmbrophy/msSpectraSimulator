#' A parallel and memory efficient way to Calculate multi-gaussian peaks 
#'
#' @param x a vector of x-coordinates from which the corresponding y-coordinates
#'   are calculated
#' @param mus a vector of means along the x-axis
#' @param sigmas a vector of standard deviations for each mu
#' @param probDensity Should the function produce a probability density function
#'   `TRUE` or a gaussian peak `FALSE` with amplitude k? default is `TRUE`.
#' @param returnComponentPks Should the function return each component peak.
#'   Default is FALSE - this can generate huge matricies.
#' @param ks Amplitude of the peak. Only used when `probDensity == FALSE`
#' @param peakIDs labels for each component peak. Only used when
#'   returnComponentPks == TRUE.
#'
#' @return a matrix
#' @export
#'
#' @examples
#' 

simulateMultiPeak_parallel_gaussian <- function(x, mus, sigmas, probDensity, ks, peakIDs){
  if(any(probDensity)){
    ks <- 1
  }
  
  #split into groups for parallel execution
  ncores <- BiocParallel::bpworkers()
  mus_list <- splitIndex(index = mus,
                         nGroups = ncores, 
                         randomize = FALSE)[[1]]
  sigmas_list <- splitIndex(index = sigmas,
                            nGroups = ncores, 
                            randomize = FALSE)[[1]]
  ks_list <- splitIndex(index = ks,
                            nGroups = ncores, 
                            randomize = FALSE)[[1]]
  
  #Do parallel
  ys <- BiocParallel::bpmapply(FUN = batch_SummedPeak,
                               mus = mus_list,
                               sigmas = sigmas_list,
                               ks = ks_list,
                               MoreArgs = list(x = x, 
                                               probDensity = probDensity))
  
  #Output
  y <- rowSums(ys)
  pks <- cbind("x" = x, "peak_sum" = y)
   
  pks
}

batch_SummedPeak <- function(x, mus, sigmas, probDensity, ks){
  #Required function
  simulatePeak_gaussian <- function(x, mu, sigma, probDensity = TRUE, k){
    if(probDensity){
      k <- (1/(sigma*sqrt(2*pi)))
    }else{
      if(missing(k)){
        stop("k not specified")
      }
    }
    
    peak <- k * exp(-0.5 * ((x-mu)/sigma)^2)
    
    peak
  }
  npeaks <- length(mus)
  length_x <- length(x)
  
  #Allocate
  y <- vector(mode = "numeric", length = length_x)
  
  i <- 1
  while(i <= npeaks){
    nth_peak <- simulatePeak_gaussian(x = x, 
                                      mu = mus[i], 
                                      sigma = sigmas[i], 
                                      probDensity = probDensity, 
                                      k = ks[i])
    
    y <- y + nth_peak
    i <- i + 1
  }
  #Return
  y
}


#' A memory efficient way to Calculate multi-gaussian peaks 
#'
#' @param x a vector of x-coordinates from which the corresponding y-coordinates
#'   are calculated
#' @param mus a vector of means along the x-axis
#' @param sigmas a vector of standard deviations for each mu
#' @param probDensity Should the function produce a probability density function
#'   `TRUE` or a gaussian peak `FALSE` with amplitude k? default is `TRUE`.
#' @param returnComponentPks Should the function return each component peak.
#'   Default is FALSE - this can generate huge matricies.
#' @param ks Amplitude of the peak. Only used when `probDensity == FALSE`
#' @param peakIDs labels for each component peak. Only used when
#'   returnComponentPks == TRUE.
#'
#' @return a matrix
#' @export
#'
#' @examples
#' 

simulateMultiPeak_gaussian <- function(x, mus, sigmas, probDensity, returnComponentPks = FALSE, ks, peakIDs){
  if(any(probDensity)){
    ks <- 1
  }
  #Loop setup
  npeaks <- length(mus)
  length_x <- length(x)
  i <- 1
  
  print(paste("Generating", npeaks,"peaks over xvector of length", length_x))
  
  #Allocate
  y <- vector(mode = "numeric", length = length_x)
  if(returnComponentPks){
    pks <- matrix(data = 0, nrow = length_x, ncol = npeaks)
  }
  
  #Memory efficient implementation - calculate each peak, then add to pre-allocated vector
  while(i <= npeaks){
    nth_peak <- simulatePeak_gaussian(x = x, 
                                      mu = mus[i], 
                                      sigma = sigmas[i], 
                                      probDensity = probDensity, 
                                      k = ks[i])
    #save component peak
    if(returnComponentPks){
      pks[,i] <- nth_peak
    }
    
    y <- y + nth_peak
    i <- i + 1
    
    #Print update
    if((i %% 100) == 0){
      print(i)
    }
  }
  
  #Output
  if(returnComponentPks){
    pks <- cbind(x, y, pks)
    colnames(pks) <- c("x", "peak_sum", paste0("peak_", as.character(peakIDs)))
  }else{
    pks <- cbind("x" = x, "peak_sum" = y)
  } 
  pks
}

#' Calculate Normal Distribution or Gaussian Peak
#'
#' Calculates either a normal distribution similar to dnorm() with an
#' integrated area of 1 or a gaussian peak with amplitude k
#'
#' @param x a vector of x-coordinates from which the corresponding y-coordinates
#'   are calculated
#' @param mu the mean
#' @param sigma the standard deviation
#' @param probDensity Should the function produce a probability density function
#'   `TRUE` or a gaussian peak `FALSE` with amplitude k? default is `TRUE`.
#' @param k Amplitude of the peak. Only used when `probDensity == FALSE`
#'
#' @return a vector of y-coordinates the same length as x
#' @export
#'
#' @examples
#' #normal distribution
#' xVec <- seq(from = 1, to = 100, by = 0.1)
#' pdensity <- func_gaussian(x = xVec, mu = 10, sigma = 1, probDensity = TRUE)
#' p1 <- plot(x = xVec, y = pdensity)
#'
#' #gaussian peak
#' gpeak <- func_gaussian(x = xVec, mu = 10, sigma = 1, probDensity = FALSE, k = 10)
#' p2 <- plot(x = xVec, y = gpeak)
#'

simulatePeak_gaussian <- function(x, mu, sigma, probDensity = TRUE, k){
  if(probDensity){
    k <- (1/(sigma*sqrt(2*pi)))
  }else{
    if(missing(k)){
      stop("k not specified")
    }
  }
  
  peak <- k * exp(-0.5 * ((x-mu)/sigma)^2)
  
  peak
}