simulateSpectrum_parallel <- function(formulaVector, chargeStates, mean, sd, ppmTol, resolution, specAbundMin = 0.0001, xResolution = 0.001, isoCalc_abundanceLimit = 0.01){
  #Calculate the charge state distribution from the formula vector
  chargeStateDist <- chargeState_distribution(formulaVector = formulaVector,
                                              chargeStates = chargeStates,
                                              mean = mean,
                                              sd = sd, 
                                              isoCalc_abundanceLimit = isoCalc_abundanceLimit)
  
  #Apply grouping algorithm
  chargeStateDist <- mzDataTable::indexMasterSpectrum(mzDt = chargeStateDist, ppmTol = ppmTol, isCentroid = TRUE)
  
  #Sum Groups
  chargeStateDist[, sumProb := sum(charge_prob), by = mzGrid_index]
  chargeStateDist[, sumProb := sumProb/max(sumProb)]
  
  spec <- unique(chargeStateDist[, list(mzGrid, mzGrid_index, sumProb, charge)])
  
  #Apply threshold and re-index
  spec <- spec[sumProb > specAbundMin,]
  data.table::setorder(x = spec, mzGrid_index)
  spec[, mzGrid_index := c(1:nrow(spec))]
  
  #Inputs for multi_gaussian
  mz_start <- floor(min(spec$mzGrid) - 10)
  mz_end <- ceiling(max(spec$mzGrid) + 10)
  x <- seq(from = mz_start, to = mz_end, by = xResolution)
  
  FWHM <- spec$mzGrid/resolution
  sigmas <- FWHM/(2*sqrt(2*log(2)))
  
  profiles <- simulateMultiPeak_parallel_gaussian(x = x, 
                                                  ks = spec$sumProb,
                                                  mus = spec$mzGrid,
                                                  sigmas = sigmas,
                                                  peakIDs =  spec$mzGrid_index,
                                                  probDensity = FALSE)
  
  #profiles <- data.table("x" = x, "y" = profiles)
  profiles <- as.data.table(profiles)
  
  list(peakProfiles = profiles, peakInfo = spec)
}




simulateSpectrum <- function(formulaVector, chargeStates, mean, sd, ppmTol, resolution, specAbundMin = 0.0001, xResolution = 0.001, returnComponentPks = FALSE, isoCalc_abundanceLimit = 0.01){
  #Calculate the charge state distribution from the formula vector
  chargeStateDist <- chargeState_distribution(formulaVector = formulaVector,
                                              chargeStates = chargeStates,
                                              mean = mean,
                                              sd = sd, 
                                              isoCalc_abundanceLimit = isoCalc_abundanceLimit)
  
  #Apply grouping algorithm
  chargeStateDist <- mzDataTable::indexMasterSpectrum(mzDt = chargeStateDist, ppmTol = ppmTol, isCentroid = TRUE)
  
  #Sum Groups
  chargeStateDist[, sumProb := sum(charge_prob), by = mzGrid_index]
  chargeStateDist[, sumProb := sumProb/max(sumProb)]
  
  spec <- unique(chargeStateDist[, list(mzGrid, mzGrid_index, sumProb, charge)])
  
  #Apply threshold and re-index
  spec <- spec[sumProb > specAbundMin,]
  data.table::setorder(x = spec, mzGrid_index)
  spec[, mzGrid_index := c(1:nrow(spec))]
  
  #Inputs for multi_gaussian
  mz_start <- floor(min(spec$mzGrid) - 10)
  mz_end <- ceiling(max(spec$mzGrid) + 10)
  x <- seq(from = mz_start, to = mz_end, by = xResolution)
  
  FWHM <- spec$mzGrid/resolution
  sigmas <- FWHM/(2*sqrt(2*log(2)))
  
  profiles <- simulateMultiPeak_gaussian(x = x, 
                                                ks = spec$sumProb,
                                                mus = spec$mzGrid,
                                                sigmas = sigmas,
                                                returnComponentPks = returnComponentPks,
                                                peakIDs =  spec$mzGrid_index,
                                                probDensity = FALSE)
  
  #profiles <- data.table("x" = x, "y" = profiles)
  profiles <- as.data.table(profiles)
  
  list(peakProfiles = profiles, peakInfo = spec)
}


#' Title
#'
#' @param formulaVector
#' @param chargeStates
#' @param mean
#' @param sd
#' @param offset
#'
#' @return
#' @export
#'
#' @examples
#'

chargeState_distribution <- function(formulaVector, chargeStates, mean, sd, isoCalc_abundanceLimit){
  #Calculate mass and natural abundance of each ion type
  ionList <- chargeState_mz(formulaVector = formulaVector,
                            chargeStates = chargeStates,
                            abundanceLimit = isoCalc_abundanceLimit)

  #multiply natural distribution by assumed charge state distribution
  chargeStateProbs <- dnorm(x = chargeStates, mean = mean, sd = sd)
  chargeStateProbs <- chargeStateProbs/max(chargeStateProbs)

  chargeStateDist <- mapply(FUN = IsoSpecR_adjustProb,
                            IsoSpecR_matrix = ionList,
                            chargeProbability = chargeStateProbs,
                            SIMPLIFY = FALSE)

  #Convert to a data.table
  chargeStateDist <- data.table::as.data.table(do.call(rbind, chargeStateDist))

  chargeStateDist
}

#' Adjust natural isotope probability by charge state distribution
#'
#' @param IsoSpecR_matrix the mass/prob/(mz) matrix returned by IsoSpecR::IsoSpecify()
#' @param chargeProbability probability of observing this chargestate - a value <= 1
#'
#' @return the same matrix with an extra column "charge_prob"
#' @export
#'
#' @examples
#'

IsoSpecR_adjustProb <- function(IsoSpecR_matrix, chargeProbability){
  charge_prob <- IsoSpecR_matrix[,"prob"] * chargeProbability

  m <- cbind(IsoSpecR_matrix, "charge_prob" = charge_prob)

  m
}

#' Calculate the m/z isotopic distributions
#'
#' @param formulaVector a named character vector
#' @param chargeStates the charge state
#' @param abundanceLimit 1-abundance limit = stopCondition for IsoSpecify()
#'
#' @return a matrix for each charge state containing mass, prob, and mz
#' @export
#'
#' @examples
#'

chargeState_mz <- function(formulaVector, chargeStates, abundanceLimit){
  #Get formula for each charge state
  ionList <- chargestate_formulas(formulaVector = formulaVector, chargeStates = chargeStates)

  #Get isotopes, calculate mz, normalize prob to 1
  stopCondition <- 1 - abundanceLimit
  massList <- mapply(FUN = IsoSpecR_call,
                     molecule = ionList, 
                     chargeState = chargeStates, 
                     MoreArgs = list(stopCondition = stopCondition), 
                     SIMPLIFY = FALSE)
  
  #Get mass list for each charge state
  #massList <- lapply(X = ionList,
  #                   FUN = IsoSpecR::IsoSpecify,
  #                   stopCondition = abundanceLimit)

  #Calculate m/z
  #massList <- mapply(FUN = IsoSpecR_mass2mz,
  #                   IsoSpecR_matrix = massList,
  #                   chargeState = chargeStates,
  #                   SIMPLIFY = FALSE)

  names(massList) <- as.character(chargeStates)
  massList
}




#' Call IsoSpecR and 1) add mz to matrix, 2) normalize prob
#'
#' @param IsoSpecR_matrix the mass/prob matrix returned by IsoSpecR::IsoSpecify()
#' @param chargeState the charge state
#'
#' @return the same matrix with an extra column "mz"
#' @export
#'
#' @examples
#'

IsoSpecR_call <- function(molecule, stopCondition, chargeState){
  IsoSpecR_matrix <- IsoSpecR::IsoSpecify(molecule = molecule, stopCondition = stopCondition)
  
  #normalize
  pmax <- max(IsoSpecR_matrix[,"prob"])
  IsoSpecR_matrix[,"prob"] <- IsoSpecR_matrix[,"prob"]/pmax
  
  #Calculate m/z
  mz <- IsoSpecR_matrix[,"mass"] / abs(chargeState)
  IsoSpecR_matrix <- cbind(IsoSpecR_matrix, "mz" = mz, "charge" = chargeState)
  
  IsoSpecR_matrix
}

#' Calculate ion formulas from neutral formulaVector
#'
#' @param formulaVector a named character vector
#' @param chargeStates an integer vector (positive or negative)
#'
#' @return a named list of formula vectors
#' @export
#'
#' @examples
#'

chargestate_formulas <- function(formulaVector, chargeStates){

  #Generate list of charge states
  if(all(chargeStates < 0)){
    #Negative
    ionList <- lapply(X = abs(chargeStates),
                      FUN = ionize_neg,
                      formulaVector = formulaVector)

  }else if(all(chargeStates > 0)){
    #Positive
    ionList <- lapply(X = chargeStates,
                      FUN = ionize_pos,
                      formulaVector = formulaVector)

  }else{
    stop("Unallowed chargeStates requested")
  }

  #Set names
  names(ionList) <- as.character(chargeStates)

  ionList
}

#' Ionize formula vector by removing hydrogen
#'
#' @param formulaVector a named character vector
#'
#' @return a named character vector
#' @export
#'
#' @examples
#'

ionize_neg <- function(formulaVector, chargeState){
  formulaVector["H"] <- formulaVector["H"] - chargeState

  formulaVector
}

#' Ionize formula vector by adding hydrogen
#'
#' @param formulaVector a named character vector
#'
#' @return a named character vector
#' @export
#'
#' @examples
#'

ionize_pos <- function(formulaVector, chargeState){
  formulaVector["H"] <- formulaVector["H"] - chargeState

  formulaVector
}
