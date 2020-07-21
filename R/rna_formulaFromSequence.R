#' Calculate the chemical formula from an RNA sequence
#'
#' @param sequence a character vector of length 1
#'
#' @return a named character vector
#' @export
#'
#' @examples
#' formula <- rna_formulaFromSequence(sequence = "AAAGCCUG")
#'
rna_formulaFromSequence <- function(sequence){
  #AMP: C10H14N5O7P
  #GMP: C10H14N5O8P
  #CMP: C9H14N3O8P
  #UMP: C9H13N2O9P
  formulas <- c("A" = "C10H14N5O7P",
                "G" = "C10H14N5O8P",
                "C" = "C9H14N3O8P",
                "U" = "C9H13N2O9P")

  sequence <- strsplit(sequence, "")[[1]]

  #Check for bad characters
  containsAllowedCharacters <- !all(sequence %in% c("A", "G", "C", "U"))
  if(containsAllowedCharacters){
    stop("sequence contains characters other than those that have been implemented")
  }

  nbases <- length(sequence)

  sequenceFormula <- formulas[sequence]
  sequenceFormula <- paste(sequenceFormula, collapse = "")

  #Translate formula string to named vector
  formulaVector <- formula_vector(sequenceFormula)

  #Remove OH
  formulaVector["O"] <- formulaVector["O"] - (nbases - 1)
  formulaVector["H"] <- formulaVector["H"] - (nbases - 1)

  formulaVector
}

#' Modify a formula vector
#'
#' @param basis_vec a formula vector to be modified
#' @param mods a list of lists containing named vectors \"add\" and \"subtract\"
#'
#' @return a list of formula vectors 
#' @export
#'
#' @examples
#' mods <- list(
#'   list(add = c("Na" = 1), 
#'        subtract = c("H" = 1)), 
#'   list(add = c("Na" = 2), 
#'        subtract = c("H" = 2)))
#' basis <- c("C" = 10, "H" = 16, "O" = 2)
#' modify_formula_vector(basis_vec = basis, mods = mods)
#' 
 
modify_formula_vector <- function(basis_vec, mods){
  unmodified <- basis_vec
  #check mods 
  if(!is.list(mods)){
    stop("mods must be a list containing lists containing vectors named \"add\" and \"subtract\"")
  }
  
  #Add 
  mods_add <- lapply(X = mods, "[[", "add")
  basis_vec <- lapply(X = mods_add, FUN = add_formula_vectors, basis_vec = basis_vec)
  
  #Subtract
  mods_subtract <- lapply(X = mods, "[[", "subtract")
  basis_vec <- mapply(FUN = subtract_formula_vectors, 
                      basis_vec = basis_vec, 
                      sub_vec = mods_subtract, 
                      SIMPLIFY = FALSE)
  
  basis_vec["basis"] <- list(unmodified)
  basis_vec
}

#' Add atoms to a formula vector
#'
#' @param basis_vec the vector to be modified
#' @param add_vec the elements to add
#'
#' @return basis_vec with additional elements in add_vec 
#' @export
#'
#' @examples
#' 

add_formula_vectors <- function(basis_vec, add_vec){
  
  elements2 <- names(add_vec)
  elements1 <- names(basis_vec)
  
  #Add empty cells to vector
  elements2add <- elements2[which(!(elements2 %in% elements1))]
  if(length(elements2add) > 0){
    basis_vec[elements2add] <- rep(x = 0, times = length(elements2add))
  }
  
  #Add values of add_vec to basis_vec by name
  basis_vec[elements2] <- basis_vec[elements2] + add_vec[elements2]
  basis_vec
}

#' Subtract atoms from a formula vector
#'
#' @param basis_vec the vector to be modified
#' @param sub_vec the elements to be removed
#'
#' @return basis_vec with elements sub_vec removed
#' @export
#'
#' @examples
#' 

subtract_formula_vectors <- function(basis_vec, sub_vec){
  
  elements2 <- names(sub_vec)
  elements1 <- names(basis_vec)
  
  #Add empty cells to vector
  elements2add <- elements2[which(!(elements2 %in% elements1))]
  if(length(elements2add) > 0){
    basis_vec[elements2add] <- rep(x = 0, times = length(elements2add))
  }
  
  #subtract values of add_vec to basis_vec by name
  basis_vec[elements2] <- basis_vec[elements2] - sub_vec[elements2]
  basis_vec
}

#' Get named elemental composition vector 
#' 
#' @details Adopted from OrgMassSpecR::ListFormula() to provide a named vector with the count of each element.
#' 
#' @param elemental.formula a molecular formula without parentheses 
#'
#' @return named vector with the count of each element
#' @export
#'

formula_vector <- function (elemental.formula) {
  elements <- c("C", "H", "N", "O", 
                "S", "P", "Br", "Cl", "F", 
                "I", "Si", "Sn", "K", "Na", 
                "Li", "Rb", "Cs")
  
  x <- checkAtoms(elemental.formula, elements = elements)

  
  result <- sapply(FUN = GetAtoms, 
                   X = elements, 
                   elemental.formula = elemental.formula)
  
  return(result[which(result > 0)])
}

#' Calculate element frequency
#'
#' @param elemental.formula the elemental formula to be converted
#' @param element a character vector of length one with a one or two letter element code
#'
#' @return
#' @export
#'
#' @examples
#' elementVector <- sapply(FUN = GetAtoms, 
#'                         X = c("C", "H", "O"), 
#'                         elemental.formula = "CH3COOH")
#'

GetAtoms <- function(elemental.formula, element) {
  reg.exp <- paste(element, "[[:digit:]]*(?![[:lower:]])", 
                   sep = "")
  x <- gregexpr(reg.exp, elemental.formula, perl = TRUE)
  if (x[[1]][1] != -1) {
    n <- vector(mode = "numeric", length = length(x[[1]]))
    for (i in 1:length(x[[1]])) {
      y <- attr(x[[1]], which = "match.length")[i]
      z <- substr(elemental.formula, x[[1]][i], x[[1]][i] + 
                    y - 1)
      number <- as.numeric(strsplit(z, split = element)[[1]][2])
      if (is.na(number)) {
        n[i] <- 1
      }
      else {
        n[i] <- number
      }
      atoms <- sum(n)
    }
  }
  else {
    atoms <- 0
  }
  return(atoms)
}

#' Internal function for formula_vector() 
#'
#' @param elemental.formula the elemental formula to be converted
#' @param elements a vector of allowed elements
#'
#' @return gregexpr() vector
#'

checkAtoms <- function(elemental.formula, elements){
  
  atoms <- gregexpr("[[:upper:]][[:lower:]]{0,1}", elemental.formula)
  
  for (i in 1:length(atoms[[1]])) {
    y <- attr(atoms[[1]], which = "match.length")[i]
    atom <- substr(x = elemental.formula, 
                   start = atoms[[1]][i], 
                   stop = atoms[[1]][i] + y - 1)
    
    
    if(!(atom %in% elements)){
      stop(paste("Elemental formula", elemental.formula, "contains element not of:", paste(elements, collapse = ", ")))
    } 
  }
  atoms
}
