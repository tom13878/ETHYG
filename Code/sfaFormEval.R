#'=================================================================================================
#' Subject:  Function to get a sfa formula which
#'           can be evaluated
#' Author:   Tom Morley
#' Contact:  tomas.morley@wur.nl
#' Output:   function for making translog formula
#'=================================================================================================

sfaFormEval <- function(modl){
  
  # get the form of the frontier equation
  form <- strsplit(as.character(formula(modl)), "~ ")[[1]][2]
  
  # combine coefs with variables
  form <- unlist(strsplit(form, " + ", fixed=TRUE))
  coef <- coef(modl)[2:(length(coef(modl))-2)]
  names(coef) <- gsub(":", "*", names(coef))
  coef <- coef[form]
  
  # paste together
  form <- paste(paste(coef,
                      form, sep = " * "), collapse=" + ")
  paste(coef(modl)[1], form, sep=" + ")
}
