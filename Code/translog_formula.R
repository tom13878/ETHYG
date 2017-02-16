#'=================================================================================================
#' Subject:  Function to make translog formula less onerous
#' Author:   Tom Morley
#' Contact:  tomas.morley@wur.nl
#' Output:   function for making translog formula
#'=================================================================================================

translog_form <- function(output="logyld", inputs=c("logN")){
  
  # get the level terms
  level_terms <- paste(inputs, collapse=" + ")
  
  # get the squared terms
  sqrd_terms <- paste(paste0("I(0.5 * ",inputs, "^2", ")"), collapse=" + ")
  
  # get the interactions
  int_terms <- t(combn(inputs, 2))
  int_terms <- paste(paste(int_terms[,1], int_terms[,2], sep="*"), collapse=" + ")
  
  # put all terms together
  all_terms <- paste(level_terms, sqrd_terms, int_terms, sep=" + ")
  
  # return a full formula
  paste0(output, " ~ ", all_terms)
}

##Example:
# mydata <- data.frame(y = rnorm(50),
# N = rnorm(50),
# L = rnorm(50),
# K = rnorm(50),
# A = rnorm(50),
# D1 = rbinom(50, 1, 0.5),
# D2 = rbinom(50, 1, 0.5),
# D3 = rbinom(50, 1, 0.5))
# inputs <- c("N", "L", "K", "A")
# output <- "y"
# form <- translog_form(output, inputs)
# lm(form, data = mydata)
# #number of interactions
# choose(4, 2) # 6
# # if we need to add more variables to the formula
# new_form <- paste(form, "D1", "D2", "D3", sep = " + ")
# lm(new_form, data = mydata)
