# Function to extract tables with variables and coefficients from class "frontier".
# Tables can be used in Rmarkdown documents

sfaTable_f <- function(model, nc = 2){
  sfaTableOut <- list()
  
  rawTable <-as.data.frame(summary(model, extraPar = F)$mleParam)
  varnames <- rownames(rawTable)
  xvar <- varnames[1:(min(grep("^Z_", varnames))-1)]
  zvar <- varnames[grep("^Z_", varnames)]
  diag <- c("sigmaSq", "gamma")
  
  rawTable <- rawTable %>%
  mutate(variable = varnames) %>%
  rename(P = `Pr(>|z|)`) %>%
  mutate(sign = ifelse(P <= 0.001, "***",
                         ifelse(P <= 0.01, "**",
                                ifelse(P <= 0.05, "*", ""))))
    
  sfaTable <- data.frame(variable = varnames, stringsAsFactors=FALSE)
  sfaTable$Coef. <- sprintf("%.3f", round(rawTable$Estimate, nc))
  sfaTable$"Std. Error" <- sprintf("%.3f", round(rawTable$"Std. Error", nc))
  sfaTable$sign <- rawTable$sign
  
  xtable <- filter(sfaTable, variable %in% xvar) %>%  setNames(c("", "Coef.", "std. Error", ""))
  ztable  <- filter(sfaTable, variable %in% zvar) %>%  setNames(c("", "Coef.", "std. Error", ""))
  diagtable <- filter(sfaTable, variable %in% diag) %>%  setNames(c("", "Coef.", "std. Error", ""))
  

  sfaTableOut$xtable <- xtable
  sfaTableOut$ytable <- ztable
  sfaTableOut$diagtable <- diagtable
  return(sfaTableOut)
}
