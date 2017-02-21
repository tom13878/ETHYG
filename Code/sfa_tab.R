


sfa_table <- function(sfa_model, name="model"){
  res <- summary(sfa_model)
  res <- round(res$mleParam, 3)
  n <- row.names(res)
  p <- cut(res[, 4], breaks=c(0, 0.01, 0.05, 1), labels=c("**", "*", ""), include.lowest = TRUE)
  res <- paste(sprintf("%.3f" ,res[,1]), " (", sprintf("%.3f" ,res[, 2]), ")", p, sep="")
  res <- data.frame(parameter=n, res)
  names(res) <- c("parameter", name)
  res
}

