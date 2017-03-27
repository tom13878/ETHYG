predict.sfa <- function(sfa.model, data){
  data <- as.data.frame(data)
  beta_names <- names(sfa.model$olsParam[-length(sfa.model$olsParam)])
  betas <- sfa.model$mleParam[beta_names]
  form <- paste("~", paste(names(betas)[-1], collapse=" + "), " ")
  M <- model.matrix(formula(form), data=data)
  as.numeric(M %*% betas)
}

# check if output is the same as in the
# frontier.fitted function
# head(fitted(sf11x9, asInData = FALSE))
# head(predict.sfa(sf11x9, db1))
