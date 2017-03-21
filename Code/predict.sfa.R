predict.sfa <- function(sfa.model, data){
  betas <- sfa.model$olsParam[-length(sfa.model$olsParam)]
  form <- paste("~", paste(names(betas)[-1], collapse=" + "), " ")
  M <- model.matrix(formula(form), data=data)
  M %*% betas
}

