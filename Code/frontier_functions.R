# remember it is just the core variables we want
# in the new x

get_mat <- function(sfa_model, inputs, var, stat){
  mat <- sfa_model$dataTable
  df <- data.frame(mat[, 4:ncol(mat)])
  names(df) <- c("Intercept", dimnames(mat)[[2]][5:ncol(mat)])
  df <- df[, names(df) %in% c("Intercept", inputs)]
  
  # get stat of each variable
  stats <- apply(df, 2, function(col) summary(col)[stat])
  new_mat <- matrix(rep(stats, nrow(mat)), nrow=nrow(mat), ncol=ncol(df), byrow=TRUE)
  new_df <- data.frame(new_mat)
  names(new_df) <- names(df)
  new_df[, var] <- df[, var]
  
  # return new_df
  new_df
}

get_prediction <- function(sfa_model, inputs, var, stat){
  
  # get the form of the frontier equation
  form <- strsplit(as.character(formula(sfa_model)), "~ ")[[1]][2]
  form <- paste("Intercept", form, sep=" + ")
  
  # combine coefs with variables
  form <- unlist(strsplit(form, " + ", fixed=TRUE))
  form <- paste(paste(coef(sfa_model)[1:length(form)],
                       form, sep = " * "), collapse=" + ")
  
  # get the matrix
  X <- get_mat(sfa_model, inputs=inputs, var=var, stat=stat)
  
  
  # evaluate the sum in the formula
  with(df, eval(parse(text=form)))
}

get_mat(sfa_CD_basic, inputs = inputs, var="logN", stat=6)
ypred <- get_prediction(sfa_CD_basic, inputs = c("logN", "logarea", "loglab", "SPEI"), var="logN", stat=6)
