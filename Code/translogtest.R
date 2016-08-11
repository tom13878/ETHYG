# Testing non-linear optimization
df <- data.frame(NX = c(0,30,60, 90, 120, 150, 0, 30, 60, 90, 120, 150),
                 Y = c(360, 2886, 3066, 3999, 3891, 3171, 1020, 2421, 3255, 3732, 3489, 3465))
df$N <- df$NX+13
df$N2 <- df$N*df$N
df$logY <- log(df$Y)
df$logN <- log(df$N)
df$logN2 <- df$logN*df$logN

plot(Y~N, data = df)
lines(df$N, predict(OLS), col = "red")


OLS <- lm(Y ~ N + N2, data =df)
summary(OLS)

TRA <- lm(logY ~logN + logN2, data = df)
summary(TRA)
df$MPP <- (coef(TRA)[2]+2*coef(TRA)[3]*df$logN)*(df$Y/(df$N))

x<-0
translog <- function(N){
  logY <-coef(TRA)[1]+coef(TRA)[2]*log(N) + coef(TRA)[3]*log(N)*log(N)
  Y = exp(logY)
  MPP = (coef(TRA)[2]+2*coef(TRA)[3]*log(N))*(Y/N)-9.6
  return(MPP)
  #return(y)
  }

translog2 <- function(N){
  logY <-coef(TRA)[1]+coef(TRA)[2]*log(N) + coef(TRA)[3]*log(N)*log(N)
  Y = -exp(logY)
  return(Y)
  #return(y)
}

df$MPP <- with(df, (coef(TRA)[2]+2*coef(TRA)[3]*log(N))*(Y/N))
exp(-7.552312)


translog2(110)


x<-10
coef(TRA)[1]*x^(coef(TRA)[2])*exp(coef(TRA)[3]*log(x))

coef(TRA)

exp(translog(log(10)))

coef(TRA)[1]+coef(TRA)[2]*log(x) + coef(TRA)[3]*(log(x)^2)


exp(coef(TRA)[1])
translog(-7.552312)
summary(TRA)

nlm(translog, 80, hessian=TRUE)
uniroot(translog, c(10, 500))

library(optimx)
optimx(90, translog2)
optimize(90, translog2)
