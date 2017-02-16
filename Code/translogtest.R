# Finding optimum for a translog function
# Example 2.1 from Jauregui and Sain (1992), p. 11
df <- data.frame(NX = c(0,30,60, 90, 120, 150, 0, 30, 60, 90, 120, 150),
                 Y = c(360, 2886, 3066, 3999, 3891, 3171, 1020, 2421, 3255, 3732, 3489, 3465))

df$N <- df$NX+13 # 13 is the N in the soil
df$N2 <- df$N*df$N
df$logY <- log(df$Y)
df$logN <- log(df$N)
df$logN2 <- df$logN*df$logN

plot(Y~N, data = df)
lines(df$N, predict(OLS), col = "red")

# OLS
OLS <- lm(Y ~ N + N2, data =df)
summary(OLS)

# Translog
TRA <- lm(logY ~logN + logN2, data = df)
summary(TRA)

# Funtion with linearlised translog and MPP (first derivative)

MPP_f <- function(N){
  logY <-coef(TRA)[1]+coef(TRA)[2]*log(N) + coef(TRA)[3]*log(N)*log(N)
  Y = exp(logY)
  MPP = (coef(TRA)[2]+2*coef(TRA)[3]*log(N))*(Y/N)-9.6 # 9.6 is the N/maize price ratio.
  return(MPP)
}

# Find point where MPP = 0
MPP <- uniroot(MPP_f, c(10, 500))

# Solution is 82.4, which is the same as 69 (82.4 -13) on p. 46
MPP[[1]]-13

# Optimum/maximum of translog function can also be found.
translog_f <- function(N){
  logY <-coef(TRA)[1]+coef(TRA)[2]*log(N) + coef(TRA)[3]*log(N)*log(N)
  Y = exp(logY)
  return(Y)
}

library(optimx)
optimx(90, translog)

