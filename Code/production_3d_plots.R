# translog production function plot

# set project root
library(rprojroot)
root <- find_root(is_rstudio_project)


# get data
db1 <- readRDS("Cache/db1.rds")
db1 <- unique(db1)
db1 <- db1[db1$lab < 3000, ]


#https://stat.ethz.ch/pipermail/r-help/2003-September/039104.html
surf.colors <- function(x, col = terrain.colors(20)) {
  
  # First we drop the 'borders' and average the facet corners
  # we need (nx - 1)(ny - 1) facet colours!
  x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
              x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
  
  # Now we construct the actual colours matrix
  colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
  
  return(colors)
}

# translog
library(frontier)
modl <- sfa(logyld ~ logN + loglab + I(logN^2) + I(loglab^2) + loglab:logN, data=db1)


N <- seq(0, 700, length.out=100) 
lab <- seq(0, 3000, length.out=100) 
f <- function(N, lab){
  logyld <- coef(modl)[1] + coef(modl)[2]*log(N + 1) + coef(modl)[3]*log(lab + 1) +
    coef(modl)[4]*I(log(N + 1)^2) + coef(modl)[5]*I(log(lab + 1)^2) +
    coef(modl)[6] * log(N + 1)*log(lab + 1)
  return(exp(logyld))
}

yld <- outer(N, lab, f)

png(file.path(root, "translog.png"))
persp(N, lab, yld, main="translog production function",
      xlab="Nitrogen", ylab="Labour", zlab="Yield",
      theta = -25, phi = 0, expand = 0.7, col=surf.colors(yld, col = heat.colors(40, alpha=0.7)))
dev.off()

# -------------------------------------
# cobb douglass 
library(frontier)
modl <- sfa(logyld ~ logN + loglab, data=db1)
f <- function(N, lab){
  logyld <- coef(modl)[1] + coef(modl)[2]*log(N + 1) + coef(modl)[3]*log(lab + 1)
  return(exp(logyld))
}

yld <- outer(N, lab, f)

png(file.path(root, "CD.png"))
persp(N, lab, yld, main="Cobb Douglass production function",
      xlab="Nitrogen", ylab="Labour", zlab="Yield",
      theta = -25, phi = 0, expand = 0.7, col=surf.colors(yld, col = heat.colors(40, alpha=0.7)))
dev.off()

