#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  relationships between community variables
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   pca output 
#'==============================================================================

# get packages
library(pacman)
p_load(char=c("rprojroot", "dplyr"), install=TRUE)

# get root path
root <- find_root(is_rstudio_project)

# read in data to run models on
db1 <- readRDS(file.path(root, "Cache/db1.rds"))
db1 <- unique(db1)
db1 <- filter(db1, lab < 2000, logN >= 0)
db1$bank <- ifelse(db1$bank %in% "NO", 0, 1)

# variables we want to cluster
pc.vars <- model.matrix(~ -1 + dist_market + popEA  +
                              cost2small_town +
                          cost2large_town + extension, data=db1) %>% unique

ir.pca <- prcomp(pc.vars[, -ncol(pc.vars)],
                 center = TRUE,
                 scale. = TRUE)

# print method
print(ir.pca)

# summary
summary(ir.pca)

# plot method - no real elbow
plot(ir.pca, type = "l")

# we can do a biplot
biplot(ir.pca)

# or a 3d plot
library(rgl)
plot3d(ir.pca$x[,1:3], col=pc.vars[, ncol(pc.vars)] + 1)
