# Section 1: Preliminary ----

rm(list=ls())
want <- c("RColorBrewer","plm","lfe","sandwich","lmtest","multiwayvcov","nlme","boot")
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)
lapply(want, function(i) require(i, character.only=TRUE))
rm(want, need)

# Section 2: Data Generation and Looping ----

# Generate new data
set.seed(12345)
r <- 10^3 # number of MC simulations
n_obs <- c(100, 200, 500, 1000, 10000) # number of observations
g_clu <- c(2, 4, 5, 10, 20, 25, 50, 100) # number of clusters
result <- list()

# Loop over simulations where data is clustered
for (n in n_obs) {
  for (g in g_clu) {
    out2 <- lapply(1:r, function(i) {
      
      print(i)
      
      # Generate data
      x1 <- rnorm(n) + rep(rnorm(g), each=n/g) 
      data <- data.frame(x1=x1)
      beta0 <- 1
      beta1 <- 1
      e1 <- rnorm(n) # unique variance for each obs
      e2 <- rep(rnorm(g), each=n/g)
      data$clu <- rep(1:g, each=n/g) # cluster id
      e <- e1 + e2
      data$y <- beta0 + beta1*data$x1 + e
      
      # Run model
      reg <- lm(y~x1, data) # OLS
      
      # Get coefficients
      b <- coef(reg)
      
      # Get SEs with different approaches
      se.canned <- sqrt(diag(vcov(reg)))
      se.clu    <- sqrt(diag(cluster.vcov(reg, data$clu))) # assumes errors are clustered
      
      # Export
      o <- c(b, se.canned, se.clu)
      names(o) <- paste(rep(c("b.ols","se.ols","se.ols.clu"), each=2), names(o))
      o
      
    })
    out2 <- do.call("rbind", out2)
    ratio <- mean(out2[,6]) / sd(out2[,2])
    result[[paste("n", n, "g", g, sep = "_")]] <- ratio  # Asked ChatGPT: How do I save my ratio calculation to the data frame
  }
}

data <- data.frame(n = rep(n_obs, each = length(g_clu)), g = rep(g_clu, times = length(n_obs)), ratio = unlist(result))

xvalue <- c(1, 2, 3, 4, 5, 6, 7, 8)
g_mapx <- c('2' = 1, '4' = 2, '5' = 3, '10' = 4, '20' = 5, '25' = 6, '50' = 7, '100' = 8) # Maps x values to maintain even spacing
data$g <- g_mapx[as.character(data$g)]
color <- c("lightgreen", "lightblue", "turquoise", "steelblue3", "purple4")

# Section 3: Plotting the Results ----

# Plot the graph of the results
png("Number_of_Clusters.png", width = 800, height = 600)
plot(xvalue, y = NULL, axes = F, pch = 19, type = "n", ylim = c(0,1),xlab = "Number of Clusters", ylab = "Std Err estimate / Std Dev of sampling distribution", main = "How many clusters do you need?")
axis(1, las=1, at=c(1, 2, 3, 4, 5, 6, 7, 8), lwd=1, labels=c("2", "4", "5", "10", "20", "25", "50", "100"))
axis(2, las=1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), lwd=1, labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1"))
for (i in 1:length(n_obs)) {
  lines(data$g[data$n == n_obs[i]], data$ratio[data$n == n_obs[i]], type = "l", col = color[i])
  points(data$g[data$n == n_obs[i]], data$ratio[data$n == n_obs[i]], pch = 19, col = color[i])
}
abline(h = 1, lty = 3)
box()
# Asked ChatGPT: How do I put a box around my legend?
legend("bottomright", legend = c("100 observations", "200 observations", "500 observations", "1000 observations", "10000 observations"), col = c("lightgreen", "lightblue", "turquoise", "steelblue3", "purple4"), pch = 19, lty = 1, lwd = 2, title = "Number of observations")

dev.off()

# The End
