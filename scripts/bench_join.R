library(xts)

iterations <- 50

ns = c(100, 1000, 10000, 100000, 1000000)
kinds = c("outer", "left", "inner")

result = matrix(0, nrow=3, ncol=length(ns))
n <- 100000
pct.overlap <- 0.2

k <- 1

for (ni in 1:length(ns)){
 n <- ns[ni]
 rng1 <- 1:n
 offset <- as.integer(n * pct.overlap)
 rng2 <- rng1 + offset
 x <- xts(matrix(rnorm(n * k), nrow=n, ncol=k),
          as.POSIXct(Sys.Date()) + rng1)
 y <- xts(matrix(rnorm(n * k), nrow=n, ncol=k),
          as.POSIXct(Sys.Date()) + rng2)
 timing <- numeric()
 for (i in 1:3) {
     kind = kinds[i]
     for(j in 1:iterations) {
       gc()  # just to be sure
       timing[j] <- system.time(merge(x,y,join=kind))[3]
     }
     #timing <- system.time(for (j in 1:iterations) merge.xts(x, y, join=kind),
     #                      gcFirst=F)
     #timing <- as.list(timing)
     result[i, ni] <- mean(timing) * 1000
     #result[i, ni] = (timing$elapsed / iterations) * 1000
   }
}

rownames(result) <- kinds
colnames(result) <- log10(ns)

mat <- matrix(rnorm(500000), nrow=100000, ncol=5)
set.seed(12345)
indexer <- sample(1:100000)

timing <- rep(0, 10)
for (i in 1:10) {
  gc()
  timing[i] = system.time(mat[indexer,])[3]
}

