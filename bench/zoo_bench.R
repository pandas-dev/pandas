library(zoo)
library(xts)

indices = rep(NA, 100000)
for (i in 1:100000)
  indices[i] <- paste(sample(letters, 10), collapse="")

timings <- numeric()

## x <- zoo(rnorm(100000), indices)
## y <- zoo(rnorm(90000), indices[sample(1:100000, 90000)])

## indices <- as.POSIXct(1:100000)

indices <- as.POSIXct(Sys.Date()) + 1:100000

x <- xts(rnorm(100000), indices)
y <- xts(rnorm(90000), indices[sample(1:100000, 90000)])

for (i in 1:10) {
  gc()
  timings[i] = system.time(x + y)[3]
}

mean(timings)
