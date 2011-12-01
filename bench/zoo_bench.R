library(zoo)
library(xts)
library(fts)
library(tseries)
library(its)
library(xtable)

## indices = rep(NA, 100000)
## for (i in 1:100000)
##   indices[i] <- paste(sample(letters, 10), collapse="")



## x <- zoo(rnorm(100000), indices)
## y <- zoo(rnorm(90000), indices[sample(1:100000, 90000)])

## indices <- as.POSIXct(1:100000)

indices <- as.POSIXct(Sys.Date()) + seq(1, 100000000, 100)

sz <- 500000

## x <- xts(rnorm(sz), sample(indices, sz))
## y <- xts(rnorm(sz), sample(indices, sz))

zoo.bench <- function(){
    x <- zoo(rnorm(sz), sample(indices, sz))
    y <- zoo(rnorm(sz), sample(indices, sz))
    timeit(function() {x + y})
}

xts.bench <- function(){
    x <- xts(rnorm(sz), sample(indices, sz))
    y <- xts(rnorm(sz), sample(indices, sz))
    timeit(function() {x + y})
}

fts.bench <- function(){
    x <- fts(rnorm(sz), sort(sample(indices, sz)))
    y <- fts(rnorm(sz), sort(sample(indices, sz))
    timeit(function() {x + y})
}

its.bench <- function(){
    x <- its(rnorm(sz), sort(sample(indices, sz)))
    y <- its(rnorm(sz), sort(sample(indices, sz)))
    timeit(function() {x + y})
}

irts.bench <- function(){
    x <- irts(sort(sample(indices, sz)), rnorm(sz))
    y <- irts(sort(sample(indices, sz)), rnorm(sz))
    timeit(function() {x + y})
}

timeit <- function(f){
  timings <- numeric()
  for (i in 1:10) {
    gc()
    timings[i] = system.time(f())[3]
  }
  mean(timings)
}

bench <- function(){
  results <- c(xts.bench(), fts.bench(), its.bench(), zoo.bench())
  names <- c("xts", "fts", "its", "zoo")
  data.frame(results, names)
}

result <- bench()
