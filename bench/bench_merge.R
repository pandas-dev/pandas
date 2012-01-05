library(plyr)
library(data.table)
N <- 10000
indices = rep(NA, N)
indices2 = rep(NA, N)
for (i in 1:N) {
  indices[i] <- paste(sample(letters, 10), collapse="")
  indices2[i] <- paste(sample(letters, 10), collapse="")
}
left <- data.frame(key=rep(indices, 10),
                   key2=rep(indices2, 10),
                   value=rnorm(100000))
right <- data.frame(key=indices,
                    key2=indices2,
                    value2=rnorm(10000))

right2 <- data.frame(key=rep(indices, 2),
                     key2=rep(indices2, 2),
                     value2=rnorm(20000))

## left <- data.frame(key=rep(indices[1:1000], 10),
##                    key2=rep(indices2[1:1000], 10),
##                    value=rnorm(100000))
## right <- data.frame(key=indices[1:1000],
##                     key2=indices2[1:1000],
##                     value2=rnorm(10000))

timeit <- function(func, niter=10) {
  timing = rep(NA, niter)
  for (i in 1:niter) {
    gc()
    timing[i] <- system.time(func())[3]
  }
  mean(timing)
}

left.join <- function(sort=FALSE) {
  result <- base::merge(left, right, all.x=TRUE, sort=sort)
}

right.join <- function(sort=FALSE) {
  result <- base::merge(left, right, all.y=TRUE, sort=sort)
}

outer.join <- function(sort=FALSE) {
  result <- base::merge(left, right, all=TRUE, sort=sort)
}

inner.join <- function(sort=FALSE) {
  result <- base::merge(left, right, sort=sort)
}

plyr.join <- function(type) {
  result <- plyr::join(left, right, by=c("key", "key2"),
                       type=type, match="first")
}

sort.options <- c(FALSE, TRUE)

results <- matrix(nrow=3, ncol=3)
colnames(results) <- c("base::merge", "plyr", "data.table")
rownames(results) <- c("inner", "outer", "left")

base.functions <- c(inner.join, outer.join, left.join)
plyr.functions <- c(function() plyr.join("inner"),
                    function() plyr.join("full"),
                    function() plyr.join("left"))
dt.functions <- c(inner.join, outer.join, left.join)
for (i in 1:3) {
  base.func <- base.functions[[i]]
  plyr.func <- plyr.functions[[i]]
  ## dt.func <- dt.functions[[i]]
  results[i, 1] <- timeit(base.func)
  results[i, 2] <- timeit(plyr.func)
}

## do.something <- function(df, f) {
##   f(df)
## }
## df <- matrix(nrow=4, ncol=2)
## functions <- c(colSums, rowSums)
## g <- functions[1]
## do.something(df, function(df) g(df))

##       dont_sort   sort
## inner    0.2297 0.2286
## outer    1.1811 1.2843
## left     0.6706 0.7766
## right    0.2995 0.3371
