N <- 10000
indices = rep(NA, N)
for (i in 1:N)
  indices[i] <- paste(sample(letters, 10), collapse="")

left <- data.frame(key=rep(indices, 10),
                   key2=sample(rep(indices, 10)),
                   value=rnorm(100000))
right <- data.frame(key=indices,
                    key2=sample(indices),
                    value2=rnorm(10000))

timeit <- function(func, niter=10) {
  timing = rep(NA, niter)
  for (i in 1:niter) {
    gc()
    timing[i] <- system.time(func())[3]
  }
  mean(timing)
}

left.join <- function(sort=TRUE) {
  result <- merge(left, right, all.x=TRUE, sort=sort)
}

right.join <- function(sort=TRUE) {
  result <- merge(left, right, all.y=TRUE, sort=sort)
}

outer.join <- function(sort=TRUE) {
  result <- merge(left, right, all=TRUE, sort=sort)
}

inner.join <- function(sort=TRUE) {
  reuslt <- merge(left, right, sort=sort)
}

sort.options <- c(FALSE, TRUE)

results <- matrix(nrow=4, ncol=2)
colnames(results) <- c("dont_sort", "sort")
rownames(results) <- c("inner", "outer", "left", "right")

join.functions <- c(inner.join, outer.join, left.join, right.join)
for (i in 1:4) {
  results[1, 1] <- timeit(function() {inner.join(sort=sort.options[1])})
  results[1, 2] <- timeit(function() {inner.join(sort=sort.options[2])})
  results[2, 1] <- timeit(function() {outer.join(sort=sort.options[1])})
  results[2, 2] <- timeit(function() {outer.join(sort=sort.options[2])})
  results[3, 1] <- timeit(function() {left.join(sort=sort.options[1])})
  results[3, 2] <- timeit(function() {left.join(sort=sort.options[2])})
  results[4, 1] <- timeit(function() {right.join(sort=sort.options[1])})
  results[4, 2] <- timeit(function() {right.join(sort=sort.options[2])})
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
