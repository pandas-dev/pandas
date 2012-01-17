library(plyr)
library(data.table)
N <- 10000
indices = rep(NA, N)
indices2 = rep(NA, N)
for (i in 1:N) {
  indices[i] <- paste(sample(letters, 10), collapse="")
  indices2[i] <- paste(sample(letters, 10), collapse="")
}
left <- data.frame(key=rep(indices[1:8000], 10),
                   key2=rep(indices2[1:8000], 10),
                   value=rnorm(80000))
right <- data.frame(key=indices[2001:10000],
                    key2=indices2[2001:10000],
                    value2=rnorm(8000))

right2 <- data.frame(key=rep(right$key, 2),
                     key2=rep(right$key2, 2),
                     value2=rnorm(16000))

left.dt <- data.table(left, key=c("key", "key2"))
right.dt <- data.table(right, key=c("key", "key2"))
right2.dt <- data.table(right2, key=c("key", "key2"))

# left.dt2 <- data.table(left)
# right.dt2 <- data.table(right)

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
  result <- base::merge(left, right, all=FALSE, sort=sort)
}

left.join.dt <- function(sort=FALSE) {
  result <- right.dt[left.dt]
}

right.join.dt <- function(sort=FALSE) {
  result <- left.dt[right.dt]
}

outer.join.dt <- function(sort=FALSE) {
  result <- merge(left.dt, right.dt, all=TRUE, sort=sort)
}

inner.join.dt <- function(sort=FALSE) {
  result <- merge(left.dt, right.dt, all=FALSE, sort=sort)
}

plyr.join <- function(type) {
  result <- plyr::join(left, right, by=c("key", "key2"),
                       type=type, match="first")
}

sort.options <- c(FALSE, TRUE)

# many-to-one

results <- matrix(nrow=4, ncol=3)
colnames(results) <- c("base::merge", "plyr", "data.table")
rownames(results) <- c("inner", "outer", "left", "right")

base.functions <- c(inner.join, outer.join, left.join, right.join)
plyr.functions <- c(function() plyr.join("inner"),
                    function() plyr.join("full"),
                    function() plyr.join("left"),
					function() plyr.join("right"))
dt.functions <- c(inner.join.dt, outer.join.dt, left.join.dt, right.join.dt)
for (i in 1:4) {
  base.func <- base.functions[[i]]
  plyr.func <- plyr.functions[[i]]
  dt.func <- dt.functions[[i]]
  results[i, 1] <- timeit(base.func)
  results[i, 2] <- timeit(plyr.func)
  results[i, 3] <- timeit(dt.func)
}


# many-to-many

left.join <- function(sort=FALSE) {
  result <- base::merge(left, right2, all.x=TRUE, sort=sort)
}

right.join <- function(sort=FALSE) {
  result <- base::merge(left, right2, all.y=TRUE, sort=sort)
}

outer.join <- function(sort=FALSE) {
  result <- base::merge(left, right2, all=TRUE, sort=sort)
}

inner.join <- function(sort=FALSE) {
  result <- base::merge(left, right2, all=FALSE, sort=sort)
}

left.join.dt <- function(sort=FALSE) {
  result <- right2.dt[left.dt]
}

right.join.dt <- function(sort=FALSE) {
  result <- left.dt[right2.dt]
}

outer.join.dt <- function(sort=FALSE) {
  result <- merge(left.dt, right2.dt, all=TRUE, sort=sort)
}

inner.join.dt <- function(sort=FALSE) {
  result <- merge(left.dt, right2.dt, all=FALSE, sort=sort)
}

sort.options <- c(FALSE, TRUE)

# many-to-one

results <- matrix(nrow=4, ncol=3)
colnames(results) <- c("base::merge", "plyr", "data.table")
rownames(results) <- c("inner", "outer", "left", "right")

base.functions <- c(inner.join, outer.join, left.join, right.join)
plyr.functions <- c(function() plyr.join("inner"),
                    function() plyr.join("full"),
                    function() plyr.join("left"),
					function() plyr.join("right"))
dt.functions <- c(inner.join.dt, outer.join.dt, left.join.dt, right.join.dt)
for (i in 1:4) {
  base.func <- base.functions[[i]]
  plyr.func <- plyr.functions[[i]]
  dt.func <- dt.functions[[i]]
  results[i, 1] <- timeit(base.func)
  results[i, 2] <- timeit(plyr.func)
  results[i, 3] <- timeit(dt.func)
}

