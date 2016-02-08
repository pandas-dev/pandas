library(microbenchmark)
library(data.table)


data.frame.subset.bench <- function (n=1e7, times=30) {
    df <- data.frame(a=rnorm(n), b=rnorm(n), c=rnorm(n))
    print(microbenchmark(subset(df, a <= b & b <= (c ^ 2 + b ^ 2 - a) & b > c),
                         times=times))
}


# data.table allows something very similar to query with an expression
# but we have chained comparisons AND we're faster BOO YAH!
data.table.subset.expression.bench <- function (n=1e7, times=30) {
    dt <- data.table(a=rnorm(n), b=rnorm(n), c=rnorm(n))
    print(microbenchmark(dt[, a <= b & b <= (c ^ 2 + b ^ 2 - a) & b > c],
                         times=times))
}


# compare against subset with data.table for good measure
data.table.subset.bench <- function (n=1e7, times=30) {
    dt <- data.table(a=rnorm(n), b=rnorm(n), c=rnorm(n))
    print(microbenchmark(subset(dt, a <= b & b <= (c ^ 2 + b ^ 2 - a) & b > c),
                         times=times))
}


data.frame.with.bench <- function (n=1e7, times=30) {
    df <- data.frame(a=rnorm(n), b=rnorm(n), c=rnorm(n))

    print(microbenchmark(with(df, a + b * (c ^ 2 + b ^ 2 - a) / (a * c) ^ 3),
                         times=times))
}


data.table.with.bench <- function (n=1e7, times=30) {
    dt <- data.table(a=rnorm(n), b=rnorm(n), c=rnorm(n))
    print(microbenchmark(with(dt, a + b * (c ^ 2 + b ^ 2 - a) / (a * c) ^ 3),
                         times=times))
}


bench <- function () {
    data.frame.subset.bench()
    data.table.subset.expression.bench()
    data.table.subset.bench()
    data.frame.with.bench()
    data.table.with.bench()
}


bench()
