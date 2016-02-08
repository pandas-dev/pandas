library(reshape2)


n <- 100000
a.size <- 5
b.size <- 5

data <- data.frame(a=sample(letters[1:a.size], n, replace=T),
                   b=sample(letters[1:b.size], n, replace=T),
                   c=rnorm(n),
                   d=rnorm(n))

timings <- numeric()

# acast(melt(data, id=c("a", "b")), a ~ b, mean)
# acast(melt(data, id=c("a", "b")), a + b ~ variable, mean)

for (i in 1:10) {
  gc()
  tim <- system.time(acast(melt(data, id=c("a", "b")), a ~ b, mean,
                           subset=.(variable=="c")))
  timings[i] = tim[3]
}

mean(timings)

acast(melt(data, id=c("a", "b")), a ~ b, mean, subset=.(variable="c"))
