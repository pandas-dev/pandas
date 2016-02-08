N <- 100000

k1 = rep(NA, N)
k2 = rep(NA, N)
for (i in 1:N){
  k1[i] <- paste(sample(letters, 1), collapse="")
  k2[i] <- paste(sample(letters, 1), collapse="")
}
df <- data.frame(a=k1, b=k2, c=rep(1:100, N / 100))
df2 <- data.frame(a=k1, b=k2)

timings <- numeric()
timings2 <- numeric()
for (i in 1:50) {
  gc()
  timings[i] = system.time(deduped <- df[!duplicated(df),])[3]
  gc()
  timings2[i] = system.time(deduped <- df[!duplicated(df[,c("a", "b")]),])[3]
}

mean(timings)
mean(timings2)
