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
