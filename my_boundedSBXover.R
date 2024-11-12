Copyright 2024 Yiwei Zhang

my_boundedSBXover <- function (parent_chromosome, lowerBounds, upperBounds, cprob, 
                            mu) 
{
  popSize = nrow(parent_chromosome)
  varNo = ncol(parent_chromosome)
  child <- parent_chromosome
  p <- 1
  for (i in 1:(popSize/2)) {
    if (runif(1) < cprob) {  # From original function: if False, two rows are skipped to consider crossover
      parent1 <- child[p,]
      parent2 <- child[p+1,]
      ### Customize: use Integer SBX instead
      c <- mosmafs::recIntSBX(list(parent1,parent2), eta = mu, p = cprob, lower = lowerBounds, upper = upperBounds)
      child[p,] <- unlist(c[1])
      child[p + 1,] <- unlist(c[2])
    }
    p <- p + 2
  }
  return(child)
}
