Copyright 2024 Yiwei Zhang


my_boundedPolyMutation <- function (parent_chromosome, lowerBounds, upperBounds, mprob, 
                                    mum) 
{
  popSize = nrow(parent_chromosome)
  varNo = ncol(parent_chromosome)
  child <- parent_chromosome
  for (i in 1:popSize) {
    for (j in 1:varNo) {
      y = child[i,j]
      if (runif(1) < mprob) { # From original function: if False, this element is skipped to consider mutation
        m = mosmafs::mutPolynomialInt(y, p = mprob, eta = mum, lower = lowerBounds[j], upper = upperBounds[j])
        child[i,j] = m
      }
    }
  }
  return(child)
}
