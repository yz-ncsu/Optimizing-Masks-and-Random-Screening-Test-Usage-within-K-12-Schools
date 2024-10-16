my_nsga2R_laptop <- function (fn, varNo, objDim, lowerBounds = rep(-Inf, varNo), 
                       upperBounds = rep(Inf, varNo), popSize = 100, tourSize = 2, 
                       generations = 20, cprob = 0.7, XoverDistIdx = 5, mprob = 0.2, 
                       MuDistIdx = 10) 
{
  cat("********** R based Nondominated Sorting Genetic Algorithm II *********")
  cat("\n")
  cat("initializing the population")
  cat("\n")
  # parent <- t(sapply(1:popSize, function(u) array(runif(length(lowerBounds), 
  #                                                       lowerBounds, upperBounds))))
  ### Customize 1: assume integer initials
  parent <- matrix(c(rdunif(varNo*(popSize-2)/2,upperBounds[1],lowerBounds[1]),rdunif(varNo/2*(popSize-2),upperBounds[varNo],lowerBounds[varNo])),
                    nrow = popSize-2, byrow = FALSE)
  ###
  ### Customize 2: specify 2 initials as no_mask_no_test and all_mask_all_test
  parent = rbind(parent, rep(0,varNo), c(rep(upperBounds[1],varNo/2),rep(upperBounds[varNo],varNo/2)))
  wb <- createWorkbook()
  addWorksheet(wb, "initial")
  writeData(wb, "initial", parent)
  #write.xlsx(parent, "~/parents in iter laptop.xlsx", sheetName = "initial", append = TRUE)
  ###
  parent <- cbind(parent, t(apply(parent, 1, fn)))
  cat("ranking the initial population")
  cat("\n")
  ranking <- fastNonDominatedSorting(parent[, (varNo + 1):(varNo + 
                                                             objDim)])
  rnkIndex <- integer(popSize)
  i <- 1
  while (i <= length(ranking)) {
    rnkIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  parent <- cbind(parent, rnkIndex)
  cat("crowding distance calculation")
  cat("\n")
  objRange <- apply(parent[, (varNo + 1):(varNo + objDim)], 
                    2, max) - apply(parent[, (varNo + 1):(varNo + objDim)], 
                                    2, min)
  cd <- crowdingDist4frnt(parent, ranking, objRange)
  parent <- cbind(parent, apply(cd, 1, sum))
  for (iter in 1:generations) {
    cat("---------------generation---------------", iter, 
        "starts")
    cat("\n")
    ### Customize 3: add a timer
    tic(iter)
    cat("\n")
    ###
    cat("tournament selection")
    cat("\n")
    matingPool <- tournamentSelection(parent, popSize, tourSize)
    cat("crossover operator")
    cat("\n")
    childAfterX <- my_boundedSBXover(matingPool[, 1:varNo], 
                                  lowerBounds, upperBounds, cprob, XoverDistIdx)
    cat("mutation operator")
    cat("\n")
    childAfterM <- my_boundedPolyMutation(childAfterX, lowerBounds, 
                                       upperBounds, mprob, MuDistIdx)
    cat("evaluate the objective fns of childAfterM")
    cat("\n")
    childAfterM <- cbind(childAfterM, t(apply(childAfterM, 
                                              1, fn)))
    cat("Rt = Pt + Qt")
    cat("\n")
    parentNext <- rbind(parent[, 1:(varNo + objDim)], childAfterM)
    cat("ranking again")
    cat("\n")
    ranking <- fastNonDominatedSorting(parentNext[, (varNo + 
                                                       1):(varNo + objDim)])
    i <- 1
    while (i <= length(ranking)) {
      rnkIndex[ranking[[i]]] <- i
      i <- i + 1
    }
    parentNext <- cbind(parentNext, rnkIndex)
    cat("crowded comparison again")
    cat("\n")
    objRange <- apply(parentNext[, (varNo + 1):(varNo + 
                                                  objDim)], 2, max) - apply(parentNext[, (varNo + 
                                                                                            1):(varNo + objDim)], 2, min)
    cd <- crowdingDist4frnt(parentNext, ranking, objRange)
    parentNext <- cbind(parentNext, apply(cd, 1, sum))
    parentNext.sort <- parentNext[order(parentNext[, varNo + 
                                                     objDim + 1], -parentNext[, varNo + objDim + 2]), 
    ]
    cat("environmental selection")
    cat("\n")
    parent <- parentNext.sort[1:popSize, ]
    ###
    addWorksheet(wb, as.character(iter))
    writeData(wb, as.character(iter), parent)
    #write.xlsx(parent, "~/parents in iter laptop.xlsx", sheetName = as.character(iter), append = TRUE)
    saveWorkbook(wb, file = "parents in iter laptop.xlsx", overwrite = TRUE)
    toc(log = TRUE)
    cat("\n")
    ###
    cat("---------------generation---------------", iter, 
        "ends")
    cat("\n")
    if (iter != generations) {
      cat("\n")
      cat("********** new iteration *********")
      cat("\n")
    }
    else {
      cat("********** stop the evolution *********")
      cat("\n")
    }
  }
  result = list(functions = fn, parameterDim = varNo, objectiveDim = objDim, 
                lowerBounds = lowerBounds, upperBounds = upperBounds, 
                popSize = popSize, tournamentSize = tourSize, generations = generations, 
                XoverProb = cprob, XoverDistIndex = XoverDistIdx, mutationProb = mprob, 
                mutationDistIndex = MuDistIdx, parameters = parent[, 
                                                                   1:varNo], objectives = parent[, (varNo + 1):(varNo + 
                                                                                                                  objDim)], paretoFrontRank = parent[, varNo + objDim + 
                                                                                                                                                       1], crowdingDistance = parent[, varNo + objDim + 
                                                                                                                                                                                       2])
  class(result) = "nsga2R"
  return(result)
}
