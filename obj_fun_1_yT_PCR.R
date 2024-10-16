# Objective Function: PCR - 1 PT + YES absence-triggering tests
obj_fun_1_yT_PCR <- function(x){
  i = 2
  t_ad = c()
  while (t <= 112){
    for (t in seq(gaps[i-1]+1,gaps[i])){
      vector <- c(S=current[t-1,1:nage],E=current[t-1,(nage+1):(2*nage)],Ia=current[t-1,(2*nage+1):(3*nage)],
                  Ips=current[t-1,(3*nage+1):(4*nage)],Is=current[t-1,(4*nage+1):(5*nage)],
                  TFNa=current[t-1,(5*nage+1):(6*nage)],TFNs=current[t-1,(6*nage+1):(7*nage)],
                  TISOa=current[t-1,(7*nage+1):(8*nage)],TISOs=current[t-1,(8*nage+1):(9*nage)],
                  UTISOa=current[t-1,(9*nage+1):(10*nage)],UTISOs=current[t-1,(10*nage+1):(11*nage)],
                  RC=current[t-1,(11*nage+1):(12*nage)],RU=current[t-1,(12*nage+1):(13*nage)],
                  D=current[t-1,(13*nage+1):(14*nage)], I=current[t-1,(14*nage+1):(15*nage)],
                  I_known=current[t-1,(15*nage+1):(16*nage)],
                  Trap_a=current[t-1,(16*nage+1):(17*nage)],Trap_s=current[t-1,(17*nage+1):(18*nage)])
      vector = unlist(vector)
      # Monday ISOLATION as indicator of the additional testing frequency
      ind_s = N_total*unname(rowSums(current[(7*i-13),c((7*nage+1):(8*nage-1),(8*nage+1):(9*nage-1))])) # indicator students
      ind_t = N_total*unname(rowSums(current[(7*i-13),c((8*nage),(9*nage))])) # indicator teachers
      t_ad[i-1] = ifelse((ind_s >= 10 | ind_t >= 2), 0.4, 0) # trigger 40% screening
      # separate weekdays and weekends
      if (t%%7 == 0 | t%%7 == 1){
        para[3] <- 0
        parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
        times <- seq(t-1, t, by = 1)
        current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
        current[t-1,(17*nage+1)]=0
      }
      if (t%%7 != 0 & t%%7 != 1){
        # mask or no mask for the week, let column Trap_s represent mask, it affects infection on next day through beta
        if (floor(x[16+i-1]) == 0){
          para[3] = beta
          # routine NO additional NO
          if (floor(x[i-1]) == 0 & t_ad[i-1] == 0){
            parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
            times <- seq(t-1, t, by = 1)
            current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
            current[t-1,(16*nage+1)] = 0
            current[t-1,(16*nage+2)] = 0
            current[t-1,(16*nage+3)] = 0
          }
          # routine NO additional YES
          if (floor(x[i-1]) == 0 & t_ad[i-1] > 0){
            # tests on day 2, it affects infections on day 3 for rapid tests
            if (t%%7 == 3){
              parameters <- c(para,  capacity_prop = test_level,Tpcr_a = test_level_add, Tpcr_s = test_level_add)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = test_level_add
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 1
            }
            # no tests on other days
            if (t%%7 != 3){
              parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = 0
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 0
            }
          }
          # routine YES additional NO
          if (floor(x[i-1]) > 0 & t_ad[i-1] == 0){
            # tests on day 3, it affects infections on day 4 for rapid tests
            if (t%%7 == 4){
              parameters <- c(para,  capacity_prop = test_level,Tpcr_a = floor(x[i-1])*test_level_rou, Tpcr_s = floor(x[i-1])*test_level_rou)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = floor(x[i-1])*test_level_rou
              current[t-1,(16*nage+2)] = 1
              current[t-1,(16*nage+3)] = 0
            }
            # no tests on other days
            if (t%%7 != 4){
              parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = 0
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 0
            }
          }
          # routine YES additional YES
          if (floor(x[i-1]) > 0 & t_ad[i-1] > 0){
            # tests on day 2 and day 3, it affects infection on day3 and day4 for rapid tests
            if (t%%7 == 3){
              parameters <- c(para, capacity_prop = test_level,Tpcr_a = test_level_add, Tpcr_s = test_level_add)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = test_level_add
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 1
            }
            if (t%%7 == 4){
              parameters <- c(para, capacity_prop = test_level,Tpcr_a = floor(x[i-1])*test_level_rou, Tpcr_s = floor(x[i-1])*test_level_rou)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = floor(x[i-1])*test_level_rou
              current[t-1,(16*nage+2)] = 1
              current[t-1,(16*nage+3)] = 0
            }
            # no tests on other days
            if (t%%7 != 3 & t%%7 != 4){
              parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = 0
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 0
            }
          }
          current[t-1,(17*nage+1)] = 0 # mask
        }
        if (floor(x[16+i-1]) == 1){
          para[3] = beta * 0.75  # Mask Effectiveness: base 25%; sensitivity analysis 10% and 40%
          # routine NO additional NO
          if (floor(x[i-1]) == 0 & t_ad[i-1] == 0){
            parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
            times <- seq(t-1, t, by = 1)
            current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
            current[t-1,(16*nage+1)] = 0
            current[t-1,(16*nage+2)] = 0
            current[t-1,(16*nage+3)] = 0
          }
          # routine NO additional YES
          if (floor(x[i-1]) == 0 & t_ad[i-1] > 0){
            # tests on day 2, it affects infections on day 3 for rapid tests
            if (t%%7 == 3){
              parameters <- c(para,  capacity_prop = test_level,Tpcr_a = test_level_add, Tpcr_s = test_level_add)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = test_level_add
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 1
            }
            # no tests on other days
            if (t%%7 != 3){
              parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = 0
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 0
            }
          }
          # routine YES additional NO
          if (floor(x[i-1]) > 0 & t_ad[i-1] == 0){
            # tests on day 3, it affects infections on day 4 for rapid tests
            if (t%%7 == 4){
              parameters <- c(para,  capacity_prop = test_level,Tpcr_a = floor(x[i-1])*test_level_rou, Tpcr_s = floor(x[i-1])*test_level_rou)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = floor(x[i-1])*test_level_rou
              current[t-1,(16*nage+2)] = 1
              current[t-1,(16*nage+3)] = 0
            }
            # no tests on other days
            if (t%%7 != 4){
              parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = 0
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 0
            }
          }
          # routine YES additional YES
          if (floor(x[i-1]) > 0 & t_ad[i-1] > 0){
            # tests on day 2 and day 3, it affects infection on day3 and day4 for rapid tests
            if (t%%7 == 3){
              parameters <- c(para,  capacity_prop = test_level,Tpcr_a = test_level_add, Tpcr_s = test_level_add)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = test_level_add
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 1
            }
            if (t%%7 == 4){
              parameters <- c(para,  capacity_prop = test_level,Tpcr_a = floor(x[i-1])*test_level_rou, Tpcr_s = floor(x[i-1])*test_level_rou)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = floor(x[i-1])*test_level_rou
              current[t-1,(16*nage+2)] = 1
              current[t-1,(16*nage+3)] = 0
            }
            # no tests on other days
            if (t%%7 != 3 & t%%7 != 4){
              parameters <- c(para, capacity_prop = 0,Tpcr_a = 0, Tpcr_s = 0)
              times <- seq(t-1, t, by = 1)
              current[t,] <- as.data.frame(ode(y = vector, times = times, func = sir_new, parms = parameters))[2,-1]
              current[t-1,(16*nage+1)] = 0
              current[t-1,(16*nage+2)] = 0
              current[t-1,(16*nage+3)] = 0
            }
          }
          current[t-1,(17*nage+1)] = 1 # mask
        }
      }
      # new infection on day 1
      if (t%%7 == 1){
        if (n_new[i-1] == 1){
          current[t,Y2[1,i-1]] = current[t,Y2[1,i-1]] - 1/N_total # S
          
          current[t,nage+Y2[1,i-1]] = current[t,nage+Y2[1,i-1]] + 1/N_total # E
        }
        if (n_new[i-1] == 2){
          current[t,Y2[2,i-1]] = current[t,Y2[2,i-1]] - 1/N_total # S
          current[t,Y2[3,i-1]] = current[t,Y2[3,i-1]] - 1/N_total # S
          
          current[t,nage+Y2[2,i-1]] = current[t,nage+Y2[2,i-1]] + 1/N_total # E
          current[t,nage+Y2[3,i-1]] = current[t,nage+Y2[3,i-1]] + 1/N_total # E
        }
        if (n_new[i-1] == 3){
          current[t,Y2[4,i-1]] = current[t,Y2[4,i-1]] - 1/N_total # S
          current[t,Y2[5,i-1]] = current[t,Y2[5,i-1]] - 1/N_total # S
          current[t,Y2[6,i-1]] = current[t,Y2[6,i-1]] - 1/N_total # S
          
          current[t,nage+Y2[4,i-1]] = current[t,nage+Y2[4,i-1]] + 1/N_total # E
          current[t,nage+Y2[5,i-1]] = current[t,nage+Y2[5,i-1]] + 1/N_total # E
          current[t,nage+Y2[6,i-1]] = current[t,nage+Y2[6,i-1]] + 1/N_total # E
        }
        if (n_new[i-1] == 4){
          current[t,Y2[7,i-1]] = current[t,Y2[7,i-1]] - 1/N_total # S$Y2[7,i-1]
          current[t,Y2[8,i-1]] = current[t,Y2[8,i-1]] - 1/N_total # S$Y2[7,i-1]
          current[t,Y2[9,i-1]] = current[t,Y2[9,i-1]] - 1/N_total # S$Y2[7,i-1]
          current[t,Y2[10,i-1]] = current[t,Y2[10,i-1]] - 1/N_total # S$Y2[7,i-1]
          current[t,nage+Y2[7,i-1]] = current[t,nage+Y2[7,i-1]] + 1/N_total # E$Y2[7,i-1]
          current[t,nage+Y2[8,i-1]] = current[t,nage+Y2[8,i-1]] + 1/N_total # E$Y2[7,i-1]
          current[t,nage+Y2[9,i-1]] = current[t,nage+Y2[9,i-1]] + 1/N_total # E$Y2[7,i-1]
          current[t,nage+Y2[10,i-1]] = current[t,nage+Y2[10,i-1]] + 1/N_total # E$Y2[7,i-1]
        }
      }
    }
    t = gaps[i]+1
    i = i + 1
    #current_all <- current_all + current
  }
  current_opt <<- current
  t_addition <<- t_ad
  # cumulative number of cases
  cost1 = N_total*unname(rowSums(current[Tn,(14*nage+1):(15*nage)]))
  # number of testing: sum(people not isolation on day t) * current[,(16*nage+1)]
  cost2 = N_total*sum((unname(rowSums(current[,1:(14*19)])-rowSums(current[,(7*19+1):(9*19)])))*current[,(16*nage+1)])
  # number of mask weeks:
  cost3 = sum(current[,(17*nage+1)])/5
  return(c(cost1,cost2,cost3))
}