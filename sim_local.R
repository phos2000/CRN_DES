source("CVD_DES_CRN.R")
# params$vN = 100000
# params_CRN = params
# params_CRN$nRN = params_CRN$vN * length(rCRN)
# params_CRN$randomNums = runif(params_CRN$nRN)

load("params_CRN.RData")
params$vN = 10000
params_CRN$vN = 10000

# need to pre-set vN for the CRN storage


# Without CRN -------------------------------------------------------------

params_PSA = params
rrStatinsASCVD_100draws = randomDraw("rrStatinsASCVD",params_PSA,100)
n_cores = detectCores() - 1

results_woCRN = mclapply(1:2, mc.cores = n_cores, function(i){
  print(paste0(i))
  params_PSA$rrStatinsASCVD = rrStatinsASCVD_100draws[i]
  
  # time_start <- Sys.time()
  # result = RIPS_CVD_run(params_PSA)
  # time_end <- Sys.time()
  # 
  # result = result %>%
  #   mutate(NHB = Effect*params_PSA$wtp - Cost,
  #          rrStatinsASCVD = rrStatinsASCVD_100draws[i],
  #          time = time_end - time_start)
  
  time_start <- Sys.time()
  patientWcea = lapply(1:10, function(j){
    RIPS_CVD_data(params_PSA)
  })
  
  strategyWcea = bind_rows(patientWcea) %>% 
    group_by(strategy) %>%
    summarize(
      time_in_model = mean(time_in_model),
      cost = mean(cost_disc), 
      QALY = mean(util_disc)) %>%
    mutate(strategy = ifelse(strategy == 1, "Statins(2013 ACC/AHA)", "Status Quo"))
  
  results_woCRN = calculate_icers(cost = strategyWcea$cost,
                                  effect = strategyWcea$QALY,
                                  strategies = strategyWcea$strategy) %>%
    left_join(strategyWcea %>% select(Strategy = strategy, LE = time_in_model))
  
  time_end <- Sys.time()
  
  results_woCRN = results_woCRN %>%
    mutate(NHB = Effect*params_PSA$wtp - Cost,
           rrStatinsASCVD = rrStatinsASCVD_100draws[i],
           time = time_end - time_start)
  
  return(result_woCRN)
})


# With CRN ----------------------------------------------------------------

params_PSA = params_CRN

results_wCRN = mclapply(1:2, mc.cores = n_cores, function(i){
  print(paste0(i))
  params_PSA$rrStatinsASCVD = rrStatinsASCVD_100draws[i]
  
  # time_start <- Sys.time()
  # result = RIPS_CVD_run(params_PSA)
  # time_end <- Sys.time()
  # 
  # result = result %>%
  #   mutate(NHB = Effect*params_PSA$wtp - Cost,
  #          rrStatinsASCVD = rrStatinsASCVD_100draws[i],
  #          time = time_end - time_start)
  
  time_start <- Sys.time()
  patientWcea = lapply(1:10, mc.cores = n_cores, function(j){
    params_PSA$nRN = params_CRN$vN * length(rCRN)
    params_PSA$randomNums = params_CRN$randomNums[((j-1)*params_PSA$nRN +1):(j*params_PSA$nRN)]
    RIPS_CVD_data(params_PSA)
  })
  
  strategyWcea = bind_rows(patientWcea) %>% 
    group_by(strategy) %>%
    summarize(
      time_in_model = mean(time_in_model),
      cost = mean(cost_disc), 
      QALY = mean(util_disc)) %>%
    mutate(strategy = ifelse(strategy == 1, "Statins(2013 ACC/AHA)", "Status Quo"))
  
  results_wCRN = calculate_icers(cost = strategyWcea$cost,
                                 effect = strategyWcea$QALY,
                                 strategies = strategyWcea$strategy) %>%
    left_join(strategyWcea %>% select(Strategy = strategy, LE = time_in_model))
  
  time_end <- Sys.time()
  
  results_wCRN = results_wCRN %>%
    mutate(NHB = Effect*params_PSA$wtp - Cost,
           rrStatinsASCVD = rrStatinsASCVD_100draws[i],
           time = time_end - time_start)
  
  return(results_wCRN)
})

save(results_woCRN, results_wCRN, file = "output/rrStatinsASCVD100.RData")
