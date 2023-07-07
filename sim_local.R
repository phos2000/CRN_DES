source("CVD_DES_CRN.R")

params$vN = 10000
params_CRN = params
params_CRN$nRN = params_CRN$vN * length(rCRN)
params_CRN$randomNums = runif(params_CRN$nRN)

# need to pre-set vN for the CRN storage


# Without CRN -------------------------------------------------------------

params_PSA = params
rrStatinsASCVD_100draws = randomDraw("rrStatinsASCVD",params_PSA,100)
n_cores = detectCores() - 1

results_woCRN = mclapply(1:100, mc.cores = n_cores, function(i){
  print(paste0(i))
  params_PSA$rrStatinsASCVD = rrStatinsASCVD_100draws[i]
  
  time_start <- Sys.time()
  result = RIPS_CVD_run(params_PSA)
  time_end <- Sys.time()
  
  result = result %>%
    mutate(NHB = Effect*params_PSA$wtp - Cost,
           rrStatinsASCVD = rrStatinsASCVD_100draws[i],
           time = time_end - time_start)
  
  return(result)
})


# With CRN ----------------------------------------------------------------

params_PSA = params_CRN

results_wCRN = mclapply(1:100, mc.cores = n_cores, function(i){
  print(paste0(i))
  params_PSA$rrStatinsASCVD = rrStatinsASCVD_100draws[i]
  
  time_start <- Sys.time()
  result = RIPS_CVD_run(params_PSA)
  time_end <- Sys.time()
  
  result = result %>%
    mutate(NHB = Effect*params_PSA$wtp - Cost,
           rrStatinsASCVD = rrStatinsASCVD_100draws[i],
           time = time_end - time_start)
  
  return(result)
})

save(results_woCRN, results_wCRN, file = "output/rrStatinsASCVD100.RData")
