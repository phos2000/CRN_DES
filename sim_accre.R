source("CVD_DES_CRN.R")
load("params_CRN.RData")
params$vN = 10000
params_CRN$vN = 10000

# Without CRN -------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
x    <- as.numeric(args[1])
set.seed(x)
n_cores = detectCores() - 1
rrStatinsASCVD_1draws = randomDraw("rrStatinsASCVD",params,1)

params_PSA = params
params_PSA$rrStatinsASCVD = rrStatinsASCVD_1draws
# time_start <- Sys.time()
# result_woCRN = RIPS_CVD_run(params_PSA)
# time_end <- Sys.time()
#   
# result_woCRN = result_woCRN %>%
#   mutate(NHB = Effect*params_PSA$wtp - Cost,
#          rrStatinsASCVD = rrStatinsASCVD_1draws,
#          time = time_end - time_start)

results_woCRN = mclapply(1:10, mc.cores = n_cores, function(i){
  
  time_start <- Sys.time()
  result = RIPS_CVD_run(params_PSA)
  time_end <- Sys.time()
  
  result = result %>%
    mutate(NHB = Effect*params_PSA$wtp - Cost,
           rrStatinsASCVD = rrStatinsASCVD_1draws,
           time = time_end - time_start)
  
  return(result)
})

# With CRN ----------------------------------------------------------------

params_PSA = params_CRN
params_PSA$rrStatinsASCVD = rrStatinsASCVD_1draws
# time_start <- Sys.time()
# result_wCRN = RIPS_CVD_run(params_PSA)
# time_end <- Sys.time()
# 
# result_wCRN = result_wCRN %>%
#   mutate(NHB = Effect*params_PSA$wtp - Cost,
#          rrStatinsASCVD = rrStatinsASCVD_1draws,
#          time = time_end - time_start)

results_wCRN = mclapply(1:10, mc.cores = n_cores, function(i){
  
  time_start <- Sys.time()
  result = RIPS_CVD_run(params_PSA)
  time_end <- Sys.time()
  
  result = result %>%
    mutate(NHB = Effect*params_PSA$wtp - Cost,
           rrStatinsASCVD = rrStatinsASCVD_1draws,
           time = time_end - time_start)
  
  return(result)
})

save(results_woCRN, results_wCRN, file = paste0("output/rrStatinsASCVD-", x, ".RData"))
