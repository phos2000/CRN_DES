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

time_start <- Sys.time()
patientWcea_woCRN = mclapply(1:10, mc.cores = n_cores, function(i){
  RIPS_CVD_data(params_PSA)
})
  
strategyWcea = bind_rows(patientWcea_woCRN) %>% 
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
         #rrStatinsASCVD = rrStatinsASCVD_1draws,
         time = time_end - time_start)


# With CRN ----------------------------------------------------------------

params_PSA = params_CRN
params_PSA$rrStatinsASCVD = rrStatinsASCVD_1draws

time_start <- Sys.time()
patientWcea_wCRN = mclapply(1:10, mc.cores = n_cores, function(i){
  params_PSA$nRN = params_CRN$vN * length(rCRN)
  params_PSA$randomNums = params_CRN$randomNums[((i-1)*params_PSA$nRN +1):(i*params_PSA$nRN)]
  RIPS_CVD_data(params_PSA)
})

strategyWcea = bind_rows(patientWcea_wCRN) %>% 
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
         #rrStatinsASCVD = rrStatinsASCVD_1draws,
         time = time_end - time_start)

save(results_woCRN, results_wCRN, file = paste0("output/rrStatinsASCVD-", x, ".RData"))
# save(patientWcea_woCRN, patientWcea_wCRN, results_woCRN, results_wCRN, file = "converge100k.RData")