source("CVD_DES_CRN.R")
load("params_CRN.RData")

# Without CRN -------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
x    <- as.numeric(args[1])
set.seed(x)
rrStatinsASCVD_1draws = randomDraw("rrStatinsASCVD",params_PSA,1)

params_PSA = params
params_PSA$rrStatinsASCVD = rrStatinsASCVD_1draws
time_start <- Sys.time()
result_woCRN = RIPS_CVD_run(params_PSA)
time_end <- Sys.time()
  
result_woCRN = result_woCRN %>%
  mutate(NHB = Effect*params_PSA$wtp - Cost,
         rrStatinsASCVD = rrStatinsASCVD_1draws,
         time = time_end - time_start)

# With CRN ----------------------------------------------------------------

params_PSA = params_CRN
params_PSA$rrStatinsASCVD = rrStatinsASCVD_1draws
time_start <- Sys.time()
result_wCRN = RIPS_CVD_run(params_PSA)
time_end <- Sys.time()

result_wCRN = result_wCRN %>%
  mutate(NHB = Effect*params_PSA$wtp - Cost,
         rrStatinsASCVD = rrStatinsASCVD_100draws[i],
         time = time_end - time_start)

save(result_woCRN, result_wCRN, file = paste0("rrStatinsASCVD-", x, ".RData"))
