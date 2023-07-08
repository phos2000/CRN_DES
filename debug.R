hist(patientWcea_wCRN$nhanesNo)
hist(patientWcea_woCRN$nhanesNo)
hist(sample(1:nrow(nhanes_pop), 200000, prob=1/nhanes_pop$WTMEC2YR, replace = TRUE)) # wrong
hist(sample(1:nrow(nhanes_pop), 200000, prob=nhanes_pop$WTMEC2YR, replace = TRUE)) # correct
NoCRN = sapply(0:99999, function(x) {params_CRN$randomNums[x*length(rCRN) + rCRN["nhanesNo"]]})
NoSEQ = sapply(NoCRN, function(x){which(nhanes_pop$cum_weight>=x)[1]})
hist(NoSEQ)