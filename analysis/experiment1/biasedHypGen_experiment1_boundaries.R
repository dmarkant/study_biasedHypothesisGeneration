setwd("~/studies/study_biasedHypothesisGeneration/analysis/experiment1")
library(arm)

DATADIR = '../../data_proc/experiment1/'
FITDIR = '../../data_proc/experiment1/fit_boundaries/'

SUBJ = c(48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
         63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76,
         77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
         91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
         104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
         115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
         126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136,
         137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147,
         148, 149, 150, 151, 154, 155, 157, 158,
         160, 161, 162, 163, 164, 165)

# test (classification) responses
for (sid in SUBJ) {

  print(sid)
  
  for (block in seq(0, 7)) {
    data = read.csv(paste(DATADIR, 'testdata/testdata_', sid, '_', block, '.csv', sep='', collapse=''))
    
    fit = bayesglm(resp ~ x1 + x2, data=data, family='binomial')
    write.csv(summary(fit)$coefficients, paste(FITDIR, 'fit_', sid, '_', block, '.csv', sep='', collapse=''))
    
    fitx1 = bayesglm(resp ~ x1, data=data, family='binomial')
    write.csv(summary(fitx1)$coefficients, paste(FITDIR, 'fitx1_', sid, '_', block, '.csv', sep='', collapse=''))
    
    fitx2 = bayesglm(resp ~ x2, data=data, family='binomial')
    write.csv(summary(fitx2)$coefficients, paste(FITDIR, 'fitx2_', sid, '_', block, '.csv', sep='', collapse=''))
    
    fitnull = bayesglm(resp ~ 1, data=data, family='binomial')
    write.csv(c(fit$aic, fitx1$aic, fitx2$aic, fitnull$aic), paste(FITDIR, 'fit_comp_', sid, '_', block, '.csv', sep='', collapse=''))
  }
}


# training (selection) data
for (sid in SUBJ) {
    print(sid)
    for (block in seq(0, 7)) {
      data = read.csv(paste(DATADIR, 'traindata/traindata_', sid, '_', block, '.csv', sep='', collapse=''))
      
      fit = bayesglm(resp ~ x1 + x2, data=data, family='binomial')
      fitx1 = bayesglm(resp ~ x1, data=data, family='binomial')
      fitx2 = bayesglm(resp ~ x2, data=data, family='binomial')
      
      write.csv(summary(fit)$coefficients, paste(FITDIR, 'fit_training_', sid, '_', block, '.csv', sep='', collapse=''))
      write.csv(c(fit$aic, fitx1$aic, fitx2$aic), paste(FITDIR, 'fit_training_comp_', sid, '_', block, '.csv', sep='', collapse=''))
      
    }
}
