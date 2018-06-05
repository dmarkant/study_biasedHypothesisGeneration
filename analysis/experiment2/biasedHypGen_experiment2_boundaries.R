setwd("~/studies/study_biasedHypothesisGeneration/analysis/experiment2")
library(arm)

DATADIR = '../../data_proc/experiment2/'
FITDIR = '../../data_proc/experiment2/fit_boundaries/'

SUBJ = c(234, 236, 238, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252,
         253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 266, 267,
         268, 269, 271, 272, 273, 274, 275, 276, 277, 278, 280, 281, 282, 283,
         284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297,
         298, 299, 300, 301, 302, 303, 305, 306, 307, 308, 309, 310, 311, 312,
         313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 325, 328, 329,
         330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 342, 343, 344, 345,
         346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 362, 363,
         364, 365, 366, 367, 373, 374, 376, 377)


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
