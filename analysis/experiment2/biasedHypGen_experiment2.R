library(lme4)
library(lmerTest)
library(plyr)
library(grid)
library(gridExtra)
require(ggplot2)
library(multcomp)

setwd("~/studies/study_biasedHypothesisGeneration/analysis/experiment2")

data = read.csv('../../data_proc/experiment2/data_by_block.csv')
data$sid = factor(data$sid)
data$cond = factor(data$cond)
data$ruletype = factor(data$ruletype)
data$dimtype = factor(data$dimtype)
data$rulerelation = factor(data$rulerelation)
#data$block = scale(data$block)
data$block = data$block - mean(data$block)
data$repMatch = ((data$ruletype=='1D') & (data$dimtype=='absolute')) | ((data$ruletype=='2D') & (data$dimtype=='relative'))



data$prev_acc = NaN
for (sid in levels(data$sid)) {
  acc = data[data$sid==sid,]$acc
  prev_acc = c(0.5, acc[1:7])
  data[data$sid==sid, 'prev_acc'] = prev_acc
}
data$change_acc = data$acc - data$prev_acc

data$prev_acc = scale(data$prev_acc)
data$change_acc = scale(data$change_acc)


data$hyp2D = as.logical(data$X2drule_aic)
data$hypMatch = ((data$ruletype=='1D') & (!data$hyp2D)) | ((data$ruletype=='2D') & (data$hyp2D))


data$prev_hypMatch = NaN
for (sid in levels(data$sid)) {
  hypm = data[data$sid==sid,]$hypMatch
  prev_hypm = c(NaN, hypm[1:7])
  data[data$sid==sid, 'prev_hypMatch'] = prev_hypm
}



# Classification accuracy ----

# the model
m1 = glmer(acc ~ dimtype*ruletype + block + (dimtype*ruletype)/rulerelation + (1|sid),
          weights=rep(32, nrow(data)),
          data=data, family=binomial, control = glmerControl(optimizer = "bobyqa"),
          nAGQ = 10)
summary(m1)

m2 = glmer(acc ~ dimtype*ruletype*block + (dimtype*ruletype)/rulerelation + (1|sid),
          weights=rep(32, nrow(data)),
          data=data, family=binomial, control = glmerControl(optimizer = "bobyqa"),
          nAGQ = 10)
summary(m2)
anova(m1,m2)

m3 = glmer(acc ~ dimtype*ruletype*block + block*((dimtype*ruletype)/rulerelation) + (1|sid),
           weights=rep(32, nrow(data)),
           data=data, family=binomial, control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 10)
summary(m3)
anova(m2,m3)

round(c(AIC(m1), AIC(m2), AIC(m3)))


m = m3


# m4 = glmer(acc ~ dimtype*ruletype*block + block*(ruletype/rulerelation) + (1|sid),
#            weights=rep(32, nrow(data)),
#            data=data, family=binomial, control = glmerControl(optimizer = "bobyqa"),
#            nAGQ = 10)
# summary(m4)
# 
# anova(m4,m3)
# 
# m = m4



# odds ratios
exp(fixef(m))

# confidence intervals
confint(m)

# 2.5 %      97.5 %
#   .sig01                                             0.96887863  1.30373385
# (Intercept)                                        2.80406236  3.99380649
# dimtyperelative                                   -2.14824950 -0.52700696
# ruletype2D                                        -2.77127464 -1.11198194
# block                                              0.83955793  1.14165926
# dimtyperelative:ruletype2D                         1.27268855  3.60756272
# dimtyperelative:block                             -0.76862718 -0.41906070
# ruletype2D:block                                  -0.79920968 -0.45508973
# dimtyperelative:ruletype2D:block                   0.38624968  0.84016996
# dimtypeabsolute:ruletype1D:rulerelation1D-F       -0.11301087  1.68078185
# dimtyperelative:ruletype1D:rulerelation1D-F       -0.09522747  1.59939805
# dimtypeabsolute:ruletype2D:rulerelation2D-N       -1.13846416  0.49617982
# dimtyperelative:ruletype2D:rulerelation2D-N       -1.56654643  0.09728136
# dimtypeabsolute:ruletype1D:block:rulerelation1D-F -0.80169339 -0.32972135
# dimtyperelative:ruletype1D:block:rulerelation1D-F -0.07737409  0.22889402
# dimtypeabsolute:ruletype2D:block:rulerelation2D-N -0.24837410 -0.02474281
# dimtyperelative:ruletype2D:block:rulerelation2D-N -0.23242508  0.06094221


# no interaction 
m1 = glmer(acc ~ dimtype*ruletype + block + ruletype/(dimtype*rulerelation) + (1|sid),
           weights=rep(32, nrow(data)),
           data=data, family=binomial, control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 10)
summary(m1)
anova(m,m1)


# contrasts
# D1_ABS_CHEM  = c(1,0,0,0,0,0,0)
# D1_ABS_FERT  = c(1,0,0,0,0,1,0)
# D1_REL_CHEM  = c(1,1,0,0,0,0,0)
# D1_REL_FERT  = c(1,1,0,0,0,1,0)
# D2_ABS_POS   = c(1,0,1,0,0,0,0)
# D2_ABS_NEG   = c(1,0,1,0,0,0,1)
# D2_REL_POS   = c(1,1,1,0,1,0,0)
# D2_REL_NEG   = c(1,1,1,0,1,0,1)
# 
# c1 = rbind("D1_ABS - D1_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D1_REL_CHEM + D1_REL_FERT)/2),
#            "D1_ABS - D2_ABS" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
#            "D1_ABS - D2_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
#            "D1_REL - D2_ABS" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
#            "D1_REL - D2_REL" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
#            "D2_ABS - D2_REL" = ((D2_ABS_POS + D2_ABS_NEG)/2 - (D2_REL_POS + D2_REL_NEG)/2))


D1_ABS_CHEM       = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_ABS_FERT       = c(1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
D1_REL_CHEM       = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_REL_FERT       = c(1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0)
D2_ABS_POS        = c(1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
D2_ABS_NEG        = c(1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0)
D2_REL_POS        = c(1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0)
D2_REL_NEG        = c(1,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0)
block_D1_ABS_CHEM = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
block_D1_ABS_FERT = c(0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0)
block_D1_REL_CHEM = c(0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0)
block_D1_REL_FERT = c(0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0)
block_D2_ABS_POS  = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0)
block_D2_ABS_NEG  = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0)
block_D2_REL_POS  = c(0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0)
block_D2_REL_NEG  = c(0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,1)


c1 = rbind("D1_ABS - D1_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D1_REL_CHEM + D1_REL_FERT)/2),
           "D1_ABS - D2_ABS" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_ABS - D2_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D1_REL - D2_ABS" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_REL - D2_REL" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D2_ABS - D2_REL" = ((D2_ABS_POS + D2_ABS_NEG)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D1_ABS_CHEM - D1_ABS_FERT" = 
           "block_D1_ABS_CHEM" = block_D1_ABS_CHEM,
           "block_D1_ABS_FERT" = block_D1_ABS_FERT,
           "block_D1_REL_CHEM" = block_D1_REL_CHEM,
           "block_D1_REL_FERT" = block_D1_REL_FERT,
           "block_D2_ABS_POS" = block_D2_ABS_POS,
           "block_D2_ABS_NEG" = block_D2_ABS_NEG,
           "block_D2_REL_POS" = block_D2_REL_POS,
           "block_D2_REL_NEG" = block_D2_REL_NEG
           )

summary(glht(m, c1), test = adjusted("bonferroni"))




y = ddply(data, c('dimtype', 'ruletype', 'rulerelation', 'sid'), function(x) {
  mn = mean(x$acc)
  data.frame(acc=mn)
})

# by condition
ddply(y, c('dimtype', 'ruletype'), function(x) {
  c(mn=mean(x$acc), sd=sd(x$acc))
})

# by rule variant
ddply(y, c('stimtype', 'ruletype', 'rulevariant'), function(x) {
  c(mn=mean(x$acc), sd=sd(x$acc))
})


# Visualize model fit

data$pred_acc = predict(m, type='response')

r = ddply(data, c('dimtype', 'ruletype', 'rulerelation', 'block'), function(x) { 
  mn = mean(x$acc) 
  mn_pr = mean(x$pred_acc)
  data.frame(acc=mn, pred=mn_pr)
})

p1 = ggplot(data=r[r$dimtype=='absolute' & r$ruletype=='1D',]) +
  geom_line(aes(x=block, y=acc, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0.5, 1)

p2 = ggplot(data=r[r$dimtype=='relative' & r$ruletype=='1D',]) +
  geom_line(aes(x=block, y=acc, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0.5, 1)

p3 = ggplot(data=r[r$dimtype=='absolute' & r$ruletype=='2D',]) +
  geom_line(aes(x=block, y=acc, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0.5, 1)

p4 = ggplot(data=r[r$dimtype=='relative' & r$ruletype=='2D',]) +
  geom_line(aes(x=block, y=acc, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0.5, 1)

grid.arrange(p1, p2, p3, p4)



# Selection distance ----

## Correlations with accuracy ----

agg = ddply(data, c('sid', 'dimtype', 'ruletype', 'rulerelation'), function(x) {
  c(mn.acc=mean(x$acc),
    mn.dist=mean(x$dist_med))
})

ddply(agg, c('dimtype', 'ruletype'), function(x) {
  res = cor.test(x$mn.dist, x$mn.acc, method='pearson')
  c(df=res$parameter, r=res$estimate, p=res$p.value)
})

## ANOVA on sample distance

m = aov(mn.dist ~ dimtype*ruletype + ruletype/rulerelation, data=agg)
summary(m)

### Effect size
r = anova(m)
for (p in 1:4) {
  print(r$'Sum Sq'[p]/(r$'Sum Sq'[p] + r$'Sum Sq'[5]))
}

TukeyHSD(m)$`stimtype:ruletype`

# M,SD
ddply(agg, c('dimtype', 'ruletype'), function(x) {
  c(mn=mean(x$mn.dist), sd=sd(x$mn.dist))
})


## New model of sample distance with block effects

m1 = lmer(dist_med ~ dimtype*ruletype + block + (dimtype*ruletype)/rulerelation + (1|sid), 
          data=data, REML=F)

m2 = lmer(dist_med ~ dimtype*ruletype*block + (dimtype*ruletype)/rulerelation + (1|sid), 
          data=data, REML=F)
anova(m1,m2)

m3 = lmer(dist_med ~ dimtype*ruletype*block + block*((dimtype*ruletype)/rulerelation) + (1|sid), 
          data=data, REML=F)
anova(m2,m3)

m4 = lmer(dist_med ~ dimtype*ruletype*block + block*((dimtype*ruletype)/rulerelation) + prev_acc + (1|sid), 
          data=data, REML=F)
anova(m3,m4)

round(c(AIC(m1), AIC(m2), AIC(m3), AIC(m4)))


m4 = lmer(dist_med ~ dimtype*ruletype*block + block*((dimtype*ruletype)/rulerelation) + prev_acc + (1|sid), 
          data=data, REML=T)
m = m4

summary(m)
confint(m)





data$pred_dist = predict(m, type='response')

r = ddply(data, c('dimtype', 'ruletype', 'rulerelation', 'block'), function(x) { 
  mn = mean(x$dist_med) 
  mn_pr = mean(x$pred_dist)
  data.frame(obs=mn, pred=mn_pr)
})

p1 = ggplot(data=r[r$dimtype=='absolute' & r$ruletype=='1D',]) +
  geom_line(aes(x=block, y=obs, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0.4, 1.1)

p2 = ggplot(data=r[r$dimtype=='relative' & r$ruletype=='1D',]) +
  geom_line(aes(x=block, y=obs, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0.4, 1.1)

p3 = ggplot(data=r[r$dimtype=='absolute' & r$ruletype=='2D',]) +
  geom_line(aes(x=block, y=obs, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0.4, 1.1)

p4 = ggplot(data=r[r$dimtype=='relative' & r$ruletype=='2D',]) +
  geom_line(aes(x=block, y=obs, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0.4, 1.1)

grid.arrange(p1, p2, p3, p4)




m = m4

D1_ABS_CHEM       = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_ABS_FERT       = c(1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
D1_REL_CHEM       = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_REL_FERT       = c(1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)
D2_ABS_POS        = c(1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D2_ABS_NEG        = c(1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
D2_REL_POS        = c(1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0)
D2_REL_NEG        = c(1,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0)
block_D1_ABS_CHEM = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
block_D1_ABS_FERT = c(0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0)
block_D1_REL_CHEM = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0)
block_D1_REL_FERT = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0)
block_D2_ABS_POS  = c(0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0)
block_D2_ABS_NEG  = c(0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0)
block_D2_REL_POS  = c(0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0)
block_D2_REL_NEG  = c(0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,1)

prev_acc = c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)

c1 = rbind("D1_ABS - D1_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D1_REL_CHEM + D1_REL_FERT)/2),
           "D1_ABS - D2_ABS" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_ABS - D2_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D1_REL - D2_ABS" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_REL - D2_REL" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D2_ABS - D2_REL" = ((D2_ABS_POS + D2_ABS_NEG)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "prev_acc" = prev_acc,
           "block_D1_ABS_CHEM" = block_D1_ABS_CHEM,
           "block_D1_ABS_FERT" = block_D1_ABS_FERT,
           "block_D1_REL_CHEM" = block_D1_REL_CHEM,
           "block_D1_REL_FERT" = block_D1_REL_FERT,
           "block_D2_ABS_POS" = block_D2_ABS_POS,
           "block_D2_ABS_NEG" = block_D2_ABS_NEG,
           "block_D2_REL_POS" = block_D2_REL_POS,
           "block_D2_REL_NEG" = block_D2_REL_NEG
)
summary(glht(m, c1), test = adjusted("bonferroni"))
confint(glht(m, c1))




# Proportion of 2D rules ----

# case of quasi-complete separation -- in the 2D-RECT-shape condition,
# all 100% 2D boundaries. So need to use Firth's bias correction
# https://web.warwick.ac.uk/statsdept/user-2011/TalkSlides/Contributed/18Aug_0950_FocusVI_3-GLM-3-Kosmidis.pdf
library(brglm)

library(rms)


ddply(data, c('ruletype', 'dimtype', 'block'), function(x) {
  c(nMatch=sum(x$hypMatch))
})


m1 = brglm(hyp2D ~ ruletype*dimtype + block + (dimtype*ruletype)/rulerelation,
           family=binomial, data=data, REML=F)
summary(m1)

m2 = brglm(hyp2D ~ ruletype*dimtype*block + ruletype/rulerelation,
           family=binomial, data=data, REML=F)
summary(m2)
anova(m1,m2)

m3 = brglm(hyp2D ~ ruletype*dimtype*block + block*((dimtype*ruletype)/rulerelation),
           family=binomial, data=data, REML=F)
summary(m3)
anova(m2,m3)



#m1 = brglm(hypMatch ~ ruletype*dimtype + block + (dimtype*ruletype)/rulerelation,
#           family=binomial, data=data, REML=F)
#summary(m1)

#m2 = brglm(hypMatch ~ ruletype*dimtype*block + ruletype/rulerelation,
#           family=binomial, data=data, REML=F)
#summary(m2)
#anova(m1,m2)

#m3 = brglm(hypMatch ~ ruletype*dimtype*block + block*((dimtype*ruletype)/rulerelation),
#           family=binomial, data=data, REML=F)
#summary(m3)
#anova(m2,m3)


# recoding based on whether dimension type matches rule
#m1 = brglm(hypMatch ~ repMatch*ruletype + block + ruletype/rulerelation,
#           family=binomial, data=data, REML=F)
#test = remove_extraneous_predictors(m1)
#summary(m1)

#m2 = brglm(hypMatch ~ repMatch*ruletype*block + ruletype/rulerelation,
#           family=binomial, data=data, REML=F)
#m2 = remove_extraneous_predictors(m2)
#summary(m2)


#m3 = brglm(hypMatch ~ repMatch*ruletype*block + block*(ruletype/rulerelation),
#           family=binomial, data=data, REML=F)
#summary(m3)


round(c(AIC(m1), AIC(m2), AIC(m3)))

m = m3


f = m
M = model.matrix(f)
# Save the "assign" and "contrasts" attributes of the model matrix.
attr.assign = attr(M, "assign")
attr.contrasts = attr(M, "contrasts")
rmvCols = c(1, which(summary(f)$aliased))
M = M[,-rmvCols] # Remove intercept column (required) and NA columns
# Reassign the "assign" attribute, removing the corresponding elements from it.
attr(M, "assign") = attr.assign[-rmvCols]
# Reassign the "contrasts" attribute to its original value
attr(M, "contrasts") = attr.contrasts

f2 = update(f, ~M)
summary(f2)
confint.default(f2)

m = f2

# visualize predictions
data$pred_hyp2D = predict(m, type='response')

r = ddply(data, c('dimtype', 'ruletype', 'block', 'rulerelation'), function(x) { 
  mn = mean(x$hyp2D) 
  mn_pr = mean(x$pred_hyp2D)
  data.frame(obs=mn, pred=mn_pr)
})

p1 = ggplot(data=r[r$dimtype=='absolute' & r$ruletype=='1D',]) +
  geom_line(aes(x=block, y=obs, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0, 1)

p2 = ggplot(data=r[r$dimtype=='relative' & r$ruletype=='1D',]) +
  geom_line(aes(x=block, y=obs, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0, 1)

p3 = ggplot(data=r[r$dimtype=='absolute' & r$ruletype=='2D',]) +
  geom_line(aes(x=block, y=obs, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0, 1)

p4 = ggplot(data=r[r$dimtype=='relative' & r$ruletype=='2D',]) +
  geom_line(aes(x=block, y=obs, col=rulerelation)) +
  geom_line(aes(x=block, y=pred, col=rulerelation), linetype=2) +
  ylim(0, 1)

grid.arrange(p1, p2, p3, p4)








D1_ABS_CHEM       = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_ABS_FERT       = c(1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
D1_REL_CHEM       = c(1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_REL_FERT       = c(1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0)
D2_ABS_POS        = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D2_ABS_NEG        = c(1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
D2_REL_POS        = c(1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0)
D2_REL_NEG        = c(1,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0)
block_D1_ABS_CHEM = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
block_D1_ABS_FERT = c(0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0)
block_D1_REL_CHEM = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0)
block_D1_REL_FERT = c(0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0)
block_D2_ABS_POS  = c(0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0)
block_D2_ABS_NEG  = c(0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0)
block_D2_REL_POS  = c(0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0)
block_D2_REL_NEG  = c(0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,1)


c1 = rbind("D1_ABS - D1_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D1_REL_CHEM + D1_REL_FERT)/2),
           "D1_ABS - D2_ABS" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_ABS - D2_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D1_REL - D2_ABS" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_REL - D2_REL" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D2_ABS - D2_REL" = ((D2_ABS_POS + D2_ABS_NEG)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D1_ABS_CHEM - D1_ABS_FERT" = D1_ABS_CHEM - D1_ABS_FERT,
           "D1_REL_CHEM - D1_REL_FERT" = D1_REL_CHEM - D1_REL_FERT,
           "D2_ABS_POS - D2_ABS_NEG" = D2_ABS_POS - D2_ABS_NEG,
           "D2_REL_POS - D2_REL_NEG" = D2_REL_POS - D2_REL_NEG,
           "block_D1_ABS_CHEM" = block_D1_ABS_CHEM,
           "block_D1_ABS_FERT" = block_D1_ABS_FERT,
           "block_D1_REL_CHEM" = block_D1_REL_CHEM,
           "block_D1_REL_FERT" = block_D1_REL_FERT,
           "block_D2_ABS_POS" = block_D2_ABS_POS,
           "block_D2_ABS_NEG" = block_D2_ABS_NEG,
           "block_D2_REL_POS" = block_D2_REL_POS,
           "block_D2_REL_NEG" = block_D2_REL_NEG
)
summary(glht(m, c1), test = adjusted("bonferroni"))
#summary(glht(m, c1), test = adjusted("none"))







data = read.csv('../../data_proc/experiment2/prop2d.csv')
data$sid = factor(data$sid)
data$cond = factor(data$cond)
data$dimtype = relevel(factor(data$dimtype), ref='absolute')
data$ruletype = relevel(factor(data$ruletype), ref='1D')
data$rulerelation = factor(data$rulerelation)
data$n2D = as.integer(data$X2drule_aic*8)




# case of quasi-complete separation -- in the 2D-RECT-shape condition,
# all 100% 2D boundaries. So need to use Firth's bias correction
# https://web.warwick.ac.uk/statsdept/user-2011/TalkSlides/Contributed/18Aug_0950_FocusVI_3-GLM-3-Kosmidis.pdf
library(brglm)

m4 = brglm(X2drule_aic ~ ruletype*dimtype + ruletype/(dimtype*rulerelation),
           weights=rep(8, nrow(data)),
           family=binomial, data=data)
summary(m4)


f = m4
M = model.matrix(f)
# Save the "assign" and "contrasts" attributes of the model matrix.
attr.assign = attr(M, "assign")
attr.contrasts = attr(M, "contrasts")
rmvCols = c(1, which(summary(f)$aliased))
M = M[,-rmvCols] # Remove intercept column (required) and NA columns
# Reassign the "assign" attribute, removing the corresponding elements from it.
attr(M, "assign") = attr.assign[-rmvCols]
# Reassign the "contrasts" attribute to its original value
attr(M, "contrasts") = attr.contrasts

f2 = update(f, ~M)
summary(f2)

confint.default(f2)


# contrasts
D1_ABS_CHEM  = c(1,0,0,0,0,0,0,0)
D1_ABS_FERT  = c(1,0,0,0,1,0,0,0)
D1_REL_CHEM  = c(1,0,1,0,0,0,0,0)
D1_REL_FERT  = c(1,0,1,0,1,0,1,0)
D2_ABS_POS   = c(1,1,0,0,0,0,0,0)
D2_ABS_NEG   = c(1,1,0,0,0,1,0,0)
D2_REL_POS   = c(1,1,1,1,0,0,0,0)
D2_REL_NEG   = c(1,1,1,1,0,1,0,1)


c1 = rbind("D1_ABS - D1_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D1_REL_CHEM + D1_REL_FERT)/2),
           "D1_ABS - D2_ABS" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_ABS - D2_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D1_REL - D2_ABS" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_REL - D2_REL" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D2_ABS - D2_REL" = ((D2_ABS_POS + D2_ABS_NEG)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D1_ABS_CHEM - D1_ABS_FERT" = D1_ABS_CHEM - D1_ABS_FERT,
           "D1_REL_CHEM - D1_REL_FERT" = D1_REL_CHEM - D1_REL_FERT,
           "D2_ABS_POS - D2_ABS_NEG" = D2_ABS_POS - D2_ABS_NEG,
           "D2_REL_POS - D2_REL_NEG" = D2_REL_POS - D2_REL_NEG
)

summary(glht(f2, c1), test = adjusted("bonferroni"))

