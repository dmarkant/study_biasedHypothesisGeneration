library(lme4)
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
data$block = scale(data$block)


# Classification accuracy ----

# the model
m = glmer(acc ~ dimtype*ruletype + block + ruletype/rulerelation + (1|sid),
          weights=rep(32, nrow(data)),
          data=data, family=binomial, control = glmerControl(optimizer = "bobyqa"),
          nAGQ = 10)
summary(m)

# odds ratios
exp(fixef(m))

# confidence intervals
confint(m)

# 2.5 %      97.5 %
#   .sig01                       0.9682921  1.30254889
# (Intercept)                  2.6476253  3.67729833
# dimtyperelative             -1.7843059 -0.57111997
# ruletype2D                  -2.3176956 -0.88011963
# block                        0.3480856  0.41844115
# dimtyperelative:ruletype2D   1.2200222  2.90207173
# ruletype1D:rulerelation1D-F  0.2866626  1.51048528
# ruletype2D:rulerelation2D-N -1.0893521  0.07511386


# no interaction 
m1 = glmer(acc ~ dimtype*ruletype + block + ruletype/(dimtype*rulerelation) + (1|sid),
           weights=rep(32, nrow(data)),
           data=data, family=binomial, control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 10)
summary(m1)
anova(m,m1)


# contrasts
D1_ABS_CHEM  = c(1,0,0,0,0,0,0)
D1_ABS_FERT  = c(1,0,0,0,0,1,0)
D1_REL_CHEM  = c(1,1,0,0,0,0,0)
D1_REL_FERT  = c(1,1,0,0,0,1,0)
D2_ABS_POS   = c(1,0,1,0,0,0,0)
D2_ABS_NEG   = c(1,0,1,0,0,0,1)
D2_REL_POS   = c(1,1,1,0,1,0,0)
D2_REL_NEG   = c(1,1,1,0,1,0,1)

c1 = rbind("D1_ABS - D1_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D1_REL_CHEM + D1_REL_FERT)/2),
           "D1_ABS - D2_ABS" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_ABS - D2_REL" = ((D1_ABS_CHEM + D1_ABS_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D1_REL - D2_ABS" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_ABS_POS + D2_ABS_NEG)/2),
           "D1_REL - D2_REL" = ((D1_REL_CHEM + D1_REL_FERT)/2 - (D2_REL_POS + D2_REL_NEG)/2),
           "D2_ABS - D2_REL" = ((D2_ABS_POS + D2_ABS_NEG)/2 - (D2_REL_POS + D2_REL_NEG)/2))

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


# Proportion of 2D rules ----
library(rms)

data = read.csv('../../data_proc/experiment2/prop2d.csv')
data$sid = factor(data$sid)
data$cond = factor(data$cond)
data$dimtype = relevel(factor(data$dimtype), ref='absolute')
data$ruletype = relevel(factor(data$ruletype), ref='1D')
data$rulerelation = factor(data$rulerelation)
data$n2D = as.integer(data$X2drule_aic*8)


# M (SD)
ddply(data, c('ruletype', 'dimtype'), function(x) {
  c(mn=mean(x$n2D/8), sd=sd(x$n2D/8))
})


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

