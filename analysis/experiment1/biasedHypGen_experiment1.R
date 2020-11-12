library(lme4)
library(lmerTest)
library(plyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(multcomp)

setwd("~/studies/study_biasedHypothesisGeneration/analysis/experiment1")


data = read.csv('../../data_proc/experiment1/data_by_block.csv')
data$sid = factor(data$sid)
data$cond = factor(data$cond)
data$stimtype = relevel(factor(data$stimtype), ref='antenna')
data$ruletype = relevel(factor(data$ruletype), ref='rb')
data$rulevariant = factor(data$rulevariant)

#data$block = scale(data$block)
data$block = data$block - mean(data$block)

data$hyp2D = as.logical(data$X2drule_aic)
data$hypMatch = ((data$ruletype=='rb') & (!data$hyp2D)) | ((data$ruletype=='ii') & (data$hyp2D))
data$repMatch = ((data$ruletype=='rb') & (data$stimtype=='antenna')) | ((data$ruletype=='ii') & (data$stimtype=='rectangle'))



data$prev_acc = NaN
for (sid in levels(data$sid)) {
  acc = data[data$sid==sid,]$acc
  prev_acc = c(0.5, acc[1:7])
  data[data$sid==sid, 'prev_acc'] = prev_acc
}
data$change_acc = data$acc - data$prev_acc

data$prev_acc = scale(data$prev_acc)
data$change_acc = scale(data$change_acc)



# Classification accuracy ----


# comparisons
m1 = glmer(acc ~ stimtype*ruletype + block + ruletype/rulevariant + (1|sid), 
          weights=rep(32, nrow(data)), data=data, family=binomial, 
          control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)

m2 = glmer(acc ~ stimtype*ruletype*block + ruletype/rulevariant + (1|sid), 
           weights=rep(32, nrow(data)), data=data, family=binomial, 
           control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
anova(m1,m2)

m3 = glmer(acc ~ stimtype*ruletype*block + block*(ruletype/rulevariant) + (1|sid), 
           weights=rep(32, nrow(data)), data=data, family=binomial, 
           control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
anova(m2,m3)
summary(m3)



round(c(AIC(m1), AIC(m2), AIC(m3)))


m = m3
summary(m)

# odds ratios
exp(fixef(m))

# confidence intervals
confint(m)

# 2.5 %      97.5 %
#   .sig01                              0.48496421  0.65821509
# (Intercept)                         3.60790583  4.54862669
# stimtyperectangle                  -3.20740272 -1.95888673
# ruletypeii                         -3.76040084 -2.62464620
# block                               0.07334379  0.32450374
# stimtyperectangle:ruletypeii        2.30636824  3.85959085
# stimtyperectangle:block            -0.23672642  0.03284893
# ruletypeii:block                   -0.23013803  0.02971146
# ruletyperb:rulevariantheight       -0.40347430  0.57731157
# ruletypeii:rulevariantneg          -0.26758294  0.59131270
# ruletyperb:rulevariantradius       -2.06937273 -0.95770916
# ruletypeii:rulevariantshape         1.33387031  2.23664463
# stimtyperectangle:ruletypeii:block -0.03988751  0.24843060
# ruletyperb:block:rulevariantheight -0.02792431  0.09022549
# ruletypeii:block:rulevariantneg    -0.10661602 -0.01579657
# ruletyperb:block:rulevariantradius -0.22134478  0.04734540
# ruletypeii:block:rulevariantshape  -0.08888605  0.05536141

ns = ddply(data, c('stimtype', 'ruletype', 'rulevariant'), 
      function(x) { c(n=nrow(x)) })


D1_DIAL_ANGLE  = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_DIAL_RADIUS = c(1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)
D1_RECT_WIDTH  = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_RECT_HEIGHT = c(1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
D2_DIAL_POS    = c(1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
D2_DIAL_NEG    = c(1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0)
D2_RECT_SIZE   = c(1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0)
D2_RECT_SHAPE  = c(1,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0)

block_D1_DIAL_ANGLE   = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
block_D1_DIAL_RADIUS  = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0)
block_D1_RECT_WIDTH   = c(0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0)
block_D1_RECT_HEIGHT  = c(0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0)
block_D2_DIAL_POS     = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0)
block_D2_DIAL_NEG     = c(0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0)
block_D2_RECT_SIZE    = c(0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0)
block_D2_RECT_SHAPE   = c(0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,1)

c1 = rbind("D1_DIAL - D1_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D1_RECT_WIDTH + D1_RECT_HEIGHT)/2),
           "D1_DIAL - D2_DIAL" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
           "D1_DIAL - D2_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D1_RECT - D2_DIAL" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
           "D1_RECT - D2_RECT" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D2_DIAL - D2_RECT" = ((D2_DIAL_POS + D2_DIAL_NEG)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D1_DIAL_ANGLE - D1_DIAL_RADIUS" = D1_DIAL_ANGLE - D1_DIAL_RADIUS,
           "D1_RECT_WIDTH - D1_RECT_HEIGHT" = D1_RECT_WIDTH - D1_RECT_HEIGHT,
           "D2_DIAL_POS - D2_DIAL_NEG" = D2_DIAL_POS - D2_DIAL_NEG,
           "D2_RECT_SIZE - D2_RECT_SHAPE" = D2_RECT_SIZE - D2_RECT_SHAPE,
           "D2_DIAL - D2_RECT_SIZE" = (D2_DIAL_POS + D2_DIAL_NEG)/2 - D2_RECT_SIZE,
           "block_D1_DIAL_ANGLE" = block_D1_DIAL_ANGLE,
           "block_D1_DIAL_RADIUS" = block_D1_DIAL_RADIUS,
           "block_D1_RECT_WIDTH" = block_D1_RECT_WIDTH,
           "block_D1_RECT_HEIGHT" = block_D1_RECT_HEIGHT,
           "block_D2_DIAL_POS" = block_D2_DIAL_POS,
           "block_D2_DIAL_NEG" = block_D2_DIAL_NEG,
           "block_D2_RECT_SIZE" = block_D2_RECT_SIZE,
           "block_D2_RECT_SHAPE" = block_D2_RECT_SHAPE
           )

mc = glht(m, c1)
summary(glht(m, c1), test = adjusted("bonferroni"))
confint(glht(m, c1))

#ctr = "stimtyperectangle = 0"
#ctr = "block + stimtyperectangle:block + ruletypeii:block + stimtyperectangle:ruletypeii:block + ruletyperb:block:rulevariantheight + ruletypeii:block:rulevariantneg + ruletyperb:block:rulevariantradius + ruletypeii:block:rulevariantshape = 0"
#g = glht(m, linfct = c(ctr))
#summary(glht(m, linfct = c(ctr)))

# D1_DIAL_ANGLE  = c(1,0,0,0,0,0,0,0,0)
# D1_DIAL_RADIUS = c(1,0,0,0,0,0,0,1,0)
# D1_RECT_WIDTH  = c(1,1,0,0,0,0,0,0,0)
# D1_RECT_HEIGHT = c(1,1,0,0,0,1,0,0,0)
# D2_DIAL_POS    = c(1,0,1,0,0,0,0,0,0)
# D2_DIAL_NEG    = c(1,0,1,0,0,0,1,0,0)
# D2_RECT_SIZE   = c(1,1,1,0,1,0,0,0,0)
# D2_RECT_SHAPE  = c(1,1,1,0,1,0,0,0,1)
# 
# 
# c1 = rbind("D1_DIAL - D1_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D1_RECT_WIDTH + D1_RECT_HEIGHT)/2),
#            "D1_DIAL - D2_DIAL" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
#            "D1_DIAL - D2_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
#            "D1_RECT - D2_DIAL" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
#            "D1_RECT - D2_RECT" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
#            "D2_DIAL - D2_RECT" = ((D2_DIAL_POS + D2_DIAL_NEG)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
#            "D1_DIAL_ANGLE - D1_DIAL_RADIUS" = D1_DIAL_ANGLE - D1_DIAL_RADIUS,
#            "D1_RECT_WIDTH - D1_RECT_HEIGHT" = D1_RECT_WIDTH - D1_RECT_HEIGHT,
#            "D2_DIAL_POS - D2_DIAL_NEG" = D2_DIAL_POS - D2_DIAL_NEG,
#            "D2_RECT_SIZE - D2_RECT_SHAPE" = D2_RECT_SIZE - D2_RECT_SHAPE,
#            "D2_DIAL_NEG - D2_RECT_SIZE" = D2_DIAL_NEG - D2_RECT_SIZE
#            )
# summary(glht(m, c1), test = adjusted("bonferroni"))
# confint(glht(m, c1))


y = ddply(data, c('stimtype', 'ruletype', 'rulevariant', 'sid'), function(x) {
  mn = mean(x$acc)
  data.frame(acc=mn)
})

# by condition
ddply(y, c('stimtype', 'ruletype'), function(x) {
  c(mn=mean(x$acc), sd=sd(x$acc))
})

# by rule variant
ddply(y, c('stimtype', 'ruletype', 'rulevariant'), function(x) {
  c(mn=mean(x$acc), sd=sd(x$acc))
})



# Visualize model fit

data$pred_acc = predict(m, type='response')

r = ddply(data, c('stimtype', 'ruletype', 'rulevariant', 'block'), function(x) { 
  mn = mean(x$acc) 
  mn_pr = mean(x$pred_acc)
  data.frame(acc=mn, pred=mn_pr)
})

p1 = ggplot(data=r[r$stimtype=='antenna' & r$ruletype=='rb',]) +
            geom_line(aes(x=block, y=acc, col=rulevariant)) +
            geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
            ylim(0.5, 1)

p2 = ggplot(data=r[r$stimtype=='rectangle' & r$ruletype=='rb',]) +
  geom_line(aes(x=block, y=acc, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0.5, 1)

p3 = ggplot(data=r[r$stimtype=='antenna' & r$ruletype=='ii',]) +
  geom_line(aes(x=block, y=acc, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0.5, 1)

p4 = ggplot(data=r[r$stimtype=='rectangle' & r$ruletype=='ii',]) +
  geom_line(aes(x=block, y=acc, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0.5, 1)

grid.arrange(p1, p2, p3, p4)


# Selection distance ----

## Correlations with accuracy ----

agg = ddply(data, c('sid', 'stimtype', 'ruletype', 'rulevariant'), function(x) {
  c(mn.acc=mean(x$acc),
    mn.dist=mean(x$dist_med))
})

ddply(agg, c('stimtype', 'ruletype'), function(x) {
  res = cor.test(x$mn.dist, x$mn.acc, method='pearson')
  c(df=res$parameter, r=res$estimate, p=res$p.value)
})


## ANOVA on sample distance

#m = aov(mn.dist ~ stimtype*ruletype + ruletype/rulevariant, data=agg)
#summary(m)

### Effect size
#r = anova(m)
#for (p in 1:4) {
#  print(r$'Sum Sq'[p]/(r$'Sum Sq'[p] + r$'Sum Sq'[5]))
#}

#TukeyHSD(m)$`stimtype:ruletype`

# M,SD
ddply(agg, c('stimtype', 'ruletype'), function(x) {
  c(mn=mean(x$mn.dist), sd=sd(x$mn.dist))
})


## New model of sample distance with block effects

m1 = lmer(dist_med ~ stimtype*ruletype + block + ruletype/rulevariant + (1|sid),
         data=data, REML=F)
summary(m1)

m2 = lmer(dist_med ~ stimtype*ruletype*block + ruletype/rulevariant + (1|sid),
          data=data, REML=F)
summary(m2)
anova(m1,m2)

m3 = lmer(dist_med ~ stimtype*ruletype*block + block*(ruletype/rulevariant) + (1|sid),
          data=data, REML=F)
summary(m3)
anova(m2,m3)

round(c(AIC(m1), AIC(m2), AIC(m3)))

m4 = lmer(dist_med ~ stimtype*ruletype*block + block*(ruletype/rulevariant) + prev_acc + (1|sid),
          data=data, REML=F)
summary(m4)
anova(m3,m4)

m4 = lmer(dist_med ~ stimtype*ruletype*block + block*(ruletype/rulevariant) + prev_acc + (1|sid),
          data=data, REML=T)
m = m4
summary(m)
confint(m)

# 2.5 %      97.5 %
#   .sig01                              0.21189653  0.29330884
# .sigma                              0.32150536  0.35485063
# (Intercept)                         0.48756370  0.83623606
# stimtyperectangle                   0.08485782  0.60609076
# ruletypeii                         -0.03603996  0.43370172
# block                              -0.09067747 -0.02388367
# prev_acc                           -0.18211377 -0.11288500
# stimtyperectangle:ruletypeii       -0.76735589 -0.08342197
# stimtyperectangle:block            -0.04388157  0.05308460
# ruletypeii:block                   -0.00517781  0.08100003
# ruletyperb:rulevariantheight       -0.10294162  0.35701781
# ruletypeii:rulevariantneg          -0.15892415  0.24887438
# ruletyperb:rulevariantradius       -0.08435683  0.35129348
# ruletypeii:rulevariantshape        -0.15081545  0.25994084
# stimtyperectangle:ruletypeii:block -0.04694519  0.08005011
# ruletyperb:block:rulevariantheight -0.06095123  0.02518902
# ruletypeii:block:rulevariantneg    -0.05060385  0.02571199
# ruletyperb:block:rulevariantradius -0.05531600  0.02613798
# ruletypeii:block:rulevariantshape  -0.10085674 -0.02467738


D1_DIAL_ANGLE  = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_DIAL_RADIUS = c(1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)
D1_RECT_WIDTH  = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_RECT_HEIGHT = c(1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
D2_DIAL_POS    = c(1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D2_DIAL_NEG    = c(1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
D2_RECT_SIZE   = c(1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0)
D2_RECT_SHAPE  = c(1,1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0)

block_D1_DIAL_ANGLE   = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
block_D1_DIAL_RADIUS  = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0)
block_D1_RECT_WIDTH   = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0)
block_D1_RECT_HEIGHT  = c(0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0)
block_D2_DIAL_POS     = c(0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0)
block_D2_DIAL_NEG     = c(0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0)
block_D2_RECT_SIZE    = c(0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0)
block_D2_RECT_SHAPE   = c(0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,0,1)

block = (block_D1_DIAL_ANGLE + block_D1_DIAL_RADIUS + 
         block_D1_RECT_HEIGHT + block_D1_RECT_WIDTH + 
         block_D2_DIAL_NEG + block_D2_DIAL_POS +
         block_D2_RECT_SHAPE + block_D2_RECT_SIZE)/8

prev_acc = c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)

c1 = rbind("D1_DIAL - D1_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D1_RECT_WIDTH + D1_RECT_HEIGHT)/2),
           "D1_DIAL - D2_DIAL" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
           "D1_DIAL - D2_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D1_RECT - D2_DIAL" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
           "D1_RECT - D2_RECT" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D2_DIAL - D2_RECT" = ((D2_DIAL_POS + D2_DIAL_NEG)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D1_DIAL_ANGLE - D1_DIAL_RADIUS" = D1_DIAL_ANGLE - D1_DIAL_RADIUS,
           "D1_RECT_WIDTH - D1_RECT_HEIGHT" = D1_RECT_WIDTH - D1_RECT_HEIGHT,
           "D2_DIAL_POS - D2_DIAL_NEG" = D2_DIAL_POS - D2_DIAL_NEG,
           "D2_RECT_SIZE - D2_RECT_SHAPE" = D2_RECT_SIZE - D2_RECT_SHAPE,
           "prev_acc" = prev_acc,
           "block_D1_DIAL_ANGLE" = block_D1_DIAL_ANGLE,
           "block_D1_DIAL_RADIUS" = block_D1_DIAL_RADIUS,
           "block_D1_RECT_WIDTH" = block_D1_RECT_WIDTH,
           "block_D1_RECT_HEIGHT" = block_D1_RECT_HEIGHT,
           "block_D2_DIAL_POS" = block_D2_DIAL_POS,
           "block_D2_DIAL_NEG" = block_D2_DIAL_NEG,
           "block_D2_RECT_SIZE" = block_D2_RECT_SIZE,
           "block_D2_RECT_SHAPE" = block_D2_RECT_SHAPE
)

summary(glht(m, c1), test = adjusted("bonferroni"))




m1 = lmer(dist_med ~ stimtype*ruletype + block + (1|sid),
         data=data)

m2 = lmer(dist_med ~ stimtype*ruletype*block + (1|sid),
         data=data)
anova(m1,m2)

m3 = lmer(dist_med ~ stimtype*ruletype + block + ruletype/(stimtype*rulevariant) + (1|sid),
          data=data)
anova(m2,m3)

m4 = lmer(dist_med ~ stimtype*ruletype*block + block*(ruletype/(stimtype*rulevariant)) + (1|sid),
          data=data)
summary(m4)
anova(m3,m4)

m = m4
D1_DIAL_ANGLE  = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_DIAL_RADIUS = c(1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)
D1_RECT_WIDTH  = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
D1_RECT_HEIGHT = c(1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
D2_DIAL_POS    = c(1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
D2_DIAL_NEG    = c(1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0)
D2_RECT_SIZE   = c(1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
D2_RECT_SHAPE  = c(1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0)

block_D1_DIAL_ANGLE   = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
block_D1_DIAL_RADIUS  = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0)
block_D1_RECT_WIDTH   = c(0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0)
block_D1_RECT_HEIGHT  = c(0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0)
block_D2_DIAL_POS     = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0)
block_D2_DIAL_NEG     = c(0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0)
block_D2_RECT_SIZE    = c(0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0)
block_D2_RECT_SHAPE   = c(0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,1)

c1 = rbind("D1_DIAL - D1_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D1_RECT_WIDTH + D1_RECT_HEIGHT)/2),
           "D1_DIAL - D2_DIAL" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
           "D1_DIAL - D2_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D1_RECT - D2_DIAL" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
           "D1_RECT - D2_RECT" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D2_DIAL - D2_RECT" = ((D2_DIAL_POS + D2_DIAL_NEG)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D1_DIAL_ANGLE - D1_DIAL_RADIUS" = D1_DIAL_ANGLE - D1_DIAL_RADIUS,
           "D1_RECT_WIDTH - D1_RECT_HEIGHT" = D1_RECT_WIDTH - D1_RECT_HEIGHT,
           "D2_DIAL_POS - D2_DIAL_NEG" = D2_DIAL_POS - D2_DIAL_NEG,
           "D2_RECT_SIZE - D2_RECT_SHAPE" = D2_RECT_SIZE - D2_RECT_SHAPE,
           "D2_DIAL_NEG - D2_RECT_SIZE" = D2_DIAL_NEG - D2_RECT_SIZE,
           "block_D1_DIAL_ANGLE - block_D1_DIAL_RADIUS" = block_D1_DIAL_ANGLE - block_D1_DIAL_RADIUS,
           "block_D1_RECT_WIDTH - block_D1_RECT_HEIGHT" = block_D1_RECT_WIDTH - block_D1_RECT_HEIGHT,
           "block_D2_DIAL_POS - block_D2_DIAL_NEG" = block_D2_DIAL_POS - block_D2_DIAL_NEG,
           "block_D2_RECT_SIZE - block_D2_RECT_SHAPE" = block_D2_RECT_SIZE - block_D2_RECT_SHAPE)


summary(glht(m, c1), test = adjusted("bonferroni"))
confint(glht(m, c1))

# Visualize model fit

data$pred_dist = predict(m4, type='response')

r = ddply(data, c('stimtype', 'ruletype', 'block', 'rulevariant'), function(x) { 
  mn = mean(x$dist_med) 
  mn_pr = mean(x$pred_dist)
  data.frame(obs=mn, pred=mn_pr)
})

p1 = ggplot(data=r[r$stimtype=='antenna' & r$ruletype=='rb',]) +
  geom_line(aes(x=block, y=obs, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0, 2)

p2 = ggplot(data=r[r$stimtype=='rectangle' & r$ruletype=='rb',]) +
  geom_line(aes(x=block, y=obs, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0, 2)

p3 = ggplot(data=r[r$stimtype=='antenna' & r$ruletype=='ii',]) +
  geom_line(aes(x=block, y=obs, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0, 2)

p4 = ggplot(data=r[r$stimtype=='rectangle' & r$ruletype=='ii',]) +
  geom_line(aes(x=block, y=obs, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0, 2)

grid.arrange(p1, p2, p3, p4)



# Proportion of 2D rules ----

# case of quasi-complete separation -- in the 2D-RECT-shape condition,
# all 100% 2D boundaries. So need to use Firth's bias correction
# https://web.warwick.ac.uk/statsdept/user-2011/TalkSlides/Contributed/18Aug_0950_FocusVI_3-GLM-3-Kosmidis.pdf
library(brglm)

ddply(data, c('ruletype', 'stimtype', 'block'), function(x) {
  c(nMatch=sum(x$hypMatch))
})



m1 = brglm(hyp2D ~ ruletype*stimtype + block + ruletype/rulevariant,
           family=binomial, data=data, REML=F)
summary(m1)

m2 = brglm(hyp2D ~ ruletype*stimtype*block + ruletype/rulevariant,
           family=binomial, data=data, REML=F)
summary(m2)
anova(m1,m2)

m3 = brglm(hyp2D ~ ruletype*stimtype*block + block*(ruletype/rulevariant),
           family=binomial, data=data, REML=F)
summary(m3)
anova(m2,m3)

#m1 = brglm(hypMatch ~ ruletype*stimtype + block + ruletype/rulevariant,
#           family=binomial, data=data, REML=F)
#summary(m1)

#m2 = brglm(hypMatch ~ ruletype*stimtype*block + ruletype/rulevariant,
#           family=binomial, data=data, REML=F)
#summary(m2)
#anova(m1,m2)

#m3 = brglm(hypMatch ~ ruletype*stimtype*block + block*(ruletype/rulevariant),
#           family=binomial, data=data, REML=F)
#summary(m3)
#anova(m2,m3)


# recoding based on whether dimension type matches rule
#m1 = brglm(hypMatch ~ repMatch*ruletype + block + ruletype/rulevariant,
#           family=binomial, data=data, REML=F)
#test = remove_extraneous_predictors(m1)
#summary(m1)

#m2 = brglm(hypMatch ~ repMatch*ruletype*block + ruletype/rulevariant,
#           family=binomial, data=data, REML=F)
#m2 = remove_extraneous_predictors(m2)
#summary(m2)


#m3 = brglm(hypMatch ~ repMatch*ruletype*block + block*(ruletype/rulevariant),
#           family=binomial, data=data, REML=F)
#summary(m3)

round(c(AIC(m1), AIC(m2), AIC(m3)))


m = m2


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




D1_DIAL_ANGLE  = c(1,0,0,0,0,0,0,0,0,0,0,0)
D1_DIAL_RADIUS = c(1,0,0,0,0,0,0,0,0,1,0,0)
D1_RECT_HEIGHT = c(1,0,1,0,0,0,0,1,0,0,0,0)
D1_RECT_WIDTH  = c(1,0,1,0,0,0,0,0,0,0,0,0)
D2_DIAL_POS    = c(1,1,0,0,0,0,0,0,0,0,0,0)
D2_DIAL_NEG    = c(1,1,0,0,0,0,0,0,1,0,0,0)
D2_RECT_SIZE   = c(1,1,1,0,1,0,0,0,0,0,0,0)
D2_RECT_SHAPE  = c(1,1,1,0,1,0,0,0,0,0,1,0)

block_D1_DIAL   = c(0,0,0,1,0,0,0,0,0,0,0,0)
block_D1_RECT   = c(0,0,0,1,0,0,1,0,0,0,0,0)
block_D2_DIAL   = c(0,0,0,1,0,1,0,0,0,0,0,0)
block_D2_RECT   = c(0,0,0,1,0,1,1,0,0,0,0,1)


c1 = rbind("D1_DIAL - D1_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D1_RECT_WIDTH + D1_RECT_HEIGHT)/2),
           "D1_DIAL - D2_DIAL" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
           "D1_DIAL - D2_RECT" = ((D1_DIAL_ANGLE + D1_DIAL_RADIUS)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D1_RECT - D2_DIAL" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_DIAL_POS + D2_DIAL_NEG)/2),
           "D1_RECT - D2_RECT" = ((D1_RECT_WIDTH + D1_RECT_HEIGHT)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D2_DIAL - D2_RECT" = ((D2_DIAL_POS + D2_DIAL_NEG)/2 - (D2_RECT_SIZE + D2_RECT_SHAPE)/2),
           "D1_DIAL_ANGLE - D1_DIAL_RADIUS" = D1_DIAL_ANGLE - D1_DIAL_RADIUS,
           "D1_RECT_WIDTH - D1_RECT_HEIGHT" = D1_RECT_WIDTH - D1_RECT_HEIGHT,
           "D2_DIAL_POS - D2_DIAL_NEG" = D2_DIAL_POS - D2_DIAL_NEG,
           "D2_RECT_SIZE - D2_RECT_SHAPE" = D2_RECT_SIZE - D2_RECT_SHAPE,
           "block_D1_DIAL" = block_D1_DIAL,
           "block_D1_RECT" = block_D1_RECT,
           "block_D2_DIAL" = block_D2_DIAL,
           "block_D2_RECT" = block_D2_RECT
)
summary(glht(m, c1), test = adjusted("bonferroni"))





# visualize predictions
data$pred_hyp2D = predict(m, type='response')

r = ddply(data, c('stimtype', 'ruletype', 'block', 'rulevariant'), function(x) { 
  mn = mean(x$hyp2D) 
  mn_pr = mean(x$pred_hyp2D)
  data.frame(obs=mn, pred=mn_pr)
})

p1 = ggplot(data=r[r$stimtype=='antenna' & r$ruletype=='rb',]) +
  geom_line(aes(x=block, y=obs, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0, 1)

p2 = ggplot(data=r[r$stimtype=='rectangle' & r$ruletype=='rb',]) +
  geom_line(aes(x=block, y=obs, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0, 1)

p3 = ggplot(data=r[r$stimtype=='antenna' & r$ruletype=='ii',]) +
  geom_line(aes(x=block, y=obs, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0, 1)

p4 = ggplot(data=r[r$stimtype=='rectangle' & r$ruletype=='ii',]) +
  geom_line(aes(x=block, y=obs, col=rulevariant)) +
  geom_line(aes(x=block, y=pred, col=rulevariant), linetype=2) +
  ylim(0, 1)

grid.arrange(p1, p2, p3, p4)



library(mbest)






m1 = glmer(hypMatch ~ stimtype*ruletype + block + ruletype/rulevariant + (1|sid), 
                      weights=rep(32, nrow(data)), data=data, family=binomial)
summary(m1)

m2 = glmer(acc ~ stimtype*ruletype*block + ruletype/rulevariant + (1|sid), 
           weights=rep(32, nrow(data)), data=data, family=binomial)
summary(m2)

m = mhglm_ml(hypMatch ~ stimtype*ruletype + block + ruletype/rulevariant + (1|sid), 
          data=data, family=binomial)
summary(m)


m = mhglm_ml(acc ~ stimtype*ruletype*block + ruletype/rulevariant + (1|sid), 
             weights=rep(32, nrow(data)), data=data, family=binomial)
summary(m)


m = glmer(hypMatch ~ ruletype*repMatch*block + ruletype/rulevariant + (1|sid), 
          data=data, family=binomial)
summary(m)



m = glm(hypMatch ~ repMatch*block + ruletype/rulevariant, 
        data=data, family = binomial,
        method = "brglm.fit")
summary(m)


m = mhglm(hypMatch ~ repMatch*block + (1|sid), 
          data=data, family=binomial,
          control=mhglm.control(fit.method="firthglm.fit"))
summary(m)


(gm.firth <- mhglm(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   data = cbpp, family = binomial,
                   control=mhglm.control(fit.method="firthglm.fit")))



library(rms)





data = read.csv('../../data_proc/experiment1/prop2d.csv')
data$sid = factor(data$sid)
data$cond = factor(data$cond)
data$stimtype = relevel(factor(data$stimtype), ref='antenna')
data$ruletype = relevel(factor(data$ruletype), ref='rb')
data$rulevariant = factor(data$rulevariant)
data$n2D = as.integer(data$X2drule_aic*8)


# M (SD)
ddply(data, c('ruletype', 'stimtype'), function(x) {
  c(mn=mean(x$n2D/8), sd=sd(x$n2D/8))
})


m1 = glm(X2drule_aic ~ ruletype + stimtype,
         weights=rep(8, nrow(data)),
         family=binomial, data=data)
summary(m1)

m2 = glm(X2drule_aic ~ ruletype*stimtype,
         weights=rep(8, nrow(data)),
         family=binomial, data=data)
summary(m2)

D1_DIAL = c(1,0,0,0)
D1_RECT = c(1,0,1,0)
D2_DIAL = c(1,1,0,0)
D2_RECT = c(1,1,1,1)

c1 = rbind("D1_DIAL - D1_RECT" = D1_DIAL - D1_RECT,
           "D1_DIAL - D2_DIAL" = D1_DIAL - D2_DIAL,
           "D1_DIAL - D2_RECT" = D1_DIAL - D2_RECT,
           "D1_RECT - D2_DIAL" = D1_RECT - D2_DIAL,
           "D1_RECT - D2_RECT" = D1_RECT - D2_RECT,
           "D2_DIAL - D2_RECT" = D2_DIAL - D2_RECT
)
summary(glht(m2, c1), test = adjusted("bonferroni"))


m3 = glm(X2drule_aic ~ ruletype*stimtype + ruletype/rulevariant,
         weights=rep(8, nrow(data)),
         family=binomial, data=data)
summary(m3)





m4 = brglm(X2drule_aic ~ ruletype*stimtype + ruletype/rulevariant,
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

summary(glht(f2, c1), test = adjusted("bonferroni"))
