library(lme4)
library(plyr)
library(grid)
library(gridExtra)
require(ggplot2)
library(multcomp)

setwd("~/studies/study_biasedHypothesisGeneration/analysis/experiment1")


# Classification accuracy ----

data = read.csv('../../data_proc/experiment1/data_by_block.csv')
data$sid = factor(data$sid)
data$cond = factor(data$cond)
data$stimtype = factor(data$stimtype)
data$ruletype = relevel(factor(data$ruletype), ref='rb')
data$rulevariant = factor(data$rulevariant)
data$block = scale(data$block)

# the model
m = glmer(acc ~ stimtype*ruletype + block + ruletype/rulevariant + (1|sid), 
          weights=rep(32, nrow(data)), data=data, family=binomial, 
          control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m)

# odds ratios
exp(fixef(m))

# confidence intervals
confint(m)

# 2.5 %     97.5 %
#   .sig01                        0.4842714  0.6572883
# (Intercept)                   3.5486993  4.4579317
# stimtyperectangle            -3.1221632 -1.8981394
# ruletypeii                   -3.6739567 -2.5646103
# block                         0.1793307  0.2458996
# stimtyperectangle:ruletypeii  2.2427218  3.7754412
# ruletyperb:rulevariantheight -0.4137913  0.5648209
# ruletypeii:rulevariantneg    -0.2572502  0.6004680
# ruletyperb:rulevariantradius -1.9890884 -0.9052187
# ruletypeii:rulevariantshape   1.3410099  2.2411377


D1_DIAL_ANGLE  = c(1,0,0,0,0,0,0,0,0)
D1_DIAL_RADIUS = c(1,0,0,0,0,0,0,1,0)
D1_RECT_WIDTH  = c(1,1,0,0,0,0,0,0,0)
D1_RECT_HEIGHT = c(1,1,0,0,0,1,0,0,0)
D2_DIAL_POS    = c(1,0,1,0,0,0,0,0,0)
D2_DIAL_NEG    = c(1,0,1,0,0,0,1,0,0)
D2_RECT_SIZE   = c(1,1,1,0,1,0,0,0,0)
D2_RECT_SHAPE  = c(1,1,1,0,1,0,0,0,1)


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
           "D2_DIAL_NEG - D2_RECT_SIZE" = D2_DIAL_NEG - D2_RECT_SIZE
           )
summary(glht(m, c1), test = adjusted("bonferroni"))
confint(glht(m, c1))


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

m = aov(mn.dist ~ stimtype*ruletype + ruletype/rulevariant, data=agg)
summary(m)

### Effect size
r = anova(m)
for (p in 1:4) {
  print(r$'Sum Sq'[p]/(r$'Sum Sq'[p] + r$'Sum Sq'[5]))
}

TukeyHSD(m)$`stimtype:ruletype`

# M,SD
ddply(agg, c('stimtype', 'ruletype'), function(x) {
  c(mn=mean(x$mn.dist), sd=sd(x$mn.dist))
})


# Proportion of 2D rules ----
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




# case of quasi-complete separation -- in the 2D-RECT-shape condition,
# all 100% 2D boundaries. So need to use Firth's bias correction
# https://web.warwick.ac.uk/statsdept/user-2011/TalkSlides/Contributed/18Aug_0950_FocusVI_3-GLM-3-Kosmidis.pdf
library(brglm)

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
