rm(list=ls())

#-1: Data preparation: Directory, Package and Custom Functions

setwd("C:/Users/Dotto/Dropbox/Università/PropensityGBM/Empirical Analysis/IncomeData")

library(haven)
library(xlsx)
library(gbm)
library(labelled)
library(caret)
library(rsample)
library(questionr)
library(mgcv)
library(sjlabelled)
library(WRTDStidal)
library(quantreg)
library(pdp)
library(pROC)
library(papeR)
library(confintr)

trim.probs = function(x){
  id1    = which(x>0.99)
  x[id1] = 0.99
  id2    = which(x<0.01)
  x[id2] =0.01
  return(x)
}

#-1.1: Load Properly formatted Data

load("data.mod.rda")
data.mod$em1contrR = data.mod$em1contr
data.mod$em1contrR = as.numeric(as.character(data.mod$em1contrR))


#-2:   Statistical Models

#-2.1: Logistic Regression

mod.base = glm(formula = em1contr~age+gen+popgrp+formalS+Urban+English
               +em1tru+incrnt+inco+hldes+fwbrelinc+rel+ems+CapInc+Rem+em1occ_isco_c
               +em1prod_c+ owncom+owncel+asacc+hometype+best_marstt+best_edu+experience
               +hhsizer+married+incchld+asacc+incfos+inccare+
               +incuif+incwc+incpfnd+incret+incint+incinh+edlitcomp
               +edlitdriv+edsaid+incretr+inclob+incgif+incloan+incsale,
               data = data.mod, family = "binomial")

pred.log = round(fitted(mod.base))
n.mod    = length(mod.base$y)
err.log  = (n.mod - sum(diag(table(pred.log,mod.base$y))))/n.mod
save(mod.base,file="LogMod.rda")

summary(mod.base)
confusionMatrix(as.factor(mod.base$y),as.factor(round(mod.base$fitted.values)))

#-2.2: GBM model

load("dat_train.rda")
load("dat_test.rda")

load("grid.rda")

i = which.min(hyper_grid$err)

gbm.mod    = gbm(formula = em1contrR~age+gen+popgrp+formalS+Urban+English+
                 +em1tru+incrnt+inco+hldes+fwbrelinc+rel+ems+CapInc+Rem+em1occ_isco_c
                 +em1prod_c+ owncom+owncel+asacc+hometype+best_marstt+best_edu+experience
                 +hhsizer+married+incchld+asacc+incfos+inccare+
                 +incuif+incwc+incpfnd+incret+incint+incinh+edlitcomp
                 +edlitdriv+edsaid+incretr+inclob+incgif+incloan+incsale,
                 bag.fraction = hyper_grid$bag.fraction[i], distribution = "bernoulli",
                 n.trees = 10000, data = dat_train,interaction.depth = 20,
                 n.minobsinnode = hyper_grid$n.minobsinnode[i], cv.folds = 10,
                 shrinkage = hyper_grid$shrinkage[i], n.cores = 4)
save(gbm.mod,file = "GBM.rda")

load("GBM.rda")
best.iter = gbm.perf(gbm.mod,plot.it = T,method = "cv" )
train.id  = as.numeric(row.names(dat_train))
test.id   = as.numeric(row.names(dat_test))

pred1    = predict(gbm.mod,n.trees = best.iter,dat_train,type="response")
pred2    = predict(gbm.mod,n.trees = best.iter,dat_test, type="response")

prob.gbm = c(pred1,pred2)

dat.gbm  = rbind(dat_train,dat_test)
dat.gbm  = as.data.frame(dat.gbm)

cl1 = ifelse(pred1<0.5,0,1)
cl2 = ifelse(pred2<0.5,0,1)

cl = c (cl1,cl2)
N.tot = nrow(dat.gbm)
err.gbm = (N.tot-sum(diag(table(dat.gbm$em1contrR,cl))))/N.tot

# 2.3: GAM

mod.gam = gam(formula = em1contr~s(age)+gen+popgrp+formalS+Urban+as.factor(English)+
              +em1tru+incrnt+inco+hldes+fwbrelinc+rel+ems+s(CapInc)+s(Rem)+em1occ_isco_c
              +em1prod_c+ owncom+owncel+asacc+hometype+best_marstt+best_edu+s(experience)
              +hhsizer+married+incchld+asacc+incfos+inccare+
              +incuif+incwc+incpfnd+incret+incint+incinh+edlitcomp
              +edlitdriv+edsaid+inclob+incgif+incloan+incsale,
              data = data.mod, family=binomial(link = "logit"),method="REML")

save(mod.gam,file="mod.gam.rda")
n.gam   = length(mod.gam$y)
err.gam = (n.gam-sum(diag(table(mod.gam$y,round(mod.gam$fitted.values)))))/n.gam
err.gam 
confusionMatrix(as.factor(mod.gam$y),as.factor(round(mod.gam$fitted.values)))
plot(roc(as.factor(mod.gam$y),((mod.gam$fitted.values))))

# 3:   Quantile Regression

# 3.1: Basic

data.mod$em1hrs[which(data.mod$em1hrs==0)]=NA
data.mod$wage.hour = data.mod$fwag_p/data.mod$em1hrs
idx.vec   = which(!is.na(data.mod$wage.hour))
reg.base  = rq(wage.hour~em1contrR,
               tau = c(0.10,0.25,0.5,0.75,0.90),
              data = data.mod[idx.vec,])
reg.base$coefficients[2,]/reg.base$coefficients[1,]

# 3.2: Logistic

log.probs   = fitted(mod.base)
log.probs   = trim.probs(log.probs)
idx.log     = as.numeric(names(mod.base$fitted.values))
data.log    = data.mod[idx.log,]
log.weights = data.log$em1contrR/log.probs+(1-data.log$em1contrR)/(1-log.probs)
data.log$em1hrs[which(data.log$em1hrs==0)]=NA
data.log$wage.hour = data.log$fwag_p/data.log$em1hrs
idx.reg     = which(!is.na(data.log$wage.hour))
reg.log     = rq(wage.hour[idx.reg]~em1contrR[idx.reg],
                 tau = c(0.10,0.25,0.5,0.75,0.90),
                 data = data.log[idx.reg,],weights = log.weights[idx.reg] )
reg.log$coefficients[2,]/reg.log$coefficients[1,]

# 3.3: GAM

gam.probs   = fitted(mod.gam)
gam.probs   = trim.probs(gam.probs)
idx.gam     = as.numeric(names(mod.base$fitted.values))
data.gam    = data.mod[idx.gam,]
gam.weights = data.gam$em1contrR/gam.probs[idx.gam]+(1-data.log$em1contrR)/(1-log.probs[idx.log])
idx.reg     = which(!is.na(data.gam$wage.hour))

reg.gam     = rq(wage.hour~em1contrR,
                 tau = c(0.10,0.25,0.5,0.75,0.90),
                 data = data.log[idx.reg,],weights = gam.weights[idx.reg] )

reg.gam$coefficients[2,]/reg.gam$coefficients[1,]

# 3.4: GBM

prob.gbm    = trim.probs(prob.gbm)
Gbm.weights = (dat.gbm$em1contrR/prob.gbm)+ ((1-dat.gbm$em1contrR)/(1-prob.gbm))
dat.gbm$em1hrs[which(dat.gbm$em1hrs==0)]=NA
dat.gbm$wage.hours = dat.gbm$fwag_p/dat.gbm$em1hrs
dat.gbm$out.var = dat.gbm$inc.tot/dat.gbm$em1hrs
idx.gbm     = which(!is.na(dat.gbm$wage.hours))
#idx.gbm     = which(!is.na(dat.gbm$wage.hours))
reg.gbm     = rq(fwag_p/em1hrs~em1contrR,
                 tau = c(0.10,0.25,0.5,0.75,0.90),
                 data = dat.gbm[idx.gbm,],weights = Gbm.weights[idx.gbm] )
reg.gbm$coefficients[2,]/reg.gbm$coefficients[1,]
