########################
#                      #
# GBM Propensity score #
#                      #
########################

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

table.rel = function(x){
  val = table(x)/sum(table(x))
  return(val)
}

trim.probs = function(x){
  id1    = which(x>0.99)
  x[id1] = 0.99
  id2    = which(x<0.01)
  x[id2] =0.01
  return(x)
}


#-1.1: Import Derived Data

hhdata             = read_dta("hhdata.dta")


#-1.3: Select relevant variables 

data.sel = hhdata[,c(1,2,c(17:21),30,54,55,56,57,239,266,270,272,273,
                     276,293,317,318,320,336,347,483,499,865,866,867,957,
                     966,1002,1003,1084,1085,1250,1271,1272,1273,1277,1286,
                     1342,1347,1356,1371,1382,1412,1589,1576,1577,1578,1579,
                     1580,1581,1582,1583,1584,1594,1601,1605,463,461,
                     1597,1598,1599,1600,1601,1602,465,467,471,473,475,477,
                     479,481,485,487,489,860,861,862,487,491,493,495,497,348,
                     277,349,284,285,326,327,335,1350)]

#-1.3.1: Assign NA 

na.ind = which(data.sel<0,arr.ind = T)

for(i in 1:nrow(na.ind)){
  
  r.ind = na.ind[i,1]
  c.ind = na.ind[i,2]
  data.sel[r.ind,c.ind] = NA
  
}

#-1.4: Check the variables and missing values

data.sel$age = 2018-data.sel$dob_y

#-1.4.1: Subset for ages

idx = which(data.sel$age>=15&data.sel$age<=65)
data.sel = data.sel[idx,]

idx.inf  = is.na(data.sel$em1contr)
data.sel = data.sel[!idx.inf,]

inc.mat = (na.omit(cbind(data.sel$fwag_p,data.sel$em1inc)))
plot(inc.mat[,1],inc.mat[,2])


summary(data.sel$fwag)
summary(data.sel$em1inc)

#-1.4.2: Create Income Variable

inc.mat = cbind(data.sel$fwag_p,data.sel$fwag_p2,data.sel$cwag_p,data.sel$swag)
inc.tot = apply(inc.mat,1,sum,na.rm=T)

data.sel$inc.tot = inc.tot
rm(inc.mat)

idx.inc  = which(data.sel$inc.tot!=0)
data.sel = data.sel[idx.inc,]
table(is.na(data.sel$inc.tot),is.na(data.sel$monthlyhours))
data.sel$wage = data.sel$inc.tot/data.sel$monthlyhours

#-1.5: Recode variables

#-1.5.1: Education: 3 levels

data.sel$best_edu[data.sel$best_edu >= 1 & data.sel$best_edu<=9]<-1
data.sel$best_edu[data.sel$best_edu >= 10 & data.sel$best_edu<=12]<-2
data.sel$best_edu[data.sel$best_edu >= 13 & data.sel$best_edu<=19]<-2
data.sel$best_edu[data.sel$best_edu >= 20 & data.sel$best_edu<=23]<-3
data.sel$best_edu[data.sel$best_edu ==25 | data.sel$best_edu==0]<-0
data.sel$best_edu[data.sel$best_edu ==24]<-NA
data.sel$best_edu = as.factor(data.sel$best_edu)

#-1.5.2: Experience

data.sel$em1strty[which(data.sel$em1strty==3333|data.sel$em1strty==5555|
                          data.sel$em1strty==9999|data.sel$em1strty==8888)]=NA
data.sel$experience = 2018-data.sel$em1strty

#-1.5.3: English Skill 1: Bad; 2:Good

table.rel(data.sel$edlitrden)
data.sel$edlitrden[data.sel$edlitrden==3|data.sel$edlitrden==4] = 1
data.sel$edlitrden[data.sel$edlitrden==1|data.sel$edlitrden==2] = 2
data.sel$edlitrden = as.factor(data.sel$edlitrden)

data.sel$edlitwrten[data.sel$edlitwrten==3|data.sel$edlitwrten==4] = 1
data.sel$edlitwrten[data.sel$edlitwrten==1|data.sel$edlitwrten==2] = 2
data.sel$edlitwrten = as.factor(data.sel$edlitwrten)

# 1.5.4: English at Home

data.sel$English = NULL
data.sel$English = ifelse(data.sel$lng==11,1,0)

#-1.5.5: Father, Mother  and primary Occupation recoded 0-1

data.sel$fthwrk_isco_c[data.sel$fthwrk_isco_c<0] = NA
data.sel$mthwrk_isco_c[data.sel$mthwrk_isco_c<0] = NA
data.sel$FISCO = NULL
data.sel$FISCO = ifelse(data.sel$fthwrk_isco_c<=4,1,0)
data.sel$MISCO = ifelse(data.sel$mthwrk_isco_c<=4,1,0)
data.sel$PISCO = data.sel$em1occ_isco_c 
data.sel$PISCO = ifelse(data.sel$em1occ_isco_c<=4,1,0)
 
# https://en.wikipedia.org/wiki/International_Standard_Classification_of_Occupations

#-1.5.7: Formal Job Search

data.sel$formalS = data.sel$em1inf
data.sel$formalS[data.sel$em1inf==1|data.sel$em1inf==2] = 1
data.sel$formalS[data.sel$em1inf==5|data.sel$em1inf<=6|data.sel$em1inf==7] = 1
data.sel$formalS[data.sel$em1inf==10] = 1

data.sel$formalS[data.sel$em1inf==3|data.sel$em1inf==4] = 2
data.sel$formalS[data.sel$em1inf==8|data.sel$em1inf==9] = 2

#-1.5.4: Urban
#da dicotomizzare

data.sel$Urban = data.sel$geo2011
data.sel$Urban[which(data.sel$geo2011==1|data.sel$geo2011==3)]=1
data.sel$Urban[which(data.sel$geo2011==2)]=2

#-1.5.6: Per capita Income

data.sel$PerInc = data.sel$hhincome/data.sel$hhsizer

#-1.5.6: Per capita Income

data.sel$em1contr_d = as.factor(data.sel$em1contr_d)
data.sel$em1occ_isco_c = as.factor(data.sel$em1occ_isco_c)
#-1.6: Create Final CodeBook

my.CodeBook = NULL

# for(j in 1:length(colnames(data.sel))){
#   vv = c(colnames(data.sel)[j], attr(data.sel[[j]],"label"))
#   my.CodeBook = rbind(my.CodeBook,vv)
# }
# row.names(my.CodeBook) = c(1:nrow(my.CodeBook))
# write.xlsx(my.CodeBook,file="CodehhDataSel.xlsx")

# 2:   Modelling the target variable formal

# 2.1: Select the variable for modelling the probability of informal

data.mod = data.sel

data.mod.tibble = data.mod


my.CodeBook.mod = NULL

for(j in 1:length(colnames(data.mod.tibble))){
  vv = c(colnames(data.mod.tibble)[j], attr(data.mod.tibble[[j]],"label"))
  my.CodeBook.mod = rbind(my.CodeBook.mod,vv)
}

row.names(my.CodeBook.mod) = c(1:nrow(my.CodeBook.mod))
write.xlsx(my.CodeBook.mod,file="CodeDataMod.xlsx")

summary(data.mod)

data.mod = as.data.frame(data.mod)
idx.sel  = is.na(data.mod$em1contr)

data.mod = data.mod[!idx.sel,]

data.mod$gen    = as.factor(data.mod$gen)
data.mod$popgrp = as.factor(data.mod$popgrp)
data.mod$PISCO  = data.mod$PISCO
data.mod$bhlive = as.factor(data.mod$bhlive)
data.mod$bhali  = as.factor(data.mod$bhali)
data.mod$FISCO  = as.factor(data.mod$FISCO)
data.mod$MISCO  = as.factor(data.mod$MISCO)
data.mod$formalS= as.factor(data.mod$formalS)
data.mod$Urban  = as.factor(data.mod$Urban)
data.mod$em1    = as.factor(data.mod$em1)

data.mod$bhlive_n = as.numeric(as.character(data.mod$bhlive_n))
data.mod$bhlive_n[which(data.mod$bhlive==2)] =0

data.mod$bhali = as.factor(data.mod$bhali)
data.mod$bhali_n[which(data.mod$bhali==2)] = 0
data.mod$em1contr = as.numeric(as.character(data.mod$em1contr))-1
data.mod$em1 = as.factor(data.mod$em1)
data.mod$em1contr = as.factor(data.mod$em1contr)
data.mod$em1tru   = as.factor(data.mod$em1tru)
data.mod$incrnt   = as.factor(data.mod$incrnt)
data.mod$inco     = as.factor(data.mod$inco)
data.mod$hldes    = as.factor(data.mod$hldes)
data.mod$fwbrelinc= as.factor(data.mod$fwbrelinc) 
data.mod$rel      = as.factor(data.mod$rel)
data.mod$owncom   = as.factor(data.mod$owncom)
data.mod$owncel   = as.factor(data.mod$owncel)
data.mod$asacc    = as.factor(data.mod$asacc)
data.mod$best_marstt = as.factor(data.mod$best_marstt)
data.mod$married  = as.factor(data.mod$married)
data.mod$formalS  = as.factor(data.mod$formalS)
data.mod$Urban    = as.factor(data.mod$Urban)
data.mod$incchld  = as.factor(data.mod$incchld)
data.mod$OINC     = rowSums(cbind(data.mod$hhgovt,data.mod$hhother,
                                  data.mod$hhinvest,data.mod$hhcapital,
                                  data.mod$hhremitt,data.mod$hhimprent),
                            na.rm = T)

data.mod$CapInc   = rowSums(cbind(data.mod$hhinvest,data.mod$hhcapital,
                                  data.mod$hhimprent),na.rm=T)

data.mod$Rem     =  rowSums(cbind(data.mod$hhremitt,data.mod$hhgovt,
                                  data.mod$hhother,data.mod$pi_hhagric),
                            na.rm=T)

data.mod$incfos   = as.factor(data.mod$incfos)
data.mod$inccare  = as.factor(data.mod$inccare)
data.mod$incwar   = as.factor(data.mod$incwar)
data.mod$incuif   = as.factor(data.mod$incuif)
data.mod$incwc    = as.factor(data.mod$incwc)
data.mod$incpfnd  = as.factor(data.mod$incpfnd)
data.mod$incret   = as.factor(data.mod$incret)
data.mod$incretp  = as.factor(data.mod$incretp)
data.mod$incint   = as.factor(data.mod$incint)
data.mod$incretr  = as.factor(data.mod$incretr) 
data.mod$incinh   = as.factor(data.mod$incinh)
data.mod$edlitcomp= as.factor(data.mod$edlitcomp)
data.mod$edlitdriv= as.factor(data.mod$edlitdriv)
data.mod$edsaid   = as.factor(data.mod$edsaid)  
data.mod$incretr  = as.factor(data.mod$incretr)
data.mod$inclob   = as.factor(data.mod$inclob)
data.mod$incgif   = as.factor(data.mod$incgif)
data.mod$incloan  = as.factor(data.mod$incloan)
data.mod$incsale  = as.factor(data.mod$incsale)
data.mod$em1prod_c = as.factor(data.mod$em1prod_c)

save(data.mod,file="data.mod.rda")

num.vec = c()

for(j in 1:ncol(data.mod)){
  num.vec[j] = is.numeric(data.mod[,j])
}
names(num.vec) = colnames(data.mod)

ll.tab = c()

for(j in 1:ncol(data.mod)){
  ll.tab[j] = length(table(data.mod[j]))
}

names(ll.tab) = colnames(data.mod)
names(ll.tab[which(ll.tab==1)])


na.lenght = c()

for(j in 1:ncol(data.mod)){
  
  check.vec = is.na(data.mod[,j])
  na.lenght[j] = length(which(check.vec==T))
  
}

names(na.lenght) = colnames(data.mod)

# https://uc-r.github.io/gbm_regression

# 2.2: Logit Model

mod.base = glm(formula = em1contr~age+gen+popgrp+formalS+Urban+English
               +em1tru+incrnt+inco+hldes+fwbrelinc+rel+ems+CapInc+Rem+em1occ_isco_c
               +em1prod_c+ owncom+owncel+asacc+hometype+best_marstt+best_edu+experience
               +hhsizer+married+monthlyhours+incchld+asacc+incfos+inccare+
               +incuif+incwc+incpfnd+incret+incint+incinh+edlitcomp
               +edlitdriv+edsaid+incretr+inclob+incgif+incloan+incsale,
                         data = data.mod, family = "binomial")

pred.log = round(fitted(mod.base))
n.mod    = length(mod.base$y)
err.log  = (n.mod - sum(diag(table(pred.log,mod.base$y))))/n.mod

summary(mod.base)
confusionMatrix(as.factor(mod.base$y),as.factor(round(mod.base$fitted.values)))

# 2.3.1: GBM Tuning NOT RUN

data.mod$em1contrR = data.mod$em1contr
data.mod$em1contrR = as.numeric(as.character(data.mod$em1contrR))

set.seed(123)

dat_split  = initial_split(data.mod, prop = .8)
dat_train  = training(dat_split)
dat_test   = testing(dat_split)

hyper_grid <- expand.grid(
  shrinkage = c(.01, .1,.001),
  n.minobsinnode = c(5,  15),
  bag.fraction = c(.75,.8, 1), 
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

for(i in 1:nrow(hyper_grid)) {
  
  set.seed(123)

  gbm.mod    = gbm(formula = em1contrR~age+gen+popgrp+formalS+Urban+English
                   +em1tru+incrnt+inco+hldes+fwbrelinc+rel+ems+CapInc+Rem+em1occ_isco_c
                   +hhsizer+married+monthlyhours+incchld+asacc+incfos+inccare+
                   +incuif+incwc+incpfnd+incret+incnt+incinh+edlitcomp
                   +edlitdriv+edsaid+incretr+inclob+incgif+incloan+incsale,
                   bag.fraction = hyper_grid$bag.fraction[i], distribution = "bernoulli", 
                   n.trees = 10000, data = dat_train,interaction.depth = 20, 
                   n.minobsinnode = hyper_grid$n.minobsinnode[i], cv.folds = 10,
                   shrinkage = hyper_grid$shrinkage[i], n.cores = 4)
  
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
  hyper_grid$err[i] =  err.gbm
  hyper_grid$best[i] = best.iter
  print(i)

  }

load("grid.rda")

i = which.min(hyper_grid$err)

gbm.mod    = gbm(formula = em1contrR~age+gen+popgrp+formalS+Urban+English+em1hrs
                 +em1tru+incrnt+inco+hldes+fwbrelinc+rel+ems+CapInc+Rem+em1occ_isco_c
                 +em1prod_c+ owncom+owncel+asacc+hometype+best_marstt+best_edu+experience
                 +hhsizer+married+monthlyhours+incchld+asacc+incfos+inccare+
                 +incuif+incwc+incpfnd+incret+incint+incinh+edlitcomp
                 +edlitdriv+edsaid+incretr+inclob+incgif+incloan+incsale,
                 bag.fraction = hyper_grid$bag.fraction[i], distribution = "bernoulli", 
                 n.trees = 10000, data = dat_train,interaction.depth = 20, 
                 n.minobsinnode = hyper_grid$n.minobsinnode[i], cv.folds = 10,
                 shrinkage = hyper_grid$shrinkage[i], n.cores = 4)

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



save(gbm.mod,file="GBM_base.rda")
save(dat_train,file="dat_train.rda")
save(dat_test,file="dat_test.rda")

plot(gbm.mod,n.trees=best.iter,i.var = 27,type="response")
confusionMatrix(as.factor(dat.gbm$em1contrR),as.factor(cl))
plot(roc(as.factor(dat.gbm$em1contrR),cl,legacy.axes=T))
roc((dat.gbm$em1contrR),cl,legacy.axes=T,plot=T,type="b")



summary(gbm.mod)
plot(gbm.mod,n.trees=best.iter,i.var = 28,type="response")
confusionMatrix(as.factor(dat.gbm$em1contrR),as.factor(cl))
roc(as.factor(dat.gbm$em1contrR),c(pred1,pred2),legacy.axes=T,plot=T)
roc(as.factor(mod.base$y),mod.base$fitted.values,legacy.axes=T,plot=T,add=T)
roc(as.factor(mod.gam$y),((mod.gam$fitted.values)),add=T,legacy.axes=T,plot=T,add=T)


# 2.4: GAM

mod.gam = gam(formula = em1contr~s(age)+gen+popgrp+formalS+Urban+English+s(em1hrs)
              +em1tru+incrnt+inco+hldes+fwbrelinc+rel+ems+s(CapInc)+s(Rem)+em1occ_isco_c
              +em1prod_c+ owncom+owncel+asacc+hometype+best_marstt+best_edu+s(experience)
              +hhsizer+married+s(monthlyhours)+incchld+asacc+incfos+inccare+
              +incuif+incwc+incpfnd+incret+incint+incinh+edlitcomp
              +edlitdriv+edsaid+incretr+inclob+incgif+incloan+incsale,
              data = data.mod, family=binomial(link = "logit"),method = "REML")
n.gam   = length(mod.gam$y)
err.gam = (n.gam-sum(diag(table(mod.gam$y,round(mod.gam$fitted.values)))))/n.gam
err.gam 
save(mod.gam,file="mod.gam.rda")
confusionMatrix(as.factor(mod.gam$y),as.factor(round(mod.gam$fitted.values)))
plot(roc(as.factor(mod.gam$y),((mod.gam$fitted.values))))

# 3:   Quantile Regression

# 3.1: Basic
data.mod$em1hrs[which(data.mod$em1hrs==0)]=NA
data.mod$
reg.base = rq(inc.tot~em1contrR,tau = c(0.10,0.25,0.5,0.75,0.90),
              data = data.mod)

# 3.2: Logistic

log.probs   = fitted(mod.base)
log.probs   = trim.probs(log.probs)
idx = as.numeric(names(mod.base$fitted.values))
log.weights = data.mod$em1contrR[idx]/log.probs[idx]+(1-data.mod$em1contrR[idx])/(1-log.probs[idx])
reg.log  = rq(inc.tot~em1contrR,tau = c(0.10,0.25,0.5,0.75,0.90),
              data = data.mod[idx,],weights = log.weights )

reg.log  = rq(inc.tot~em1contrR,tau = c(0.10,0.25,0.5,0.75,0.90),
              data = data.mod[idx,],weights = log.weights )

# 3.3: GAM

log.gam   = fitted(mod.gam)
log.gam   = trim.probs(log.gam)
idx       = as.numeric(names(mod.base$fitted.values))
gam.weights = data.mod$em1contrR[idx]/log.gam[idx]+(1-data.mod$em1contrR[idx])/(1-log.gam[idx])

reg.gam  = rq(inc.tot~em1contrR,tau = c(0.10,0.25,0.5,0.75,0.90),
              data = data.mod[idx,],weights = gam.weights )

reg.gam  = rq(wage~em1contrR,tau = c(0.10,0.25,0.5,0.75,0.90),
              data = data.mod[idx,],weights = gam.weights )


# 3.4: GBM
prob.gbm = trim.probs(prob.gbm)
Gbm.weights = (dat.gbm$em1contrR/prob.gbm)+ ((1-dat.gbm$em1contrR)/(1-prob.gbm))

reg.gbm = rq(wage~em1contrR,tau = c(0.10,0.25,0.5,0.75,0.90),
             data = dat.gbm,weights = Gbm.weights )

reg.gbm1 = rq(inc.tot~em1contrR,tau = c(0.10,0.25,0.5,0.6,0.7,0.75,0.8,0.85,0.90,0.95,0.99),
             data = dat.gbm,weights = Gbm.weights )

dat.gbm$em1hrs[which(dat.gbm$em1hrs==0)]=NA
idx.reg = which(!is.na(dat.gbm$em1hrs))
out.var = dat.gbm$inc.tot/dat.gbm$em1hrs
summary(out.var[idx.reg])
reg.gbm1 = rq(out.var[idx.reg]~em1contrR[idx.reg],
              tau = c(0.10,0.25,0.5,0.6,0.7,0.75,0.8,0.85,0.90,0.95,0.99),
              data = dat.gbm,weights = Gbm.weights[idx.reg] )
reg.gbm1$coefficients[2,]/reg.gbm1$coefficients[1,]

reg.gbm2 = rq(log(out.var[idx.reg])~em1contrR[idx.reg],
              tau = c(0.10,0.25,0.5,0.6,0.7,0.75,0.8,0.85,0.90,0.95,0.99),
              data = dat.gbm,weights = Gbm.weights[idx.reg] )
reg.gbm2$coefficients[2,]/reg.gbm2$coefficients[1,]




