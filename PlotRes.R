## Plot results ##

# -0: Directory, Packages and functions and data

setwd("C:/Users/Dotto/Dropbox/Università/PropensityGBM/Empirical Analysis/IncomeData")

load("GBM_base.rda")
load("dat_test.rda")
load("dat_train.rda")
load("mod.gam.rda")
load("data.mod.rda")
load("LogMod.rda")

#load("EmpiricalAnalysis.RData")

dat.gbm  = rbind(dat_train,dat_test)
dat.gbm  = as.data.frame(dat.gbm)

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


library(plotrix)
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
library(stargazer)
library(xtable)
library(dismo)
library(ggplot2)
library(tidygam)
library(tidymv)
library(itsadug)
library(xtable)
library(xlsx)
library(dplyr)
library(papeR)

line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}

# 1. Classification Output

# 1.1: Roc Curve
X11()
par(mfrow=c(1,2),pty = "s",mai=c(0,0,0,0))
ROC.GBM = roc(as.factor(dat.gbm$em1contrR),c(pred1,pred2),legacy.axes=T,plot=F)
plot(ROC.GBM,main = "ROC Curves",col="blue")
ROC.LOG = roc(as.factor(mod.base$y),mod.base$fitted.values,legacy.axes=T)
plot(ROC.LOG,add=T,col="green",lty=2)
ROC.GAM = roc(as.factor(mod.gam$y),((mod.gam$fitted.values)),legacy.axes=T)
plot(ROC.GAM,add=T,col="red",lty=3)
legend('bottomright',c("GBM","Logistic","GAM"),col=c("blue","green","red"),bty="n",lty=c(1,2,3))

# 1.1: Confusion Matrix

ConfMatGBM = confusionMatrix(as.factor(dat.gbm$em1contrR),as.factor(cl))
ConfMatGAM = confusionMatrix(as.factor(mod.gam$y),as.factor(round(mod.gam$fitted.values)))
ConfMatLOG = confusionMatrix(as.factor(mod.base$y),as.factor(round(mod.base$fitted.values)))

Class.mat = cbind(c(ConfMatGBM$byClass,ConfMatGBM$overall[1],c(ROC.GBM$auc)),
                  c(ConfMatGAM$byClass,ConfMatGAM$overall[1],ROC.GAM$auc),
                  c(ConfMatLOG$byClass,ConfMatLOG$overall[1],ROC.LOG$auc))
row.names(Class.mat)[nrow(Class.mat)] = "AUC"
colnames(Class.mat) = c("GBM","GAM","Logistic")
Class.mat = round(Class.mat,2)
Class.mat = data.frame(row.names(Class.mat),Class.mat)
colnames(Class.mat) = c("Metric","GBM","GAM","Logistic")
Class.mat = Class.mat[c(1,2,11,12,13),]

plot(2, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '',xlim=c(0,10),ylim=c(0,7))
title("Classification Metrics",line=-1.5,cex=2.5,cex.main=1.5)
addtable2plot(1,3,Class.mat, 
              xpad=1, ypad=1,
              bty='o',
              display.rownames = FALSE, 
              hlines = TRUE,
              vlines = TRUE, 
              title = "",
              cex = 1.5)



# 1.2 Variable Coefficients

summary.gam = summary(mod.gam)
names(summary.gam)
summary.gam$p.coeff
summary.gam$family
summary(mod.base)
summary(mod.gam)

summary.log = summary(mod.base)
print(xtable(summary.log$coefficients[,c(1,4)],digits = c(0,3,2)))
print(xtable(cbind(summary.gam$p.coeff,summary.gam$p.pv),digits = c(0,3,2)))
print(xtable(summary.gam$s.table))

summary.log =as.data.frame(summary.log$coefficients)

# 1.3 Variable Importance

var.imp = as.data.frame(summary.gbm(gbm.mod,n.trees = best.iter))
n.var   = gbm.mod$var.names

X11()
par(mar=c(5,6,3,0.4))
barplot(var.imp[order(var.imp[,2]),2],
        horiz=T,col="blue",las=1,xlab="Relative Importance",
        main="Variable Importance",cex.names = .75,
        names=var.imp$var[order(var.imp[,2])])


X11()
#par(mfrow=c(3,2),mar=c(4.5,4,4.5,1))
layout(matrix(c(1,1,2,2,3,4,5,6),byrow=T, ncol=2))
# (1): Productive sector

boxplot(prob.gbm~dat.gbm$em1prod_c,outline=F,main="Productive sector",
        col="green",xlab="",ylab="P(T=1)",xaxt = "n")
names=c("Private \n households","Agriculture \n Fishing","\n Mining"," Manufacturing",
        "Electricity \n supply","Construction","Trade and \n Restaurants","Transport",
        "Financial \n intermediation","Community and \n social services")

axis(side=1,at=c(1:length(names)),labels=names,tick=F,cex.axis=.7)

# (3): Occupation

boxplot(prob.gbm~dat.gbm$em1occ_isco_c,outline=F,xlab="",
        ylab="P(T=1)",col="green",xaxt="n",main="Occupation")
names = c("Armed forces","Managers","Professionals","Technicians",
          "Clerical \n support","Service and \n sales workers",
          "Agricultural \n and fishery workers","Craft",
          "Plant and \n machine operators","Elementary \n occupations")
axis(side=1,at=c(1:length(names)),labels=names,tick=F,cex.axis=.8)



# (2): Bank Account

boxplot(prob.gbm~dat.gbm$asacc,outline=F,main="Bank Account",
        col="green",xlab="",ylab="P(T=1)",xaxt="n")
axis(side=1,at=c(1,2),labels=c("Yes","No"),tick=F,cex.axis=.8)


# (4) Trade Union Membership

boxplot(prob.gbm~dat.gbm$em1tru,outline=F,col="green",
        main ="Trade union memberhsip",names = c("Yes","No"),
        xlab="",ylab="P(T=1)",xaxt="n")
axis(side=1,at=c(1,2),labels=c("Yes","No"),tick=F,cex.axis=.8)

# (5): Experience

plotdata <- plot(gbm.mod, return.grid=TRUE,
                 i.var=which(n.var==var.imp[5,1]),
                 type="response",n.trees=best.iter)
plot(plotdata[c(1:100),1],plotdata[c(1:100),2], 
     type="l",col="green",
     main="Experience",
     ylab="P(T=1)",xlab="",)


# (6): Capital Income

plotdata <- plot(gbm.mod, return.grid=TRUE,
                 i.var=which(n.var==var.imp[6,1]),
                 type="response",n.trees=best.iter)
plot(plotdata[c(1:100),1],plotdata[c(1:100),2], 
     type="l",col="green",
     main="Capital Income",
     ylab="P(T=1)",xlab="",)

X11()
par(mfrow=c(2,2))


# (7): Capital Income

plotdata <- plot(gbm.mod, return.grid=TRUE,
                 i.var=which(n.var==var.imp[7,1]),
                 type="response",n.trees=best.iter)
plot(plotdata[c(1:100),1],plotdata[c(1:100),2], 
     type="l",col="green",
     main="Age",
     ylab="P(T=1)",xlab="",)


# (8): Perceived Income

boxplot(prob.gbm~dat.gbm$fwbrelinc,outline=F,col="green",
        main ="Perceived Income",
        xlab="",ylab="P(T=1)",xaxt="n")
names = c("Much above \n average"," Above \n average",
          "Average","Below \n average","Much below \n average")
axis(side=1,at=c(1:5),labels=names,tick=F,cex.axis=.8)

# (9): Self perceived income

plotdata <- plot(gbm.mod, return.grid=TRUE,
                 i.var=which(n.var==var.imp[9,1]),
                 type="response",n.trees=best.iter)
plot(plotdata[c(1:100),1],plotdata[c(1:100),2], 
     type="l",col="green",
     main="Remittances",
     ylab="P(T=1)",xlab="",)

# (10): Education
idx.vec = which(dat.gbm$best_edu==0|dat.gbm$best_edu==1|dat.gbm$best_edu==2|dat.gbm$best_edu==3)
dat.gbm1 = dat.gbm[idx.vec,]
dat.gbm1$best_edu = as.character(dat.gbm1$best_edu)
boxplot(prob.gbm[idx.vec]~dat.gbm1$best_edu,
        outline=F,col="green", main ="Education",
        xlab="",ylab="P(T=1)",xaxt="n")
names = c("No ducation"," Low education",
          "Medium education","High education")
axis(side=1,at=c(1:4),labels=names,tick=F,cex.axis=.8)

# (10): Population Group

X11()
idx.pop = which(dat.gbm$popgrp ==1|dat.gbm$popgrp ==2|dat.gbm$popgrp ==3|dat.gbm$popgrp ==4)
dat.gbm$popgrpN = as.character(dat.gbm$popgrp)
boxplot(prob.gbm[idx.pop]~dat.gbm$popgrpN[idx.pop],
        outline=F,col="green",
        main ="Population Group",
        xlab="",ylab="P(T=1)",xaxt="n")
names = c("African"," Coloured",
          "Asian/Indian","White")
axis(side=1,at=c(1:4),labels=names,tick=F,cex.axis=.8)


# Quantile regression

base.summary = summary(reg.base,se="boot")
gam.summary  = summary(reg.gam,se="boot")
gbm.summary  = summary(reg.gbm,se="boot")
log.summary  = summary(reg.log,se="boot")

save.image("ExPostAnalysis.RData")

tab.res = NULL

for(j in 1:length(base.summary)){
  
  v1 =  round(base.summary[[j]]$coefficients[2,1]/base.summary[[j]]$coefficients[1,1],2)
  v2 =  round(log.summary[[j]]$coefficients[2,1]/log.summary[[j]]$coefficients[1,1],2)
  v3 =  round(gam.summary[[j]]$coefficients[2,1]/gam.summary[[j]]$coefficients[1,1],2)
  v4 =  round(gbm.summary[[j]]$coefficients[2,1]/gbm.summary[[j]]$coefficients[1,1],2)
  
  v1.vec = c(v1,v2,v3,v4)
  
  v1.s = paste0("(",round(base.summary[[j]]$coefficients[2,2],2),")")
  v2.s = paste0("(",round(log.summary[[j]]$coefficients[2,2],2),")")
  v3.s = paste0("(",round(gam.summary[[j]]$coefficients[2,2],2),")")
  v4.s = paste0("(",round(gbm.summary[[j]]$coefficients[2,2],2),")")
  
  res.v =rbind(rbind(v1,v1.s),rbind(v2,v2.s),rbind(v3,v3.s),rbind(v4,v4.s))
  
  tab.res = cbind(tab.res,res.v)
  
}

xtable((tab.res))

# (4): Descriptive Statistics of the variables
ind.vec = NULL
for(j in 1:length(n.var)){
  ind = which(colnames(data.mod)==n.var[j])
  ind.vec =c(ind,ind.vec)
}
ind.vec = c(ind.vec,which(colnames(data.mod)=="em1hrs"))

data.mod$ems = as.factor(data.mod$ems)
data.mod$hometype = as.factor(data.mod$hometype)
data.mod$English = as.factor(data.mod$English)
data.mod$incretr = as.factor(data.mod$incretr)
str(data.mod)

xtable(summarize(data.mod[,c(ind.vec)],type = "numeric"))
xtable(summarize(data.mod[,c(ind.vec)],type = "factor"))


names.vec = c("q10","","q25","","q50","","q75","","q90","")
tab.res = cbind(names.vec,tab.res)
colnames(tab.res) = c("Quantile","Empirical","Logistic","GAM","GBM")
print(xtable(tab.res),include.rownames=F)

setwd("C:/Users/Dotto/Dropbox/Università/PropensityGBM/Tex")
tab.var = read.xlsx("TabVar.xlsx",sheetIndex = 1)
print(xtable(tab.var),include.rownames=F)
