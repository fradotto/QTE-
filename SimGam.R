
library(gbm)
library(mgcv)
library(quantreg)
library(mclust)
library(caret)
library(rsample)
library(questionr)
library(qte)
library(dplyr)
library(Sieve)
library(plotrix)


trim.probs = function(x){
  id1    = which(x>0.99)
  x[id1] = 0.99
  id2    = which(x<0.01)
  x[id2] =0.01
  return(x)
}

# Jittered Linear

n.sim = 500

res.mat = NULL

for(l in 1:n.sim){
  
  set.seed(l)
  
  # Data Gen
  
  N.tot  = 5000
  int    = -.2
  X1     = rnorm(N.tot,3,)
  X2     = rnorm(N.tot,4,2)
  X3     = rnorm(N.tot,10,5)
  X4     = rbinom(n=N.tot,prob=.7,size=1)+1
  X5     = runif(N.tot,-3,3)
  X6     = ifelse(X3>5,0,1)
  X7     = sample(x=c(0,1,2,3,4),size = N.tot,prob = c(0.1,0.25,0.2,0.05,0.4),replace = T)
  X8     = rbinom(n=N.tot,prob=.3,size=1)
  X9     = rbinom(n=N.tot,prob=.5,size=1)
  X10    = runif(N.tot,0,1)
  X11    = runif(N.tot,-1,1)
  X12    = runif(N.tot,0,1)
  
  
  probs.vec = 1/(1+exp(-(int-log(abs(X10))*X4*X11+0.4*X6
                         +(I(X7==1)*0.4+I(X7==2)*0.3+I(X7==3)*(-0.6)+I(X7==4)*0.5)*X5
                         -.5*exp(-(X5))+0.9*X12)))
  probs.vec[which(probs.vec<0.01)]=0.01
  Y0 = c(); Y1=c(); Y=c()
  mu0 = 10; mu1=25; sigma1=8; sigma0=3
  
  t.var = c()
  
  for(i in 1:length(X1)){
    t.var[i] = rbinom(n=1,size=1,prob = probs.vec[i])
  }
  
  for(i in 1:length(X1)){
    
    Y0[i] = X1[i]+X2[i]
    Y1[i] = X1[i]+X3[i]
    Y[i]  = (t.var[i]*Y1[i])+(1-t.var[i])*Y0[i]
    
     }
  
  unf.qte = quantile(Y1,probs =c(0.10,0.25,0.50,0.75,0.90))-quantile(Y0,probs =c(0.10,0.25,0.50,0.75,0.90))
  emp.qte = quantile(Y[which(t.var==1)],probs =c(0.10,0.25,0.50,0.75,0.90))-quantile(Y[which(t.var==0)],probs =c(0.10,0.25,0.50,0.75,0.90))
  true.qte= qnorm(c(0.10,0.25,0.50,0.75,0.90),13,sqrt(26))-qnorm(c(0.10,0.25,0.50,0.75,0.90),7,sqrt(5))
  
  dat = as.data.frame(cbind(Y,X1,X2,X3,X4,X5,X6,X7,
                            X8,X9,X10,X11,X12,t.var)) 
  
  dat$X4 = as.factor(dat$X4)
  dat$X6 = as.factor(dat$X6)
  dat$X7 = as.factor(dat$X7)
  dat$X8 = as.factor(dat$X8)
  dat$X9 = as.factor(dat$X9)
  
  for(i in 1:nrow(dat)){
    for(j in c(2:5)){
      mis = sample(x=c(0,1),size=1,prob=c(1,0))
      if(mis==1){
        dat[i,j] = NA
      }
    }
  }
  
  true.weights = dat$t.var/probs.vec+ ((1-dat$t.var)/(1-probs.vec))
  
  RegOracle    = rq(Y~t.var,weights = true.weights,
                    tau = c(0.10,0.25,0.5,0.75,0.90),
                    data=dat)
  Oracle.Sum   = try(summary(RegOracle))
  
  cov.vec.Or   = c()
  
  if(!inherits(Oracle.Sum,"try-error")){  
    
    for(j in 1:length(Oracle.Sum)){
      int.cov      = c(Oracle.Sum[[j]]$coefficients[2,1]-1.96*Oracle.Sum[[j]]$coefficients[2,2],
                       Oracle.Sum[[j]]$coefficients[2,1]+1.96*Oracle.Sum[[j]]$coefficients[2,2])
      cov.vec.Or[j]= between(true.qte[j],int.cov[1],int.cov[2])
      
    }
  }else{cov.vec.Or=c(rep(NA,5))}
  
  na.mat   = is.na(dat)
  na.check = apply(na.mat,1,sum)
  idx = which(na.check==0)
  
  log.reg      = glm(t.var~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12,
                     data=dat, family=binomial(link = "logit"))
  log.probs    = fitted(log.reg)
  log.probs    = trim.probs(log.probs)
  log.weights  = dat$t.var[idx]/log.probs+ ((1-dat$t.var[idx])/(1-log.probs))  
  log.pred     = ifelse(log.probs<0.5,0,1)
  
  Reg.Log      = rq(Y~t.var, weights = log.weights, 
                    tau = c(0.10,0.25,0.5,0.75,0.90), 
                    data = dat[idx,])
  
  Log.Sum      = try(summary(Reg.Log))
  cov.vec.Log  = c()
  
  if(!inherits(Oracle.Sum,"try-error")){
    
    for(j in 1:length(Log.Sum)){
      int.cov       = c(Log.Sum[[j]]$coefficients[2,1]-1.96*Log.Sum[[j]]$coefficients[2,2],
                        Log.Sum[[j]]$coefficients[2,1]+1.96*Log.Sum[[j]]$coefficients[2,2])
      cov.vec.Log[j]= between(true.qte[j],int.cov[1],int.cov[2])
    }
    
  }else{cov.vec.Log=c(rep(NA,5))}
  
  log.gam      = gam(t.var~s(X1)+s(X2)+s(X3)+X4+s(X5)+X6+X7
                     +X8+X9+s(X10)+s(X11)+s(X12),
                     data=dat,family=binomial(link = "logit"))
  gam.pobs     = fitted(log.gam)
  gam.pobs     = trim.probs(gam.pobs)
  gam.weights  = dat$t.var[idx]/gam.pobs+ ((1-dat$t.var[idx])/(1-gam.pobs))
  
  Reg.Gam      = rq(Y~t.var, weights = gam.weights, 
                    tau = c(0.10,0.25,0.5,0.75,0.90),
                    data = dat[idx,])
  Gam.Sum      = try(summary(Reg.Gam))
  
  cov.vec.Gam = c()
  
  if(!inherits(Gam.Sum,"try-error")){
    
    for(j in 1:length(Gam.Sum)){
      int.cov       = c(Gam.Sum[[j]]$coefficients[2,1]-1.96*Gam.Sum[[j]]$coefficients[2,2],
                        Gam.Sum[[j]]$coefficients[2,1]+1.96*Gam.Sum[[j]]$coefficients[2,2])
      cov.vec.Gam[j]= between(true.qte[j],int.cov[1],int.cov[2])
    }
    
  }else{cov.vec.Gam=c(rep(NA,5))}
  
  # GBM
  
  dat_split  = initial_split(dat, prop = .8)
  dat_train  = training(dat_split)
  dat_test   = testing(dat_split)
  
  # Launch Model
  
  gbm.fit   = gbm(formula = t.var~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12, 
                  distribution = "bernoulli", data = dat_train, n.trees = 6000, 
                  interaction.depth = 6,n.minobsinnode = 10,
                  shrinkage = .01,cv.folds = 5, n.cores = NULL,  verbose = F)  
  best.iter = gbm.perf(gbm.fit,plot.it = T,method = "cv" )
  train.id  = as.numeric(row.names(dat_train))
  test.id   = as.numeric(row.names(dat_test))
  
  pred1    = predict(gbm.fit,n.trees = best.iter,dat_train,type="response")
  pred2    = predict(gbm.fit,n.trees = best.iter,dat_test, type="response")
  
  prob.gbm = c(pred1,pred2)
  prob.gbm = trim.probs(prob.gbm)
  dat.gbm  = rbind(dat_train,dat_test)
  dat.gbm  = as.data.frame(dat.gbm)
  colnames(dat.gbm) = c("Y","X1","X2","X3","X4","X5","X6",
                        "X7","X8","X9","X10","X11","X12","t.var")
  dat.gbm  = as.data.frame(dat.gbm) 
  
  cl1 = ifelse(pred1<0.5,0,1)
  cl2 = ifelse(pred2<0.5,0,1)
  
  cl = c (cl1,cl2)
  err.gbm = (N.tot-sum(diag(table(dat.gbm$t.var,cl))))/N.tot
  
  # GBM Regression
  
  Gbm.weights = (dat.gbm$t.var/prob.gbm)+ ((1-dat.gbm$t.var)/(1-prob.gbm))
  GbmMod      = rq(Y~t.var,weights = Gbm.weights,
                   tau = c(0.10,0.25,0.5,0.75,0.90),
                   data=dat.gbm)
  
  Gbm.sum      = try(summary(GbmMod))
  
  cov.vec.gbm = c()
  
  if(!inherits(Gbm.sum,"try-error")){
    
    for(j in 1:length(Gbm.sum)){
      int.cov       = c(Gbm.sum[[j]]$coefficients[2,1]-1.96*Gbm.sum[[j]]$coefficients[2,2],
                        Gbm.sum[[j]]$coefficients[2,1]+1.96*Gbm.sum[[j]]$coefficients[2,2])
      cov.vec.gbm[j]= between(true.qte[j],int.cov[1],int.cov[2])
    }
    
  }else{cov.vec.Gbm=c(rep(NA,5))}
  
  res =  rbind(true.qte,unf.qte,emp.qte, 
               c(RegOracle$coefficients[2,],cov.vec.Or),
               c(Reg.Log$coefficients[2,],cov.vec.Log),
               c(Reg.Gam$coefficients[2,],cov.vec.Gam),
               c(GbmMod$coefficients[2,],cov.vec.gbm))
  
  row.names(res) = c("True","Unfeasible","Empirical","Oracle","Log","Gam","Gbm")
  
  res.mat = rbind(res.mat,res)
  
  print(l)
  
  if(l %% 100==0){
    save(res.mat,file="ResMatJLin.rda")
  }
  
}

save(res.mat,file="res.rda")

Unf.row = which(row.names(res.mat)=="Unfeasible")
Orac.row= which(row.names(res.mat)=="Oracle")
GBM.row = which(row.names(res.mat)=="Gbm")
Log.row = which(row.names(res.mat)=="Log")
Gam.row = which(row.names(res.mat)=="Gam")
Emp.row = which(row.names(res.mat)=="Empirical")
true.row= which(row.names(res.mat)=="True")

X11()
par(mfrow=c(3,2))
boxplot(res.mat[Unf.row,1],res.mat[Orac.row,1],res.mat[Emp.row,1],
        res.mat[Log.row,1], res.mat[Gam.row,1],res.mat[GBM.row,1],
        names=c("Unfeasible","Oracle","Empirical","Log","GAM","GBM"),
        main="q=0.10",col="blue",outline = T)
abline(h=res.mat[1,1],col="red")

boxplot(res.mat[Unf.row,2],res.mat[Orac.row,2],res.mat[Emp.row,2],
        res.mat[Log.row,2], res.mat[Gam.row,2],res.mat[GBM.row,2],
        names=c("Unfeasible","Oracle","Empirical","Log","GAM","GBM"),
        main="q=0.25",col="blue",outline = T)
abline(h=res.mat[1,2],col="red")

boxplot(res.mat[Unf.row,3],res.mat[Orac.row,3],res.mat[Emp.row,3],
        res.mat[Log.row,3], res.mat[Gam.row,3],res.mat[GBM.row,3],
        names=c("Unfeasible","Oracle","Empirical","Log","GAM","GBM"),
        main="q=0.50",col="blue",outline = T)
abline(h=res.mat[1,3],col="red")

boxplot(res.mat[Unf.row,4],res.mat[Orac.row,4],res.mat[Emp.row,4],
        res.mat[Log.row,4], res.mat[Gam.row,4],res.mat[GBM.row,4],
        names=c("Unfeasible","Oracle","Empirical","Log","GAM","GBM"),
        main="q=0.75",col="blue",outline = T)
abline(h=res.mat[1,4],col="red")

boxplot(res.mat[Unf.row,5],res.mat[Orac.row,5],res.mat[Emp.row,5],
        res.mat[Log.row,5], res.mat[Gam.row,5],res.mat[GBM.row,5],
        names=c("Unfeasible","Oracle","Empirical","Log","GAM","GBM"),
        main="q=0.90",col="blue",outline = T)
abline(h=res.mat[1,5],col="red")


true.vec = res.mat[1,1:5]

MSE.Log = sweep(res.mat[Log.row,c(1:5)],FUN="-",MARGIN = 2,STATS = true.vec)
MSE.Log2= MSE.Log^2
MSE.Log.vec = apply(MSE.Log2,2,mean)

MSE.GAM = sweep(res.mat[Gam.row,c(1:5)],FUN="-",MARGIN = 2,STATS = true.vec)
MSE.GAM2= MSE.GAM^2
MSE.GAM.vec = apply(MSE.GAM2,2,mean)

MSE.GBM = sweep(res.mat[GBM.row,c(1:5)],FUN="-",MARGIN = 2,STATS = true.vec)
MSE.GBM2= MSE.GBM^2
MSE.GBM.vec = apply(MSE.GBM2,2,mean)

MSE.mat = round(rbind(MSE.Log.vec,MSE.GAM.vec,MSE.GBM.vec),2)
Model      = row.names(MSE.mat) = c("Log","GAM","GBM")
MSE.mat = cbind(Model,MSE.mat)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '',
     main="Mean Squared Error",xlim=c(0,10),ylim=c(0,10))
addtable2plot(-0.5,0,MSE.mat, 
              xpad=1, ypad=1,
              bty='o',
              display.rownames = FALSE, 
              hlines = TRUE,
              vlines = TRUE, 
              title = "")

save(res.mat,file="ResCat.rda")
