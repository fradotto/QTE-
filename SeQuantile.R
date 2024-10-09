vcov.rq <- function(object, ...) summary(object, se = "boot", covariance = TRUE)$cov

quant.vec = c(0.10,0.25,0.5,0.75,0.90)

mat.coef = NULL

for(i in 1:length(quant.vec)){

   # Base
  
   reg.base  = rq(wage.hour~em1contrR,tau = quant.vec[i],
                           data = data.mod[idx.vec,])
   b = reg.base$coefficients
   se.base <- deltamethod(g = ~ x2/x1, b, vcov(reg.base))  
   
  
   # Logistic
   
   reg.log     = rq(wage.hour[idx.reg]~em1contrR[idx.reg],
                    tau  = quant.vec[i],
                    data = data.log[idx.reg,],weights = log.weights[idx.reg] )
   b = reg.log$coefficients
   se.log <- deltamethod(g = ~ x2/x1, b, vcov(reg.log)) 
   
   # GAM
   
   reg.gam     = rq(wage.hour~em1contrR,
                    tau  = quant.vec[i],
                    data = data.log[idx.reg,],weights = gam.weights[idx.reg] )
   b = reg.gam$coefficients
   se.gam <- deltamethod(g = ~ x2/x1, b, vcov(reg.gam)) 
   
   # GBM
   
   reg.gbm     = rq(fwag_p/em1hrs~em1contrR,
                    tau = quant.vec[i],
                    data = dat.gbm[idx.gbm,],weights = Gbm.weights[idx.gbm] )
   b = reg.gbm$coefficients
   se.gbm <- deltamethod(g = ~ x2/x1, b, vcov(reg.gbm)) 
   
   c1 = round(reg.base$coefficients[2]/reg.base$coefficients[1],3)
   s1 = round(se.base,3)
   s1 = paste0("(",s1,")")
   
   c2 = round(reg.log$coefficients[2]/reg.log$coefficients[1],3)
   s2 = round(se.log,3)
   s2 = paste0("(",s2,")")
   
   c3 = round(reg.gam$coefficients[2]/reg.gam$coefficients[1],3)
   s3 = round(se.gam,3)
   s3 = paste0("(",s3,")")
   
   c4 = round(reg.gbm$coefficients[2]/reg.gbm$coefficients[1],3)
   s4 = round(se.gbm,3)
   s4 = paste0("(",s4,")")
   col = c(c1,s1,c2,s2,c3,s3,c4,s4)
   
   mat.coef = cbind(mat.coef,col)
   
   
}

row.names(mat.coef) = c("Empirical","","Log","","GAM","","GBM","")
xtable(mat.coef)
