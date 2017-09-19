#################################################################################################
##
##            Sample size of 24 and 32
##
##
################################################################################################

## accrual -> exponential distribution
## visits -> at fixed time

library(markovchain)

Nsim = 2000               # number of simulations

J = 6                        # number of cycles
n = 40                       # number of patients
obswin = 6                 # days
rate = 6                # rate for accrual

K=5 # number of doses


#######################################################################################
# likelihood function PROBIT


skeleton<-function(delta, nu, K , target, a=1)
{
  d<- matrix(NA, nrow=K, ncol=2)
  s<- matrix(NA, nrow=K, ncol=2)
  d[nu,1]<-nu
  d[nu,2]<--(a+qnorm(1-target))
  for (l in (nu+1):K){
    d[l,1]<-l
    d[l,2]<-((qnorm(1-target-delta)+a)*d[l-1,2])/(qnorm(1-target+delta)+a)
  }
  for (l in 0:(nu-2)){
    d[nu-(l+1),1]<-(nu-(l+1))
    d[nu-(l+1),2]<-((qnorm(1-target+delta)+a)*d[nu-l,2])/(qnorm(1-target-delta)+a)
  }
  
  
  s[,1]<-d[,1]
  s[,2]<-round(1-pnorm(-(a+d[,2])), 3)
  d<-round(d,5)
  #cat ("k", "\t", "Doses", "\t", "Skeleton", "\n")
  #for (i in 1:K) {
  #cat(s[i,1], "\t" , d[i,2], "\t", s[i,2], "\n")
  #}
  return(list(d=d, s=s))
}

skel1 <- skeleton(0.04,2,5,0.25)$s[,2]
doselabel <- skeleton(0.04,2,5,0.25)$d[,2]

ltite3bis<-function(b,x1p,z1p,w1p,w2p,alpha0){
  v<-0
  for (i in (1:length(x1p))) {
    v<- v +  ((z1p[i]==1)*log( max(1-w1p[i]*pnorm(alpha0+b[1]*x1p[i]), 2^(-1074)   )) 
              +  (z1p[i]==2)*log( max(w1p[i]*pnorm(alpha0+b[1]*x1p[i]) -  w2p[i]*pnorm(alpha0+b[1]*x1p[i]-b[2]), 2^(-1074) ))
              +  (z1p[i]==3)*log( max(w2p[i]*pnorm(alpha0+b[1]*x1p[i]-b[2]), 2^(-1074) )))
  }
  return(-v)
}



titecrmmc <- function(x,doselabel,y,follow,alpha0,Tmax,target1,target2){
  # x dose level
  # doselabel    numbers used in regression
  # y grade
  # weight
  # alpha0 intercept
  
  x1p <- doselabel[x]
  w1p <- w2p <- follow/Tmax
  w1p[y==2] <- rep(1, sum(y==2) )
  w2p[y==3] <- rep(1, sum(y==3) )
  #est <- optim(c(1,1), ltite3, x1p=x1p, z1p=y, w1p=w1p, w2p=w2p, alpha0=alpha0)$par
  est <- optim(c(1,1), ltite3bis, method = "L-BFGS-B", lower = rep(0.1, 2),
               upper = rep(Inf,2), x1p=x1p, z1p=y, w1p=w1p, w2p=w2p, alpha0=alpha0)$par
  
  p1tox<-pnorm(alpha0+est[1]*doselabel)
  p2tox<-pnorm(alpha0+est[1]*doselabel-est[2])
  cur1<-which(abs(p1tox-target1)==min(abs(p1tox-target1)))
  cur2<-which(abs(p2tox-target2)==min(abs(p2tox-target2)))
  cur<-min(cur1, cur2)
  
  list(newdose=cur, p1tox=p1tox, p2tox=p2tox)
  
}

######################################################################################################
# run simulation for probit T fixed
library(dfcrm)
# J cycles 
# K doses
# n sample size


#invlogit <- function(x) exp(x)/(1-exp(x))

cohort = 3
cohort2 = 1

prior = skel1 #c(0.05,0.1,0.2,0.3)
alpha0=1
target1 = 0.5
target2 = 0.25


TR = Nsim

MTD_titecrmmc24 <- NULL
MTD_titecrmmc32 <- NULL

ntite = resultsTfixprobscen1$ntite
doses_titecrmmc = resultsTfixprobscen1$doses_titecrmmc
tox_titecrmmc = resultsTfixprobscen1$tox_titecrmmc
#- list(MTD_titecrmmc = MTD_titecrmmc, doses_titecrmmc = doses_titecrmmc, tox_titecrmmc = tox_titecrmmc,
                      #       ntite=ntite)


for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  y=tox_titecrmmc[tr,1:24]
  x=doses_titecrmmc[tr,1:24]
  
  cond1 <- any(y==3) & any(y==2) & any(y==1)
  
  cond2 <- (any(y==3) & any(y==2) &! any(y==1)) |
    (any(y==3) & any(y==1) &! any(y==2))  
  
  cond3 <- (any(y==2) & any(y==1) &! any(y==3))    
  
  
  if (cond1) {
    results <- titecrmmc(x,doselabel,y=y,follow=rep(6,24),alpha0,Tmax=J,target1,target2)
    MTD <- min(results$newdose, max(x)+1, K)
    MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
  } else {
      if (cond2) {
        grad = 3
        target_reg = target2
      } else {
        grad = 2
        target_reg = target1
      }
      # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
      y <- rep(0, length(x)) 
      postox <- which(x==grad )
      y[postox] <- rep(1, length(postox))
      
      results <- titecrm(prior, target_reg, y, x, followup=rep(J,24), obswin=J)
      MTD <- min(results$mtd, max(x)+1, K)
      MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
    } 
    
  
  y=tox_titecrmmc[tr,1:32]
  x=doses_titecrmmc[tr,1:32]
  
  cond1 <- any(y==3) & any(y==2) & any(y==1)
  
  cond2 <- (any(y==3) & any(y==2) &! any(y==1)) |
    (any(y==3) & any(y==1) &! any(y==2))  
  
  cond3 <- (any(y==2) & any(y==1) &! any(y==3))     
  
  
  if (cond1) {
    results <- titecrmmc(x,doselabel,y=tox_titecrmmc[tr,1:32],follow=rep(6,32),alpha0,Tmax=J,target1,target2)
    MTD <- min(results$newdose, max(x)+1, K)
    MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
  } else {
    if (cond2) {
      grad = 3
      target_reg = target2
    } else {
      grad = 2
      target_reg = target1
    }
    # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
    y <- rep(0, length(x)) 
    postox <- which(x==grad )
    y[postox] <- rep(1, length(postox))
    
    results <- titecrm(prior, target_reg, y, x, followup=rep(J,32), obswin=J)
    MTD <- min(results$mtd, max(x)+1, K)
    MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
  } 
  
}
  
  resultsTfixprobscen1$MTD_titecrmmc24 = MTD_titecrmmc24
  resultsTfixprobscen1$MTD_titecrmmc32 = MTD_titecrmmc32
  write.table(MTD_titecrmmc24, file="scenario1fixprob24.txt")
  write.table(MTD_titecrmmc32, file="scenario1fixprob34.txt")
    
save.image()
###############################################################################################################
# for poisson accrual 

library(dfcrm)


cohort = 3
cohort2 = 1

prior = skel1 
alpha0=1
target1 = 0.5
target2 = 0.25


TR = Nsim

MTD_titecrmmc24 <- NULL
MTD_titecrmmc32 <- NULL

ntite = resultsTpoisprobscen1$ntite
doses_titecrmmc = resultsTpoisprobscen1$doses_titecrmmc
tox_titecrmmc = resultsTpoisprobscen1$tox_titecrmmc
#- list(MTD_titecrmmc = MTD_titecrmmc, doses_titecrmmc = doses_titecrmmc, tox_titecrmmc = tox_titecrmmc,
#       ntite=ntite)


for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  y=tox_titecrmmc[tr,1:24]
  x=doses_titecrmmc[tr,1:24]
  
  cond1 <- any(y==3) & any(y==2) & any(y==1)
  
  cond2 <- (any(y==3) & any(y==2) &! any(y==1)) |
    (any(y==3) & any(y==1) &! any(y==2))  
  
  cond3 <- (any(y==2) & any(y==1) &! any(y==3))    
  
  
  if (cond1) {
    results <- titecrmmc(x,doselabel,y=y,follow=rep(6,24),alpha0,Tmax=J,target1,target2)
    MTD <- min(results$newdose, max(x)+1, K)
    MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
  } else {
    if (cond2) {
      grad = 3
      target_reg = target2
    } else {
      grad = 2
      target_reg = target1
    }
    # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
    y <- rep(0, length(x)) 
    postox <- which(x==grad )
    y[postox] <- rep(1, length(postox))
    
    results <- titecrm(prior, target_reg, y, x, followup=rep(J,24), obswin=J)
    MTD <- min(results$mtd, max(x)+1, K)
    MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
  } 
  
  
  y=tox_titecrmmc[tr,1:32]
  x=doses_titecrmmc[tr,1:32]
  
  cond1 <- any(y==3) & any(y==2) & any(y==1)
  
  cond2 <- (any(y==3) & any(y==2) &! any(y==1)) |
    (any(y==3) & any(y==1) &! any(y==2))  
  
  cond3 <- (any(y==2) & any(y==1) &! any(y==3))     
  
  
  if (cond1) {
    results <- titecrmmc(x,doselabel,y=tox_titecrmmc[tr,1:32],follow=rep(6,32),alpha0,Tmax=J,target1,target2)
    MTD <- min(results$newdose, max(x)+1, K)
    MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
  } else {
    if (cond2) {
      grad = 3
      target_reg = target2
    } else {
      grad = 2
      target_reg = target1
    }
    # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
    y <- rep(0, length(x)) 
    postox <- which(x==grad )
    y[postox] <- rep(1, length(postox))
    
    results <- titecrm(prior, target_reg, y, x, followup=rep(J,32), obswin=J)
    MTD <- min(results$mtd, max(x)+1, K)
    MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
  } 
  
}

resultsTpoisprobscen1$MTD_titecrmmc24 = MTD_titecrmmc24
resultsTpoisprobscen1$MTD_titecrmmc32 = MTD_titecrmmc32
write.table(MTD_titecrmmc24, file="scenario1poisprob24.txt")
write.table(MTD_titecrmmc32, file="scenario1poisprob34.txt")


save.image()
#######################################################################################################
#######################################################################################################
#######################################################################################################
# Likelihood for empiric model

#the skeleton assumes that F(d,beta)=d^(beta))

skeleton_empiric<-function(delta, nu, K , target, a=1)
{
  d<- matrix(NA, nrow=K, ncol=2)
  s<- matrix(NA, nrow=K, ncol=2)
  d[nu,1]<-nu
  d[nu,2]<-target #-(a+qnorm(1-target))
  for (l in (nu+1):K){
    d[l,1]<-l
    d[l,2]<-exp((log(target+delta)*log(d[l-1,2]))/log(target-delta))
  }
  for (l in 0:(nu-2)){
    d[nu-(l+1),1]<-(nu-(l+1))
    d[nu-(l+1),2]<-exp((log(target-delta)*log(d[nu-l,2]))/log(target+delta))
  }
  
  
  s[,1]<-d[,1]
  s[,2]<-d[,2]#round(1-pnorm(-(a+d[,2])), 3)
  d<-round(d,5)
  #cat ("k", "\t", "Doses", "\t", "Skeleton", "\n")
  #for (i in 1:K) {
  #cat(s[i,1], "\t" , d[i,2], "\t", s[i,2], "\n")
  #}
  return(list(d=d, s=s))
}

skel1 <- doselabel <- skeleton_empiric(0.06,2,5,0.25)$d[,2]

################################) Emipiric working model



ltite3empiric<-function(b,x1p,z1p,w1p,w2p){
  v<-0
  for (i in (1:length(x1p))) {
    v<- v +  ((z1p[i]==1)*log( max(1-w1p[i]*x1p[i]^b[1], 2^(-1074)   )) 
              +  (z1p[i]==2)*log( max(w1p[i]*x1p[i]^b[1] -  w2p[i]*x1p[i]^(b[1]+b[2]), 2^(-1074) ))
              +  (z1p[i]==3)*log( max(w2p[i]*x1p[i]^(b[1]+b[2]), 2^(-1074) )))
  }
  return(-v)
}



titecrmmc <- function(x,doselabel,y,follow,alpha0,Tmax,target1,target2){
  # x dose level
  # doselabel    numbers used in regression
  # y grade
  # weight
  # alpha0 intercept
  
  x1p <- doselabel[x]
  w1p <- w2p <- follow/Tmax
  w1p[y==2] <- rep(1, sum(y==2) )
  w2p[y==3] <- rep(1, sum(y==3) )
  #est <- optim(c(1,1), ltite3, x1p=x1p, z1p=y, w1p=w1p, w2p=w2p, alpha0=alpha0)$par
  est <- optim(c(1,1), ltite3empiric, method = "L-BFGS-B", lower = rep(0.1, 2),
               upper = rep(Inf,2), x1p=x1p, z1p=y, w1p=w1p, w2p=w2p)$par
  
  p1tox<-doselabel^(est[1])
  p2tox<-doselabel^(est[1] + est[2])
  cur1<-which(abs(p1tox-target1)==min(abs(p1tox-target1)))
  cur2<-which(abs(p2tox-target2)==min(abs(p2tox-target2)))
  cur<-min(cur1, cur2)
  
  list(newdose=cur, p1tox=p1tox, p2tox=p2tox)
  
}

####################################################################################################
### starting simulations T fixed


library(dfcrm)
# J cycles 
# K doses
# n sample size


#invlogit <- function(x) exp(x)/(1-exp(x))

cohort = 3
cohort2 = 1

prior = skel1 #c(0.05,0.1,0.2,0.3)
alpha0=1
target1 = 0.5
target2 = 0.25


TR = Nsim

MTD_titecrmmc24 <- NULL
MTD_titecrmmc32 <- NULL

ntite = resultsTfixempscen1$ntite
doses_titecrmmc = resultsTfixempscen1$doses_titecrmmc
tox_titecrmmc = resultsTfixempscen1$tox_titecrmmc
#- list(MTD_titecrmmc = MTD_titecrmmc, doses_titecrmmc = doses_titecrmmc, tox_titecrmmc = tox_titecrmmc,
#       ntite=ntite)


for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  y=tox_titecrmmc[tr,1:24]
  x=doses_titecrmmc[tr,1:24]
  
  cond1 <- any(y==3) & any(y==2) & any(y==1)
  
  cond2 <- (any(y==3) & any(y==2) &! any(y==1)) |
    (any(y==3) & any(y==1) &! any(y==2))  
  
  cond3 <- (any(y==2) & any(y==1) &! any(y==3))    
  
  
  if (cond1) {
    results <- titecrmmc(x,doselabel,y=y,follow=rep(6,24),alpha0,Tmax=J,target1,target2)
    MTD <- min(results$newdose, max(x)+1, K)
    MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
  } else {
    if (cond2) {
      grad = 3
      target_reg = target2
    } else {
      grad = 2
      target_reg = target1
    }
    # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
    y <- rep(0, length(x)) 
    postox <- which(x==grad )
    y[postox] <- rep(1, length(postox))
    
    results <- titecrm(prior, target_reg, y, x, followup=rep(J,24), obswin=J)
    MTD <- min(results$mtd, max(x)+1, K)
    MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
  } 
  
  
  y=tox_titecrmmc[tr,1:32]
  x=doses_titecrmmc[tr,1:32]
  
  cond1 <- any(y==3) & any(y==2) & any(y==1)
  
  cond2 <- (any(y==3) & any(y==2) &! any(y==1)) |
    (any(y==3) & any(y==1) &! any(y==2))  
  
  cond3 <- (any(y==2) & any(y==1) &! any(y==3))     
  
  
  if (cond1) {
    results <- titecrmmc(x,doselabel,y=tox_titecrmmc[tr,1:32],follow=rep(6,32),alpha0,Tmax=J,target1,target2)
    MTD <- min(results$newdose, max(x)+1, K)
    MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
  } else {
    if (cond2) {
      grad = 3
      target_reg = target2
    } else {
      grad = 2
      target_reg = target1
    }
    # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
    y <- rep(0, length(x)) 
    postox <- which(x==grad )
    y[postox] <- rep(1, length(postox))
    
    results <- titecrm(prior, target_reg, y, x, followup=rep(J,32), obswin=J)
    MTD <- min(results$mtd, max(x)+1, K)
    MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
  } 
  
}

resultsTfixempscen1$MTD_titecrmmc24 = MTD_titecrmmc24
resultsTfixempscen1$MTD_titecrmmc32 = MTD_titecrmmc32
write.table(MTD_titecrmmc24, file="scenario1fixemp24.txt")
write.table(MTD_titecrmmc32, file="scenario1fixemp34.txt")

save.image()

######################## pois accrual #####################################################################
library(dfcrm)
# J cycles 
# K doses
# n sample size


#invlogit <- function(x) exp(x)/(1-exp(x))

cohort = 3
cohort2 = 1

prior = skel1 #c(0.05,0.1,0.2,0.3)
alpha0=1
target1 = 0.5
target2 = 0.25


TR = Nsim

MTD_titecrmmc24 <- NULL
MTD_titecrmmc32 <- NULL

ntite = resultsTpoisempscen1$ntite
doses_titecrmmc = resultsTpoisempscen1$doses_titecrmmc
tox_titecrmmc = resultsTpoisempscen1$tox_titecrmmc
#- list(MTD_titecrmmc = MTD_titecrmmc, doses_titecrmmc = doses_titecrmmc, tox_titecrmmc = tox_titecrmmc,
#       ntite=ntite)


for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  y=tox_titecrmmc[tr,1:24]
  x=doses_titecrmmc[tr,1:24]
  
  cond1 <- any(y==3) & any(y==2) & any(y==1)
  
  cond2 <- (any(y==3) & any(y==2) &! any(y==1)) |
    (any(y==3) & any(y==1) &! any(y==2))  
  
  cond3 <- (any(y==2) & any(y==1) &! any(y==3))    
  
  
  if (cond1) {
    results <- titecrmmc(x,doselabel,y=y,follow=rep(6,24),alpha0,Tmax=J,target1,target2)
    MTD <- min(results$newdose, max(x)+1, K)
    MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
  } else {
    if (cond2) {
      grad = 3
      target_reg = target2
    } else {
      grad = 2
      target_reg = target1
    }
    # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
    y <- rep(0, length(x)) 
    postox <- which(x==grad )
    y[postox] <- rep(1, length(postox))
    
    results <- titecrm(prior, target_reg, y, x, followup=rep(J,24), obswin=J)
    MTD <- min(results$mtd, max(x)+1, K)
    MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
  } 
  
  
  y=tox_titecrmmc[tr,1:32]
  x=doses_titecrmmc[tr,1:32]
  
  cond1 <- any(y==3) & any(y==2) & any(y==1)
  
  cond2 <- (any(y==3) & any(y==2) &! any(y==1)) |
    (any(y==3) & any(y==1) &! any(y==2))  
  
  cond3 <- (any(y==2) & any(y==1) &! any(y==3))     
  
  
  if (cond1) {
    results <- titecrmmc(x,doselabel,y=tox_titecrmmc[tr,1:32],follow=rep(6,32),alpha0,Tmax=J,target1,target2)
    MTD <- min(results$newdose, max(x)+1, K)
    MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
  } else {
    if (cond2) {
      grad = 3
      target_reg = target2
    } else {
      grad = 2
      target_reg = target1
    }
    # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
    y <- rep(0, length(x)) 
    postox <- which(x==grad )
    y[postox] <- rep(1, length(postox))
    
    results <- titecrm(prior, target_reg, y, x, followup=rep(J,32), obswin=J)
    MTD <- min(results$mtd, max(x)+1, K)
    MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
  } 
  
}

resultsTpoisempscen1$MTD_titecrmmc24 = MTD_titecrmmc24
resultsTpoisempscen1$MTD_titecrmmc32 = MTD_titecrmmc32
write.table(MTD_titecrmmc24, file="scenario1poisemp24.txt")
write.table(MTD_titecrmmc32, file="scenario1poisemp34.txt")

save.image()









