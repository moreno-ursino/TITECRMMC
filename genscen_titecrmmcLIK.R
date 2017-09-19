#################################################################################################
##
##            First part: generating scenario's clinical trials
##            Second part: run clinical trial using TITECRMMC-likelihood (probit and empirical)
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

PI <-matrix( c( 0.892, 0.078, 0.03, 0.86, 0.09, 0.05, 0.83, 0.10, 0.07, 0.79, 0.12, 0.09, 0.74, 0.14, 0.12,
                 0, 0.872, 0.128,    0, 0.87, 0.13,    0, 0.84, 0.16,    0, 0.82, 0.18,   0, 0.80, 0.20,
                 0,    0,    1,    0,    0,    1,    0,    0,    1,    0,    0,    1, 0,    0,    1),
              byrow=T, nrow=3 )


accrual = "fixed"
#accrual = "poisson"

SED = 32
# generating dataset
for (tr in 1:Nsim){
  
  data_complete = NULL
  
  set.seed(SED)
  if (accrual=="fixed") {
    T_entrance = c(0,rep(obswin/rate,n-1))
  } else T_entrance = c(0,rexp(n-1,rate/obswin))          # time of accrual for each patients
  T_entrance = cumsum(T_entrance)#round(cumsum(T_entrance))  
  
  for (m in 1:n){                                       # for each patient
    
    times <- cumsum(c(T_entrance[m]+1,rep(1,J-1)))
    datap <- cbind(rep(m,J),1:J,times)                        # number of patients and cycle number
    
    for (d in 1:K){
      patient <- new("markovchain", states = as.character(1:3),
                     transitionMatrix = PI[,(3*d-2):(3*d)],
                     name = "Patmarkov")
      grades=rmarkovchain(n = J, object = patient, t0 = "1")
      datap= cbind(datap,as.numeric(grades))
    }
    data_complete= rbind(data_complete,datap)
  }
  dimnames(data_complete)[[2]] <- c("patient","cycle","time",paste('dose_',1:K,sep=''))
  SED = SED + tr
  eval(parse(text = paste("data",tr, " <- list(data_complete=data_complete, T_entrance=T_entrance)", sep="")))
}
save.image()



SED = 400
# generating dataset
accrual="poisson"
for (tr in 1:Nsim){
  
  
  set.seed(SED)
  if (accrual=="fixed") {
    T_entrance = c(0,rep(obswin/rate,n-1))
  } else T_entrance = c(0,rexp(n-1,rate/obswin))          # time of accrual for each patients
  T_entrance = cumsum(T_entrance)#round(cumsum(T_entrance))  
  
  SED = SED + tr
  eval(parse(text = paste("data",tr, "$T_entrance2 <- T_entrance", sep="")))
}
save.image()

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

MTD_titecrmmc <- NULL
doses_titecrmmc <- NULL
tox_titecrmmc <- NULL
xtox_titecrmmc <- NULL
ntite <- NULL

for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  eval(parse(text = paste("tox <- data",tr,"$data_complete", sep="")))
  eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance", sep="")))
  
  
  tox <- data.frame(tox)
  T_entrance <- c(T_entrance,Inf)
  MTD=-1
  ntitetr=0
  
  x <- rep(1,cohort)                     # first dose to the first cohort
  M = n/cohort                           # number of total cohorts
  
  # data_reg = data_complete[1:(cohort*J),1:4]
  # y <- tox[cbind(1:length(x),x)]  
  
  
  
  time = T_entrance[cohort+1]   # time of entrance of the first patient of second cohort
  n_fin = 0                     # number of patient with complete follow up
  
  data_reg_c <- NULL              # building the dataset for TITE: complete, with patients who ended the followup
  data_reg_i <- NULL              # during follow up
  
  for (l in 1:cohort){
    pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
    pos <- max(pos,0.1)
    if (pos == J) {
      follow = J
      grade = tox[l*J, (x[l]+3)]
      n_fin = l
      data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
    }
    else {
      follow = pos
      grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
      data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
    }
  }
  
  #dimnames(data_reg_c)[[2]] <- c("patient", "follow", "grade")
  #dimnames(data_reg_i)[[2]] <- c("patient", "follow", "grade")
  
  data_reg <- rbind(data_reg_c,data_reg_i)
  # check if we can apply TITECRM or not
  if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
      any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
      any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
    stage1 <- stage2 <- FALSE
  } else {
    if ( length(unique(data_reg[,3]))==1 ) {
      stage1 <- TRUE
      stage2 <- FALSE
    } else {
      stage1 <- FALSE
      stage2 <- TRUE
    } 
  }
  
  pat=cohort
  
  if (stage1) {
    for(i in 2:M) {
      
      if ( sum(data_reg[,3]==3)==length(x) ) {
        MTD=0
        break
      } else {
        if ( sum(data_reg[,3]==2)==length(x) ) {
          x <- c(x,rep(x[length(x)], cohort)) 
        } else {
          x <- c(x,rep(min(max(x)+1,K), cohort))
        }
      }
      
      # dose for the cohort
      time = T_entrance[i*cohort+1]                   # time of entrance of the first of the next cohort
      pat=pat+cohort
      
      data_reg_i <- NULL                               # reset the temporary dataset
      
      for (l in (n_fin+1):(i*cohort)){
        pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
        pos <- max(pos,0.1)
        if (pos == J) {
          follow = J
          grade = tox[l*J, (x[l]+3)]
          n_fin = l
          data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
        }
        else {
          follow = pos
          grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
          data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
        }
      }
      
      data_reg <- rbind(data_reg_c,data_reg_i)
      
      if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
          any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
          any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
        stage2 <- FALSE
        ntitetr = length(x)
        break
      } else {
        if ( (any(data_reg[,3]==2) & any(data_reg[,3]==3)) | 
             (any(data_reg[,3]==1) & any(data_reg[,3]==2)) |
             (any(data_reg[,3]==1) & any(data_reg[,3]==3)) ) {
          stage2 <- TRUE
          break
        }
      }
    }
  }
  
  if (MTD==0) {
    MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
    doses_titecrmmc <- rbind(doses_titecrmmc,rep(0,40))
    tox_titecrmmc <- rbind(tox_titecrmmc,rep(1,40))
    ntite <- c(ntite,0)
  } else {
    
    for (i in seq(pat+1,n,1)){ 
      
      if (stage2) {  
        if (any(data_reg[,3]==3)) {
          grad = 3
          target_reg = target2
        } else {
          grad = 2
          target_reg = target1
        }
        # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
        y <- rep(0, dim(data_reg)[1]) 
        postox <- which(data_reg[,3]==grad )
        y[postox] <- rep(1, length(postox))
        
        results <- titecrm(prior, target_reg, y, x, followup=data_reg[,2], obswin=J)
        
        newdose <- min(results$mtd, max(x)+1, K)
        # check on the skipping dose
        x <- c(x,rep(newdose,cohort2))
        
        time = T_entrance[i*cohort2+1]  # first patient of the next cohort
        
        data_reg_i <- NULL  
        for (l in (n_fin+1):(i*cohort2)){
          pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
          pos <- max(pos,0.1)
          if (pos == J) {
            follow = J
            grade = tox[l*J, (x[l]+3)]
            n_fin = l
            data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
          }
          else {
            follow = pos
            grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
            data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
          }
        }
        
        data_reg <- rbind(data_reg_c,data_reg_i)
        
        if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
            any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
            any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
          stage2 <- FALSE
          ntitetr = length(x)
        }
        
      } else {
        ##### here stage 3
        
        
        results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
        
        newdose <- min(results$newdose, max(x)+1, K)
        # check on the skipping dose
        x <- c(x,rep(newdose,cohort2))
        
        time = T_entrance[i*cohort2+1]  # first patient of the next cohort
        
        data_reg_i <- NULL  
        for (l in (n_fin+1):(i*cohort2)){
          pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
          pos <- max(pos,0.1)
          if (pos == J) {
            follow = J
            grade = tox[l*J, (x[l]+3)]
            n_fin = l
            data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
          }
          else {
            follow = pos
            grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
            data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
          }
        }
        
        data_reg <- rbind(data_reg_c,data_reg_i)
        
        
      }
    }
    
    ##### analysis of complete data
    if (stage2) {
      if (any(data_reg[,3]==3)) {
        grad = 3
        target_reg = target2
      } else {
        grad = 2
        target_reg = target1
      }
      # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
      y <- rep(0, dim(data_reg)[1]) 
      postox <- which(data_reg[,3]==grad )
      y[postox] <- rep(1, length(postox))
      
      results <- titecrm(prior, target_reg, y, x, followup=data_reg[,2], obswin=J)
      MTD <- min(results$mtd, max(x)+1, K)
      ntitetr=0
    } else {
      results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
      MTD <- min(results$newdose, max(x)+1, K)
    }
    
    MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
    doses_titecrmmc <- rbind(doses_titecrmmc,x)
    tox_titecrmmc <- rbind(tox_titecrmmc,data_reg[,3])
    ntite <- c(ntite,ntitetr)
  }
}

library(xtable)

resultsTfixprobscen1 <- list(MTD_titecrmmc = MTD_titecrmmc, doses_titecrmmc = doses_titecrmmc, tox_titecrmmc = tox_titecrmmc,
                         ntite=ntite)
write.table(MTD_titecrmmc, file="scenario1fixprob.txt")
write.table(ntite, file="scenario1fixprobntite.txt")


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

MTD_titecrmmc <- NULL
doses_titecrmmc <- NULL
tox_titecrmmc <- NULL
xtox_titecrmmc <- NULL
ntite <- NULL

for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  eval(parse(text = paste("tox <- data",tr,"$data_complete", sep="")))
  # eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance", sep="")))
  
  ################################## for pois
  eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance2", sep="")))
  tpois <- NULL
  for (poiss in 1:length(T_entrance) ){
    times <- cumsum(c(T_entrance[poiss]+1,rep(1,J-1)))
    tpois <- c(tpois,times)
  }
  tox <- data.frame(tox)
  tox$time <- tpois
  #########################################################################
  
  
  tox <- data.frame(tox)
  T_entrance <- c(T_entrance,Inf)
  MTD=-1
  ntitetr=0
  
  x <- rep(1,cohort)                     # first dose to the first cohort
  M = n/cohort                           # number of total cohorts
  
  # data_reg = data_complete[1:(cohort*J),1:4]
  # y <- tox[cbind(1:length(x),x)]  
  
  
  
  time = T_entrance[cohort+1]   # time of entrance of the first patient of second cohort
  n_fin = 0                     # number of patient with complete follow up
  
  data_reg_c <- NULL              # building the dataset for TITE: complete, with patients who ended the followup
  data_reg_i <- NULL              # during follow up
  
  for (l in 1:cohort){
    pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
    pos <- max(pos,0.1)
    if (pos == J) {
      follow = J
      grade = tox[l*J, (x[l]+3)]
      n_fin = l
      data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
    }
    else {
      follow = pos
      grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
      data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
    }
  }
  
  #dimnames(data_reg_c)[[2]] <- c("patient", "follow", "grade")
  #dimnames(data_reg_i)[[2]] <- c("patient", "follow", "grade")
  
  data_reg <- rbind(data_reg_c,data_reg_i)
  # check if we can apply TITECRM or not
  if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
      any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
      any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
    stage1 <- stage2 <- FALSE
  } else {
    if ( length(unique(data_reg[,3]))==1 ) {
      stage1 <- TRUE
      stage2 <- FALSE
    } else {
      stage1 <- FALSE
      stage2 <- TRUE
    } 
  }
  
  pat=cohort
  
  if (stage1) {
    for(i in 2:M) {
      
      if ( sum(data_reg[,3]==3)==length(x) ) {
        MTD=0
        break
      } else {
        if ( sum(data_reg[,3]==2)==length(x) ) {
          x <- c(x,rep(x[length(x)], cohort)) 
        } else {
          x <- c(x,rep(min(max(x)+1,K), cohort))
        }
      }
      
      # dose for the cohort
      time = T_entrance[i*cohort+1]                   # time of entrance of the first of the next cohort
      pat=pat+cohort
      
      data_reg_i <- NULL                               # reset the temporary dataset
      
      for (l in (n_fin+1):(i*cohort)){
        pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
        pos <- max(pos,0.1)
        if (pos == J) {
          follow = J
          grade = tox[l*J, (x[l]+3)]
          n_fin = l
          data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
        }
        else {
          follow = pos
          grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
          data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
        }
      }
      
      data_reg <- rbind(data_reg_c,data_reg_i)
      
      if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
          any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
          any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
        stage2 <- FALSE
        ntitetr = length(x)
        break
      } else {
        if ( (any(data_reg[,3]==2) & any(data_reg[,3]==3)) | 
             (any(data_reg[,3]==1) & any(data_reg[,3]==2)) |
             (any(data_reg[,3]==1) & any(data_reg[,3]==3)) ) {
          stage2 <- TRUE
          break
        }
      }
    }
  }
  
  if (MTD==0) {
    MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
    doses_titecrmmc <- rbind(doses_titecrmmc,rep(0,40))
    tox_titecrmmc <- rbind(tox_titecrmmc,rep(1,40))
    ntite <- c(ntite,0)
  } else {
    
    for (i in seq(pat+1,n,1)){ 
      
      if (stage2) {  
        if (any(data_reg[,3]==3)) {
          grad = 3
          target_reg = target2
        } else {
          grad = 2
          target_reg = target1
        }
        # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
        y <- rep(0, dim(data_reg)[1]) 
        postox <- which(data_reg[,3]==grad )
        y[postox] <- rep(1, length(postox))
        
        results <- titecrm(prior, target_reg, y, x, followup=data_reg[,2], obswin=J)
        
        newdose <- min(results$mtd, max(x)+1, K)
        # check on the skipping dose
        x <- c(x,rep(newdose,cohort2))
        
        time = T_entrance[i*cohort2+1]  # first patient of the next cohort
        
        data_reg_i <- NULL  
        for (l in (n_fin+1):(i*cohort2)){
          pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
          pos <- max(pos,0.1)
          if (pos == J) {
            follow = J
            grade = tox[l*J, (x[l]+3)]
            n_fin = l
            data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
          }
          else {
            follow = pos
            grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
            data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
          }
        }
        
        data_reg <- rbind(data_reg_c,data_reg_i)
        
        if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
            any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
            any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
          stage2 <- FALSE
          ntitetr = length(x)
        }
        
      } else {
        ##### here stage 3
        
        
        results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
        
        newdose <- min(results$newdose, max(x)+1, K)
        # check on the skipping dose
        x <- c(x,rep(newdose,cohort2))
        
        time = T_entrance[i*cohort2+1]  # first patient of the next cohort
        
        data_reg_i <- NULL  
        for (l in (n_fin+1):(i*cohort2)){
          pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
          pos <- max(pos,0.1)
          if (pos == J) {
            follow = J
            grade = tox[l*J, (x[l]+3)]
            n_fin = l
            data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
          }
          else {
            follow = pos
            grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
            data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
          }
        }
        
        data_reg <- rbind(data_reg_c,data_reg_i)
        
        
      }
    }
    
    ##### analysis of complete data
    if (stage2) {
      if (any(data_reg[,3]==3)) {
        grad = 3
        target_reg = target2
      } else {
        grad = 2
        target_reg = target1
      }
      # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
      y <- rep(0, dim(data_reg)[1]) 
      postox <- which(data_reg[,3]==grad )
      y[postox] <- rep(1, length(postox))
      
      results <- titecrm(prior, target_reg, y, x, followup=data_reg[,2], obswin=J)
      MTD <- min(results$mtd, max(x)+1, K)
      ntitetr=0
    } else {
      results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
      MTD <- min(results$newdose, max(x)+1, K)
    }
    
    MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
    doses_titecrmmc <- rbind(doses_titecrmmc,x)
    tox_titecrmmc <- rbind(tox_titecrmmc,data_reg[,3])
    ntite <- c(ntite,ntitetr)
  }
}

library(xtable)

resultsTpoisprobscen1 <- list(MTD_titecrmmc = MTD_titecrmmc, doses_titecrmmc = doses_titecrmmc, tox_titecrmmc = tox_titecrmmc,
                         ntite=ntite)
write.table(MTD_titecrmmc, file="scenario1poisprob.txt")
write.table(ntite, file="scenario1poisprobntite.txt")


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

MTD_titecrmmc <- NULL
doses_titecrmmc <- NULL
tox_titecrmmc <- NULL
xtox_titecrmmc <- NULL
ntite <- NULL

for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  eval(parse(text = paste("tox <- data",tr,"$data_complete", sep="")))
  eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance", sep="")))

  
  tox <- data.frame(tox)
  T_entrance <- c(T_entrance,Inf)
  MTD=-1
  ntitetr=0
  
  x <- rep(1,cohort)                     # first dose to the first cohort
  M = n/cohort                           # number of total cohorts
  
  # data_reg = data_complete[1:(cohort*J),1:4]
  # y <- tox[cbind(1:length(x),x)]  
  
  
  
  time = T_entrance[cohort+1]   # time of entrance of the first patient of second cohort
  n_fin = 0                     # number of patient with complete follow up
  
  data_reg_c <- NULL              # building the dataset for TITE: complete, with patients who ended the followup
  data_reg_i <- NULL              # during follow up
  
  for (l in 1:cohort){
    pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
    pos <- max(pos,0.1)
    if (pos == J) {
      follow = J
      grade = tox[l*J, (x[l]+3)]
      n_fin = l
      data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
    }
    else {
      follow = pos
      grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
      data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
    }
  }
  
  #dimnames(data_reg_c)[[2]] <- c("patient", "follow", "grade")
  #dimnames(data_reg_i)[[2]] <- c("patient", "follow", "grade")
  
  data_reg <- rbind(data_reg_c,data_reg_i)
  # check if we can apply TITECRM or not
  if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
      any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
      any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
    stage1 <- stage2 <- FALSE
  } else {
    if ( length(unique(data_reg[,3]))==1 ) {
      stage1 <- TRUE
      stage2 <- FALSE
    } else {
      stage1 <- FALSE
      stage2 <- TRUE
    } 
  }
  
  pat=cohort
  
  if (stage1) {
    for(i in 2:M) {
      
      if ( sum(data_reg[,3]==3)==length(x) ) {
        MTD=0
        break
      } else {
        if ( sum(data_reg[,3]==2)==length(x) ) {
          x <- c(x,rep(x[length(x)], cohort)) 
        } else {
          x <- c(x,rep(min(max(x)+1,K), cohort))
        }
      }
      
      # dose for the cohort
      time = T_entrance[i*cohort+1]                   # time of entrance of the first of the next cohort
      pat=pat+cohort
      
      data_reg_i <- NULL                               # reset the temporary dataset
      
      for (l in (n_fin+1):(i*cohort)){
        pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
        pos <- max(pos,0.1)
        if (pos == J) {
          follow = J
          grade = tox[l*J, (x[l]+3)]
          n_fin = l
          data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
        }
        else {
          follow = pos
          grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
          data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
        }
      }
      
      data_reg <- rbind(data_reg_c,data_reg_i)
      
      if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
          any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
          any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
        stage2 <- FALSE
        ntitetr = length(x)
        break
      } else {
        if ( (any(data_reg[,3]==2) & any(data_reg[,3]==3)) | 
             (any(data_reg[,3]==1) & any(data_reg[,3]==2)) |
             (any(data_reg[,3]==1) & any(data_reg[,3]==3)) ) {
          stage2 <- TRUE
          break
        }
      }
    }
  }
  
  if (MTD==0) {
    MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
    doses_titecrmmc <- rbind(doses_titecrmmc,rep(0,40))
    tox_titecrmmc <- rbind(tox_titecrmmc,rep(1,40))
    ntite <- c(ntite,0)
  } else {
    
    for (i in seq(pat+1,n,1)){ 
      
      if (stage2) {  
        if (any(data_reg[,3]==3)) {
          grad = 3
          target_reg = target2
        } else {
          grad = 2
          target_reg = target1
        }
        # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
        y <- rep(0, dim(data_reg)[1]) 
        postox <- which(data_reg[,3]==grad )
        y[postox] <- rep(1, length(postox))
        
        results <- titecrm(prior, target_reg, y, x, followup=data_reg[,2], obswin=J)
        
        newdose <- min(results$mtd, max(x)+1, K)
        # check on the skipping dose
        x <- c(x,rep(newdose,cohort2))
        
        time = T_entrance[i*cohort2+1]  # first patient of the next cohort
        
        data_reg_i <- NULL  
        for (l in (n_fin+1):(i*cohort2)){
          pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
          pos <- max(pos,0.1)
          if (pos == J) {
            follow = J
            grade = tox[l*J, (x[l]+3)]
            n_fin = l
            data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
          }
          else {
            follow = pos
            grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
            data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
          }
        }
        
        data_reg <- rbind(data_reg_c,data_reg_i)
        
        if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
            any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
            any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
          stage2 <- FALSE
          ntitetr = length(x)
        }
        
      } else {
        ##### here stage 3
        
        
        results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
        
        newdose <- min(results$newdose, max(x)+1, K)
        # check on the skipping dose
        x <- c(x,rep(newdose,cohort2))
        
        time = T_entrance[i*cohort2+1]  # first patient of the next cohort
        
        data_reg_i <- NULL  
        for (l in (n_fin+1):(i*cohort2)){
          pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
          pos <- max(pos,0.1)
          if (pos == J) {
            follow = J
            grade = tox[l*J, (x[l]+3)]
            n_fin = l
            data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
          }
          else {
            follow = pos
            grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
            data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
          }
        }
        
        data_reg <- rbind(data_reg_c,data_reg_i)
        
        
      }
    }
    
    ##### analysis of complete data
    if (stage2) {
      if (any(data_reg[,3]==3)) {
        grad = 3
        target_reg = target2
      } else {
        grad = 2
        target_reg = target1
      }
      # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
      y <- rep(0, dim(data_reg)[1]) 
      postox <- which(data_reg[,3]==grad )
      y[postox] <- rep(1, length(postox))
      
      results <- titecrm(prior, target_reg, y, x, followup=data_reg[,2], obswin=J)
      MTD <- min(results$mtd, max(x)+1, K)
      ntitetr=0
    } else {
      results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
      MTD <- min(results$newdose, max(x)+1, K)
    }
    
    MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
    doses_titecrmmc <- rbind(doses_titecrmmc,x)
    tox_titecrmmc <- rbind(tox_titecrmmc,data_reg[,3])
    ntite <- c(ntite,ntitetr)
  }
}

library(xtable)

resultsTfixempscen1 <- list(MTD_titecrmmc = MTD_titecrmmc, doses_titecrmmc = doses_titecrmmc, tox_titecrmmc = tox_titecrmmc,
                         ntite=ntite)
write.table(MTD_titecrmmc, file="scenario1fixemp.txt")
write.table(ntite, file="scenario1fixempntite.txt")


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

MTD_titecrmmc <- NULL
doses_titecrmmc <- NULL
tox_titecrmmc <- NULL
xtox_titecrmmc <- NULL
ntite <- NULL

for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  eval(parse(text = paste("tox <- data",tr,"$data_complete", sep="")))
  # eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance", sep="")))
  
  ################################## for pois
  eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance2", sep="")))
  tpois <- NULL
  for (poiss in 1:length(T_entrance) ){
    times <- cumsum(c(T_entrance[poiss]+1,rep(1,J-1)))
    tpois <- c(tpois,times)
  }
  tox <- data.frame(tox)
  tox$time <- tpois
  #########################################################################
  
  
  tox <- data.frame(tox)
  T_entrance <- c(T_entrance,Inf)
  MTD=-1
  ntitetr=0
  
  x <- rep(1,cohort)                     # first dose to the first cohort
  M = n/cohort                           # number of total cohorts
  
  # data_reg = data_complete[1:(cohort*J),1:4]
  # y <- tox[cbind(1:length(x),x)]  
  
  
  
  time = T_entrance[cohort+1]   # time of entrance of the first patient of second cohort
  n_fin = 0                     # number of patient with complete follow up
  
  data_reg_c <- NULL              # building the dataset for TITE: complete, with patients who ended the followup
  data_reg_i <- NULL              # during follow up
  
  for (l in 1:cohort){
    pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
    pos <- max(pos,0.1)
    if (pos == J) {
      follow = J
      grade = tox[l*J, (x[l]+3)]
      n_fin = l
      data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
    }
    else {
      follow = pos
      grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
      data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
    }
  }
  
  #dimnames(data_reg_c)[[2]] <- c("patient", "follow", "grade")
  #dimnames(data_reg_i)[[2]] <- c("patient", "follow", "grade")
  
  data_reg <- rbind(data_reg_c,data_reg_i)
  # check if we can apply TITECRM or not
  if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
      any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
      any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
    stage1 <- stage2 <- FALSE
  } else {
    if ( length(unique(data_reg[,3]))==1 ) {
      stage1 <- TRUE
      stage2 <- FALSE
    } else {
      stage1 <- FALSE
      stage2 <- TRUE
    } 
  }
  
  pat=cohort
  
  if (stage1) {
    for(i in 2:M) {
      
      if ( sum(data_reg[,3]==3)==length(x) ) {
        MTD=0
        break
      } else {
        if ( sum(data_reg[,3]==2)==length(x) ) {
          x <- c(x,rep(x[length(x)], cohort)) 
        } else {
          x <- c(x,rep(min(max(x)+1,K), cohort))
        }
      }
      
      # dose for the cohort
      time = T_entrance[i*cohort+1]                   # time of entrance of the first of the next cohort
      pat=pat+cohort
      
      data_reg_i <- NULL                               # reset the temporary dataset
      
      for (l in (n_fin+1):(i*cohort)){
        pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
        pos <- max(pos,0.1)
        if (pos == J) {
          follow = J
          grade = tox[l*J, (x[l]+3)]
          n_fin = l
          data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
        }
        else {
          follow = pos
          grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
          data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
        }
      }
      
      data_reg <- rbind(data_reg_c,data_reg_i)
      
      if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
          any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
          any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
        stage2 <- FALSE
        ntitetr = length(x)
        break
      } else {
        if ( (any(data_reg[,3]==2) & any(data_reg[,3]==3)) | 
             (any(data_reg[,3]==1) & any(data_reg[,3]==2)) |
             (any(data_reg[,3]==1) & any(data_reg[,3]==3)) ) {
          stage2 <- TRUE
          break
        }
      }
    }
  }
  
  if (MTD==0) {
    MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
    doses_titecrmmc <- rbind(doses_titecrmmc,rep(0,40))
    tox_titecrmmc <- rbind(tox_titecrmmc,rep(1,40))
    ntite <- c(ntite,0)
  } else {
    
    for (i in seq(pat+1,n,1)){ 
      
      if (stage2) {  
        if (any(data_reg[,3]==3)) {
          grad = 3
          target_reg = target2
        } else {
          grad = 2
          target_reg = target1
        }
        # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
        y <- rep(0, dim(data_reg)[1]) 
        postox <- which(data_reg[,3]==grad )
        y[postox] <- rep(1, length(postox))
        
        results <- titecrm(prior, target_reg, y, x, followup=data_reg[,2], obswin=J)
        
        newdose <- min(results$mtd, max(x)+1, K)
        # check on the skipping dose
        x <- c(x,rep(newdose,cohort2))
        
        time = T_entrance[i*cohort2+1]  # first patient of the next cohort
        
        data_reg_i <- NULL  
        for (l in (n_fin+1):(i*cohort2)){
          pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
          pos <- max(pos,0.1)
          if (pos == J) {
            follow = J
            grade = tox[l*J, (x[l]+3)]
            n_fin = l
            data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
          }
          else {
            follow = pos
            grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
            data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
          }
        }
        
        data_reg <- rbind(data_reg_c,data_reg_i)
        
        if (any( data_reg[,3]==1) & any(data_reg[which(data_reg[,3]==1),2]==J) & 
            any(data_reg[,3]==2) & any(data_reg[which(data_reg[,3]==2),2]==J) & 
            any(data_reg[,3]==3) &  any(data_reg[which(data_reg[,3]==3),2]==J)) {
          stage2 <- FALSE
          ntitetr = length(x)
        }
        
      } else {
        ##### here stage 3
        
        
        results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
        
        newdose <- min(results$newdose, max(x)+1, K)
        # check on the skipping dose
        x <- c(x,rep(newdose,cohort2))
        
        time = T_entrance[i*cohort2+1]  # first patient of the next cohort
        
        data_reg_i <- NULL  
        for (l in (n_fin+1):(i*cohort2)){
          pos <- which( tox$time[((l-1)*J+1):(l*J)] <= time )
          pos <- max(pos,0.1)
          if (pos == J) {
            follow = J
            grade = tox[l*J, (x[l]+3)]
            n_fin = l
            data_reg_c <- rbind( data_reg_c, c(l,follow,grade))
          }
          else {
            follow = pos
            grade = tox[((l-1)*J+ceiling(pos)), (x[l]+3)]
            data_reg_i <- rbind( data_reg_i, c(l,follow,grade))
          }
        }
        
        data_reg <- rbind(data_reg_c,data_reg_i)
        
        
      }
    }
    
    ##### analysis of complete data
    if (stage2) {
      if (any(data_reg[,3]==3)) {
        grad = 3
        target_reg = target2
      } else {
        grad = 2
        target_reg = target1
      }
      # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
      y <- rep(0, dim(data_reg)[1]) 
      postox <- which(data_reg[,3]==grad )
      y[postox] <- rep(1, length(postox))
      
      results <- titecrm(prior, target_reg, y, x, followup=data_reg[,2], obswin=J)
      MTD <- min(results$mtd, max(x)+1, K)
      ntitetr=0
    } else {
      results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
      MTD <- min(results$newdose, max(x)+1, K)
    }
    
    MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
    doses_titecrmmc <- rbind(doses_titecrmmc,x)
    tox_titecrmmc <- rbind(tox_titecrmmc,data_reg[,3])
    ntite <- c(ntite,ntitetr)
  }
}

library(xtable)

resultsTpoisempscen1 <- list(MTD_titecrmmc = MTD_titecrmmc, doses_titecrmmc = doses_titecrmmc, tox_titecrmmc = tox_titecrmmc,
                         ntite=ntite)
write.table(MTD_titecrmmc, file="scenario1poisemp.txt")
write.table(ntite, file="scenario1poisempntite.txt")


save.image()









