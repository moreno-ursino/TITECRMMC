

cohort = 1
cohort2 = 1

prior = skel1 #c(0.05,0.1,0.2,0.3)
alpha0=1
target1 = 0.5
target2 = 0.25


TR = 2000



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

library(rstan)
sm <- stan_model(file="probit_stan2.stan", model_name='prob2',verbose=FALSE)


titecrmmc <- function(x,doselabel,y,follow,alpha0,Tmax,target1,target2){
  # x dose level
  # doselabel    numbers used in regression
  # y grade
  # weight
  # alpha0 intercept
  
  N=length(x)
  K=length(doselabel)
  #x1p <- doselabel[x]
  w1p <- w2p <- follow/Tmax
  w1p[y==2] <- rep(1, sum(y==2))
  w2p[y==3] <- rep(1, sum(y==3))
  
  data_p = list( a=alpha0, N=N, K=K, x=x, y=y, w1=w1p, w2=w2p, d=doselabel)
  
  fit <- sampling(sm, data=data_p ,iter=4000, chains=4, control = list(adapt_delta = 0.8))
  #coeff=get_posterior_mean(fit)
  
  est <- extract(fit, pars="theta")
  theta <- apply(est$theta,2,median)
  theta <- min(theta)
  
  est2 <- extract(fit, pars=c("beta", "gam"))
  bet <- exp(median(est2$beta))
  gam <- exp(median(est2$gam))
  
  p1tox<-pnorm(alpha0+bet*doselabel)
  p2tox<-pnorm(alpha0+bet*doselabel-gam)
  #cur1<-which(abs(est[3:(2+K),5]-target1)==min(abs(est[3:(2+K),5]-target1)))
  #cur2<-which(abs(est[(3+K):(2+2*K),5]-target2)==min(abs(est[(3+K):(2+2*K),5]-target2)))
  cur<-order(abs(doselabel-theta))[1]
  
  list(newdose=cur, p1tox=p1tox, p2tox=p2tox)
  
}






MTD_titecrmmc <- NULL
MTD_titecrmmc24 <- NULL
MTD_titecrmmc32 <- NULL
doses_titecrmmc <- NULL
tox_titecrmmc <- NULL
xtox_titecrmmc <- NULL


for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  eval(parse(text = paste("tox <- data",tr,"$data_complete", sep="")))
  #eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance", sep="")))
  
  #   ################################## for pois
    eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance2", sep="")))
    tpois <- NULL
    for (poiss in 1:length(T_entrance) ){
      times <- cumsum(c(T_entrance[poiss]+1,rep(1,J-1)))
      tpois <- c(tpois,times)
    }
    tox <- data.frame(tox)
    tox$time <- tpois
  #   #########################################################################
  
  
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
  
  
  for(i in 2:M) {
    if (length(x)==1){
      results <- titecrmmc(array(x, dim = 1),doselabel,y=array(data_reg[,3], dim = 1),
                           follow=array(data_reg[,2], dim = 1),alpha0,Tmax=J,target1,target2)
    } else {
      
      results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
    }
    closeAllConnections()
    
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
  
  # complete data analysis
  results <- titecrmmc(x,doselabel,y=data_reg[,3],follow=data_reg[,2],alpha0,Tmax=J,target1,target2)
  MTD <- min(results$newdose, max(x)+1, K)
  
  MTD_titecrmmc <- c(MTD_titecrmmc,MTD)
  doses_titecrmmc <- rbind(doses_titecrmmc,x)
  tox_titecrmmc <- rbind(tox_titecrmmc,data_reg[,3])
  
  results <- titecrmmc(x[1:24],doselabel,y=data_reg[1:24,3],follow=rep(J,24),alpha0,Tmax=J,target1,target2)
  MTD <- min(results$newdose, max(x)+1, K)
  MTD_titecrmmc24 <- c(MTD_titecrmmc24, MTD)
  
  results <- titecrmmc(x[1:32],doselabel,y=data_reg[1:32,3],follow=rep(J,32),alpha0,Tmax=J,target1,target2)
  MTD <- min(results$newdose, max(x)+1, K)
  MTD_titecrmmc32 <- c(MTD_titecrmmc32, MTD)
}

library(xtable)


write.table(MTD_titecrmmc, file="scenario1poisprobbayes.txt")
write.table(MTD_titecrmmc24, file="scenario1poisprob24bayes.txt")
write.table(MTD_titecrmmc32, file="scenario1poisprob32bayes.txt")

write.table(doses_titecrmmc, file="scenario1poisprobbayes_doses.txt")
write.table(tox_titecrmmc, file="scenario1poisprobbayes_tox.txt")


#save.image()

