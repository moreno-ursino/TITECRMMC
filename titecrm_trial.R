#### TTICRMCM with TITE as medium stage


library(dfcrm)
# J cycles 
# K doses
# n sample size


#invlogit <- function(x) exp(x)/(1-exp(x))

cohort = 1
cohort2 = 1

prior = skel1 #c(0.05,0.1,0.2,0.3)
alpha0=1
#target1 = 0.6
target2 = 0.25


TR = Nsim

MTD_titecrm <- NULL
doses_titecrm <- NULL
tox_titecrm <- NULL
xtox_titecrm <- NULL


for (tr in 1:TR){
  # x doses assigned
  # y toxicities
  
  eval(parse(text = paste("tox <- data",tr,"$data_complete", sep="")))
  eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance", sep="")))
  
#   ################################## for pois
#   eval(parse(text = paste("T_entrance <- data",tr,"$T_entrance2", sep="")))
#   tpois <- NULL
#   for (poiss in 1:length(T_entrance) ){
#     times <- cumsum(c(T_entrance[poiss]+1,rep(1,J-1)))
#     tpois <- c(tpois,times)
#   }
#   tox <- data.frame(tox)
#   tox$time <- tpois
#   #########################################################################
  
  
  tox <- data.frame(tox)
  T_entrance <- c(T_entrance,Inf)
  MTD=-1

  
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

  # check if we can apply TITECRM or not
  data_reg <- rbind(data_reg_c,data_reg_i)
  # check if we can apply TITECRM or not
  if (any(data_reg[,3]==3)) {
    stage1=FALSE
  } else stage1=FALSE # put true if you want a two stage design
  
  
  pat=cohort
  
  if (stage1) {
  for(i in 2:M) {
      x <- c(x,rep(min(max(x)+1,K), cohort))

      
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
      
      if (any(data_reg[,3]==3)) {
        stage1=FALSE
        break
      }
    }
  }  
      for (i in seq(pat+1,n,1)){ 
        
          # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
          y <- rep(0, dim(data_reg)[1]) 
          postox <- which(data_reg[,3]==3 )
          y[postox] <- rep(1, length(postox))
          
          results <- titecrm(prior, target2, y, x, followup=data_reg[,2], obswin=J)
          
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
        
    }
  
  y <- rep(0, dim(data_reg)[1]) 
  postox <- which(data_reg[,3]==3 )
  y[postox] <- rep(1, length(postox))
  results <- titecrm(prior, target2, y, x, followup=data_reg[,2], obswin=J)  
  MTD <- min(results$mtd, max(x)+1, K)
  
  MTD_titecrm <- c(MTD_titecrm,MTD)
  doses_titecrm <- rbind(doses_titecrm,x)
  tox_titecrm <- rbind(tox_titecrm,y)
}

library(xtable)

resultsTfixedtite <- list(MTD_titecrm = MTD_titecrm, doses_titecrm = doses_titecrm, tox_titecrm = tox_titecrm)
write.table(MTD_titecrmmc, file="scenario1fixedtite.txt")

save.image()

####################################################################################################


MTD_titecrm <- NULL
doses_titecrm <- NULL
tox_titecrm <- NULL
xtox_titecrm <- NULL


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
  
  # check if we can apply TITECRM or not
  data_reg <- rbind(data_reg_c,data_reg_i)
  # check if we can apply TITECRM or not
  if (any(data_reg[,3]==3)) {
    stage1=FALSE
  } else stage1=FALSE # put true if you want a two stage design
  
  
  pat=cohort
  
  if (stage1) {
    for(i in 2:M) {
      x <- c(x,rep(min(max(x)+1,K), cohort))
      
      
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
      
      if (any(data_reg[,3]==3)) {
        stage1=FALSE
        break
      }
    }
  }  
  for (i in seq(pat+1,n,1)){ 
    
    # data_reg <- rbind(data_reg_c,data_reg_i)    # preparation of data for regression
    y <- rep(0, dim(data_reg)[1]) 
    postox <- which(data_reg[,3]==3 )
    y[postox] <- rep(1, length(postox))
    
    results <- titecrm(prior, target2, y, x, followup=data_reg[,2], obswin=J)
    
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
    
  }
  
  y <- rep(0, dim(data_reg)[1]) 
  postox <- which(data_reg[,3]==3 )
  y[postox] <- rep(1, length(postox))
  results <- titecrm(prior, target2, y, x, followup=data_reg[,2], obswin=J)  
  MTD <- min(results$mtd, max(x)+1, K)
  
  MTD_titecrm <- c(MTD_titecrm,MTD)
  doses_titecrm <- rbind(doses_titecrm,x)
  tox_titecrm <- rbind(tox_titecrm,y)
}

library(xtable)

resultsTpoistite <- list(MTD_titecrm = MTD_titecrm, doses_titecrm = doses_titecrm, tox_titecrm = tox_titecrm)
write.table(MTD_titecrmmc, file="scenario1poistite.txt")

save.image()

