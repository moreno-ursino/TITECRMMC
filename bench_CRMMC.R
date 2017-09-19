### Run simulation for benchmark for the CRMMC under the 11 scenarios - original code from Ken Cheung.


##  Number of simulation replicates
nsim = 10000  # 10000 replicates used in the manuscript

############ DO NOT MODIFY BELOW
############
## DESIGN SETUP
n=40
#n = 32
#n = 24
tL = 1
tH = 2
thresh = c(tL, tH)
targetL = 0.50
targetH = 0.25
target = c(targetL, targetH)
J = length(target)

PI1 = matrix(c(100,100,100,100,100,50,67,76,84,90,25,37,45,56,64)/100, byrow = F,ncol=3) #scen1
PI2 = matrix(c(100,100,100,100,100,35,50,67,76,84,13,25,37,45,56)/100, byrow = F,ncol=3) #scen2
PI3 = matrix(c(100,100,100,100,100,22,35,50,67,76,7,13,25,37,45)/100, byrow = F,ncol=3) #scen3
PI4 = matrix(c(100,100,100,100,100,11,22,35,50,67,1,7,13,25,37)/100, byrow = F,ncol=3) #scen4
PI5 = matrix(c(100,100,100,100,100,6,11,22,35,50,0,1,7,13,25)/100, byrow = F,ncol=3) #scen5
PI6 = matrix(c(100,100,100,100,100,17,31,39,50,67,8,14,24,38,55)/100, byrow = F,ncol=3) #scen6
PI7 = matrix(c(100,100,100,100,100,50,67,76,84,90,11,25,37,46,56)/100, byrow = F,ncol=3) #scen7
PI8 = matrix(c(100,100,100,100,100,50,67,76,84,90,7,11,25,37,46)/100, byrow = F,ncol=3) #scen8
PI9 = matrix(c(100,100,100,100,100,50,67,76,84,90,3,8,15,25,37)/100, byrow = F,ncol=3) #scen9
PI10 = matrix(c(100,100,100,100,100,35,50,67,76,84,7,11,25,37,46)/100, byrow = F,ncol=3) #scen10
PI11 = matrix(c(100,100,100,100,100,17,26,35,50,67,0,1,2,14,26)/100, byrow = F,ncol=3) #scen11


toprint <- NULL
for (scen in 1:11){
  eval(parse(text = paste("PI <- PI",scen, sep="")))
dimnames(PI)[[2]] <- c(0,1,2)
  
L = ncol(PI) - 1
K = nrow(PI)
RAT = PI
for (j in 2:(L+1)) {
	RAT[,j] = PI[,j] / PI[,(j-1)]
}
RAT = RAT[,(2:(L+1))]

x = as.numeric(colnames(PI))

set.seed(2022)
mtd = matrix(rep(NA,J*nsim),nrow=J)

for (r in 1:nsim) {  
	
	### Step 3 in Algorithm: get ycomp, the complete profile
	ycomp = matrix(rep(NA,n*K),nrow=n)
	for (i in 1:n) {
		u = runif(L)
		for (k in 1:K) {
			rat = RAT[k,]
			pos  = which(u > rat)
			if (length(pos)==0) ycomp[i,k] = x[L+1]
			else ycomp[i,k] = x[min(pos)]
		}
	}

	### Step 4 in Algorithm: Estimate PI (only at the two constraint threshold values)
	pihat = matrix(rep(NA,J*K),nrow=J)
	for (j in 1:J) {
		for (k in 1:K) {
			pihat[j,k] = length(which(ycomp[,k] >= thresh[j]))/n
		}
	}
	
	### Step 5 in Algorithm: Estimate target dose d(PI)
	# for (j in 1:J) {
	# 	if (all(pihat[j,] <= target[j])) mtd[j,r] = K
	# 	else if (all(pihat[j,] >= target[j])) mtd[j,r] = 1
	# 	else {
	# 		dist = abs(pihat[j,] - target[j])
	# 		lpos = max(which(pihat[j,] <= target[j]))
	# 		upos = min(which(pihat[j,] >= target[j]))
	# 		if (dist[lpos] <= dist[upos]) mtd[j,r] = lpos
	# 		else mtd[j,r] = upos
	# 	}
	# }
	
	for (j in 1:J) {
	  mtd[j,r] = min(which(abs(pihat[j,]-target[j])==min(abs(pihat[j,]-target[j]))))
	}
	
	
	
}
MTD = apply(mtd,2,min)
rec = rep(NA,K)
for (k in 1:K) { rec[k] = length(which(MTD==k))/nsim; }
obj = t(cbind(1:K,100*rec))
toprint <- rbind(toprint, round(obj[2,]))
}
rownames(toprint) = c("scenario 1", "scenario 2", "scenario 3","scenario 4","scenario 5","scenario 6",
                      "scenario 7","scenario 8","scenario 9","scenario 10","scenario 11")
colnames(toprint) = 1:K

library(xtable)
xtable(toprint)


# cat("Benchmark design for CRMMC\n")
# cat("nsim:",nsim,"\n")
# cat("n:", n,"\n")
# cat("Toxicity threshold:", thresh, "\n")
# cat("Target rate for threshold:",target,"\n")
# obj = t(cbind(1:K,100*rec))
# rownames(obj) = c("Dose level","Percent selected")
# print(round(obj))

