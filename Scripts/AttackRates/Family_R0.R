#Country = "Kenya"
#Country = "Italy"
Country = "China"
#HHsize <- 3.9   # Kenya
#HHsize <- 2.40  # Italy
HHsize <- 2.85  # Wuhan
R0 <- 2.5
#expectedVal <- 0.158
expectedVal <- 0.30

DefaultDir = paste("/home/parisia/Projects/Infect/Programs/CoV/Public/Data/", Country, "/", sep='');

KKhome <- read.table( paste(DefaultDir, "Contacts/", Country, "ContactMatrix_home.csv", sep=''), header=FALSE )
KKwork <- read.table( paste(DefaultDir, "Contacts/", Country, "ContactMatrix_work.csv", sep=''), header=FALSE )
KKschool <- read.table( paste(DefaultDir, "Contacts/", Country, "ContactMatrix_school.csv", sep=''), header=FALSE )
KKother  <- read.table( paste(DefaultDir, "Contacts/", Country, "ContactMatrix_other.csv", sep=''), header=FALSE )

KK = KKwork + KKschool + KKother
KKhome = rbind(KKhome, KKhome[16,] )
KKhome = cbind(KKhome, KKhome[,16] )
KK = rbind(KK, KK[16,])
KK = cbind(KK, KK[,16])

NN <- read.table( paste(DefaultDir, "Setup/", Country, "_5km_stats.dat", sep=''), header=FALSE )
NN <- NN[,3]
nAgeGroups <- length(NN);


print(paste("Age groups:", nAgeGroups))
#print(NN)

FF <- KKhome
cols.sums <- apply(FF, 2, sum)
for (ii in 1:nAgeGroups)  {
	for (jj in 1:nAgeGroups)  {
		FF[ii,jj] <- FF[ii,jj]/cols.sums[jj]
	}
}

#print(KKhome)
# Now build LL households
LL <- 1000
MM <- 1000
lambdaCH <- 16.7448
tau <- 4
R0_wm <- 1
R0_HH <- 1
betamul <- 1

attackrateMul <- 1.0
attackrate <- 0.158
beta <- (R0/tau)/lambdaCH
betagen <- beta*betamul
for (zz in 1:10)  {
	if (expectedVal > 0)  {
		attackrateMul <- attackrateMul * expectedVal/attackrate
	}
	prob <- beta*betamul*attackrateMul
	print( paste("Attack rate multiplier:", attackrateMul) )
	print( paste("Beta multiplier:",betamul) )

	pb <- txtProgressBar(min = 0, max = MM, style = 3)
	result <- c()
	for (trial in 1:MM)  {
		size <- rep(0, LL)
		for (jj in 1:LL)  {
			while (size[jj] == 0)  {
				size[jj] = rpois(1, 3.9)
			}
		}

		counter = 0
		for (kk in 1:length(size))  {
			pop <- rep(0, size[kk])
			pop[size[kk]] <- 1
			for (ii in 1:tau)  {
				for (jj in 1:size[kk])  {
					if (pop[jj] == 0 && runif(1) < prob)  {
						pop[jj] <- 1
					}
				}
			}
			if (length(pop) > 1)  {
				counter <- counter + (sum(pop)-1)/(size[kk]-1)
			}
		}
		result <- c( result, counter/LL )
		setTxtProgressBar(pb, trial)
	}
	attackrate <- mean(result)
	print( paste("Household attack rate:", attackrate, '\u00B1', 1.96*sd(result)/sqrt(MM) ) )

	min.attackrate = attackrate - 1.96*sd(result)/sqrt(MM)
	max.attackrate = attackrate + 1.96*sd(result)/sqrt(MM)

	lambdaKKhome <- abs(eigen(KKhome)[[1]])[1]

	eigenValsFull <- eigen(KK + KKhome)[1]
	eigenValsHH   <- eigen(betamul*KK + (attackrate*HHsize/lambdaKKhome/tau) * KKhome)[1]
	min.eigenValsHH <- eigen(betamul*KK + (min.attackrate*HHsize/lambdaKKhome/tau) * KKhome)[1]
	max.eigenValsHH <- eigen(betamul*KK + (max.attackrate*HHsize/lambdaKKhome/tau) * KKhome)[1]
	R0_wm <- max(unlist(lapply(eigenValsFull, Mod)))
	R0_HH <- max(unlist(lapply(eigenValsHH, Mod)))
	R0_HH.min <- max(unlist(lapply(min.eigenValsHH, Mod)))
	R0_HH.max <- max(unlist(lapply(max.eigenValsHH, Mod)))
	print(paste("R0/Reff:",R0_wm/R0_HH, "C.I. [", R0_wm/R0_HH.max, ",", R0_wm/R0_HH.min, "]") )
	print(paste("betamul:",betamul*R0_wm/R0_HH, "C.I. [", betamul*R0_wm/R0_HH.max, ",", betamul*R0_wm/R0_HH.min, "]") )
	if (expectedVal != -1)  {
		print(paste("HH bmul:",attackrateMul * expectedVal/attackrate, "C.I. [", attackrateMul * expectedVal/max.attackrate, ",", attackrateMul * expectedVal/min.attackrate, "]") )
	}
	betamul <- betamul*R0_wm/R0_HH
}


