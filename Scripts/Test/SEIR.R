#if not available
library(deSolve)

model <- function( t, x, params )  {
	with (as.list(x, params),  {
		ds <- -beta*ss*ii/NN
		de <- +beta*ss*ii/NN - sigma*ee
		di <- +sigma*ee - gamma*ii
		list( c(ds,de,di) )
	})
}


NN <- 47983657
R0 <- 2.5
sigma <- 1/3.0
gamma <- 1/4.0

beta <- R0*gamma
times <- seq( from = 0, to = 366, by = 0.1)
E0 <-  0
I0 <-  50
S0 <-  NN-E0-I0
params <- c( NN=NN, beta=beta, sigma=sigma, gamma=gamma )

res <- ode( c( ss=S0, ee=E0, ii=I0 ), times, model, params, method=rk4 )

dat <- read.table("DeryaSE-summary-Daily.dat", header=FALSE)
dat$Cases <- dat$V2 + dat$V3 + dat$V4 + dat$V5
dat$Cumul <- cumsum(dat$Cases)

plot(  res[,1], res[,4], type="l", col="red", lwd=2, ylim=c(0,NN), xlab="Time", ylab="Incidence in population" )
lines( res[,1], NN-res[,2], col="red", lt=2 )
lines( dat$V1, dat$Cases, col="blue", lt=1 )
lines( dat$V1, dat$Cumul, col="blue", lt=2 )
#lines( res[,1], res[,2], col=2, lwd=2, ylim=c(0,NN) )
legend( "topright", legend=c("I"), bg=rgb(1,1,1), lwd=2, col=c(3,2) )

