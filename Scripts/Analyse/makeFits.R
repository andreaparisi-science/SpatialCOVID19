dat <- read.table("Fits/output.dat", header=TRUE)
#edit(dat)


timer.max <- max(dat$Timer)
plot(dat[dat$Timer == timer.max,])

ndata <- length(colnames(dat))-2
nrows <- floor((ndata-1) / 3)+1
par(mfrow=c(nrows, 3))
for (kk in seq(1,ndata))  {
	plot(density(dat[,kk+2]), main=colnames(dat)[kk+2])
}
