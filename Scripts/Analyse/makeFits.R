args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2)  {
	print( "Rscript makeFits  COUNTRY   STARTDATE")
	print( "    - COUNTRY: Country to be analized")
	print( "    - STARTDATE: the date corresponding to day zero")
	print( "                 for instance '2020-03-16'")
	quit('no')
}

previous <- 0
ts.previous <- 0
country <- args[1]
startDate <- args[2]
if (length(args) > 2)  {
        previous <- as.numeric(args[3])
		ts.previous <- previous
}
if (length(args) > 3)  {
        ts.previous <- as.numeric(args[4])
}

dat <- read.table("Fits/output.dat", header=TRUE)
#edit(dat)


timer.max <- max(dat$Timer) - previous
#plot(dat[dat$Timer == timer.max,])
dat <- dat[dat$Timer == timer.max,]
plot(dat[3:nrow(dat),3:ncol(dat)])

dat <- dat[,3:ncol(dat)]
ndata <- length(colnames(dat))
nrows <- floor((ndata-1) / 3)+1
par(mfrow=c(nrows, 3))
for (kk in seq(1,ndata))  {
	plot(density(dat[,kk]), main=colnames(dat)[kk])
}
par(mfrow=c(nrows, 3))
for (kk in seq(1,ndata))  {
	hist(dat[,kk], main=colnames(dat)[kk], n=20)
}

for (jj in 1:ncol(dat))  {
	min    = quantile(dat[,jj], 0.025)
	max    = quantile(dat[,jj], 0.975)
	median = median(dat[,jj])
	print( paste( "Param ",jj,":", median, "CI [",min,"-",max,"]; exp() -> ",exp(median), "CI [",exp(min),"-",exp(max),"]" ) )
}



#ts <- read.table( paste("../../Data/", country, "/", country, "_timeseries.dat", sep='') )
ts <- read.table( paste(country, "_timeseries.dat", sep=''), header=FALSE )
#print(ts)
directory <- paste("Fits/Generation-", (timer.max-ts.previous), sep='')
#fullname  <- paste("directory", "successfulRun-", sep='/')
fullnames <- list.files( path=directory, pattern='successfulRun*', full.names=TRUE )
kk <- 1
particles <- data.frame( t=as.numeric(), cases=as.numeric())
max.upper <- -1
min.lower <- 1000000
for (file in fullnames)  {
	readData <- read.table( file, header=TRUE, blank.lines.skip=FALSE )
	seps <- which(is.na(readData[,1]))
	upper <- seps-1
	lower <- as.vector( c(0,(seps+1))[1:(length(seps))] )
	for (jj in 1:length(upper))  {
		section <- readData[lower[jj]:upper[jj],]
		#particles <- rbind(particles, data.frame( t=c(section[,1]), cases=c(section[,2]) ))
		particles <- merge(particles, data.frame( t=c(section[,1]), cases=c(section[,2]) ), by = 't', all.y=TRUE)
		kk <- kk+1
	}
}
particles[is.na(particles)] <- 0
particles <- particles[1:(nrow(particles)-1),]
times    <- unlist(particles[,1]) + as.Date(startDate)
particles <- as.matrix(particles[,2:ncol(particles)])
p.min    <- apply(particles, 1, quantile, 0.025)
p.low    <- apply(particles, 1, quantile, 0.250)
p.high   <- apply(particles, 1, quantile, 0.750)
p.max    <- apply(particles, 1, quantile, 0.975)
p.median <- apply(particles, 1, median)
p.mean   <- apply(particles, 1, mean)

pdf("fitTimeseries.pdf", width=7, height=5)
size <- 1.6
par(mfrow=c(1,1))
x.min <- min(times)
x.max <- max(times)
y.max <- max(p.max)
plot( c(x.min, x.max), c(0, y.max), ty='n', cex.axis=1.6, cex.lab=1.6, cex.main=1.6, xlab="Time", ylab="Deaths" )
abline( v=as.Date(startDate), lty=2, col="black", lwd=1.6 )
polygon( c(times, rev(times)), c(p.min, rev(p.max)),  col=rgb(0.7, 0.7, 0.9) )
polygon( c(times, rev(times)), c(p.low, rev(p.high)), col=rgb(0.5, 0.5, 0.9) )
lines(times, p.median, type='l', col=rgb(0.2,0.2,1), lwd=1.6)
#lines(times, p.mean,   type='l', col=rgb(0,0,1))
lines(times, p.min, type='l', col=rgb(0.3,0.3,0.8))
lines(times, p.max, type='l', col=rgb(0.3,0.3,0.8))
lines(times, p.min, type='l', col=rgb(0.2,0.2,0.8))
lines(times, p.max, type='l', col=rgb(0.2,0.2,0.8))

times <- seq(1,length(ts[,2])) + as.Date(startDate)
lines(x=times, y=ts[,2], type='l', lty='dashed', col=rgb(1,0,0), lwd=1.6)
legend((x.max-x.min)*0.02+x.min, y.max*0.95, legend=c("Data", "Median", "50% PI", "95% PI"), lty=1, lwd=c(1.6, 1.6, 6.4, 6.4), col=c(rgb(1,0,0), rgb(0.2,0.2,1), rgb(0.5,0.5,0.9), rgb(0.7,0.7,0.9)), cex=1.2)

