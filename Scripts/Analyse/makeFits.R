dat <- read.table("Fits/output.dat", header=TRUE)
#edit(dat)


timer.max <- max(dat$Timer)
plot(dat[dat$Timer == timer.max,])

