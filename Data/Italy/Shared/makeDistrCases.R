library(dplyr)
dat <- read.table("Data_2020_05_01.dat", header=FALSE)
#maxVal <- max(dat$V3)
#dat$V3 <- dat$V3/maxVal
dat <- select(dat, V1, V3)
write.table(dat, file = "distrCases_byId.dat", col.names = FALSE, row.names = FALSE)

