library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
COUNTRY <- args[1]
START_DATE <- args[2]

dat <- read.csv( "../../Data/Other/Private/Google/Global_Mobility_Report_2020_05_26.csv", header = TRUE )
dat$date <- as.Date(dat$date)
dat <- dat[dat$country_region == COUNTRY & dat$sub_region_1 == "" & dat$date >= as.Date(START_DATE),]
dat$day <- as.numeric(dat$date - as.Date(START_DATE))
dat$week <- floor(dat$day / 7)
dat <- summarise_all( group_by( dat, week ), mean )
dat$day <- dat$week*7
names <- colnames(dat)
names <- strsplit(names, "_")
names <- lapply(names, function(x) {return(x[1])})
colnames(dat) <- names

dat$date <- as.Date(START_DATE)+dat$week*7  # Overwrite with start-of-week day
print(dat[c(1,6:ncol(dat))], n=20, width=Inf)


