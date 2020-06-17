args = commandArgs(trailingOnly=TRUE)
country   = args[1]
startDate = args[2]

cases <- read.csv(paste("../../Data/Other/Private/JHSSE/time_series_covid19_confirmed_global_", startDate, ".csv", sep=''), header=TRUE)
linelist <- cases[cases$Country.Region == country,]
#linelist[56:length(linelist)] # From 13.03.2020, day 0
cases <- as.vector(unlist(linelist[56:length(linelist)]))

deaths <- read.csv(paste("../../Data/Other/Private/JHSSE/time_series_covid19_deaths_global_", startDate, ".csv", sep=''), header=TRUE)
linelist <- deaths[deaths$Country.Region == country,]
deaths <- as.vector(unlist(linelist[56:length(linelist)]))

dat <- data.frame( cases=cases, deaths=deaths )
write.table(dat, file=paste("../../Data/",country,"/", country, "_timeseries.dat", sep=''), row.names = FALSE, col.names = FALSE)

