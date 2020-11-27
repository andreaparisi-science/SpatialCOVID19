args = commandArgs(trailingOnly=TRUE)
country   = args[1]
startDate = args[2]
if (length(args) > 2)  {
	countryWrite = args[3]
} else {
	countryWrite = country
}

cases <- read.csv(paste("../../Data/Other/Private/JHSSE/time_series_covid19_confirmed_global.csv", sep=''), header=TRUE)
linelist <- cases[cases$Country.Region == country & cases$Province.State == "",]
#linelist[56:length(linelist)] # From 13.03.2020, day 0
startDate <- as.Date(startDate, "%Y-%m-%d")
ndays = startDate - as.Date( "2020-01-22", "%Y-%m-%d" )
cases <- as.vector(unlist(linelist[(5+ndays):length(linelist)]))
edit(cases)
deaths <- read.csv(paste("../../Data/Other/Private/JHSSE/time_series_covid19_deaths_global.csv", sep=''), header=TRUE)
linelist <- deaths[deaths$Country.Region == country & deaths$Province.State == "",]
deaths <- as.vector(unlist(linelist[(5+ndays):length(linelist)]))

dat <- data.frame( cases=cases, deaths=deaths )
country <- gsub( "[^A-Za-z]", "_", country )
write.table(dat, file=paste("../../Data/",countryWrite,"/", countryWrite, "_timeseries.dat", sep=''), row.names = FALSE, col.names = FALSE)

