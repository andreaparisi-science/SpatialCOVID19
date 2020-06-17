library(dplyr)

cases  <- read.csv( "../../Data/Other/Private/JHSSE/time_series_covid19_confirmed_global.csv", header=TRUE )
deaths <- read.csv( "../../Data/Other/Private/JHSSE/time_series_covid19_deaths_global.csv", header=TRUE )

ndays <- as.numeric( as.Date("2020-03-13") - as.Date("2020-01-22") )
cases <- unlist(cases[cases$Country.Region == "Kenya",][(5+ndays):ncol(cases)])
deaths <- unlist(deaths[deaths$Country.Region == "Kenya",][(5+ndays):ncol(deaths)])

out <- data.frame( cases=cases, deaths= deaths )
write.table( out, file="../../Data/Kenya/Kenya_timeseries_JHSSE.dat", col.names = FALSE, row.names=FALSE )


