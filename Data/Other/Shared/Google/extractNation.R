args <- commandArgs(trailingOnly=TRUE)
country <- args[1]
if (country == "UK")  {
	country = "GB"
}

dat <- read.csv("Global_Mobility_Report.csv", header=TRUE, stringsAsFactors=FALSE, na.strings='')
#print(unique(dat$country_region_code))
dat <- dat[ dat$country_region_code == country, ]
dat <- dat[ is.na(dat$sub_region_1), ]
print( paste('[', unique(dat$country_region), ']', sep='') )
dat$date <- as.Date(dat$date, "%Y-%m-%d")
dat <- data.frame( date = dat$date, 
	workplaces = dat$workplaces_percent_change_from_baseline, 
	retail_and_recreation = dat$retail_and_recreation_percent_change_from_baseline, 
	transit=dat$transit_stations_percent_change_from_baseline
)

write.table( dat, file=paste(args[1], '_google_mobility.dat', sep=''), row.names=FALSE ) 
