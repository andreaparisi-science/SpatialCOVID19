library(dplyr)
library(sf)
setwd("/home/parisia/Projects/Infect/Programs/CoV/PPP/Data/UK/Shared")


# Load data on population size of underlying gridded map.
stats <- read.table("../Maps/Counties/2019/stats.dat", header=TRUE, stringsAsFactors = FALSE)
stats <- select( stats, Name, SortId, Population )

# Names from UK_2019 data set. This is obtained from the 2019 LA files,
# hand edited to remove NI (currently not included)
england.names <- st_read("../Maps/Shapefiles/UK_2019_no_NIx.shp")
england.names <- data.frame( Local.Authority = england.names$LAD19NM, code = england.names$LAD19CD )
england.names <- england.names[ unlist(lapply(england.names$code, function(x) {return(substring(x,1,1))})) == "E", ]
# Get normalized names
england.names$Local.Authority <- gsub("[^A-Za-z]", "_", england.names$Local.Authority)

# List of cases in England from NHS website, by Local Authority
england.cases <- read.csv("coronavirus-cases_latest.csv", header=TRUE, stringsAsFactors = FALSE)
#edit(england.cases)
england.cases <- england.cases[ england.cases$Area.type == "Lower tier local authority", ]
england.cases$date <- as.Date( england.cases$Specimen.date, "%Y-%m-%d" )
england.cases$Area.name <- gsub("[^A-Za-z]", "_", england.cases$Area.name)
england.cases <- england.cases[ england.cases$date <= as.Date("2020-03-06", "%Y-%m-%d"), ]

# We group and sum as I thought there was some error in the data cumulative, but I am
# almost surely wrong. Anyway, the code is here now, so let us use it.
g.england.cases <- group_by( england.cases, Area.name )
s.england.cases <- summarise( g.england.cases, cumul.cases=sum(Daily.lab.confirmed.cases))
s.england.cases <- merge( s.england.cases, stats, by.x = "Area.name", by.y="Name", all.x=TRUE )
s.england.cases <- merge( s.england.cases, england.names, by.x = "Area.name", by.y="Local.Authority", all.x=TRUE )
s.england.cases <- select( s.england.cases, Local.Authority = Area.name, code, SortId, Population, 
				 cases=cumul.cases )
# Output is total cases per 100,000 individuals
s.england.cases$cases <- s.england.cases$cases


# List of cases in Wales. This is age stratified, thus we are forced to group and sum.
# NOTE: actually they are not age stratified.  What was going on when I first did this?
wales.cases <- read.csv("Wales_Rapid COVID-19 surveillance data.csv", header=TRUE, stringsAsFactors = FALSE)
wales.cases$date <- as.Date(wales.cases$Specimen.date, "%d/%m/%Y")
wales.cases <- wales.cases[ wales.cases$date <= as.Date("2020-03-06", "%Y-%m-%d"), ]
# We sum up all cases up to the date, as the cumulative misses data for the given date
g.wales.cases <- group_by(wales.cases, Local.Authority )
s.wales.cases <- summarise(g.wales.cases, cases=sum(Cases..new.) )
s.wales.cases <- s.wales.cases[ s.wales.cases$Local.Authority != "Outside Wales", ]
s.wales.cases <- s.wales.cases[ s.wales.cases$Local.Authority != "Unknown", ]
s.wales.cases[ is.na(s.wales.cases) ] <- 0
#edit(wales.cases)

# Names and codes should be available from shapefiles, but I had already implemented this way:
# if it works, do not break it!
wales.names <- list( 
	list("Blaenau Gwent", "W06000019"),
	list("Bridgend", "W06000013"),
	list("Caerphilly", "W06000018"),
	list("Cardiff", "W06000015"),
	list("Carmarthenshire", "W06000010"),
	list("Ceredigion", "W06000008"),
	list("Conwy", "W06000003"),
	list("Denbighshire", "W06000004"),
	list("Flintshire", "W06000005"),
	list("Gwynedd", "W06000002"),
	list("Isle of Anglesey", "W06000001"),
	list("Merthyr Tydfil", "W06000024"),
	list("Monmouthshire", "W06000021"),
	list("Neath Port Talbot", "W06000012"),
	list("Newport", "W06000022"),
	list("Pembrokeshire", "W06000009"),
	list("Powys", "W06000023"),
	list("Rhondda Cynon Taf", "W06000016"),
	list("Swansea", "W06000011"),
	list("Torfaen", "W06000020"),
	list("Vale of Glamorgan", "W06000014"),
	list("Wrexham", "W06000006")
)
wales.names <- data.frame( name=sapply(wales.names,function(x) x[[1]]), code=sapply(wales.names,function(x) x[[2]]), stringsAsFactors = FALSE )

# Group and su,
s.wales.cases <- merge( s.wales.cases, wales.names, by.x = "Local.Authority", by.y = "name", all.x=TRUE )
#s.wales.cases[ s.wales.cases$Local.Authority == "Vale of Glamorgan", ]$Local.Authority <- "The Vale of Glamorgan"
s.wales.cases$Local.Authority <- gsub("[^A-Za-z]", "_", s.wales.cases$Local.Authority)
s.wales.cases <- merge( s.wales.cases, stats, by.x = "Local.Authority", by.y = "Name", all.x=TRUE )
s.wales.cases <- select( s.wales.cases, Local.Authority, code, SortId, Population, cases )
# Cases per 100,000 individuals
s.wales.cases$cases <- s.wales.cases$cases

# Merge to unique list
cases <- rbind(s.england.cases, s.wales.cases)

# Cases at 6-03-2020 for Scotland. Not age stratified and not Local Authority distributed
scotland.total <- 21

		
scotland.names <- st_read("../Maps/Shapefiles/UK_2019.shp")
scotland.names <- data.frame( Local.Authority = scotland.names$LAD19NM, code = scotland.names$LAD19CD )
scotland.names <- scotland.names[ unlist(lapply(scotland.names$code, function(x) {return(substring(x,1,1))})) == "S", ]
scotland.names$Local.Authority <- gsub("[^A-Za-z]", "_", scotland.names$Local.Authority)
scotland.cases <- merge(scotland.names, stats, by.x = "Local.Authority", by.y = "Name", all.x = TRUE)
scotland.cases$cases <- scotland.total*scotland.cases$Population/sum(scotland.cases$Population)

cases <- rbind(cases, scotland.cases)
final.cases <- merge( cases, stats, by.x = "Local.Authority", by.y = "Name", all = TRUE)
final.cases <- data.frame( SortId=final.cases$SortId.y, cases=final.cases$cases )
final.cases[ is.na(final.cases) ] <- 0
output.cases <- final.cases[ order(final.cases$SortId), ]
write.table(output.cases, file = "distrCases_byId.dat", col.names = FALSE, row.names = FALSE)



