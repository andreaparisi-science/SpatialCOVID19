library(dplyr)
#setwd("/home/parisia/Projects/Infect/Data/Mobility/LegacyItaly/AnalysisCommutersItaly/MATRICE PENDOLARISMO 2011")
dat <- read.table("matrix_pendo2011_10112014.txt", header = FALSE, 
		   stringsAsFactors = FALSE, 
		   col.names = c("record.type", "residence.type", "province.from", "comune.from", "sex", "reason", "destination.type", "province.to", "comune.to", "state.to", "transport", "time.out", "time.long", "estimate", "count"),
		   colClasses = c("character","character","character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character"))
#Select the simplified records
dat <- dat[dat$record.type=="S",]
head(dat)
dat$count <- as.numeric(dat$count)
# Remove ppl commuting abroad
dat <- dat[which(dat$province.to != "000"),]
# Group by orig and target locations
gg <- group_by(dat, province.from, province.to)
comm <- summarise(gg, count=sum(count))

# Add name of provinces
loc <- read.table("provences.dat", header=TRUE)
loc <- loc[c(2,1,3,4,5)]
#sardinia <- c("cagliari","nuoro","sassari","oristano","carbonia","medio-campidano","ogliastra","olbia")
#loc <- loc[which(!(loc$provence %in% sardinia)),]
comm$province.to   = as.numeric(comm$province.to)
comm$province.from = as.numeric(comm$province.from)
comm <- merge(comm, select(loc, name.to=provence, code), by.x="province.to", by.y = "code")
head(comm)
comm <- merge(comm, select(loc, name.from=provence, code), by.x="province.from", by.y = "code")
head(comm)
comm<-comm[c(1,5,2,4,3)]
head(comm)


write.table(loc,  "locations.csv", sep=",", col.names=FALSE, row.names=FALSE)
write.table(comm, "commuting.csv", sep=",", col.names=FALSE, row.names=FALSE)

#edit(comm)
maxId <- max( max(comm[,1]), max(comm[,3]) )
comm <- comm[c(1,3,5)]
for (ii in 1:maxId)  {
	for (jj in 1:maxId)  {
		if (nrow(comm[ comm$province.from == ii & comm$province.to == jj, ]) == 0)  {
			comm <- rbind( comm, data.frame( province.from = ii, province.to = jj, count = 0 ) )
#			print( c(ii, jj, comm[ comm$province.from == ii & comm$province.to == jj, ]) )
		} else if (ii == jj)  {
			comm[ comm$province.from == ii & comm$province.to == jj, ]$count <- 0
#			comm <- rbind( comm, data.frame( province.from = ii, province.to = jj, count = 0 ) )
		}
		comm[ comm$province.from == ii & comm$province.to == jj, ]$count = comm[ comm$province.from == ii & comm$province.to == jj, ]$count / loc[ loc$code == ii,]$Population
	}
}
warnings()
comm <- comm[ order( comm$province.from, comm$province.to ), ]
matc <- matrix( comm$count, ncol = maxId, byrow=TRUE )

write.table(matc, "ItalyCommuting.csv", sep="\t", col.names=FALSE, row.names=FALSE)

