library("foreign")
# READ the Stata file, as this is easier to be read properly (NA vs spaces)
dat <- data.frame()
pb <- txtProgressBar(0, 8, style=3)
kk <- 0
for (year in 2018:2019)  {
	for (jj in 1:4)  {
		setTxtProgressBar( pb, kk )
	#	name = paste("STATA/EPA_2019T", jj, ".dta", sep='')
	#	contents <- read.dta(name)
		name = paste("CSV/EPA_", year, "T", jj, ".csv", sep='')
		contents <- read.table(name, header=TRUE, stringsAsFactors=FALSE)
		dat <- rbind(dat, contents)
		kk <- kk+1
	}
}
print( "" )
#edit(dat)
print(nrow(dat))
print(colnames(dat))
dat[ is.na(dat) ] <- 0
work <- dat[ dat$TRAREM == 1, ]
work$Hours <- floor(work$HORASE / 100) * 100
work$Fract <- floor(work$HORASE %% 100) * 100 / 60
print(sum((work$Hours+work$Fract)*work$FACTOREL)/sum(work$FACTOREL)/100.0)

#edit(dat$HORASE)
#sys.exit(0)


dat$PROV <- as.integer(dat$PROV)
dat$PROEST <- as.integer(dat$PROEST)
#provcode <- sort(c(unique(dat$PROV), unique(dat$PROEST)))
#nProvinces <- provcode[length(provcode)]
#print(nProvinces)

#popdata <- read.table("../../Counties/stats.dat", header=TRUE, stringsAsFactors=FALSE)
#popdata <- popdata[ order(popdata$SortId), ]


provnames <- read.csv("dr_EPA_2005.csv", header=TRUE)
provnames$Void <- NULL
colnames(provnames) <- c("Code", "Name")
provnames$Name <- gsub("[^A-Za-z]", "_", provnames$Name)
provnames$target <- c(
		5,	#01,Araba/Álava,
		3,	#02,Albacete,
		2,	#03,Alicante/Alacant,
		4,	#04,Almería,
		48,	#05,Ávila,
		7,	#06,Badajoz,
		24,	#07,"Balears, Illes",
		8,	#08,Barcelona,
		10,	#09,Burgos,
		11,	#10,Cáceres,
		12,	#11,Cádiz,
		15,	#12,Castellón /Castelló,
		16,	#13,Ciudad Real,
		13,	#14,Córdoba,
		1,	#15,"Coruña, A",
		17,	#16,Cuenca,
		19,	#17,Girona,
		20,	#18,Granada,
		21,	#19,Guadalajara,
		18,	#20,Gipuzkoa,
		22,	#21,Huelva,
		23,	#22,Huesca,
		25,	#23,Jaén,
		27,	#24,León,
		28,	#25,Lleida,
		26,	#26,"Rioja, La",
		29,	#27,Lugo,
		31,	#28,Madrid,
		30,	#29,Málaga,
		32,	#30,Murcia,
		33,	#31,Navarra,
		34,	#32,Ourense,
		6,	#33,Asturias,
		35,	#34,Palencia,
		-1,	#35,"Palmas, Las",
		36,	#36,Pontevedra,
		37,	#37,Salamanca,
		-1,	#38,Santa Cruz de Tenerife,
		14,	#39,Cantabria,
		38,	#40,Segovia,
		39,	#41,Sevilla,
		40,	#42,Soria,
		41,	#43,Tarragona,
		42,	#44,Teruel,
		43,	#45,Toledo,
		44,	#46,Valencia/València,
		45,	#47,Valladolid,
		9,	#48,Bizkaia,
		46,	#49,Zamora,
		47,	#50,Zaragoza,
		-1,	#51,Ceuta,
		-1	#52,Melilla,
)
#print(provnames)
#print(popdata)
nProvinces <- rev(sort(provnames$target))[1]

comm <- matrix( 0, nrow=nProvinces, ncol=nProvinces )
base <- rep(0, nProvinces)

pb <- txtProgressBar(0, nrow(dat), style=3)
for (jj in 1:nrow(dat))  {
	setTxtProgressBar( pb, jj )
	elem <- dat[jj,]

	if (elem$PROV == 0)  next
	orig = provnames[elem$PROV,]$target
	if (orig == -1) next
	# We need to count all contributions
	base[orig] <- base[orig] + elem$FACTOREL
#	base[orig] <- base[orig] + 1

	# Exclude only from commuting data
	if (elem$PROEST == 0)  next
	dest = provnames[elem$PROEST,]$target
	if (dest == -1) next
#	print(elem$PROV)
#	print(orig)
#	print(elem$PROEST)
#	print(dest)
#	print("---")

	comm[orig, dest] <- comm[orig, dest] + elem$FACTOREL
#	comm[orig, dest] <- comm[orig, dest] + 1
}
#print(comm)
#print(popdata)

for (orig in 1:nProvinces)  {
	for (dest in 1:nProvinces)  {
		if (base[orig] > 0)  {
			comm[orig, dest] <- comm[orig, dest] / base[orig]
		}
	}
}
print("")

print(comm)
write.table(comm, file="Commuting.dat", row.names = FALSE, col.names=FALSE)


