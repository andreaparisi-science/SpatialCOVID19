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

provnames <- read.csv("dr_EPA_2005.csv", header=TRUE)
provnames$Void <- NULL
colnames(provnames) <- c("Code", "Name")
provnames$Name <- gsub("[^A-Za-z]", "_", provnames$Name)
provnames$target <- c(5, 3, 2, 4, 7, 8, 25, 9, 11, # A, B
	12, 13, 15, 16, 17, 1, 18, 
	20, 21, 22, 19, 23, 24, 26,
	28, 29, 27, 30, 31, 32, 33,
	34, 35, 6, 36, -1, 37, 38,	# =>| Salamanca (missing Las Palmas)
	-1, 14, 39, 40, 41, 42, 43, # =>| Teruel (missing Santa Crus de Tenerife)
	44, 45, 46, 10, 47, 48, -1, -1) # =>| Ceuta, Melilla (excluded)
print(provnames)
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

#	comm[orig, dest] <- comm[orig, dest] + elem$FACTOREL
	comm[orig, dest] <- comm[orig, dest] + 1
}

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


