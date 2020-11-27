library(dplyr)
dat <- read.csv("CasesOrig.csv",header=TRUE)
names <- read.table("EditedCounties.dat",header=TRUE)
names <- names[c(2,4)]
ages <- data.frame(agecat="0-4", age=0, stringsAsFactors=FALSE)
ages <- rbind(ages, c("4-9", 0))
ages <- rbind(ages, c("10-14", 1))
ages <- rbind(ages, c("15-19", 1))
ages <- rbind(ages, c("20-24", 2))
ages <- rbind(ages, c("25-29", 2))
ages <- rbind(ages, c("30-34", 3))
ages <- rbind(ages, c("35-39", 3))
ages <- rbind(ages, c("40-44", 4))
ages <- rbind(ages, c("45-49", 4))
ages <- rbind(ages, c("50-54", 5))
ages <- rbind(ages, c("55-59", 5))
ages <- rbind(ages, c("60-64", 6))
ages <- rbind(ages, c("65-69", 6))
ages <- rbind(ages, c("70-74", 7))
ages <- rbind(ages, c("75-79", 7))
ages <- rbind(ages, c("80+", 8))
dat <- merge(dat, names, by.x="county_where_the_case_was_diagonised", by.y="Name", all=TRUE)
dat <- merge(dat, ages, by.x="agecat_model", by.y="agecat", all=TRUE)
edit(dat)
head(dat)
dat$age <- as.numeric(dat$age)
datOrdering <- order(dat$SortId, dat$age)
dat[datOrdering,]
dat.grouped <- group_by(dat, age)
dat.summarised <-summarise(dat.grouped, count = sum(n))
write.table(dat.summarised, "CaseByAge.dat", col.names = FALSE, row.names = FALSE)

dat1 <- read.csv("GlobalCases.csv", header=TRUE)
dat2 <- read.csv("GlobalDeaths.csv", header=TRUE)
vals1 <- dat1[dat1$Country.Region == "Kenya",]
vals1 <- as.vector(unlist(vals1[-(1:4)]))
vals1
vals2 <- dat2[dat2$Country.Region == "Kenya",]
vals2 <- as.vector(unlist(vals2[-(1:4)]))
vals2
dat <- data.frame(cases=vals1, deaths=vals2)[-(1:51),]
write.table(dat, file="Kenya_timeseries.dat", row.names = FALSE, col.names = FALSE)
dat
class(dat)
