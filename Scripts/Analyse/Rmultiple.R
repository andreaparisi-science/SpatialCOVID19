library(dplyr)

args = commandArgs(trailingOnly=TRUE)
startDate = as.numeric(args[1])

directory <- (list.files(".", pattern="Run-\\d{3}"))
directory <- sort(directory)
if (length(args) == 2)  {
	limit = as.numeric(args[2])
} else {
	limit = length(directory)
}

kk <- 1
cumul <- data.frame()
for (jj in 1:length(directory))  {
	if (jj > limit)  break
	locdata <- read.table(paste(directory[jj],"DeryaSE-summary-Daily.dat", sep="/"), header=FALSE)
	header <- c("Day", "Hosp.cases", "ICU.cases", "Home.cases", "Asympt.cases", "Recovs", "Deaths", "Cases.cumul", "Asympt.cumul", "Hosp.occupancy", "ICU.occupancy", "Home.occupancy", "Ctrl1", "Ctrl2", "Deaths.cumul.0", "Deaths.cumul.5", "Deaths.cumul.10", "Deaths.cumul.15", "Deaths.cumul.20", "Deaths.cumul.25", "Deaths.cumul.30", "Deaths.cumul.35", "Deaths.cumu.40", "Deaths.cumul.45", "Deaths.cumul.50", "Deaths.cumul.55", "Deaths.cumul.60", "Deaths.cumul.65", "Deaths.cumul.70", "Deaths.cumul.75", "Deaths.cumul.80", "Deaths.cumul")
	colnames(locdata) <- header
	if (tail(locdata,1)$Asympt.cumul < 1000)  {
		print( paste("Early extinction for [", directory[jj], "]", sep='') )
		next
	}
	print(directory[jj])

	locdata <- locdata[-1,]
	locdata <- locdata[locdata$Day > startDate,]
	locdata$Day <- locdata$Day-(startDate)
	if (kk == 1)  {
		fulldata <- locdata
		cumul <- data.frame( tail(locdata, 1)[15:31] )
	} else  {
		fulldata <- rbind(fulldata, locdata)
		cumul <- rbind(cumul, data.frame(tail(locdata, 1)[15:31]) )
	}
	locdata$Week = floor((locdata$Day-1)/7)+1
	glocdata <- group_by(locdata, Week)
	locdata <- summarise(glocdata, Hosp.occupancy=max(Hosp.occupancy), ICU.occupancy=max(ICU.occupancy))
	if (kk == 1)  {
		data <- locdata
		stats <- data.frame( Hosp.occupancy = max(locdata$Hosp.occupancy), ICU.occupancy = max(locdata$ICU.occupancy) )
	} else  {
		data <- rbind(data, locdata)
		stats <- rbind( stats, data.frame( Hosp.occupancy = max(locdata$Hosp.occupancy), ICU.occupancy = max(locdata$ICU.occupancy) ) )
	}
	kk <- kk+1
}

par(mar=c(5,5,4,3))


# OCCUPANCY
#
pdf("occupancyWeekly.pdf", width=8, height=5)
max.Hosp <- max( data$Hosp.occupancy )
max.ICUs <- max( data$ICU.occupancy )
max.vals <- max( max.Hosp, max.ICUs )

boxplot( Hosp.occupancy ~ Week, data = data, cex=0.7, col="darkolivegreen3", border="darkblue", cex.axis=1.2, cex.lab=1.2, xlab="Day", ylab="Counts", outline=FALSE,
	  main="Maximal weekly bed occupancy", cex.main=1.4)
boxplot( ICU.occupancy ~ Week, data=data, col="orange", border="darkred", outline=FALSE, xaxt="n", yaxt="n", add=TRUE)
legend(33, 0.9*max.vals, c("Hospital occupancy", "ICU occupancy"), fill=c("darkolivegreen4", "orange"), border=c("darkblue","darkred"), cex=1.2)



# AGE DISTRIBUTION OF DEATHS
#
#cumul <- unlist(tail(fulldata, 1)[15:31])
cumul  <- as.matrix(cumul)
mean   <- apply(cumul, 2, mean)
stddev <- apply(cumul, 2, sd)
pdf("deathsByAge.pdf", width=8, height=5)
max.y  <- max(mean+2.05*stddev)  # This avoid thinning of highest error bar
centers <- barplot(mean, names.arg=c("0-4","5-9","10-14","15-19","20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+"),
			main="Expected number of deaths by age", cex.axis = 1.2, cex.lab=1.2, cex.main=1.4, ylim = c(0,max.y))
#segments(centers, mean-1.96*stddev, centers, mean+1.96*stddev, lwd=1.5)
arrows(centers, mean, centers, mean-1.96*stddev, angle=90, lwd=1.5, length=0.08)
arrows(centers, mean, centers, mean+1.96*stddev, angle=90, lwd=1.5, length=0.08)
write( c(sum(mean), sum(mean-1.96*stddev), sum(mean+1.96*stddev)), file="deathCount.pdf" )

# DAILY INCIDENCE
#
fulldata$cases.val <- fulldata$ICU.cases + fulldata$Hosp.cases + fulldata$Home.cases
gfulldata <- group_by(fulldata, Day)
fulldata <- summarise(gfulldata, cases.max=quantile(cases.val, prob=0.975), cases.min=quantile(cases.val, prob=0.025), cases.median=quantile(cases.val, prob=0.5))

pdf("incidenceDaily.pdf", width=8, height=5)

t.max = max(fulldata$Day)
y.max = max(fulldata$cases.max)
plot( c(0,t.max), c(1,y.max), type="n", log="", xlab="Time", ylab="Individuals", main="Daily incidence" )
polygon( 
	c(fulldata$Day, rev(fulldata$Day)),
	c(fulldata$cases.max, rev(fulldata$cases.min)),
	border=NA,
	col=rgb(.8,.8,1)
)
lines(fulldata$Day, fulldata$cases.min, col=rgb(.6,.6,1), lwd=1)
lines(fulldata$Day, fulldata$cases.max, col=rgb(.6,.6,1), lwd=1)
lines(fulldata$Day, fulldata$cases.median, col="blue", lwd=1)
# lines( data$Day, data$ICU.occupancy, type="o", cex=0.7, col="red" )
# legend(220, 35000, c("Hospital occupancy", "ICU occupancy"), col=c("blue", "red"), lty=1, lwd=2, cex=1.2)
# 
# pdf("incidence.pdf", width=8, height=5)
# plot( data$Day, data$Hosp.case+data$ICU.cases+data$Home.cases, type="o", cex=0.7, col="black", cex.axis=1.2, cex.lab=1.2, xlab="Day", ylab="Incidence",
# 	main="Daily number of cases", cex.main=1.4)
# #legend(220, 20000, c("Total cases"), col=c("black"), lty=1, lwd=2, cex=1.2)
 


