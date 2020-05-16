print( "#if not available\n" )
print( "library(deSolve)\n" )
print( "\n" )
print( "model <- function( t, x, params )  {\n" )
print( "\twith (as.list(x, params),  {\n" )

for xx in range(0,17):
	print( "\t\tff_" + str(xx) + " <- beta/eigen*(", end='' )
	for yy in range(0,17):
		if yy > 0:
			print( "+", end='' )
		print( "KK_" +str(xx) + "_" + str(yy) + "*ss_"+str(xx)+"*(ii_"+str(yy)+"+aa_"+str(yy)+")/NN_"+str(yy), end='' )
	print( ")" )

	print( "\t\tds_" + str(xx) + " <- -ff_" + str(xx) )
	print( "\t\tde_" + str(xx) + " <- +ff_" + str(xx) + " - sigma*ee_"+str(xx) )

	print( "\t\tdi_" + str(xx) + " <- +zz_"+str(xx)+"*sigma*ee_"+str(xx) + " - gamma*ii_" + str(xx) )
	print( "\t\tda_" + str(xx) + " <- +(1.0-zz_"+str(xx)+")*sigma*ee_"+str(xx) + " - gamma*aa_" + str(xx) )

print( "\t\tlist( c(", end='' )
for xx in range(0,17):
	if xx > 0:
		print( ", ", end='' )
	print( "ds_" + str(xx) + ", de_" + str(xx) + ", di_" + str(xx) + ", da_" + str(xx), end='' )
print( ") )" )
print( "\t})" )
print( "}" )

print( "dat <- read.table( '../../../Test/Contacts/KenyaContactMatrix.csv', header=FALSE )" )
print( "KK <- data.matrix(dat)" )
print( "KK <- cbind(KK, KK[,16])" )
print( "KK <- rbind(KK, KK[16,])" )

print( "zz = c( 0.0014432192790456596, 0.0009357685327922587, 0.0010761342759786922, 0.0013458086273749326, " )
print( "     0.013695268824959014, 0.016889531485125703, 0.03167807081661833, 0.024290277337257693, 0.028326981798033054, " )
print( "     0.040739451216136126, 0.07430476961279056, 0.13743204603117717, 0.17923146887955166, 0.3315015659109177, " )
print( "     0.2923192122329743, 0.375935907058367, 1.00)" )
for xx in range(0,17):
	print( "zz_"+str(xx)+" = zz[" + str(xx+1) + "]" )

print( "dat <- read.table( '../../../Test/Setup/Test_5000km_stats.dat', header=FALSE)" )
print( "NN <- 47983657" )
for xx in range(0,17):
	print("NN_"+str(xx)+" <- dat$V3["+str(xx+1)+"]" )
print( "R0 <- 2.5" )
print( "sigma <- 1/3.0" )
print( "gamma <- 1/4.0" )
print( "eigen <- 16.7448" )
print()
print( "beta <- R0*gamma" )
#print( "KK = array(1, c(17,17))" )
for xx in range(0,17):
	for yy in range(0,17):
		print( "KK_"+str(xx)+"_"+str(yy)+"=KK["+str(xx+1)+", "+str(yy+1)+"]" )

print( "times <- seq( from = 0, to = 366, by = 0.1)" )
for xx in range(0,17):
	print( "S"+ str(xx) + " <- NN_"+str(xx)+"-2" )
	print( "E"+ str(xx) + " <- 0" )
	print( "I"+ str(xx) + " <- 2" )
	print( "A"+ str(xx) + " <- 0" )

print( "params <- c( ", end='' )
for xx in range(0,17):
	if xx > 0:
		print( ", ", end='' )
	print( "NN_" + str(xx) + "=NN_"+str(xx), end='' )
print( ", beta=beta, sigma=sigma, gamma=gamma, eigen=eigen", end='' )
for xx in range(0,17):
	for yy in range(0,17):
		print( ", KK_" + str(xx) + "_" + str(yy) + "=KK["+str(xx+1)+", "+str(yy+1)+"]", end='' )
for xx in range(0,17):
	print( ", zz_"+str(xx) + " = zz_"+str(xx), end='' )
print( ")" )

print( "res <- ode( c( ", end='' )
for xx in range(0,17):
	if xx > 0:
		print( ", ", end='' )
	print( "ss_" +str(xx) + "=S" +str(xx), end='' )
	print( ", ee_" +str(xx) + "=E" +str(xx), end='' )
	print( ", ii_" +str(xx) + "=I" +str(xx), end='' )
	print( ", aa_" +str(xx) + "=A" +str(xx), end='' )
	#res <- ode( c( ss=S0, ee=E0, ii=I0 ), times, model, params, method=rk4 )
print( "), times, model, params, method=rk4 )" )

print()
print( "tt <- res[,1]" )
print( "cases <- ", end='')
for xx in range(0,17):
	if xx > 0:
		print( "+", end='' )
	print( "res[," + str(4+4*xx) + "]+res[,"+str(5+4*xx)+"]", end='' )
print()
print( "cases <- cases[tt==floor(tt)]*gamma" )
print( "tt <- tt[tt==floor(tt)]+1" )
print( "cumul <- cumsum(cases)" )

print( "dat <- read.table(\"DeryaSE-summary-Daily.dat\", header=FALSE)" )
print( "dat$Cases <- dat$V2 + dat$V3 + dat$V4 + dat$V5" )
print( "dat$Cumul <- cumsum(dat$Cases)" )
print()
print( "plot(  tt, cases, type='l', col='red', lwd=2, ylim=c(0,NN), xlab='Time', ylab='Incidence in population' )" )
print( "lines( tt, cumul, col='red', lt=2 )" )
print( "lines( dat$V1, dat$Cases, col='blue', lt=1 )" )
print( "lines( dat$V1, dat$Cumul, col='blue', lt=2 )" )
#lines( res[,1], res[,2], col=2, lwd=2, ylim=c(0,NN) )
print( "legend( 'topright', legend=c('I'), bg=rgb(1,1,1), lwd=2, col=c(3,2) )" )

