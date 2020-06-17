import math
import sys

PROCESSES=12
MAXWEEK=18
MAXDATA = 1
#COLUMN  = 3  # 5% sampling of cases
COLUMN  = 1  # Cases
TYPE    = 1


if len(sys.argv) > 1:
	country = sys.argv[1]
	MAXDATA = int(sys.argv[2])
	if MAXDATA < 0:
		MAXDATA = -MAXDATA
		directory = "Fits/Replica-{:03d}/"
		TYPE = -1
	else:
		directory = "Run-{:03d}/"
		TYPE = 1

handler = open((directory+"DeryaSE-map-Weekly-000000-000.dat").format(1), "r")
lines = handler.readlines()
vec = lines[0].split(" ")
ncols = int(vec[0])
nrows = int(vec[1])
handler.close()



pop  = [ [-1 for x in range(0, nrows)] for y in range(0, ncols) ]
mmax = [ [-1 for x in range(0, nrows)] for y in range(0, ncols) ]
tmax = [ [-1 for x in range(0, nrows)] for y in range(0, ncols) ]
tmin = [ [-1 for x in range(0, nrows)] for y in range(0, ncols) ]
ids  = [ [-1 for x in range(0, nrows)] for y in range(0, ncols) ]


def mixcolours( col1, col2, x ):
	col = col1
	col[0] = (1-x)*col1[0] + x*col2[0]
	col[1] = (1-x)*col1[1] + x*col2[1]
	col[2] = (1-x)*col1[2] + x*col2[2]
	return  col


def convertToColour( val, method ):
	black = [0.0, 0.0, 0.0]
	red = [1.0,0.0,0.0]
	orange = [1.0,0.5,0.0]
	yellow = [1.0,1.0,0.0]
	green  = [0.0,1.0,0.0]
	cyan   = [0.0,1.0,1.0]
	blue   = [0.0,0.0,1.0]
	white  = [1.0,1.0,1.0]

	ret = red
	if method == 1:
		col = mixcolours(white, red, val)
	elif method == 2:
		if val < 0.25:
			col = mixcolours(red, orange, val/0.25)
		elif val < 0.5:
			col = mixcolours(orange, yellow, (val-0.25)/0.25)
		elif val < 0.75:
			col = mixcolours(yellow, green, (val-0.50)/0.25)
		else:
			col = mixcolours(green, blue, (val-0.75)/0.25)
	#	col = mixcolours(red, yellow, val)
	ret[0] = int(255*col[0])
	ret[1] = int(255*col[1])
	ret[2] = int(255*col[2])
	return ret


# Load counties ids
handler = open("../../Data/"+country+"/Maps/"+country+"_5km_ids.asc")
lines = handler.readlines()
lines = lines[6:]
maxid = -1
for yy in range(len(lines)):
	line = lines[yy].split()
	for xx in range(len(line)):
		ids[xx][nrows-yy-1] = int(line[xx])
		maxid = max(maxid, ids[xx][nrows-yy-1])





# Read frames and elaborate
handler = open("legend-tmax.ppm", "w")
handler.write("P3 ")
handler.write(str(20)+" "+str(nrows)+ " 255\n")
for yy in range(0,nrows):
	for xx in range(0, 20):
		cols = convertToColour( (1.0*yy)/nrows, 2 )
		handler.write( str(cols[0]) + " " + str(cols[1]) + " " + str(cols[2]) + " " )
handler.close()




countyStore = [ [ [0 for z in range(MAXWEEK+1)] for x in range(MAXDATA)] for y in range(maxid+1)]
for kk in range(0,MAXWEEK+1):
	print( "Handling frame [{}/{}]\r".format(kk, MAXWEEK) )
	mm = [ [-1.0 for x in range(0, nrows)] for y in range(0, ncols) ]
	for jj in range(1,(MAXDATA+1)):
		for qq in range(0,PROCESSES):
			if TYPE == -1:
				filename = (directory+"DeryaSE-map-Weekly-{:06d}-{:03d}.dat").format(jj-1, kk, qq)
			else:
				filename = (directory+"DeryaSE-map-Weekly-{:06d}-{:03d}.dat").format(jj, kk, qq)
			handler = open(filename, "r")
			lines = handler.readlines()
			for zz in range(2, len(lines)):
				data = lines[zz].split()
				xx = int(data[0])
				yy = int(data[1])
				pp = int(data[2])
				val = float(data[2+COLUMN])
				if pop[xx][yy] == -1:
					pop[xx][yy] = pp

				if (mm[xx][yy] == -1.0):
					mm[xx][yy] = 0.0
				mm[xx][yy] = mm[xx][yy] + val/(1.0*MAXDATA)
				ww = ids[xx][yy]
				if ww > 0:
					countyStore[ww][jj-1][kk] = countyStore[ww][jj-1][kk] + val
#					countyStore[1][1]=1

	#				print( str(xx)+" "+str(yy)+" "+str(val)+"\n" )
			handler.close()


	handler = open("out-{:03d}.ppm".format(kk), "w")
	handler.write("P3 ")
	handler.write(str(ncols)+" "+str(nrows)+ " 255\n")
	minVal = 0.0
	maxVal = 0.0
	for xx in range(0,ncols):
		for yy in range(0,nrows):
			maxVal = max(mm[xx][yy], maxVal)

	for yy in range(nrows-1,-1,-1):
		for xx in range(0,ncols):
			if (ids[xx][yy] == 0):
				handler.write( "127 255 255 ")
			elif (mm[xx][yy] == -1):
				handler.write( "255 255 255 ")
			else:
	#			level = int(255*mm[xx][yy]/(maxVal))
				if (maxVal > 0):
					level = math.log(mm[xx][yy]+1)/math.log(maxVal+1)
				else:
					level = 1
				cols = convertToColour(level,1)
	#			print(math.log(mm[xx][yy]+1)/math.log(maxVal+1))
	#			level = int(255*math.log(mm[xx][yy]/(maxVal)+1))
				handler.write( str(cols[0]) + " " + str(cols[1]) + " " + str(cols[2]) + " " )
			if mm[xx][yy] > mmax[xx][yy]:
				mmax[xx][yy] = mm[xx][yy]
				tmax[xx][yy] = kk
			if mm[xx][yy] > 0:
				if pop[xx][yy] > 0:
					if tmin[xx][yy] == -1:
						tmin[xx][yy] = kk
					else:
						tmin[xx][yy] = min(tmin[xx][yy], kk)
		handler.write( "\n" )
	handler.write( "\n" )



# Determine early extinctions
extinct = []
for run in range(0, MAXDATA):
	cases = 0
	for county in range(1, maxid+1):
		for week in range(0, MAXWEEK+1):
			cases = cases + countyStore[county][run][week]
	if cases < 10000000:
		extinct.append(run)




for yy in range(nrows-1,-1,-1):
	for xx in range(0,ncols):
		pass
		if pop[xx][yy] > 0:
			if tmin[xx][yy] == -1:
				tmin[xx][yy] = 0



maxVal = 0.0
minVal = 1000000
for xx in range(0,ncols):
	for yy in range(0,nrows):
		maxVal = max(tmax[xx][yy], maxVal)
		if tmax[xx][yy] > 0:
			minVal = min(tmax[xx][yy], minVal)
print("minVal [" + str(minVal) + "] - maxVal [" + str(maxVal) + "]")
minVal=6
maxVal=26 # Max number of weeks in plots
minVal=0
maxVal=18 # Max number of weeks in plots

handler = open("out-tmax.ppm".format(kk), "w")
handler.write("P3 ")
handler.write(str(ncols)+" "+str(nrows)+ " 255\n")
for yy in range(nrows-1,-1,-1):
	for xx in range(0,ncols):
		if (tmax[xx][yy] == -1 and ids[xx][yy] == 0):
			handler.write( "192 255 255 ")
		elif (tmax[xx][yy] == -1):
			handler.write( "0 0 0 " )
		elif (tmax[xx][yy] == 0):
			handler.write( "0 0 0 " )
		else:
			level = (1.0*tmax[xx][yy]-minVal)/(maxVal-minVal)
			if level > 1.0:
#				cols = [143, 0, 255]
				cols = [0, 0, 255]
			else:
				cols = convertToColour(level, 2)
#			level = int(255*math.log(tmax[xx][yy]+1)/math.log(maxVal+1))
#			print(math.log(mm[xx][yy]+1)/math.log(maxVal+1))
#			level = int(255*math.log(mm[xx][yy]/(maxVal)+1))
			handler.write( str(cols[0]) + " " + str(cols[1]) + " " + str(cols[2]) + " " )
	handler.write( "\n" )
handler.write( "\n" )


countyTime = [ -1 for ii in range(maxid+1)]
countyNums = [  0 for ii in range(maxid+1)]
countyMaxt = [ [ -1 for jj in range(MAXDATA)] for ii in range(maxid+1)]
countyMaxw = [ [ -1 for jj in range(MAXDATA)] for ii in range(maxid+1)]
handler = open("out-tmax2.ppm".format(kk), "w")
handler.write("P3 ")
handler.write(str(ncols)+" "+str(nrows)+ " 255\n")
for jj in range(maxid+1):
	for qq in range(0, MAXDATA):
		for zz in range(0, MAXWEEK+1):
			if countyMaxw[jj][qq] < countyStore[jj][qq][zz]:
				countyMaxw[jj][qq] = countyStore[jj][qq][zz]
				countyMaxt[jj][qq] = zz

for jj in range(maxid+1):
	for qq in range(0, MAXDATA):
		if countyMaxt[jj][qq] > 0:
			countyTime[jj] = countyTime[jj] + countyMaxt[jj][qq]*1.0
			countyNums[jj] = countyNums[jj] + 1

for jj in range(maxid+1):
	if countyNums[jj] > 0:
		countyTime[jj] = countyTime[jj]*1.0/countyNums[jj]


#for yy in range(nrows-1,-1,-1):
#	for xx in range(0,ncols):
#		jj = ids[xx][yy]
#		if jj > 0:
#			if tmax[xx][yy] > 0:
#				countyTime[jj] = min( tmax[xx][yy], countyTime[jj] )

for yy in range(nrows-1,-1,-1):
	for xx in range(0,ncols):
		jj = ids[xx][yy]
		val = countyTime[jj]
#		print(val)
		if (val == 1000000 or val == -1):
			handler.write( "192 255 255 ")
		elif (val == 0):
			handler.write( "0 0 0 " )
		else:
			level = (1.0*val-minVal)/(maxVal-minVal)
			if level > 1.0:
#				cols = [143, 0, 255]
				cols = [0, 0, 255]
			else:
				cols = convertToColour(level, 2)
#			level = int(255*math.log(tmax[xx][yy]+1)/math.log(maxVal+1))
#			print(math.log(mm[xx][yy]+1)/math.log(maxVal+1))
#			level = int(255*math.log(mm[xx][yy]/(maxVal)+1))
			handler.write( str(cols[0]) + " " + str(cols[1]) + " " + str(cols[2]) + " " )
	handler.write( "\n" )
handler.write( "\n" )



for jj in range(maxid+1):
	handler = open("county-ts-{:02d}.dat".format(jj), "w")
	t0 = 0
	for zz in range(0, MAXWEEK+1):
		val = 0
		if countyNums[jj] > 0:
			for qq in range(0, MAXDATA):
				if len(extinct) > 0 and qq == extinct[0]:
					extinct = extinct[1:]
					continue

				val = val + countyStore[jj][qq][zz]*1.0/countyNums[jj]
		handler.write(str(t0) + " " + str(val) + "\n")
		t0 = t0+7









maxVal = 0.0
minVal = 1000000
for xx in range(0,ncols):
	for yy in range(0,nrows):
		maxVal = max(tmin[xx][yy], maxVal)
		if tmin[xx][yy] > 0:
			minVal = min(tmin[xx][yy], minVal)
print("minVal [" + str(minVal) + "] - maxVal [" + str(maxVal) + "]")
#minVal= 

maxVal=26 # Max number of weeks in plots
maxVal=18 # Max number of weeks in plots
handler = open("out-tmin.ppm".format(kk), "w")
handler.write("P3 ")
handler.write(str(ncols)+" "+str(nrows)+ " 255\n")
for yy in range(nrows-1,-1,-1):
	for xx in range(0,ncols):
		if (tmin[xx][yy] == -1 and ids[xx][yy] == 0):
			handler.write( "192 255 255 ")
		elif (tmin[xx][yy] == -1):
			handler.write( "0 0 0 " )
		elif (tmin[xx][yy] == 0):
			handler.write( "0 0 0 " )
		else:
			level = (1.0*tmin[xx][yy]-minVal)/(maxVal-minVal)
			if level > 1.0:
#				cols = [143, 0, 255]
				cols = [0, 0, 255]
			else:
				cols = convertToColour(level, 2)
#			level = int(255*math.log(tmax[xx][yy]+1)/math.log(maxVal+1))
#			print(math.log(mm[xx][yy]+1)/math.log(maxVal+1))
#			level = int(255*math.log(mm[xx][yy]/(maxVal)+1))
			handler.write( str(cols[0]) + " " + str(cols[1]) + " " + str(cols[2]) + " " )
	handler.write( "\n" )
handler.write( "\n" )



