import os.path
import sys

linelist = [140, 141, 142, 143, 144]
shorten =  [ 60,  60,  62,  63,  63]


if os.path.isfile( "UK_5km_base.asc" ):
	print("Operation already performed")
	sys.exit()

header = []
for el in ["_0", "_1", "_2", "_3", "_4", "_5", "_6", "_7", "_8", ""]:
	filein_name    = "UK_5km" + str(el) + ".asc"
	filein_handler = open(filein_name, "r")
	fileout_name    = "UK_5km" + str(el) + "_base.asc"
	fileout_handler = open(fileout_name, "w")
	for line in filein_handler:
		fileout_handler.write( line )
	filein_handler.close()
	fileout_handler.close()

	filein_name    = "UK_5km" + str(el) + "_base.asc"
	filein_handler = open(filein_name, "r")
	fileout_name    = "UK_5km" + str(el) + ".asc"
	fileout_handler = open(fileout_name, "w")
	lcount = 0
	for line in filein_handler:
		vec = line.split()
		if lcount < 6:
			pass
		elif lcount >= 140 and lcount <= 175:
			if lcount >= 140 and lcount <= 144:
				lim = shorten[lcount-140]
			else:
				lim = 79
			vec = ["0" for x in range(0, lim+1)] + vec[lim+1:]
		lcount += 1
		fileout_handler.write( " ".join(vec) + "\n" )
	filein_handler.close()
	fileout_handler.close()



