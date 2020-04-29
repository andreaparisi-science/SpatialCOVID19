import sys
from generate_base import *

country = sys.argv[1]

cores = simList[country]['cores']
mobility = simList[country]['mobility']
simLen   = simList[country]['simLen']
gridres  = simList[country]['gridres']

if len(sys.argv) == 2:
	readinit = "\tReadinit		Read\n"
	simlen   = "\tSimlength\t\t10\n"
else:
	readinit = "\tReadinit		" + sys.argv[2] + "\n"
	if sys.argv[2] == "Write":
		simlen = "\tSimlength\t\t10\n"
	elif sys.argv[2] == "Read":
		simlen = "\tSimlength\t\t" + str(simLen) + "\n"
	else:
		print("Unrecognized parameter value.\n")
		sys.exit(1)

nAgeGroups = 17
TOTCLASSES = (10)*nAgeGroups

outputLse  = "GENERAL:\n"
outputLse += "\tPartfile\t\t\"../../../" + country + "/Setup/" + country + "_" + str(gridres) + "km-{:03d}.ppm\"\n".format(cores)
#outputLse += "\tPartfile\t\t\"../Simple.txt\"\n"
outputLse += "\tLocations\t\t50\n"
outputLse += readinit
outputLse += "\tStorefile\t\t\"../storage\"\n"
outputLse += "\tReduceFactor\t1\n"
outputLse += simlen
outputLse += "\tIndivPrefLocs\t0\n"
outputLse += "\tIndivHasData\t1\n"
outputLse += "\tAccessCycle\t\t1\n\n"
outputLse += "POPULATIONS:\n"
outputLse += "\tDefault	\"../../../"+country+"/Setup/"+country+"_"+str(gridres)+"km.asc\"\thandleMobility\tRadiationNorm\t"+str(mobility)+"\n\n"
#outputLse += "\tDefault	\"../../../"+country+"/Setup/"+country+"_"+str(gridres)+"km.asc\"\tNone\tRadiationNorm\t"+str(mobility)+"\n\n"

outputClasses = "CLASSES:\n"
for agegr in range(0,nAgeGroups):
	outputClasses += "\tSus_" + str(agegr) + "\t1\n"
	outputClasses += "\tEsp_" + str(agegr) + "\t2\n"
	outputClasses += "\tInf_" + str(agegr) + "\t2\n"
	outputClasses += "\tAsy_" + str(agegr) + "\t2\n"
	outputClasses += "\tHos_" + str(agegr) + "\t1\n"
	outputClasses += "\tIcu_" + str(agegr) + "\t1\n"
	outputClasses += "\tHom_" + str(agegr) + "\t1\n"
	outputClasses += "\tRec_" + str(agegr) + "\t1\n"
	outputClasses += "\tDed_" + str(agegr) + "\t1\n"
	outputClasses += "\tAsyrec_" + str(agegr) + "\t1\n"

sympt_to_hosp    = [3.8, 3.8, 2.6, 2.6, 2.8, 2.8, 2.7, 2.7, 5.4, 5.4, 12.6, 12.6, 19.7, 19.7, 28.7, 28.7, 27.3]
sympt_to_critcal = [0.01, 0.01, 0.02, 0.02, 0.08, 0.08, 0.18, 0.18, 0.34, 0.34, 1.5, 1.5, 5.4, 5.4, 12.4, 12.4, 19.3]
hosp_to_death    = [3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 4.0, 4.5, 5.6, 7.8, 11.3, 16.9, 23.2, 29.1, 34.8, 53.5]
critic_to_death  = 0.33  # Age independent?

outputParams = "PARAMETERS:\n"
outputParams += "\tRestart\t\tNone\t0\n"
outputParams += "\ttau\t\tNone\t1.0\n"
outputParams += "\tzmax\t\tNone\t1.0\n"
outputParams += "\tt0\t\tNone\t1000000.0\n"
outputParams += "\teta\t\tNone\t1.0\n"
outputParams += "\tR0\t\tNone\t2.5\n"
outputParams += "\tgamma\tNone\t1.0/2.76\n"
outputParams += "\tsigma\tNone\t1.0/5.2\n"
outputParams += "\tbeta\tNone\tR0*gamma\n"
outputParams += "\t#eigen\tNone\t23.1567\n"
outputParams += "\t# Kenya's eigenvalue is 23.1567\n"
outputParams += "\t# Italy eigenvalue is 17.0450\n"
outputParams += "\t# Using the Chinese eigenvalue (see Sam's idea)\n"
outputParams += "\teigen\tNone\t16.7448\n"
#outputParams += "\tbeta\tNone\tR0*gamma/20.0891\n"
outputParams += "\tmuYoung\tNone\t1.0/14.0\n"
outputParams += "\tmuOld\tNone\t1.0/20.0\n"
outputParams += "\ttheta\tNone\t1.0\n"
outputParams += "\tzeta\tNone\t1.0/7.0\n"
outputParams += "\thome\tNone\t1.0\n"
outputParams += "\twork\tNone\t1.0\n"
outputParams += "\tschool\tNone\t1.0\n"
outputParams += "\tother\tNone\t1.0\n"
outputParams += "\tmobility\tNone\t1.0\n"
for xx in range(0, nAgeGroups):
	for yy in range(0, nAgeGroups):
		outputParams += "\tKK_" + str(xx) + "_" + str(yy) + "\tNone\t1.0\n"
		outputParams += "\tKK_home_" + str(xx) + "_" + str(yy) + "\tNone\t1.0\n"
		outputParams += "\tKK_work_" + str(xx) + "_" + str(yy) + "\tNone\t1.0\n"
		outputParams += "\tKK_school_" + str(xx) + "_" + str(yy) + "\tNone\t1.0\n"
		outputParams += "\tKK_other_" + str(xx) + "_" + str(yy) + "\tNone\t1.0\n"
for xx in range(0, nAgeGroups):
	outputParams += "\tzz_" + str(xx) + "\tNone\t1.0\n"
	outputParams += "\tsy2ho_" + str(xx) + "\tNone\t"+str(sympt_to_hosp[xx]/100.0)+"\n"
	outputParams += "\tsy2ic_" + str(xx) + "\tNone\t"+str(sympt_to_critcal[xx]/100.0)+"\n"
	outputParams += "\tho2de_" + str(xx) + "\tNone\t"+str(hosp_to_death[xx]/100.0)+"\n"
outputParams += "\tcr2de\tNone\t" + str(critic_to_death)+"\n"

outputTrans  = "TRANSITIONS:\n"
for xx in range(0, nAgeGroups):
	for yy in range(0, nAgeGroups):
		outputTrans += "\tSus_" + str(xx) + "\tInf_" + str(yy) + "\tEsp_" + str(xx) + "(First)\tinfectother\tNone\t(beta/eigen)*KK_"+str(xx)+"_"+str(yy)
		outputTrans += "*@Default/(@Sus_"+str(yy)+"+@Esp_"+str(yy)+"+@Inf_"+str(yy)+"+@Asy_"+str(yy)+"+@Hos_"+str(yy)+"+@Icu_"+str(yy)+"+@Hom_"+str(yy)+"+@Ded_"+str(yy)+"+@Asyrec_"+str(yy)+"+@Rec_"+str(yy)+")"
		outputTrans += "\n"
		outputTrans += "\tSus_" + str(xx) + "\tAsy_" + str(yy) + "\tEsp_" + str(xx) + "(First)\tinfectother\tNone\t(beta/eigen)*KK_"+str(xx)+"_"+str(yy)
		outputTrans += "*@Default/(@Sus_"+str(yy)+"+@Esp_"+str(yy)+"+@Inf_"+str(yy)+"+@Asy_"+str(yy)+"+@Hos_"+str(yy)+"+@Icu_"+str(yy)+"+@Hom_"+str(yy)+"+@Ded_"+str(yy)+"+@Asyrec_"+str(yy)+"+@Rec_"+str(yy)+")"
		outputTrans += "\n"
		outputTrans += "\tSus_" + str(xx) + "\tInf_" + str(yy) + "\tSus_" + str(xx) + "(First)\tinfectnone\tNone\t(1-beta/eigen)*KK_"+str(xx)+"_"+str(yy)
		outputTrans += "*@Default/(@Sus_"+str(yy)+"+@Esp_"+str(yy)+"+@Inf_"+str(yy)+"+@Asy_"+str(yy)+"+@Hos_"+str(yy)+"+@Icu_"+str(yy)+"+@Hom_"+str(yy)+"+@Ded_"+str(yy)+"+@Asyrec_"+str(yy)+"+@Rec_"+str(yy)+")"
		outputTrans += "\n"

		outputTrans += "\tSus_" + str(xx) + "\tInf_" + str(yy) + "\tEsp_" + str(xx) + "(First)\tinfectatwork\tNone\t(beta/eigen)*KK_work_"+str(xx)+"_"+str(yy)
		outputTrans += "*@Default/(@Sus_"+str(yy)+"+@Esp_"+str(yy)+"+@Inf_"+str(yy)+"+@Asy_"+str(yy)+"+@Hos_"+str(yy)+"+@Icu_"+str(yy)+"+@Hom_"+str(yy)+"+@Ded_"+str(yy)+"+@Asyrec_"+str(yy)+"+@Rec_"+str(yy)+")"
		outputTrans += "\n"
		outputTrans += "\tSus_" + str(xx) + "\tAsy_" + str(yy) + "\tEsp_" + str(xx) + "(First)\tinfectatwork\tNone\t(beta/eigen)*KK_work_"+str(xx)+"_"+str(yy)
		outputTrans += "*@Default/(@Sus_"+str(yy)+"+@Esp_"+str(yy)+"+@Inf_"+str(yy)+"+@Asy_"+str(yy)+"+@Hos_"+str(yy)+"+@Icu_"+str(yy)+"+@Hom_"+str(yy)+"+@Ded_"+str(yy)+"+@Asyrec_"+str(yy)+"+@Rec_"+str(yy)+")"
		outputTrans += "\n"
		outputTrans += "\tSus_" + str(xx) + "\tInf_" + str(yy) + "\tSus_" + str(xx) + "(First)\tinfectnone\tNone\t(1-beta/eigen)*KK_work_"+str(xx)+"_"+str(yy)
		outputTrans += "*@Default/(@Sus_"+str(yy)+"+@Esp_"+str(yy)+"+@Inf_"+str(yy)+"+@Asy_"+str(yy)+"+@Hos_"+str(yy)+"+@Icu_"+str(yy)+"+@Hom_"+str(yy)+"+@Ded_"+str(yy)+"+@Asyrec_"+str(yy)+"+@Rec_"+str(yy)+")"
		outputTrans += "\n"
		#outputTrans += "\tSus_" + str(xx) + "\tAsy_" + str(yy) + "\tSus_" + str(xx) + "\tinfect\tNone\t(1-beta/eigen)*KK_"+str(xx)+"_"+str(yy)
		#outputTrans += "*@Default/(eigen*(@Sus_"+str(yy)+"+@Esp_"+str(yy)+"+@Inf_"+str(yy)+"+@Asy_"+str(yy)+"+@Rec_"+str(yy)+"))"
		#outputTrans += "\n"

for xx in range(0, nAgeGroups):
	outputTrans += "\tEsp_" + str(xx) + "\tNone\tEsp_" + str(xx) + "(Next)\tNone\tNone\t$Esp_"+str(xx)+"*sigma\n"
	outputTrans += "\tEsp_" + str(xx) + "(Last)\tNone\tInf_" + str(xx) + "(First)\tmoveToInfect\tNone\t$Esp_"+str(xx)+"*sigma*zz_"+str(xx)+"\n"
	outputTrans += "\tEsp_" + str(xx) + "(Last)\tNone\tAsy_" + str(xx) + "(First)\tNone\tNone\t$Esp_"+str(xx)+"*sigma*(1-zz_"+str(xx)+")\n"
	outputTrans += "\tInf_" + str(xx) + "\tNone\tInf_" + str(xx) + "(Next)\tNone\tNone\t$Inf_"+str(xx)+"*gamma\n"
	outputTrans += "\tAsy_" + str(xx) + "\tNone\tAsy_" + str(xx) + "(Next)\tNone\tNone\t$Asy_"+str(xx)+"*gamma\n"
	outputTrans += "\tInf_" + str(xx) + "(Last)\tNone\tHos_" + str(xx) + "\tmoveToHos\tCaseHos\t$Inf_"+str(xx)+"*gamma*sy2ho_"+str(xx)+"\n"
	outputTrans += "\tInf_" + str(xx) + "(Last)\tNone\tIcu_" + str(xx) + "\tmoveToIcu\tCaseIcu\t$Inf_"+str(xx)+"*gamma*sy2ic_"+str(xx)+"\n"
	outputTrans += "\tInf_" + str(xx) + "(Last)\tNone\tHom_" + str(xx) + "\tmoveToHom\tCaseHom\t$Inf_"+str(xx)+"*gamma*(1.0-sy2ho_"+str(xx)+"-sy2ic_"+str(xx)+")\n"
	outputTrans += "\tAsy_" + str(xx) + "(Last)\tNone\tAsyrec_" + str(xx) + "\trecoveryHidden\tHidden\t$Asy_"+str(xx)+"*gamma\n"
	if xx >= 14:
		outputTrans += "\tHos_" + str(xx) + "\tNone\tDed_" + str(xx) + "\tdedHos\tDeaths\tmuOld*theta*ho2de_"+str(xx)+"\n"
		outputTrans += "\tHos_" + str(xx) + "\tNone\tRec_" + str(xx) + "\trecHos\tRecovs\tmuOld*theta*(1.0-ho2de_"+str(xx)+")\n"
		outputTrans += "\tIcu_" + str(xx) + "\tNone\tDed_" + str(xx) + "\tdedIcu\tDeaths\tmuOld*theta*cr2de\n"
		outputTrans += "\tIcu_" + str(xx) + "\tNone\tRec_" + str(xx) + "\trecIcu\tRecovs\tmuOld*theta*(1.0-cr2de)\n"
	else:
		outputTrans += "\tHos_" + str(xx) + "\tNone\tDed_" + str(xx) + "\tdedHos\tDeaths\tmuYoung*theta*ho2de_"+str(xx)+"\n"
		outputTrans += "\tHos_" + str(xx) + "\tNone\tRec_" + str(xx) + "\trecHos\tRecovs\tmuYoung*theta*(1.0-ho2de_"+str(xx)+")\n"
		outputTrans += "\tIcu_" + str(xx) + "\tNone\tDed_" + str(xx) + "\tdedIcu\tDeaths\tmuYoung*theta*cr2de\n"
		outputTrans += "\tIcu_" + str(xx) + "\tNone\tRec_" + str(xx) + "\trecIcu\tRecovs\tmuYoung*theta*(1.0-cr2de)\n"
	outputTrans += "\tHom_" + str(xx) + "\tNone\tRec_" + str(xx) + "\trecHom\tRecovs\tzeta\n"

#for xx in range(0, nAgeGroups-1):
#	outputTrans += "\tSus_" + str(xx) + "\tNone\tSus_" + str(xx+1) + "\tNone\tNone\tageing_"+str(xx)+"\n"
#	outputTrans += "\tEsp_" + str(xx) + "\tNone\tEsp_" + str(xx+1) + "\tNone\tNone\tageing_"+str(xx)+"\n"
#	outputTrans += "\tInf_" + str(xx) + "\tNone\tInf_" + str(xx+1) + "\tNone\tNone\tageing_"+str(xx)+"\n"
#	outputTrans += "\tRec_" + str(xx) + "\tNone\tRec_" + str(xx+1) + "\tNone\tNone\tageing_"+str(xx)+"\n"

outputInit  = "INIT:\n"
#for xx in range(0, nAgeGroups):
outputInit += "\tSus_0\tinit\t@Default\tNone\tNone\tNone\n"
#outputInit += "\tInf_0\tinit\t1\tNone\tNone\tNone\n\n"

outputLine  = "OUTPUTLINE:\n"
outputLine += "\tDaily\tprintextra\t@@CaseHos\t@@CaseIcu\t@@CaseHom\t@@Hidden\t@@Recovs\t@@Deaths\n"
outputLine += "\tWeekly\tNone\t@@CaseHos\t@@CaseIcu\t@@CaseHom\t@@Hidden\t@@Recovs\t@@Deaths"
#for xx in range(0, nAgeGroups):
#	if xx > 0:
#		outputLine += "+"
#	outputLine += "@@Inf_" + str(xx)


outputMap  = "OUTPUTMAP:\n"
outputMap += "\tWeekly\tprintmapextra\t@CaseHos+@CaseIcu+@CaseHom\t@Hidden\t@Recovs\t@Deaths"
#for xx in range(0, nAgeGroups):
#	if xx > 0:
#		outputMap += "+"
#	outputMap += "@Inf_" + str(xx)


print(outputLse + "\n")
print(outputClasses + "\n")
print(outputParams + "\n")
print(outputTrans + "\n")
print(outputInit + "\n")
print(outputLine + "\n")
print(outputMap + "\n")


contactTypes = ["_home", "_school", "_work", "_other"]
outputContact  = "void  updateContactMatrix()  {\n"
for ii in range(0, nAgeGroups):
	for jj in range(0, nAgeGroups):
		outputContact += "\tparams.KK_" + str(ii) + "_" + str(jj) + " = "
		outputContact += "params.home*params.KK_home_" + str(ii) + "_" + str(jj) + "+"
		#outputContact += "params.work*params.KK_work_" + str(ii) + "_" + str(jj) + "+"
		outputContact += "params.school*params.KK_school_" + str(ii) + "_" + str(jj) + "+"
		outputContact += "params.other*params.KK_other_" + str(ii) + "_" + str(jj) + ";\n"
#outputContact += "\tarma::mat Kmat(" + str(nAgeGroups) + ", " + str(nAgeGroups) + ", arma::fill::zeros);\n"
#for ii in range(0, nAgeGroups):
#	for jj in range(0, nAgeGroups):
#		outputContact += "\tKmat(" + str(ii) + ", " + str(jj) + ") = params.KK_" + str(ii) + "_" + str(jj) + ";\n"
#outputContact += "\tarma::cx_vec eigenvals = arma::eig_gen(Kmat);\n"
#outputContact += "\tparams.eigen = std::abs( eigenvals.max() );\n"
#outputContact += "\tstd::cout << \"EIGENVAL : \" << params.eigen;\n"
outputContact += "}\n\n"


outputContact += "void  initContactMatrix()  {\n"
outputContact += "\tint nAgeGroups = groups[ groups.size()-1 ] + 1;\n"
outputContact += "\tdouble  totContacts;\n"
outputContact += "\tifstream  handler;\n"

for el in contactTypes:
	outputContact += "\thandler.open( contactMatrixFile + \"" + el + ".csv\" );\n"
	outputContact += "\tif (!handler.good())  {\n"
	outputContact += "\t\tsimStatus.exit(\"Contact matrix file [\" + contactMatrixFile + \"" + el + ".csv] not found\" );\n"
	outputContact += "\t}\n"
	#outputContact += "std::cout << sizes[0] << \"\t\" << params.KK_home_0_0 << \"\t\" << simStatus.getTotalPopulationSize() << \"\\n\" << std::flush;\n"
	for ii in range(0, 16):
		for jj in range(0, 16):
			outputContact += "\thandler >> params.KK" + el + "_" + str(ii) + "_" + str(jj) + ";\n"
			outputContact += "\ttotContacts += sizes[" +str(ii) + "]*params.KK" + el + "_"+str(ii)+"_"+str(jj)+";\n"
	#	outputContact += "\tparams.KK_" + str(ii) + "_16 = params.KK_" + str(ii) + "_15;\n"
	#for jj in range(0, nAgeGroups-1):
	#	outputContact += "\tparams.KK_16_" + str(jj) + " = params.KK_15_" + str(jj) + ";\n"
	for ii in range(0, 16):
		outputContact += "\tparams.KK" + el + "_" + str(ii) + "_16 = params.KK" + el + "_" + str(ii) + "_15;\n"
	for ii in range(0, nAgeGroups):
		outputContact += "\tparams.KK" + el + "_16_" + str(ii) + " = params.KK" + el + "_15_" + str(ii) + ";\n"
	outputContact += "\tstd::cout << totContacts/simStatus.getTotalPopulationSize() << \" \" << params.R0/(totContacts/simStatus.getTotalPopulationSize()) << \"\\n\";\n"
	outputContact += "\thandler.close();\n\n"

outputContact += "\tlocalKK_home.resize( " + str(nAgeGroups) + " );\n";
outputContact += "\tfor (int qq = 0; qq < " + str(nAgeGroups) + "; qq++)  {\n";
outputContact += "\t\tlocalKK_home[qq].resize( " + str(nAgeGroups) + " );\n";
outputContact += "\t}\n";
for ii in range(0, nAgeGroups):
	for jj in range(0, nAgeGroups):
		outputContact += "\tlocalKK_home[" + str(ii) + "][" + str(jj) + "] = params.KK_home_" + str(ii) + "_" + str(jj) + ";\n"
outputContact += "}\n\n"


handler = open( "lse-userdefined-base.cpp", "r" )
contents = handler.readlines()
handler.close()

handler = open( "lse-userdefined.cpp", "w" )
#handler.writelines( "int constexpr TOTCLASSES = " + str(TOTCLASSES) + ";\n" )
handler.writelines( "static std::vector< std::vector<double> >  localKK_home;\n" )
handler.writelines( contents )
handler.write( outputContact )
handler.close()

