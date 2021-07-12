# CDmetaPOP.py
# Author: Erin L Landguth
# Created: February 2008
# v 1.0 Release: MARCH 2014
# ----------------------------------------------------------------------------
# General CDmetaPOP information
appName = "CDmetaPOP"
appVers = "version 1.70"
appRele = "2021.7.12-14:54:01MST"
authorNames = "Erin L Landguth, Casey Day, Andrew Bearlin"

# ---------------
# Global symbols
#----------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False
# File absolute paths for importing functions
SRC_PATH =  "../src/"

# ------------------------------------------
# Import Modules with Except/Try statements
# ------------------------------------------
# Python specific functions
import datetime,time,pdb,os,sys,shutil,gc,multiprocessing,warnings
	
# Numpy functions
try:
	import numpy as np                    
except ImportError as eMsg:
	print("ImportError (%s) Numpy required."%(eMsg))
	sys.exit(-1)

#Import the package specific folders
CDPOP_folder = os.path.dirname(os.path.abspath(SRC_PATH+"CDmetaPOP"))

if CDPOP_folder not in sys.path:
     sys.path.insert(0, CDPOP_folder)

# CDmetaPOP functions
try:
	from CDmetaPOP_Modules import * 
	from CDmetaPOP_PostProcess import *
	from CDmetaPOP_PreProcess import *
	from CDmetaPOP_Mate import *
	from CDmetaPOP_Emigration import *
	from CDmetaPOP_Immigration import *
	from CDmetaPOP_Offspring import *
	from CDmetaPOP_Mortality import *
except ImportError:
	raise ImportError, "CDmetaPOP Modules required."	
				
#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 
if __name__ == '__main__':
	warnings.filterwarnings("ignore")	
	# ------------------------------------------------------	
	# Start timer, get script arguments, create log writeout
	# ------------------------------------------------------
	# Timing events: start
	start_time = datetime.datetime.now()
	foldertime = int(time.time())
	
	if len(sys.argv) >= 4:
		datadir = sys.argv[1]+'/'
		fileans = datadir+sys.argv[2]
		outdir = datadir+sys.argv[3]+str(foldertime)+'/'
		if len(sys.argv) == 5:
			noproc = int(sys.argv[4])
		else: # assume one processor
			noproc = 1
	
	# If user did not specify .rip file
	else:
		print "User must specify data directory, input file name, and output file directory (e.g., at command line type CDmetaPOP.py ../CDmetaPOP_data/ PopVars.csv exampleout_foldername)."
		sys.exit(-1)	
	
	# If .ip file does not exist
	if not os.path.exists(fileans):
		print("Cannot find or open runtime inputs file(%s)"%(fileans))
		sys.exit(-1)
	
	# If entered more processors than machine
	if noproc >= multiprocessing.cpu_count():
		print('Warning: Specified more CPUs than on local machine. Using one less than '+str(multiprocessing.cpu_count()))
	
	# Create output file directory - will automatically put in the data directory
	os.mkdir(outdir)
	
	# This properly names log file
	logSessionPath = outdir+"CDmetaPOP.log"
	logfHndl =open(logSessionPath,'w')
	
	msgVerbose = True
	logMsg(logfHndl,"\n%s Release %s Version %s\n"%(appName,appRele,appVers))
	logMsg(logfHndl,"Author(s): %s"%(authorNames)+'\n')
	logMsg(logfHndl,"Session runtime inputs from: %s"%(fileans)+'\n\n')    
	msgVerbose = False
	
	# ------------------------------------	
	# Call DoUserInput()
	# ------------------------------------
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Call function and store inputvariables
	batchVars,batchVarsIndex,nSimulations = loadFile(fileans,1,',',True)
	
	# Print to log
	stringout = 'DoUserInput(): '+str(datetime.datetime.now() -start_time1) + ''
	logMsg(logfHndl,stringout)
	print 'DoUserInput(): ',str(datetime.datetime.now() -start_time1),''

	# ----------------------------------------	
	# Begin Batch Looping - assign processors
	# ----------------------------------------
	# This loop is defined by the number of rows in inputvariables.csv
	for ibatch in xrange(nSimulations):
		
		# Timing events: start
		start_timeB = datetime.datetime.now()
		
		# Store all information and the type of each, also do some error checks 
		xyfilename = datadir+batchVars['xyfilename'][ibatch]
		#agefilename = datadir+batchVars['agefilename'][ibatch]
		constMortans = batchVars['constMortans'][ibatch]
		mcruns = int(batchVars['mcruns'][ibatch])
		looptime = int(batchVars['runtime'][ibatch])
		nthfile_out = batchVars['output_years'][ibatch]
		outputans = str(batchVars['summaryOutput'][ibatch])
		cdclimgentimelist = batchVars['cdclimgentime'][ibatch]
		matecdmatfile = batchVars['mate_cdmat'][ibatch]
		dispOutcdmatfile = batchVars['migrateout_cdmat'][ibatch]
		dispBackcdmatfile = batchVars['migrateback_cdmat'][ibatch]
		straycdmatfile = batchVars['stray_cdmat'][ibatch]
		dispLocalcdmatfile = batchVars['disperseLocal_cdmat'][ibatch]
		matemoveno = batchVars['matemoveno'][ibatch]
		matemoveparA = batchVars['matemoveparA'][ibatch]
		matemoveparB = batchVars['matemoveparB'][ibatch]
		matemoveparC = batchVars['matemoveparC'][ibatch]
		matemovethreshval = batchVars['matemovethresh'][ibatch]
		freplace = batchVars['Freplace'][ibatch]
		mreplace = batchVars['Mreplace'][ibatch]
		selfing = batchVars['selfans'][ibatch]
		sexans = batchVars['sexans'][ibatch]
		assortmateModel_pass = batchVars['AssortativeMate_Model'][ibatch]
		assortmateC_pass = batchVars['AssortativeMate_Factor'][ibatch]
		dispmoveOutno = batchVars['migratemoveOutno'][ibatch]
		dispmoveOutparA = batchVars['migratemoveOutparA'][ibatch]
		dispmoveOutparB = batchVars['migratemoveOutparB'][ibatch]
		dispmoveOutparC = batchVars['migratemoveOutparC'][ibatch]
		dispmoveOutthreshval = batchVars['migratemoveOutthresh'][ibatch]
		dispmoveBackno = batchVars['migratemoveBackno'][ibatch]
		dispmoveBackparA = batchVars['migratemoveBackparA'][ibatch]
		dispmoveBackparB = batchVars['migratemoveBackparB'][ibatch]
		dispmoveBackparC = batchVars['migratemoveBackparC'][ibatch]
		dispmoveBackthreshval = batchVars['migratemoveBackthresh'][ibatch]
		StrBackno = batchVars['StrayBackno'][ibatch]
		StrBackparA = batchVars['StrayBackparA'][ibatch]
		StrBackparB = batchVars['StrayBackparB'][ibatch]
		StrBackparC = batchVars['StrayBackparC'][ibatch]
		StrBackthreshval = batchVars['StrayBackthresh'][ibatch]
		dispLocalno = batchVars['disperseLocalno'][ibatch]
		dispLocalparA = batchVars['disperseLocalparA'][ibatch]
		dispLocalparB = batchVars['disperseLocalparB'][ibatch]
		dispLocalparC = batchVars['disperseLocalparC'][ibatch]
		dispLocalthreshval = batchVars['disperseLocalthresh'][ibatch]
		homeattempt = batchVars['HomeAttempt'][ibatch]
		offno = batchVars['offno'][ibatch]
		inheritans_classfiles = batchVars['offans_InheritClassVars'][ibatch]
		equalClutch = batchVars['equalClutchSize'][ibatch]
		eggFreq = float(batchVars['eggFrequency'][ibatch])
		gridformat = batchVars['gridformat'][ibatch]
		gridsample = batchVars['gridsampling'][ibatch]
		muterate = float(batchVars['muterate'][ibatch])
		mutationans = batchVars['mutationtype'][ibatch]
		loci = int(batchVars['loci'][ibatch])
		alleles = batchVars['alleles'][ibatch]
		mtdna = batchVars['mtdna'][ibatch]
		geneswap = int(batchVars['startGenes'][ibatch])
		cdevolveans = batchVars['cdevolveans'][ibatch]
		burningen_cdevolve = int(batchVars['startSelection'][ibatch])
		timecdevolve = batchVars['implementSelection'][ibatch]
		
		plasticans = batchVars['plasticgeneans'][ibatch]
		burningen_plastic = int(batchVars['startPlasticgene'][ibatch])
		timeplastic = batchVars['implementPlasticgene'][ibatch]
		
		cdinfect = batchVars['cdinfect'][ibatch]
		transmissionprob = float(batchVars['transmissionprob'][ibatch])
		growans = batchVars['growth_option'][ibatch]
		sizeLoo = batchVars['growth_Loo'][ibatch] # Check sex ratios
		sizeR0 = batchVars['growth_R0'][ibatch]# Check sex ratios
		size_eqn_1 = batchVars['growth_temp_max'][ibatch]# Check sex ratios
		size_eqn_2 = batchVars['growth_temp_CV'][ibatch]# Check sex ratios
		size_eqn_3 = batchVars['growth_temp_t0'][ibatch]# Check sex ratios
		sizeans = batchVars['sizecontrol'][ibatch]
		mat_set = batchVars['mature_length_set'][ibatch]
		mat_slope = batchVars['mature_eqn_slope'][ibatch]
		mat_int = batchVars['mature_eqn_int'][ibatch]
		egg_mean_ans = batchVars['Egg_Mean_ans'][ibatch]
		egg_mean_1 = float(batchVars['Egg_Mean_par1'][ibatch])
		egg_mean_2 = float(batchVars['Egg_Mean_par2'][ibatch])
		egg_percmort_mu = float(batchVars['Egg_Mortality'][ibatch])
		egg_percmort_sd = float(batchVars['Egg_Mortality_StDev'][ibatch])
		Femalepercent_egg = batchVars['Egg_FemalePercent'][ibatch]
		packans = batchVars['popmodel'][ibatch]
		packpar1 = float(batchVars['popmodel_par1'][ibatch])
		cor_mat_ans = batchVars['correlation_matrix'][ibatch]
		defaultAgeMature = batchVars['mature_defaultAge'][ibatch]
		subpopmort_pass = batchVars['subpopmort_file'][ibatch]
			
		# -------------------------------
		# Distill some vars
		# -------------------------------
		# Grab the nthfile list range specific to user input, list or sequence
		if not isinstance(nthfile_out, (list,tuple)):
			nthfile_out = int(nthfile_out)
			if nthfile_out != 0:
				nthfile = range(0,looptime+nthfile_out,nthfile_out)
				del(nthfile[-1]) # Delete the last value 0, looptime - 1
			else:
				nthfile = [0]
		# If specified years with |
		else:
			nthfile = []
			# Split up list, removing space values, and appending to nthfile
			for inum in xrange(len(nthfile_out)):
				# Error check here if | at the end
				if len(nthfile_out[inum]) != 0:
					nthfile.append(int(nthfile_out[inum]))
		
		# Error check on nthfile, must be 1 less than looptime for indexing
		if max(nthfile) >= looptime:
			print('nthfile selection maximum value must be less than to looptime.')
			sys.exit(-1)
				
		# Store cdmat file information - header file (loadFile()) passes tuple or string if only 1
		if not isinstance(cdclimgentimelist, (list,tuple)):
			cdclimgentime = [cdclimgentimelist]
		else: 
			cdclimgentime = cdclimgentimelist
			
		# Create allele array
		if len(alleles.split(':')) == 1:
			alleles = int(batchVars['alleles'][ibatch])*np.ones(loci,int)
		else:
			alleles = np.asarray(alleles.split(':'),dtype = int)
		
		# ----------------------------
		# For Sex ratios option splits
		# ----------------------------
		# Deterministic mature set value either age or size
		if len(mat_set.split('~')) == 1:
			Fmat_set = mat_set.split('~')[0]
			Mmat_set = mat_set.split('~')[0]
			YYmat_set = mat_set.split('~')[0]
		elif len(mat_set.split('~')) == 2:
			Fmat_set = mat_set.split('~')[0]
			Mmat_set = mat_set.split('~')[1]
			YYmat_set = mat_set.split('~')[1]
		elif len(mat_set.split('~')) == 3:
			Fmat_set = mat_set.split('~')[0]
			Mmat_set = mat_set.split('~')[1]
			YYmat_set = mat_set.split('~')[2]
		else:
			print('mature_length_set must be 1 value for all sex classes or separated by :')
			sys.exit(-1)
		
		# Logistic equation for maturation as a function of size
		if len(mat_slope.split('~')) == 1:
			Fmat_slope = float(mat_slope.split('~')[0])
			Mmat_slope = float(mat_slope.split('~')[0])
			YYmat_slope = float(mat_slope.split('~')[0])
		elif len(mat_slope.split('~')) == 2:
			Fmat_slope = float(mat_slope.split('~')[0])
			Mmat_slope = float(mat_slope.split('~')[1])
			YYmat_slope = float(mat_slope.split('~')[1])
		elif len(mat_slope.split('~')) == 3:
			Fmat_slope = float(mat_slope.split('~')[0])
			Mmat_slope = float(mat_slope.split('~')[1])
			YYmat_slope = float(mat_slope.split('~')[2])
		else:
			print('logistic maturation equation parameter values must be 1 value for all sex classes or separated by :')
			sys.exit(-1)
		if len(mat_int.split('~')) == 1:
			Fmat_int = float(mat_int.split('~')[0])
			Mmat_int = float(mat_int.split('~')[0])
			YYmat_int = float(mat_int.split('~')[0])
		elif len(mat_int.split('~')) == 2:
			Fmat_int = float(mat_int.split('~')[0])
			Mmat_int = float(mat_int.split('~')[1])
			YYmat_int = float(mat_int.split('~')[1])
		elif len(mat_int.split('~')) == 3:
			Fmat_int = float(mat_int.split('~')[0])
			Mmat_int = float(mat_int.split('~')[1])
			YYmat_int = float(mat_int.split('~')[2])
		else:
			print('logistic maturation equation parameter values must be 1 value for all sex classes or separated by :')
			sys.exit(-1)
			
		# ---------------------------------
		# Some Error checking
		# ---------------------------------
		
		# DoEmigration() skipped and warning for if selection on
		if cdevolveans != 'N':
			if timecdevolve.find('Out') != -1:
				if dispOutcdmatfile == 'N':
					print('Warning: DoEmigration module skipped and spatial selection during this time frame specified, which will also be skipped.')
		# 	
		# Constant mortality checks
		if not (constMortans == '1' or constMortans == '2'):
			print('Constant mortalities are compounded using option 1 or 2 specifiy correct values. If no constant mortalities are entered, then enter 1.')
			sys.exit(-1)
		
		# Check on cdevolve answer input
		if not (cdevolveans == '1' or cdevolveans == '2' or cdevolveans == '1_mat' or cdevolveans == '2_mat' or cdevolveans == 'N' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == 'stray' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans.split('_')[0] == 'Hindex'):
			print('CDEVOLVE answer either N, 1, 2, M, G, MG_ind, MG_link, 1_mat, 2_mat, stray, 1_G_ind, 1_G_link, Hindex or Plastic.')
			sys.exit(-1)
			
		# For mature and size ans
		if (cdevolveans == 'M' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == 'G' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and sizeans == 'N':
			print('CDEVOLVE answer is M or G and size answer must be Y.')
			sys.exit(-1)
		
		# For Hindex answer
		if cdevolveans.split('_')[0] == 'Hindex':
			# Split for Gaussian
			if cdevolveans.split('_')[1] == 'Gauss':
				if len(cdevolveans.split('_')[2].split(':')) != 6:
					print('CDEVOLVE answer is Hindex and 6 parameters for the Gaussian function must be specified, see user manual and example files.')
					sys.exit(-1)
			elif cdevolveans.split('_')[1] == 'Para':
				if len(cdevolveans.split('_')[2].split(':')) != 3:
					print('CDEOLVE answer is Hindex and 3 parameters for the Parabolic function must be specified, see user manual and example files.')
					sys.exit(-1)
			elif cdevolveans.split('_')[1] == 'Step':
				if len(cdevolveans.split('_')[2].split(':')) != 3:
					print('CDEOLVE answer is Hindex and 3 parameters for the Step function must be specified, see user manual and example files.')
					sys.exit(-1)
			else:
				print('CDEVOLVE and Hindex parameter not entered correctly, check user manual and example files.')
				sys.exit(-1)
		
		# If cdevolve is turned on must have 2 alleles
		if cdevolveans != 'N' and alleles[0] != 2:
			print('Warning: More than 2 alleles per locus specified. CDEVOLVE only considers first 2 alleles in selection models (except Hindex scenario).')
		if plasticans != 'N' and alleles[0] != 2:
			print('Warning: More than 2 alleles per locus specified. Plastic gene turned on and only considers first 2 alleles in this model.')
		
		# For Plastic answer
		if plasticans != 'N':
			# Split for temp
			if (plasticans.split('_')[0] != 'Temp') and (plasticans.split('_')[0] != 'Hab'):
				print('Plastic parameter not entered corectly, check user manual and example files.')
				sys.exit(-1)
		if plasticans != 'N':	
			if ((timeplastic.find('Out') == -1) and (timeplastic.find('Back') == -1)):
				print('Plastic timing must be specified (e.g., Out or Back).')
				sys.exit(-1)
			
		# Must have more than 1 loci
		if loci <= 1:
			print('Currently, CDmetaPOP needs more than 1 locus to run.')
			sys.exit(-1)
		if cdevolveans == '1' or cdevolveans == '2' or cdevolveans == '1_mat' or cdevolveans == '2_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans.split('_')[0] == 'Hindex':
			if ((timecdevolve.find('Eggs') == -1) and (timecdevolve.find('Out') == -1) and (timecdevolve.find('Back') == -1) and timecdevolve.find('packing') == -1):
				print('CDEVOLVE timing must be specified (e.g., Out, Back or Eggs).')
				sys.exit(-1)
			
		# Error check on forward mutation in A and backward mutation in B
		#	Can only happen if cdevolve == 2.
		if mutationans == 'forwardAbackwardBrandomN' and (cdevolveans != '2' or cdevolveans == '2_mat'):
			print('This special case of mutation is for AAbb ancestors and 2-locus selection.')
			sys.exit(-1)		
					
		# grid format
		if (gridformat == 'cdpop' or gridformat == 'general' or gridformat == 'genalex' or gridformat == 'genepop' or gridformat == 'structure') == False:
			print('Grid format parameter not an option.')
			sys.exit(-1)
		
		# If genepop, some conditions
		if gridformat == 'genepop' and (alleles[0] > 99 or loci > 99):
			print('GENEPOP format requires less than 99 alleles and 99 loci.')
			sys.exit(-1)
		
		# Check burn in times
		if cdevolveans != 'N' and burningen_cdevolve < geneswap:
			print('Warning: Selection burnin time < time at which genetic exchange is to initialize, setting burnin time = start genetic exchange time.')
			burningen_cdevolve = geneswap
		if plasticans != 'N' and burningen_plastic < geneswap:
			print('Warning: Selection burnin time < time at which genetic exchange is to initialize, setting burnin time = start genetic exchange time.')
			burningen_plastic = geneswap
			
		# Egg frequency
		if eggFreq > 1:
			print('Egg frequency must be less than or equal to 1.')
			sys.exit(-1)
			
		# Inherit answer can only be:
		if not (inheritans_classfiles == 'random' or inheritans_classfiles == 'Hindex' or inheritans_classfiles == 'mother'):
			print('Inherit answer for multiple class files is not correct: enter either random or Hindex.')
			sys.exit(-1)
			
		# If inherit answer uses Hindex, mutation can't be on
		if muterate != 0.0 and (inheritans_classfiles == 'Hindex' or inheritans_classfiles == 'mother'):
			print('Currently, mutation is not operating with Hindex inheritance options.')
			sys.exit(-1)
				
		# ---------------------------------------------	
		# Begin Monte-Carlo Looping
		# ---------------------------------------------
		
		# xrange(mcruns) is typically 10 - 50...and it takes a long time.
		for ithmcrun in xrange(mcruns):	
			
			# Timing events: start
			start_timeMC = datetime.datetime.now()
		
			# -----------------------------------------
			# Create storage variables
			# ------------------------------------------	
			# These variables will be stored in output.csv at the end of the simulation						
			
			# GetMetrics()
			Track_p1 = []
			Track_p2 = []
			Track_q1 = []
			Track_q2 = []
			Track_Alleles = []
			Track_He = []
			Track_Ho = []
			Track_N_Init_pop = []
			Track_N_Init_age = []
			Track_N_Init_class = []
			Track_K = []
			Track_CaptureCount_Out = []
			Track_CaptureCount_ClassOut = []
			Track_CaptureCount_Back = []
			Track_CaptureCount_ClassBack = []
									
			# DoMate()
			Track_FAvgMate = []
			Track_MAvgMate = []
			Track_FSDMate = []
			Track_MSDMate = []
			Track_MateDistCD = []
			Track_MateDistCDstd = []
			Track_BreedEvents = []
			Track_AAaaMates = []
			Track_AAAAMates =[]
			Track_aaaaMates = []
			Track_AAAaMates = []
			Track_aaAaMates = []
			Track_AaAaMates = []
			Track_YYsAdded = []
			Track_BreedFemales = []
			Track_BreedMales = []
			Track_BreedYYMales = []
			Track_MatureCount = []
			Track_ImmatureCount = []
			Track_ToTFemales = []
			Track_ToTMales = []
			Track_ToTYYMales = []
			
			# DoOffspring
			Track_Births = []
			Track_EggDeaths = []
			Track_BirthsYY = []
			
			# DoUpdate
			Track_N_back_age = []
			Track_N_out_age = []
			
			# Emigration()
			N_Emigration_pop = []
			N_Emigration_age = []
			subpopemigration = []			
			F_EmiDist = []
			M_EmiDist = []						
			F_EmiDist_sd = []
			M_EmiDist_sd = []
			SelectionDeathsEmi = []	
			DisperseDeathsEmi = []
			PackingDeathsEmi = []
			PackingDeathsEmiAge = []
			MgSuccess = []
			AdultNoMg = []
			Track_YYSelectionPackDeathsEmi = []
			Track_WildSelectionPackDeathsEmi = []
			
			# Mortlity after Emigration
			N_EmiMortality = []
			PopDeathsOUT = []
			AgeDeathsOUT = []
			SizeDeathsOUT = []
			
			# Immigration
			N_beforePack_Immi_pop = []
			N_beforePack_Immi_age = []
			N_Immigration_pop = []
			N_Immigration_age = []
			subpopimmigration = []
			F_HomeDist = []
			M_HomeDist = []					
			F_HomeDist_sd = []
			M_HomeDist_sd = []
			F_StrayDist = []
			M_StrayDist = []					
			F_StrayDist_sd = []
			M_StrayDist_sd = []
			F_ZtrayDist = []
			M_ZtrayDist = []					
			F_ZtrayDist_sd = []
			M_ZtrayDist_sd = []
			SelectionDeathsImm = []
			SelectionDeathsImm_Age0s = []
			DisperseDeathsImm = []
			PackingDeathsImmAge = []
			PackingDeathsImm = []
			StrSuccess = []
			Track_YYSelectionPackDeathsImmi = []
			Track_WildSelectionPackDeathsImmi = []

			# Mortality after immigration
			N_ImmiMortality = []
			PopDeathsIN = []
			AgeDeathsIN = []
			SizeDeathsIN = []
			
			# DoOutput()
			Infected = []
			Residors = []
			Strayers1 = []
			Strayers2 = []
			Immigrators = []
			IDispersers = []
			RDispersers = []
			PopSizes_Mean = []
			PopSizes_Std = []
			AgeSizes_Mean = []
			AgeSizes_Std = []
			ClassSizes_Mean = []
			ClassSizes_Std = []
					
			# ------------------------------------	
			# Call DoPreProcess()
			# ------------------------------------
			
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			# Call function
			tupPreProcess = DoPreProcess(outdir,datadir,ibatch,ithmcrun,\
			xyfilename,loci,alleles,\
			0,logfHndl,cdevolveans,cdinfect,\
			subpopemigration,subpopimmigration,sizeans,eggFreq,Fmat_set,Mmat_set,Fmat_int,Fmat_slope,Mmat_int,Mmat_slope,burningen_cdevolve,cor_mat_ans,inheritans_classfiles,sexans,YYmat_set,YYmat_slope,YYmat_int,defaultAgeMature)
			
			ithmcrundir = tupPreProcess[0]			
			fitvals_pass = tupPreProcess[1] # sex ratio check throughout
			allelst = tupPreProcess[2]
			subpopemigration = tupPreProcess[3]
			subpopimmigration = tupPreProcess[4]
			age_percmort_out_mu = tupPreProcess[5]
			age_percmort_back_mu = tupPreProcess[6]
			age_Mg = tupPreProcess[7]# sex ratio check throughout
			age_S = tupPreProcess[8]# sex ratio check throughout			
			age_mu = tupPreProcess[9]
			age_size_mean = tupPreProcess[10]
			age_size_std = tupPreProcess[11]			
			xgridpop = tupPreProcess[12]
			ygridpop = tupPreProcess[13]			
			SubpopIN_init = tupPreProcess[14]
			N0 = tupPreProcess[15]
			K_mu = tupPreProcess[16]
			dtype = tupPreProcess[17]
			outsizevals_pass = tupPreProcess[18]
			backsizevals_pass = tupPreProcess[19]
			popmort_out_pass = tupPreProcess[20]
			popmort_back_pass = tupPreProcess[21]
			Mg_pass = tupPreProcess[22]
			Str_pass = tupPreProcess[23]
			eggmort_pass = tupPreProcess[24]
			setmigrate = tupPreProcess[25]
			#M_mature = tupPreProcess[26]
			#F_mature = tupPreProcess[27]
			age_mature = tupPreProcess[26]# sex ratio check throughout
			age_sigma = tupPreProcess[27]
			outgrowdays_pass = tupPreProcess[28]
			backgrowdays_pass = tupPreProcess[29]
			Kmu_pass = tupPreProcess[30]
			age_capture_out = tupPreProcess[31] # sex ratio check throughout
			age_capture_back = tupPreProcess[32] # sex ratio check throughout
			Kstd_pass = tupPreProcess[33]
			K_std = tupPreProcess[34]
			popmort_out_sd_pass = tupPreProcess[35]
			popmort_back_sd_pass = tupPreProcess[36]
			eggmort_sd_pass = tupPreProcess[37]
			outsizevals_sd_pass = tupPreProcess[38]
			backsizevals_sd_pass = tupPreProcess[39]
			outgrowdays_sd_pass = tupPreProcess[40]
			backgrowdays_sd_pass = tupPreProcess[41]
			size_percmort_out_mu = tupPreProcess[42] 
			size_percmort_back_mu = tupPreProcess[43] 
			age_percmort_out_sd = tupPreProcess[44] 
			age_percmort_back_sd = tupPreProcess[45] 
			size_percmort_out_sd = tupPreProcess[46] 
			size_percmort_back_sd = tupPreProcess[47] 
			pop_capture_back_pass = tupPreProcess[48]
			pop_capture_out_pass = tupPreProcess[49]
			pop_capture_back = tupPreProcess[50]
			natal = tupPreProcess[51]
			cor_mat = tupPreProcess[52]
			migrate = tupPreProcess[53]
			N0_pass = tupPreProcess[54]
			allefreqfiles_pass = tupPreProcess[55]
			classvarsfiles_pass = tupPreProcess[56]
			PopTag = tupPreProcess[57] # To pass into AddIndividuals, Emigration, Immigration
			outhabvals_pass = tupPreProcess[58]
			backhabvals_pass = tupPreProcess[59]
			
			# Grab first one only
			K = K_mu # Initialize K with mu
						
			# Print to log
			stringout = 'DoPreProcess(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print 'DoPreProcess(): ',str(datetime.datetime.now() -start_time1),''
			
			# ---------------------------------
			# Call GetMetrics()
			# ---------------------------------
			
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			GetMetrics(SubpopIN_init,K,Track_N_Init_pop,Track_K,loci,alleles,0,Track_Ho,Track_Alleles,Track_He,Track_p1,Track_p2,Track_q1,Track_q2,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,Track_N_Init_age,sizeans,age_size_mean,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,sexans,packans,RDispersers,IDispersers)
			
			# Print to log
			stringout = 'GetMetrics() Initial: '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print 'GetMetrics() Initial: ',str(datetime.datetime.now() -start_time1),''
			
			# ---------------------------------
			# Error statements
			# ---------------------------------			
			# Error statement here in case no females or males, then break
			#if Track_ToTFemales[0][0]==0 or (Track_ToTMales[0][0] + Track_ToTYYMales[0][0])==0:
			if Track_N_Init_pop[0][0] == 0:
				print('There are no individuals to begin time loop.\n')
				break
				
			# ----------------------------------------------------
			# Call DoUpdate() - output initial file here ind-1.csv
			# ----------------------------------------------------
							
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			DoUpdate(packans,SubpopIN_init,K,xgridpop,ygridpop,-1,nthfile,ithmcrundir,loci,alleles,logfHndl,'Initial','N','N',defaultAgeMature,[],burningen_cdevolve,[],[],[],[],[],[])
			
			# Print to log
			stringout = 'DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print 'First DoUpdate(): ',str(datetime.datetime.now() -start_time1),''	
			
			# -------------------------------------------
			# Start Generation Looping 
			# -------------------------------------------
			# Begin generation loop
			for gen in xrange(looptime):
				
				# Timing events: start
				start_timeGen = datetime.datetime.now()
									
				# If initial generation - update with initial populations
				if gen == 0:
					SubpopIN = SubpopIN_init
					del SubpopIN_init
					# Use NatalPop 
					sourcePop = 'NatalPop'		
				else: # This was for versions previous v1.37
					#sourcePop = 'ImmiPop'
					sourcePop = 'NatalPop'
				
				# Exit the system if population is 0 or 1
				if sum(N0) <= 1:
					stringout = 'Population went extinct, program ended.'
					logMsg(logfHndl,stringout)
					print('Population went extinct after generation '+str(gen-1)+'.\n')
					break
				
				# ---------------------------------
				# Call CDClimate()
				# ---------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				# Check gen time equal to cdclimgentime
				for icdtime in xrange(len(cdclimgentime)): 
					if gen == int(cdclimgentime[icdtime]):
						tupClimate = DoCDClimate(datadir,icdtime,cdclimgentime,matecdmatfile,dispOutcdmatfile,\
						dispBackcdmatfile,straycdmatfile,matemoveno,dispmoveOutno,dispmoveBackno,StrBackno,matemovethreshval,dispmoveOutthreshval,dispmoveBackthreshval,StrBackthreshval,matemoveparA,matemoveparB,matemoveparC,dispmoveOutparA,dispmoveOutparB,dispmoveOutparC,dispmoveBackparA,dispmoveBackparB,dispmoveBackparC,StrBackparA,StrBackparB,StrBackparC,Mg_pass,Str_pass,Kmu_pass,outsizevals_pass,backsizevals_pass,outgrowdays_pass,backgrowdays_pass,fitvals_pass,popmort_back_pass,popmort_out_pass,eggmort_pass,Kstd_pass,popmort_back_sd_pass,popmort_out_sd_pass,eggmort_sd_pass,outsizevals_sd_pass,backsizevals_sd_pass,outgrowdays_sd_pass,backgrowdays_sd_pass,pop_capture_back_pass,pop_capture_out_pass,cdevolveans,N0_pass,allefreqfiles_pass,classvarsfiles_pass,assortmateModel_pass,assortmateC_pass,subpopmort_pass,PopTag,dispLocalcdmatfile,dispLocalno,dispLocalparA,dispLocalparB,dispLocalparC,dispLocalthreshval,outhabvals_pass,backhabvals_pass)	

						cdmatrix_mate = tupClimate[0]
						cdmatrix_FOut = tupClimate[1]
						cdmatrix_MOut = tupClimate[2]
						cdmatrix_FBack = tupClimate[3]
						cdmatrix_MBack = tupClimate[4]
						cdmatrix_StrBack = tupClimate[5]				
						thresh_mate = tupClimate[6]
						thresh_FOut = tupClimate[7]
						thresh_MOut = tupClimate[8]				
						thresh_FBack = tupClimate[9]
						thresh_MBack = tupClimate[10]
						thresh_Str = tupClimate[11]
						Mg = tupClimate[12]
						Str = tupClimate[13]
						Str_ScaleMin = tupClimate[14]
						Str_ScaleMax = tupClimate[15]
						FdispBack_ScaleMin = tupClimate[16]
						FdispBack_ScaleMax = tupClimate[17]
						MdispBack_ScaleMin = tupClimate[18]
						MdispBack_ScaleMax = tupClimate[19]
						FdispOut_ScaleMin = tupClimate[20]
						FdispOut_ScaleMax = tupClimate[21]
						MdispOut_ScaleMin = tupClimate[22]
						MdispOut_ScaleMax = tupClimate[23]
						mate_ScaleMin = tupClimate[24]
						mate_ScaleMax = tupClimate[25]
						outsizevals_mu = tupClimate[26]
						backsizevals_mu = tupClimate[27]
						outgrowdays_mu = tupClimate[28]
						backgrowdays_mu = tupClimate[29]
						fitvals = tupClimate[30]
						K_mu = tupClimate[31]
						popmort_back_mu = tupClimate[32]
						popmort_out_mu = tupClimate[33]
						eggmort_mu = tupClimate[34]
						K_std = tupClimate[35]
						popmort_back_sd = tupClimate[36]
						popmort_out_sd = tupClimate[37]
						eggmort_sd = tupClimate[38]
						outsizevals_sd = tupClimate[39]
						backsizevals_sd = tupClimate[40]
						outgrowdays_sd = tupClimate[41]
						backgrowdays_sd = tupClimate[42]
						pop_capture_back = tupClimate[43]
						pop_capture_out = tupClimate[44]
						mateno = tupClimate[45]
						FdispOutno = tupClimate[46]
						MdispOutno = tupClimate[47]
						FdispBackno = tupClimate[48]
						MdispBackno = tupClimate[49]
						Strno = tupClimate[50]
						tempN0 = tupClimate[51]
						tempAllelefile = tupClimate[52]
						tempClassVarsfile = tupClimate[53]
						assortmateModel = tupClimate[54]
						assortmateC = tupClimate[55]
						subpopmort_mat = tupClimate[56]
						FdispmoveOutparA = tupClimate[57]
						MdispmoveOutparA = tupClimate[58]
						FdispmoveOutparB = tupClimate[59]
						MdispmoveOutparB = tupClimate[60]
						FdispmoveOutparC = tupClimate[61]
						MdispmoveOutparC = tupClimate[62]
						FdispmoveBackparA = tupClimate[63]
						MdispmoveBackparA = tupClimate[64]
						FdispmoveBackparB = tupClimate[65]
						MdispmoveBackparB = tupClimate[66]
						FdispmoveBackparC = tupClimate[67]
						MdispmoveBackparC = tupClimate[68]
						cdmatrix_dispLocal = tupClimate[69]
						dispLocalparA = tupClimate[70]
						dispLocalparB = tupClimate[71]
						dispLocalparC = tupClimate[72]
						thresh_dispLocal = tupClimate[73]
						dispLocal_ScaleMin = tupClimate[74]
						dispLocal_ScaleMax = tupClimate[75]
						outhabvals = tupClimate[76]
						backhabvals = tupClimate[77]
						
						# ----------------------------------------
						# Introduce new individuals
						# ----------------------------------------
						if (gen != 0 and len(N0_pass[0].split('|')) > 1):
							SubpopIN = AddIndividuals(SubpopIN,tempN0,tempAllelefile,tempClassVarsfile,datadir,loci,alleles,sizeans,cdinfect,cdevolveans,burningen_cdevolve,fitvals,eggFreq,Fmat_set,Mmat_set,Fmat_int,Fmat_slope,Mmat_int,Mmat_slope,dtype,N0,natal,gen,PopTag,sexans,YYmat_set,YYmat_slope,YYmat_int,defaultAgeMature)			
				
				# -------------------------------------------
				# Update stochastic parameters each year here
				# -------------------------------------------
				tupStoch = DoStochasticUpdate(K_mu,K_std,popmort_back_mu,popmort_back_sd,popmort_out_mu,popmort_out_sd,eggmort_mu,eggmort_sd,outsizevals_mu,outsizevals_sd,backsizevals_mu,backsizevals_sd,outgrowdays_mu,outgrowdays_sd,backgrowdays_mu,backgrowdays_sd,age_percmort_out_mu,age_percmort_out_sd,age_percmort_back_mu,age_percmort_back_sd,size_percmort_out_mu,size_percmort_out_sd,size_percmort_back_mu,size_percmort_back_sd,egg_percmort_mu,egg_percmort_sd,cor_mat)
				K = tupStoch[0]
				popmort_back = tupStoch[1]
				popmort_out = tupStoch[2]
				eggmort_patch = tupStoch[3]
				outsizevals = tupStoch[4]
				backsizevals = tupStoch[5]
				outgrowdays = tupStoch[6]
				backgrowdays = tupStoch[7]
				age_percmort_out = tupStoch[8] # sex ratio check throughout
				age_percmort_back = tupStoch[9] # sex ratio check throughout
				size_percmort_out = tupStoch[10] # sex ratio check throughout
				size_percmort_back = tupStoch[11] # sex ratio check throughout
				eggmort_pop = tupStoch[12]
				
				# Print to log
				stringout = 'DoCDClimate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoCDClimate(): ',str(datetime.datetime.now() -start_time1),''								
				# ---------------------------------------
				# Call DoMate()
				# ---------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()				
				
				Bearpairs = DoMate(SubpopIN,K,\
				freplace,mreplace,mateno,thresh_mate,\
				cdmatrix_mate,Track_MateDistCD,xgridpop,\
				ygridpop,Track_MateDistCDstd,Track_FAvgMate,Track_MAvgMate,Track_FSDMate,Track_MSDMate,Track_BreedEvents,gen,sourcePop,mate_ScaleMax,mate_ScaleMin,matemoveparA,matemoveparB,matemoveparC,Femalepercent_egg,eggFreq,sexans,selfing,assortmateC,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,assortmateModel,subpopmort_mat,Track_BreedFemales,Track_BreedMales,Track_BreedYYMales,Track_MatureCount, Track_ImmatureCount,Track_ToTFemales,Track_ToTMales,Track_ToTYYMales)
				
				# Print to log
				stringout = 'DoMate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoMate(): ',str(datetime.datetime.now() -start_time1),''
				
				# Error statement to stop if no females or males
				if Track_ToTFemales[gen][0]==0 or (Track_ToTMales[gen][0] + Track_ToTYYMales[gen][0])==0:			
					print('There are no more females or males left in population after year '+str(gen)+'.\n')
					break
						
				# ---------------------------------------
				# Call DoOffspring()
				# ---------------------------------------
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()			
				
				noOffspring,Bearpairs = DoOffspring(offno,Bearpairs,\
				Track_Births,transmissionprob,gen,K,sourcePop,\
				age_mu,age_sigma,sizeans,\
				egg_mean_1,egg_mean_2,egg_mean_ans,equalClutch,dtype,eggmort_patch,Track_EggDeaths,eggmort_pop,Track_BirthsYY)
							
				# Print to log
				stringout = 'DoOffspring(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoOffspring(): ',str(datetime.datetime.now() -start_time1),''				
							
				# ----------------------------------------------------------------
				# Call 2nd DoUpdate() - grow, age/mature (selection option),egglay,capture, output ind.csv file;no Age0s; ind.csv
				# ----------------------------------------------------------------
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				SubpopIN = DoUpdate(packans,SubpopIN,K,xgridpop,ygridpop,gen,nthfile,ithmcrundir,loci,alleles,logfHndl,'Middle',growans,cdevolveans,defaultAgeMature,fitvals,burningen_cdevolve,age_capture_back,pop_capture_back,Track_CaptureCount_Back,Track_CaptureCount_ClassBack,sizeans,age_size_mean,Track_N_back_age,eggFreq,backsizevals,sizeLoo,sizeR0,size_eqn_1,size_eqn_2,size_eqn_3,backgrowdays,sourcePop,plasticans,burningen_plastic,timeplastic,geneswap,backhabvals)
												
				# Print to log
				stringout = 'Second DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'Second DoUpdate(): ',str(datetime.datetime.now() -start_time1),''
				
				# ------------------------------------------
				# Call DoEmigration()
				# ------------------------------------------
				#pdb.set_trace()				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				SubpopIN = DoEmigration(SubpopIN,K,FdispOutno,\
				MdispOutno,cdmatrix_FOut,cdmatrix_MOut,gen,xgridpop,ygridpop,F_EmiDist,M_EmiDist,cdevolveans,fitvals,F_EmiDist_sd,M_EmiDist_sd,subpopemigration,\
				SelectionDeathsEmi,DisperseDeathsEmi,burningen_cdevolve,Mg,\
				MgSuccess,AdultNoMg,sum(alleles),age_Mg,thresh_FOut,thresh_MOut,N_Emigration_pop,sourcePop,dtype,setmigrate,sizeans,age_size_mean,PackingDeathsEmi,N_Emigration_age,loci,muterate,mtdna,mutationans,FdispOut_ScaleMax,FdispOut_ScaleMin,MdispOut_ScaleMax,MdispOut_ScaleMin,FdispmoveOutparA,FdispmoveOutparB,FdispmoveOutparC,MdispmoveOutparA,MdispmoveOutparB,MdispmoveOutparC,packans,PackingDeathsEmiAge,ithmcrundir,packpar1,timecdevolve,age_percmort_out,migrate,outsizevals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,plasticans,burningen_plastic,timeplastic)
				
				# Print to log
				stringout = 'DoEmigration(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoEmigration(): ',str(datetime.datetime.now() -start_time1),''					
						
				# ----------------------------------------
				# Call DoMortality() - when 'Out'
				# ----------------------------------------			
				#pdb.set_trace()
				start_time1 = datetime.datetime.now() # Timing events: start
				
				SubpopIN = DoMortality(SubpopIN,K,PopDeathsOUT,\
				popmort_out,age_percmort_out,\
				gen,N_EmiMortality,AgeDeathsOUT,sizeans,age_size_mean,size_percmort_out,SizeDeathsOUT,constMortans,packans,'OUT')
				
				# Print to log
				stringout = 'DoOutMortality(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoOutMortality(): ',str(datetime.datetime.now() -start_time1),''	
				
				# ----------------------------------------------------
				# Call DoUpdate() - grow, capture, and optional output indSample.csv
				# ----------------------------------------------------
				
				start_time1 = datetime.datetime.now() # Timing events: start
				
				SubpopIN = DoUpdate(packans,SubpopIN,K,xgridpop,ygridpop,gen,nthfile,ithmcrundir,loci,alleles,logfHndl,gridsample,growans,cdevolveans,defaultAgeMature,fitvals,burningen_cdevolve,age_capture_out,pop_capture_out,Track_CaptureCount_Out,Track_CaptureCount_ClassOut,sizeans,age_size_mean,Track_N_out_age,eggFreq,outsizevals,sizeLoo,sizeR0,size_eqn_1,size_eqn_2,size_eqn_3,outgrowdays,'EmiPop',plasticans,burningen_plastic,timeplastic,geneswap,outhabvals,age_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,YYmat_int,YYmat_slope,YYmat_set)
					
				# Print to log
				stringout = 'Third DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'Third DoUpdate(): ',str(datetime.datetime.now() -start_time1),''								
				#pdb.set_trace()
				# ------------------------------------------
				# Call DoImmigration()
				# ------------------------------------------			
				#pdb.set_trace()
				start_time1 = datetime.datetime.now() # Timing events: start
						
				SubpopIN = DoImmigration(SubpopIN,K,N0,natal,FdispBackno,\
				MdispBackno,cdmatrix_FBack,cdmatrix_MBack,gen,\
				xgridpop,ygridpop,cdevolveans,fitvals,subpopimmigration,\
				SelectionDeathsImm,DisperseDeathsImm,burningen_cdevolve,Str,\
				StrSuccess,\
				Strno,cdmatrix_StrBack,age_S,thresh_FBack,thresh_MBack,thresh_Str,N_Immigration_pop,dtype,sizeans,age_size_mean,PackingDeathsImm,N_Immigration_age,FdispBack_ScaleMax,FdispBack_ScaleMin,MdispBack_ScaleMax,MdispBack_ScaleMin,FdispmoveBackparA,FdispmoveBackparB,FdispmoveBackparC,MdispmoveBackparA,MdispmoveBackparB,MdispmoveBackparC,Str_ScaleMax,Str_ScaleMin,StrBackparA,StrBackparB,StrBackparC,packans,PackingDeathsImmAge,ithmcrundir,packpar1,noOffspring,Bearpairs,age_size_std,Femalepercent_egg,sourcePop,transmissionprob,age_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,loci,muterate,mtdna,mutationans,geneswap,allelst,homeattempt,timecdevolve,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,backsizevals,assortmateModel,inheritans_classfiles,PopTag,subpopmort_mat,eggFreq,sexans,YYmat_slope,YYmat_int,YYmat_set,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,cdmatrix_dispLocal,dispLocalparA,dispLocalparB,dispLocalparC,thresh_dispLocal,dispLocal_ScaleMin,dispLocal_ScaleMax,dispLocalno,alleles,plasticans,burningen_plastic,timeplastic)
				del(Bearpairs)
				# Print to log
				stringout = 'DoImmigration(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoImmigration(): ',str(datetime.datetime.now() -start_time1),''					
				#pdb.set_trace()
				# ------------------------------------------
				# Call DoMortality() - when 'Back'
				# ------------------------------------------
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				SubpopIN = DoMortality(SubpopIN,K,PopDeathsIN,\
				popmort_back,age_percmort_back,\
				gen,N_ImmiMortality,AgeDeathsIN,sizeans,age_size_mean,size_percmort_back,SizeDeathsIN,
				constMortans,packans,'BACK')
				
				# Print to log
				stringout = 'DoInMortality(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoInMortality(): ',str(datetime.datetime.now() -start_time1),''
				
				# ---------------------------------
				# Call GetMetrics()
				# ---------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				GetMetrics(SubpopIN,K,Track_N_Init_pop,Track_K,loci,alleles,gen+1,Track_Ho,Track_Alleles,Track_He,Track_p1,Track_p2,Track_q1,Track_q2,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,Track_N_Init_age,sizeans,age_size_mean,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,sexans,packans,RDispersers,IDispersers)
				
				# Print to log
				stringout = 'GetMetrics(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'GetMetrics(): ',str(datetime.datetime.now() -start_time1),''

				# Error statement here in case no females or males, then break
				if Track_N_Init_pop[gen+1][0] == 0:
					print('Population went extinct after year'+str(gen)+'.\n')
					break
								
				# Print to log
				stringout = 'End Generation/Year Loop'+str(gen)+': '+str(datetime.datetime.now() -start_timeGen) + '\n'
				logMsg(logfHndl,stringout)
				print 'End Generation/Year Loop',str(gen),': ',str(datetime.datetime.now() -start_timeGen),'\n'
			# End::generation loop
						
			# ------------------------------------------
			# Call DoPostProcess()
			# ------------------------------------------
			
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			DoPostProcess(ithmcrundir,loci,alleles,looptime,\
			Track_ToTFemales,Track_ToTMales,Track_BreedFemales,Track_BreedMales,Track_Births,PopDeathsIN,\
			PopDeathsOUT,Track_Alleles,Track_He,Track_Ho,Track_MateDistCD,Track_MateDistCDstd,nthfile,logfHndl,\
			Track_p1,Track_p2,Track_q1,Track_q2,subpopemigration,\
			subpopimmigration,Track_FAvgMate,Track_MAvgMate,Track_FSDMate,Track_MSDMate,\
			SelectionDeathsEmi,SelectionDeathsImm,\
			DisperseDeathsEmi,DisperseDeathsImm,\
			Track_BreedEvents,gridformat,\
			MgSuccess,AdultNoMg,StrSuccess,\
			Track_EggDeaths,Track_K,Track_N_Init_pop,N_Emigration_pop,N_EmiMortality,N_Immigration_pop,N_ImmiMortality,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,PackingDeathsEmi,PackingDeathsImm,Track_N_Init_age,N_Emigration_age,N_Immigration_age,AgeDeathsOUT,AgeDeathsIN,PackingDeathsEmiAge,PackingDeathsImmAge,Track_MatureCount,Track_ImmatureCount,Track_N_back_age,Track_N_out_age,outputans,gen,Track_CaptureCount_Back,Track_CaptureCount_ClassBack,Track_CaptureCount_Out,Track_CaptureCount_ClassOut,age_size_mean,sizeans,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,SizeDeathsOUT,SizeDeathsIN,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,F_EmiDist,M_EmiDist,F_EmiDist_sd,M_EmiDist_sd,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,Track_ToTYYMales,Track_BreedYYMales,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,RDispersers,IDispersers,Track_BirthsYY)
			
			# Print to log
			stringout = 'DoPostProcess(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print 'DoPostProcess(): ',str(datetime.datetime.now() -start_time1),''
				
			# Print to log
			stringout = 'End Monte Carlo Loop'+str(ithmcrun)+': '+str(datetime.datetime.now() -start_timeMC) + '\n'
			logMsg(logfHndl,stringout)
			print 'End Monte Carlo Loop',str(ithmcrun),': ',str(datetime.datetime.now() -start_timeMC),'\n'
			
		# End::Monte Carlo Loop
		
		# Print to log
		stringout = 'End Batch Loop'+str(ibatch)+': '+str(datetime.datetime.now() -start_timeB) + '\n'
		logMsg(logfHndl,stringout)
		print 'End Batch Loop',str(ibatch),': ',str(datetime.datetime.now() -start_timeB),'\n'
	#End::Batch Loop
	
# End::Main Loop	
# Print to log
stringout = 'Total CDmetaPOP Simulation Time: '+str(datetime.datetime.now() -start_time) + ''
logMsg(logfHndl,stringout)
logfHndl.close()
print 'Total CDmetaPOP Simulation Time: ',str(datetime.datetime.now() -start_time),''