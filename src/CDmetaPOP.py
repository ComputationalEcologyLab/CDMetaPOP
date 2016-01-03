# CDmetaPOP.py
# Author: Erin L Landguth
# Created: February 2008
# v 1.0 Release: MARCH 2014
# ----------------------------------------------------------------------------
# General CDmetaPOP information
appName = "CDmetaPOP"
appVers = "version 0.99.07"
appRele = "2015.06.10-09:19:01MDT"
authorNames = "Erin L Landguth"

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
import datetime,time,pdb,os,sys,shutil,gc,multiprocessing
	
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
except ImportError:
	raise ImportError, "CDmetaPOP_Modules required."	
				
#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 
if __name__ == '__main__':
		
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
		print "User must specify data directory, input file name, and output file directory (e.g., at command line type CDmetaPOP.py ../CDmetaPOP_data/ inputvariables.csv exampleout_foldername)."
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
		dispOutcdmatfile = batchVars['dispout_cdmat'][ibatch]
		dispBackcdmatfile = batchVars['dispback_cdmat'][ibatch]
		straycdmatfile = batchVars['stray_cdmat'][ibatch]
		matemoveno = batchVars['matemoveno'][ibatch]
		matemoveparA = batchVars['matemoveparA'][ibatch]
		matemoveparB = batchVars['matemoveparB'][ibatch]
		matemoveparC = batchVars['matemoveparC'][ibatch]
		matemovethreshval = batchVars['matemovethresh'][ibatch]
		freplace = batchVars['Freplace'][ibatch]
		mreplace = batchVars['Mreplace'][ibatch]
		selfing = batchVars['selfans'][ibatch]
		sexans = batchVars['sexans'][ibatch]
		FdispmoveOutno = batchVars['FdispmoveOutno'][ibatch]
		FdispmoveOutparA = batchVars['FdispmoveOutparA'][ibatch]
		FdispmoveOutparB = batchVars['FdispmoveOutparB'][ibatch]
		FdispmoveOutparC = batchVars['FdispmoveOutparC'][ibatch]
		FdispmoveOutthreshval = batchVars['FdispmoveOutthresh'][ibatch]
		MdispmoveOutno = batchVars['MdispmoveOutno'][ibatch]
		MdispmoveOutparA = batchVars['MdispmoveOutparA'][ibatch]
		MdispmoveOutparB = batchVars['MdispmoveOutparB'][ibatch]
		MdispmoveOutparC = batchVars['MdispmoveOutparC'][ibatch]
		MdispmoveOutthreshval = batchVars['MdispmoveOutthresh'][ibatch]
		FdispmoveBackno = batchVars['FdispmoveBackno'][ibatch]
		FdispmoveBackparA = batchVars['FdispmoveBackparA'][ibatch]
		FdispmoveBackparB = batchVars['FdispmoveBackparB'][ibatch]
		FdispmoveBackparC = batchVars['FdispmoveBackparC'][ibatch]
		FdispmoveBackthreshval = batchVars['FdispmoveBackthresh'][ibatch]
		MdispmoveBackno = batchVars['MdispmoveBackno'][ibatch]
		MdispmoveBackparA = batchVars['MdispmoveBackparA'][ibatch]
		MdispmoveBackparB = batchVars['MdispmoveBackparB'][ibatch]
		MdispmoveBackparC = batchVars['MdispmoveBackparC'][ibatch]
		MdispmoveBackthreshval = batchVars['MdispmoveBackthresh'][ibatch]
		StrBackno = batchVars['StrayBackno'][ibatch]
		StrBackparA = batchVars['StrayBackparA'][ibatch]
		StrBackparB = batchVars['StrayBackparB'][ibatch]
		StrBackparC = batchVars['StrayBackparC'][ibatch]
		StrBackthreshval = batchVars['StrayBackthresh'][ibatch]
		homeattempt = batchVars['HomeAttempt'][ibatch]
		offno = batchVars['offno'][ibatch]
		equalClutch = batchVars['equalClutchSize'][ibatch]
		eggFreq = float(batchVars['eggFrequency'][ibatch])
		gendmatans = batchVars['gendmatans'][ibatch]
		gridformat = batchVars['gridformat'][ibatch]
		gridsample = batchVars['gridsampling'][ibatch]
		muterate = float(batchVars['muterate'][ibatch])
		mutationans = batchVars['mutationtype'][ibatch]
		loci = int(batchVars['loci'][ibatch])
		alleles = int(batchVars['alleles'][ibatch])*np.ones(loci,int)
		mtdna = batchVars['mtdna'][ibatch]
		geneswap = int(batchVars['startGenes'][ibatch])
		cdevolveans = batchVars['cdevolveans'][ibatch]
		burningen = int(batchVars['startSelection'][ibatch])
		timecdevolve = batchVars['implementSelection'][ibatch]
		cdinfect = batchVars['cdinfect'][ibatch]
		transmissionprob = float(batchVars['transmissionprob'][ibatch])
		growans = batchVars['growth_option'][ibatch]
		sizeLoo = float(batchVars['growth_Loo'][ibatch])
		sizeR0 = float(batchVars['growth_R0'][ibatch])
		size_eqn_1 = float(batchVars['growth_temp_max'][ibatch])
		size_eqn_2 = float(batchVars['growth_temp_CV'][ibatch])
		size_eqn_3 = float(batchVars['growth_temp_t0'][ibatch])
		sizeans = batchVars['sizecontrol'][ibatch]
		Mmat_set = batchVars['mature_length_male'][ibatch]
		Fmat_set = batchVars['mature_length_female'][ibatch]
		Mmat_slope = float(batchVars['mature_slope_male'][ibatch])
		Mmat_int = float(batchVars['mature_int_male'][ibatch])
		Fmat_slope = float(batchVars['mature_slope_female'][ibatch])
		Fmat_int = float(batchVars['mature_int_female'][ibatch])
		egg_mean_ans = batchVars['Egg_Mean_ans'][ibatch]
		egg_mean_1 = float(batchVars['Egg_Mean_par1'][ibatch])
		egg_mean_2 = float(batchVars['Egg_Mean_par2'][ibatch])
		egg_percmort_mu = float(batchVars['Egg_Mortality'][ibatch])
		egg_percmort_sd = float(batchVars['Egg_Mortality_StDev'][ibatch])
		Femalepercent_egg = batchVars['Egg_FemalePercent'][ibatch]
		packans = batchVars['popmodel'][ibatch]
		packpar1 = float(batchVars['popmodel_par1'][ibatch])
		cor_mat_ans = batchVars['correlation_matrix'][ibatch]
				
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
		
		# ---------------------------------
		# Some Error checking
		# ---------------------------------
				
		# Constant mortality checks
		if not (constMortans == '1' or constMortans == '2'):
			print('Constant mortalities are compounded using option 1 or 2 specifiy correct values. If no constant mortalities are entered, then enter 1.')
			sys.exit(-1)
		
		# Check on cdevolve answer input
		if not (cdevolveans == '1' or cdevolveans == '2' or cdevolveans == '1_mat' or cdevolveans == '2_mat' or cdevolveans == 'N' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link'):
			print('CDEVOLVE answer either N, 1, 2, M, G, MG_ind, MG_link, 1_mat, or 2_mat.')
			sys.exit(-1)
			
		# For mature and size ans
		if (cdevolveans == 'M' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link') and sizeans == 'N':
			print('CDEVOLVE answer is M and size answer must be Y.')
			sys.exit(-1)
		
		# If cdevolve is turned on must have 2 alleles
		if cdevolveans != 'N' and alleles[0] != 2:
			print('Warning: More than 2 alleles per locus specified. CDEVOLVE only considers first 2 alleles in selection models.')
		# Must have more than 1 loci
		if loci <= 1:
			print('Currently, CDmetaPOP needs more than 1 locus to run.')
			sys.exit(-1)
		if cdevolveans == '1' or cdevolveans == '2' or cdevolveans == '1_mat' or cdevolveans == '2_mat':
			if not (timecdevolve != 'Out' or timecdevolve != 'Back' or timecdevolve != 'Both'):
				print('CDEVOLVE timing must be specified (e.g., Out, Back or Both).')
				sys.exit(-1)
			
		# Error check on forward mutation in A and backward mutation in B
		#	Can only happen if cdevolve == 2.
		if mutationans == 'forwardAbackwardBrandomN' and (cdevolveans != '2' or cdevolveans == '2_mat'):
			print('This special case of mutation is for AAbb ancestors and 2-locus selection.')
			sys.exit(-1)		
			
		# Gd matrix answers
		if (gendmatans == 'N' or gendmatans == 'Dps' or gendmatans == 'braycurtis' or gendmatans == 'Da') == False:
			print('Genetic distance matrix output not an option.')
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
		if cdevolveans != 'N' and burningen < geneswap:
			print('Warning: Burnin time < time at which genetic exchange is to initialize, setting burnin time = start genetic exchange time.')
			burningen = geneswap
			
		# Egg frequency
		if eggFreq > 1:
			print('Egg frequency must be less than or equal to 1.')
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
			Track_MatureCount = []
			Track_ImmatureCount = []
			
			# DoMate()
			Track_FAvgMate = []
			Track_MAvgMate = []
			Track_FSDMate = []
			Track_MSDMate = []
			Track_MateDistCD = []
			Track_MateDistCDstd = []
			Track_ToTFemales = []
			Track_ToTMales = []
			Track_BreedFemales = []
			Track_BreedMales = []
			Track_BreedEvents = []
			
			# DoOffspring
			Track_Births = []
			Track_EggDeaths = []
			
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
			subpopemigration,subpopimmigration,sizeans,geneswap,eggFreq,Fmat_set,Mmat_set,Fmat_int,Fmat_slope,Mmat_int,Mmat_slope,burningen,cor_mat_ans)
			
			ithmcrundir = tupPreProcess[0]			
			fitvals_pass = tupPreProcess[1]
			allelst = tupPreProcess[2]
			subpopemigration = tupPreProcess[3]
			subpopimmigration = tupPreProcess[4]
			age_percmort_out_mu = tupPreProcess[5]
			age_percmort_back_mu = tupPreProcess[6]
			age_Mg = tupPreProcess[7]
			age_S = tupPreProcess[8]			
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
			M_mature = tupPreProcess[26]
			F_mature = tupPreProcess[27]
			age_sigma = tupPreProcess[28]
			outgrowdays_pass = tupPreProcess[29]
			backgrowdays_pass = tupPreProcess[30]
			Kmu_pass = tupPreProcess[31]
			age_capture_out = tupPreProcess[32]
			age_capture_back = tupPreProcess[33]
			Kstd_pass = tupPreProcess[34]
			K_std = tupPreProcess[35]
			popmort_out_sd_pass = tupPreProcess[36]
			popmort_back_sd_pass = tupPreProcess[37]
			eggmort_sd_pass = tupPreProcess[38]
			outsizevals_sd_pass = tupPreProcess[39]
			backsizevals_sd_pass = tupPreProcess[40]
			outgrowdays_sd_pass = tupPreProcess[41]
			backgrowdays_sd_pass = tupPreProcess[42]
			size_percmort_out_mu = tupPreProcess[43]
			size_percmort_back_mu = tupPreProcess[44]
			age_percmort_out_sd = tupPreProcess[45]
			age_percmort_back_sd = tupPreProcess[46]
			size_percmort_out_sd = tupPreProcess[47]
			size_percmort_back_sd = tupPreProcess[48]
			pop_capture_back_pass = tupPreProcess[49]
			pop_capture_out_pass = tupPreProcess[50]
			pop_capture_back = tupPreProcess[51]
			natal = tupPreProcess[52]
			cor_mat = tupPreProcess[53]
			migrate = tupPreProcess[54]
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
			
			GetMetrics(SubpopIN_init,K,Track_N_Init_pop,Track_K,loci,alleles,0,Track_Ho,Track_Alleles,Track_He,Track_p1,Track_p2,Track_q1,Track_q2,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,Track_ToTMales,Track_ToTFemales,Track_BreedMales,Track_BreedFemales,Track_N_Init_age,Track_MatureCount,Track_ImmatureCount,sizeans,age_size_mean,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,sexans)
			
			# Print to log
			stringout = 'GetMetrics() Initial: '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print 'GetMetrics() Initial: ',str(datetime.datetime.now() -start_time1),''
			
			# ---------------------------------
			# Error statements
			# ---------------------------------			
			# Error statement here in case no females or males, then break
			if Track_ToTFemales[0][0]==0 or Track_ToTMales[0][0]==0:						
				print('There are no females or males to begin time loop.\n')
				break
				
			# ------------------------------------------
			# Call DoUpdate() - output initial file here ind-1.csv
			# ------------------------------------------	
							
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			DoUpdate(SubpopIN_init,K,xgridpop,ygridpop,-1,nthfile,ithmcrundir,loci,alleles,logfHndl,'Initial','N','N',[],burningen,[],[],[],[],[],[])
			
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
				else:
					sourcePop = 'ImmiPop'				
				
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
						dispBackcdmatfile,straycdmatfile,matemoveno,FdispmoveOutno,MdispmoveOutno,FdispmoveBackno,MdispmoveBackno,StrBackno,matemovethreshval,FdispmoveOutthreshval,MdispmoveOutthreshval,FdispmoveBackthreshval,MdispmoveBackthreshval,StrBackthreshval,\
						matemoveparA,matemoveparB,matemoveparC,FdispmoveOutparA,\
						FdispmoveOutparB,FdispmoveOutparC,MdispmoveOutparA,MdispmoveOutparB,MdispmoveOutparC,\
						FdispmoveBackparA,\
						FdispmoveBackparB,FdispmoveBackparC,MdispmoveBackparA,MdispmoveBackparB,MdispmoveBackparC,StrBackparA,StrBackparB,StrBackparC,Mg_pass,Str_pass,Kmu_pass,outsizevals_pass,backsizevals_pass,outgrowdays_pass,backgrowdays_pass,fitvals_pass,popmort_back_pass,popmort_out_pass,eggmort_pass,Kstd_pass,popmort_back_sd_pass,popmort_out_sd_pass,eggmort_sd_pass,outsizevals_sd_pass,backsizevals_sd_pass,outgrowdays_sd_pass,backgrowdays_sd_pass,pop_capture_back_pass,pop_capture_out_pass,cdevolveans)	

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
				age_percmort_out = tupStoch[8]
				age_percmort_back = tupStoch[9]
				size_percmort_out = tupStoch[10]
				size_percmort_back = tupStoch[11]
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
				ygridpop,Track_MateDistCDstd,Track_FAvgMate,Track_MAvgMate,Track_FSDMate,Track_MSDMate,Track_BreedEvents,gen,sourcePop,dtype,mate_ScaleMax,mate_ScaleMin,matemoveparA,matemoveparB,matemoveparC,Femalepercent_egg,eggFreq,sexans,selfing)
				
				# Print to log
				stringout = 'DoMate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoMate(): ',str(datetime.datetime.now() -start_time1),''							
				
				# ---------------------------------------
				# Call DoOffspring()
				# ---------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()			
				
				noOffspring,Bearpairs = DoOffspring(offno,Bearpairs,\
				Track_Births,transmissionprob,gen,K,sourcePop,\
				age_mu,age_sigma,sizeans,\
				egg_mean_1,egg_mean_2,egg_mean_ans,equalClutch,dtype,Mmat_set,Fmat_set,eggmort_patch,Track_EggDeaths,eggmort_pop)
							
				# Print to log
				stringout = 'DoOffspring(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoOffspring(): ',str(datetime.datetime.now() -start_time1),''				
								
				# ----------------------------------------------------------------
				# Call 2nd DoUpdate() - grow, age/mature (selection option),egglay,capture, output ind.csv file;no Age0s; ind.csv
				# ----------------------------------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				SubpopIN = DoUpdate(SubpopIN,K,xgridpop,ygridpop,gen,nthfile,ithmcrundir,loci,alleles,logfHndl,'Middle',growans,cdevolveans,fitvals,burningen,age_capture_back,pop_capture_back,Track_CaptureCount_Back,Track_CaptureCount_ClassBack,sizeans,age_size_mean,Track_N_back_age,eggFreq,backsizevals,sizeLoo,sizeR0,size_eqn_1,size_eqn_2,size_eqn_3,backgrowdays,sourcePop,sizeans,M_mature,F_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set)
												
				# Print to log
				stringout = 'Second DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'Second DoUpdate(): ',str(datetime.datetime.now() -start_time1),''					
				
				# ------------------------------------------
				# Call DoEmigration()
				# ------------------------------------------		
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				SubpopIN = DoEmigration(SubpopIN,K,FdispOutno,\
				MdispOutno,cdmatrix_FOut,cdmatrix_MOut,gen,xgridpop,ygridpop,F_EmiDist,M_EmiDist,cdevolveans,fitvals,F_EmiDist_sd,M_EmiDist_sd,subpopemigration,\
				SelectionDeathsEmi,DisperseDeathsEmi,\
				burningen,Mg,\
				MgSuccess,AdultNoMg,sum(alleles),age_Mg,thresh_FOut,thresh_MOut,N_Emigration_pop,sourcePop,dtype,setmigrate,sizeans,age_size_mean,PackingDeathsEmi,N_Emigration_age,loci,muterate,mtdna,mutationans,FdispOut_ScaleMax,FdispOut_ScaleMin,MdispOut_ScaleMax,MdispOut_ScaleMin,FdispmoveOutparA,FdispmoveOutparB,FdispmoveOutparC,MdispmoveOutparA,MdispmoveOutparB,MdispmoveOutparC,packans,PackingDeathsEmiAge,ithmcrundir,packpar1,timecdevolve,age_percmort_out,migrate)
				
				# Print to log
				stringout = 'DoEmigration(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoEmigration(): ',str(datetime.datetime.now() -start_time1),''					
								
				# ----------------------------------------
				# Call DoMortality()
				# ----------------------------------------			
				start_time1 = datetime.datetime.now() # Timing events: start
				
				SubpopIN = DoMortality(SubpopIN,K,PopDeathsOUT,\
				popmort_out,age_percmort_out,\
				gen,N_EmiMortality,AgeDeathsOUT,sizeans,age_size_mean,size_percmort_out,SizeDeathsOUT,constMortans,packans)
				
				# Print to log
				stringout = 'DoOutMortality(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoOutMortality(): ',str(datetime.datetime.now() -start_time1),''	
				
				# ----------------------------------------------------
				# Call DoUpdate() - grow, capture, and optional output indSample.csv
				# ----------------------------------------------------
				start_time1 = datetime.datetime.now() # Timing events: start
				
				SubpopIN = DoUpdate(SubpopIN,K,xgridpop,ygridpop,gen,nthfile,ithmcrundir,loci,alleles,logfHndl,gridsample,growans,'N',[],burningen,age_capture_out,pop_capture_out,Track_CaptureCount_Out,Track_CaptureCount_ClassOut,sizeans,age_size_mean,Track_N_out_age,eggFreq,outsizevals,sizeLoo,sizeR0,size_eqn_1,size_eqn_2,size_eqn_3,outgrowdays,'EmiPop')
				
				# Print to log
				stringout = 'Third DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'Third DoUpdate(): ',str(datetime.datetime.now() -start_time1),''								
				
				# ------------------------------------------
				# Call DoImmigration()
				# ------------------------------------------			
				start_time1 = datetime.datetime.now() # Timing events: start
						
				SubpopIN = DoImmigration(SubpopIN,K,N0,natal,FdispBackno,\
				MdispBackno,cdmatrix_FBack,cdmatrix_MBack,gen,\
				xgridpop,ygridpop,cdevolveans,fitvals,subpopimmigration,\
				SelectionDeathsImm,DisperseDeathsImm,burningen,Str,\
				StrSuccess,\
				Strno,cdmatrix_StrBack,age_S,thresh_FBack,thresh_MBack,thresh_Str,N_Immigration_pop,dtype,sizeans,age_size_mean,PackingDeathsImm,N_Immigration_age,FdispBack_ScaleMax,FdispBack_ScaleMin,MdispBack_ScaleMax,MdispBack_ScaleMin,FdispmoveBackparA,FdispmoveBackparB,FdispmoveBackparC,MdispmoveBackparA,MdispmoveBackparB,MdispmoveBackparC,Str_ScaleMax,Str_ScaleMin,StrBackparA,StrBackparB,StrBackparC,packans,PackingDeathsImmAge,ithmcrundir,packpar1,noOffspring,Bearpairs,age_size_std,Femalepercent_egg,sourcePop,transmissionprob,M_mature,F_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,loci,muterate,mtdna,mutationans,geneswap,allelst,homeattempt,timecdevolve,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd)
				del(Bearpairs)
				# Print to log
				stringout = 'DoImmigration(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoImmigration(): ',str(datetime.datetime.now() -start_time1),''					
				
				# ------------------------------------------
				# Call DoMortality()
				# ------------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				SubpopIN = DoMortality(SubpopIN,K,PopDeathsIN,\
				popmort_back,age_percmort_back,\
				gen,N_ImmiMortality,AgeDeathsIN,sizeans,age_size_mean,size_percmort_back,SizeDeathsIN,
				constMortans,packans)
				
				# Print to log
				stringout = 'DoInMortality(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoInMortality(): ',str(datetime.datetime.now() -start_time1),''
				
				# ---------------------------------
				# Call GetMetrics()
				# ---------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				GetMetrics(SubpopIN,K,Track_N_Init_pop,Track_K,loci,alleles,gen+1,Track_Ho,Track_Alleles,Track_He,Track_p1,Track_p2,Track_q1,Track_q2,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,Track_ToTMales,Track_ToTFemales,Track_BreedMales,Track_BreedFemales,Track_N_Init_age,Track_MatureCount,Track_ImmatureCount,sizeans,age_size_mean,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,sexans)
				
				# Print to log
				stringout = 'GetMetrics(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'GetMetrics(): ',str(datetime.datetime.now() -start_time1),''

				# Error statement here in case no females or males, then break
				if Track_N_Init_pop[gen+1][0] == 0:
					print('Population went extinct after year'+str(gen)+'.\n')
					break
				if Track_ToTFemales[gen+1][0]==0 or Track_ToTMales[gen+1][0]==0:			
					print('There are no more females or males left in population after year '+str(gen)+'.\n')
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
			
			DoPostProcess(ithmcrundir,gendmatans,loci,alleles,looptime,\
			Track_ToTFemales,Track_ToTMales,Track_BreedFemales,Track_BreedMales,Track_Births,PopDeathsIN,\
			PopDeathsOUT,Track_Alleles,Track_He,Track_Ho,Track_MateDistCD,Track_MateDistCDstd,nthfile,logfHndl,\
			Track_p1,Track_p2,Track_q1,Track_q2,subpopemigration,\
			subpopimmigration,Track_FAvgMate,Track_MAvgMate,Track_FSDMate,Track_MSDMate,\
			SelectionDeathsEmi,SelectionDeathsImm,\
			DisperseDeathsEmi,DisperseDeathsImm,\
			Track_BreedEvents,gridformat,\
			MgSuccess,AdultNoMg,StrSuccess,\
			Track_EggDeaths,Track_K,Track_N_Init_pop,N_Emigration_pop,N_EmiMortality,N_Immigration_pop,N_ImmiMortality,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,PackingDeathsEmi,PackingDeathsImm,Track_N_Init_age,N_Emigration_age,N_Immigration_age,AgeDeathsOUT,AgeDeathsIN,PackingDeathsEmiAge,PackingDeathsImmAge,Track_MatureCount,Track_ImmatureCount,Track_N_back_age,Track_N_out_age,outputans,gen,Track_CaptureCount_Back,Track_CaptureCount_ClassBack,Track_CaptureCount_Out,Track_CaptureCount_ClassOut,age_size_mean,sizeans,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,SizeDeathsOUT,SizeDeathsIN,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,F_EmiDist,M_EmiDist,F_EmiDist_sd,M_EmiDist_sd)
			
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