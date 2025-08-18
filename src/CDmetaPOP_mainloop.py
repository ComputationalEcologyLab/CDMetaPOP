# -------------------------------------------------------------------------------------------------
# CDmetaPOP_mainloop.py
# Author: Erin L Landguth, Casey Day
# Created: October 2020
# Description: This is the function/module file for target execution of main loop in parallel runs
# --------------------------------------------------------------------------------------------------

# Python specific functions
import numpy as np 
import pdb, random, copy, os, sys,datetime,signal
from ast import literal_eval 
from CDmetaPOP_Modules import * 
from CDmetaPOP_PostProcess import *
from CDmetaPOP_PreProcess import *
from CDmetaPOP_Mate import *
from CDmetaPOP_Emigration import *
from CDmetaPOP_Immigration import *
from CDmetaPOP_Offspring2 import *
from CDmetaPOP_Mortality import *

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

# --------------------------------------------------------------------------
def logMsg(outf,msg):
	'''
	logMsg() --log file message handler.
	Inputs:
	outf - open file handle
	msg -- string containing formatted message
	--always outputs to log file by default.
	--using msgVerbose, can be set to "Tee" output to stdout as well
	'''
	outf.write(msg+ '\n')
	if msgVerbose:
		print(("%s"%(msg)))
		
	# End::logMsg()

# --------------------------------------------------------------------------------------------------------------------
def main_loop(spcNO,fileans,irun,datadir,sizeans,constMortans,mcruns,looptime,nthfile_out,gridformat,gridsample,outputans,cdclimgentimelist,outdir,startcomp,implementcomp,passlogfHndl,XQs, nspecies, extinctQ, global_extinctQ,current_system_pid):
	'''
	Main loop here
	'''
	#if multiprocessing.current_process().name == "S1":
	#	pdb.set_trace()
	# Open the logfHndl file
	#passlogfHndl = outdir+"CDmetaPOP"+str(spcNO)+".log"
	logfHndl = open(passlogfHndl,'a')
	
	# Call function and store inputvariables
	batchVars,batchVarsIndex,nSimulations = loadFile(fileans,1,',',True)

	# ----------------------------------------	
	# Begin Batch Looping 
	# ----------------------------------------
	# This loop is defined by the number of rows in PopVars.csv
	for ibatch in range(nSimulations):
	
		# Store all information and the type of each, also do some error checks 
		xyfilename = datadir+batchVars['xyfilename'][ibatch]
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
		sexchromo = int(batchVars['sex_chromo'][ibatch])
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
		muterate_pass = batchVars['muterate'][ibatch]
		mutationans = batchVars['mutationtype'][ibatch]
		loci = int(batchVars['loci'][ibatch])
		alleles = batchVars['alleles'][ibatch]
		mtdna = batchVars['mtdna'][ibatch]
		geneswap = int(batchVars['startGenes'][ibatch])
		cdevolveans = batchVars['cdevolveans'][ibatch]
		burningen_cdevolve = int(batchVars['startSelection'][ibatch])
		timecdevolve = batchVars['implementSelection'][ibatch]
		betaFile_selection = batchVars['betaFile_selection'][ibatch]
		
		plasticans = batchVars['plasticgeneans'][ibatch]
		plastic_signalresp_pass = batchVars['plasticSignalResponse'][ibatch]
		plastic_behaviorresp_pass = batchVars['plasticBehavioralResponse'][ibatch]
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
		#mat_set = batchVars['mature_length_set'][ibatch]
		mat_slope = batchVars['mature_eqn_slope'][ibatch]
		mat_int = batchVars['mature_eqn_int'][ibatch]
		eggFreq_mu = float(batchVars['Egg_Freq_Mean'][ibatch])
		eggFreq_sd = float(batchVars['Egg_Freq_StDev'][ibatch])
		egg_mean_ans = batchVars['Egg_Mean_ans'][ibatch]
		egg_mean_1 = float(batchVars['Egg_Mean_par1'][ibatch])
		egg_mean_2 = float(batchVars['Egg_Mean_par2'][ibatch])
		egg_percmort_mu = batchVars['Egg_Mortality'][ibatch]
		egg_percmort_sd = batchVars['Egg_Mortality_StDev'][ibatch]
		Femalepercent_egg = batchVars['Egg_FemaleProb'][ibatch]
		packans = batchVars['popmodel'][ibatch]
		packpar1 = float(batchVars['popmodel_par1'][ibatch])
		cor_mat_ans = batchVars['correlation_matrix'][ibatch]
		defaultMature = batchVars['mature_default'][ibatch]
		subpopmort_pass = batchVars['subpopmort_file'][ibatch]
		egg_delay = int(batchVars['egg_delay'][ibatch])
		egg_add = batchVars['egg_add'][ibatch]
			
		# -------------------------------
		# Distill some vars
		# -------------------------------
		# Grab the nthfile list range specific to user input, list or sequence
		#pdb.set_trace()
		if not isinstance(nthfile_out, (list,tuple)):
			nthfile_out = int(nthfile_out)
			if nthfile_out != 0:
				nthfile = list(range(0,looptime+nthfile_out,nthfile_out))
				del(nthfile[-1]) # Delete the last value 0, looptime - 1
			else:
				nthfile = [0]
		# If specified years with |
		else:
			nthfile = []
			# Split up list, removing space values, and appending to nthfile
			for inum in range(len(nthfile_out)):
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
		if sexchromo not in [2,3,4]:
			print('Number of sex chromosome options must be 2,3, or 4.')
			sys.exit()
		# Check Deterministic mature set value either age or size
		tupVal = sexsplit(defaultMature, sexchromo)
		FXXmat_set = tupVal[0]
		MXYmat_set = tupVal[1]
		MYYmat_set = tupVal[2]
		FYYmat_set = tupVal[3]
		
		# Logistic equation for maturation as a function of size - slope
		tupVal = sexsplit(mat_slope, sexchromo)
		FXXmat_slope = tupVal[0]
		MXYmat_slope = tupVal[1]
		MYYmat_slope = tupVal[2]
		FYYmat_slope = tupVal[3]
		
		# Logistic equation for maturation as a function of size - intercept
		tupVal = sexsplit(mat_int, sexchromo)
		FXXmat_int= tupVal[0]
		MXYmat_int= tupVal[1]
		MYYmat_int= tupVal[2]
		FYYmat_int= tupVal[3]
		
		# ---------------------------------
		# Some Error checking
		# ---------------------------------
		
		# DoEmigration() skipped and warning for if selection on
		if cdevolveans != 'N':
			if timecdevolve.find('Out') != -1:
				if dispOutcdmatfile == 'N':
					stringout = 'Warning: DoEmigration module skipped and spatial selection during this time frame specified, which will also be skipped.'
					logMsg(logfHndl,stringout)
						# 	
		# Constant mortality checks
		if not (constMortans == '1' or constMortans == '2'):
			print('Constant mortalities are compounded using option 1 or 2 specifiy correct values. If no constant mortalities are entered, then enter 1.')
			sys.exit(-1)
			
		# Warning check on multiplicative mortlaity	
		if constMortans == '2':
<<<<<<< Updated upstream
			print('Warning: Using multiplicative mortality option, double check for patch and class values when 0%')
=======
			logMsg(logfHndl,'Warning: Using multiplicative mortality option, double check for patch and class values when 0%')
				
		# Check on cdevolve answer input
		valid_values = ['1', '2', '1_mat', '2_mat', 'N', 'M', 'G', 'MG_ind', 'MG_link', 'stray', '1_G_ind', '1_G_link', 'Hindex', 'F', 'Plastic', 'Multilocus','runtiming']
		validate(cdevolveans not in valid_values, 'Check CDEvolve answer options.')
>>>>>>> Stashed changes
		
		# Check on cdevolve answer input
		if not (cdevolveans == '1' or cdevolveans == '2' or cdevolveans == '1_mat' or cdevolveans == '2_mat' or cdevolveans == 'N' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == 'stray' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans.split('_')[0] == 'F' or cdevolveans.split('_')[0] == 'Hindex' or cdevolveans.split('_')[0] == 'P' or cdevolveans.split('_')[0] == 'FHindex'):
			print('CDEVOLVE answer either N, 1, 2, M, G, MG_ind, MG_link, 1_mat, 2_mat, stray, 1_G_ind, 1_G_link, Hindex, F, Plastic, or Multilocus.')
			sys.exit(-1)
			
		# For mature and size ans
<<<<<<< Updated upstream
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
					print('CDEVOLVE answer is Hindex and 3 parameters for the Parabolic function must be specified, see user manual and example files.')
					sys.exit(-1)
			elif cdevolveans.split('_')[1] == 'Step':
				if len(cdevolveans.split('_')[2].split(':')) != 3:
					print('CDEVOLVE answer is Hindex and 3 parameters for the Step function must be specified, see user manual and example files.')
					sys.exit(-1)
			elif cdevolveans.split('_')[1] == 'Linear':
				if len(cdevolveans.split('_')[2].split(':')) != 2:
					print('CDEVOLVE answer is Hindex and 2 parameters for the Linear function must be specified, see user manual and example files.')
					sys.exit(-1)					
			else:
				print('CDEVOLVE and Hindex parameter not entered correctly, check user manual and example files.')
				sys.exit(-1)
		
=======
		valid_values = ['M', 'MG_ind', 'MG_link', 'G', '1_G_ind', '1_G_link']
		validate(cdevolveans in valid_values and sizeans != 'Y','CDEVOLVE answer is M or G and size answer must be Y.')
				
		# For Hindex answer and each function
		validate('Hindex' in cdevolveans and 'Gauss' in cdevolveans and len(cdevolveans.split('_')[2].split(':')) != 6, 'CDEVOLVE answer is Hindex and 6 parameters for the Gaussian function must be specified, see user manual and example files.')
		validate('Hindex' in cdevolveans and 'Para' in cdevolveans and len(cdevolveans.split('_')[2].split(':')) != 3, 'CDEVOLVE answer is Hindex and 3 parameters for the Parabolic function must be specified, see user manual and example files.')
		validate('Hindex' in cdevolveans and 'Step' in cdevolveans and len(cdevolveans.split('_')[2].split(':')) != 3, 'CDEVOLVE answer is Hindex and 3 parameters for the Step function must be specified, see user manual and example files.')
		validate('Hindex' in cdevolveans and 'Linear' in cdevolveans and len(cdevolveans.split('_')[2].split(':')) != 3, 'CDEVOLVE answer is Hindex and 2 parameters for the Linear function must be specified, see user manual and example files.')
				
>>>>>>> Stashed changes
		# If cdevolve is turned on must have 2 alleles
		if cdevolveans != 'N' and alleles[0] != 2:
			print('Warning: More than 2 alleles per locus specified. CDEVOLVE only considers first 2 alleles in selection models (except Hindex scenario).')
		if plasticans != 'N' and alleles[0] != 2:
			print('Warning: More than 2 alleles per locus specified. Plastic gene turned on and only considers first 2 alleles in this model.')
		
		# For Plastic answer
		if plasticans != 'N':
			# Split for temp
			if ((plasticans.split('_')[0] != 'Temp') and (plasticans.split('_')[0] != 'Hab')):
				print('Plastic type (Temp/Hab) not entered corectly, check user manual and example files.')
				sys.exit(-1)
				
		if plasticans != 'N':
			# Split for temp
			if ((plasticans.split('_')[1] != 'dom') and (plasticans.split('_')[1] != 'rec') and (plasticans.split('_')[1] != 'codom')):
				print('Plastic type (dom/codom/rec) not entered corectly, check user manual and example files.')
				sys.exit(-1)
				
		if plasticans != 'N':	
			if ((timeplastic.find('Out') == -1) and (timeplastic.find('Back') == -1)):
				print('Plastic timing must be specified (e.g., Out or Back).')
				sys.exit(-1)
				
		# Must have more than 1 loci
		if loci <= 1:
			print('Currently, CDmetaPOP needs more than 1 locus to run.')
			sys.exit(-1)
		if cdevolveans == '1' or cdevolveans == '2' or cdevolveans == '1_mat' or cdevolveans == '2_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans.split('_')[0] == 'Hindex' or cdevolveans.split('_')[0] == 'FHindex':
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
			stringout = 'Warning: Selection burnin time < time at which genetic exchange is to initialize, setting burnin time = start genetic exchange time.'
			logMsg(logfHndl,stringout)
			
			burningen_cdevolve = geneswap
		if plasticans != 'N' and burningen_plastic < geneswap:
			stringout = 'Warning: Selection burnin time < time at which genetic exchange is to initialize, setting burnin time = start genetic exchange time.'
			logMsg(logfHndl,stringout)
			burningen_plastic = geneswap
			
		# Inherit answer can only be:
		if not (inheritans_classfiles == 'random' or inheritans_classfiles == 'Hindex' or inheritans_classfiles == 'mother'):
			print('Inherit answer for multiple class files is not correct: enter either random or Hindex.')
			sys.exit(-1)
		
		# If inherit answer uses Hindex, mutation can't be on
		if len(muterate_pass) > 1:
			if sum(np.asarray(muterate_pass,dtype=float)) != 0.0 and (inheritans_classfiles == 'Hindex' or inheritans_classfiles == 'mother'):
				print('Mutation is not operating with Hindex inheritance options in this version.')
				sys.exit(-1)
		else:
			if sum(np.asarray([muterate_pass],dtype=float)) != 0.0 and (inheritans_classfiles == 'Hindex' or inheritans_classfiles == 'mother'):
				print('Mutation is not operating with Hindex inheritance options in this version.')
				sys.exit(-1)
			
		# If egg_delay is gretter than 1 v2.68 turned off for further testing
<<<<<<< Updated upstream
		if egg_delay > 0:
			print('Currently, egg delay is not operating beyond 0 year/time unit. Testing sourcePop and gen = 0 issue.')
			sys.exit(-1)
		
		# Reproduction answers checks - check with Erin		
		if (sexans == 'N' or sexans == 'Y' or sexans == 'H') == False:
			print('Reproduction choices either N, Y or H, check user manual.')
			sys.exit(-1)
		if sexans == 'H' and (selfing == 'N' or selfing == 'Y'):
			print('Hermaphroditic mating structure specified - H - then must specify the selfing probability.')
			sys.exit(-1)
		'''if sexans == 'H':
			try:
				float(selfing)
			except ValueError:
				# Currently this causes a crash but I couldn't get try/except to work right.
				raise ValueError('Sexans "H" for hermaphroditic system must be paired with a probability for selfing (selfans).') from None
				print('Sexans "H" for hermaphroditic system must be paired with a probability for selfing (selfans).')
				sys.exit(-1)'''
				
		# Spot to add age 0s / eggs / pups /liter, etc
		if (egg_add == 'mating' or egg_add == 'nonmating') == False:
			print('Egg add choices either mating or nonmating, check user manual.')
			sys.exit(-1)
		if egg_add == 'nonmating' and packans != 'logistic':
			print('Logistic model should be used with nonmating add age 0 locations. See user manual')
			sys.exit(-1)
			
=======
		validate(egg_delay > 0,'Currently, egg delay is not operating beyond 0 year/time unit. Testing sourcePop and gen = 0 issue.')
					
		# Reproduction answers checks 	
		validate((sexans == 'N' or sexans == 'Y' or sexans == 'H') == False,'Reproduction choices either N, Y or H, check user manual.')
		validate(sexans == 'H' and (selfing == 'N' or selfing == 'Y'),'Hermaphroditic mating structure specified - H - then must specify the selfing probability.')
									
		# Spot to add age 0s / eggs / pups /liter, etc
		validate(egg_add not in ['mating', 'nonmating'], 'Egg add choice incorrect.')
		#validate(egg_add not in ['mating', 'nonmating'] and packans.split('_')[0] not in 'logistic', 'Nonmating option should use logistic model.')
		validate(egg_add == 'nonmating' and packans != 'logistic_back', 'Nonmating option must be used with logistic in the back location.')
					
		valid_values = ['N','Both','Back','Out']
		validate(implementdisease not in valid_values, 'Implement disease value incorrect.')
		validate(implementdisease!='N' and alleles[0] != 2,'More than 2 alleles per locus specified with disease defense options.')
		
		valid_values = ['N','Both','Back','Out']
		validate(implementcomp not in valid_values, 'Implement disease value incorrect.')
		'''
		# Error check here for runtiming and 4 mats
		if len(dispBackcdmatfile) > 1:
			validate(cdevolveans == 'runtiming' and len(dispBackcdmatfile[0].split('~')[0].split(';')) != 4, 'runtiming CDEvolve answer specified, 4 migrate back cdmats required.')
		else:
			validate(cdevolveans == 'runtiming' and len(dispBackcdmatfile.split('~')[0].split(';')) != 4, 'runtiming CDEvolve answer specified, 4 migrate back cdmats required.')
		'''	
>>>>>>> Stashed changes
		# ---------------------------------------------	
		# Begin Monte-Carlo Looping
		# ---------------------------------------------
		
		# range(mcruns) is typically 10 - 50...and it takes a long time.
		for ithmcrun in range(mcruns):	
			
			# Timing events: start
			start_timeMC = datetime.datetime.now()
			# Keep track so extinction message are printed only once
			temp_extinct = 0
		
			# -----------------------------------------
			# Create storage variables
			# ------------------------------------------	
			# These variables will be stored in output.csv at the end of the simulation						
			
			# GetMetrics()
			Track_p1, Track_p2, Track_q1, Track_q2, Track_Alleles, Track_He, Track_Ho, Track_N_Init_pop, Track_N_Init_age, Track_N_Init_class, Track_K, Track_CaptureCount_Out, Track_CaptureCount_ClassOut, Track_CaptureCount_Back, Track_CaptureCount_ClassBack, maxfit, minfit = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
									
			# DoMate()
			Track_FAvgMate, Track_MAvgMate, Track_FSDMate, Track_MSDMate, Track_MateDistCD, Track_MateDistCDstd, Track_BreedEvents, Track_AAaaMates, Track_AAAAMates, Track_aaaaMates, Track_AAAaMates, Track_aaAaMates, Track_AaAaMates, Track_BreedFemales, Track_BreedMales, Track_BreedYYMales, Track_BreedYYFemales, Track_MatureCount, Track_ImmatureCount, Track_ToTFemales, Track_ToTMales, Track_ToTYYMales, Track_ToTYYFemales = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
			
			# DoOffspring
			Track_Births, Track_EggDeaths, Track_BirthsMYY, Track_BirthsFYY = [], [], [], []
			
			# DoUpdate
			Track_N_back_age, Track_N_out_age = [], []
			
			# Emigration()
			N_Emigration_pop, N_Emigration_age, subpopemigration, F_EmiDist, M_EmiDist, F_EmiDist_sd, M_EmiDist_sd, SelectionDeathsEmi, DisperseDeathsEmi, PackingDeathsEmi, PackingDeathsEmiAge, MgSuccess, AdultNoMg, Track_YYSelectionPackDeathsEmi, Track_WildSelectionPackDeathsEmi, SelectionDeaths_Age0s, N_beforePack_pop, N_beforePack_age, Track_KadjEmi = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
			
			# Mortlity after Emigration
			N_EmiMortality, PopDeathsOUT, AgeDeathsOUT, SizeDeathsOUT = [], [], [], []
			
			# Immigration
			N_Immigration_pop, N_Immigration_age, subpopimmigration, F_HomeDist, M_HomeDist, F_HomeDist_sd, M_HomeDist_sd, F_StrayDist, M_StrayDist, F_StrayDist_sd, M_StrayDist_sd, F_ZtrayDist, M_ZtrayDist, F_ZtrayDist_sd, M_ZtrayDist_sd, SelectionDeathsImm, DisperseDeathsImm, PackingDeathsImmAge, PackingDeathsImm, StrSuccess, Track_YYSelectionPackDeathsImmi, Track_WildSelectionPackDeathsImmi, Track_KadjImmi = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

			# Mortality after immigration
			N_ImmiMortality, PopDeathsIN, AgeDeathsIN, SizeDeathsIN = [], [], [], []
			
			# DoOutput()
			Infected, Residors, Strayers1, Strayers2, Immigrators, IDispersers, RDispersers, PopSizes_Mean, PopSizes_Std, AgeSizes_Mean, AgeSizes_Std, ClassSizes_Mean, ClassSizes_Std = [], [], [], [], [], [], [], [], [], [], [], [], []
			
			# Non-tracking variables - create empty 2-D list
			noOffspring_temp = [np.asarray([]),np.asarray([])] 
			Bearpairs_temp = [[[-9999,-9999]],[[-9999,-9999]]]
					
			# ------------------------------------	
			# Call DoPreProcess()
			# ------------------------------------
			
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			# Call function
			tupPreProcess = DoPreProcess(outdir,datadir,irun,ithmcrun,\
			xyfilename,loci,alleles,0,logfHndl,cdevolveans,cdinfect,\
			subpopemigration,subpopimmigration,sizeans,burningen_cdevolve,cor_mat_ans,inheritans_classfiles,sexans,spcNO,ibatch,betaFile_selection,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,eggFreq_mu,eggFreq_sd)
			ithmcrundir = tupPreProcess[0]			
			fitvals_pass = tupPreProcess[1] # sex ratio check throughout
			allelst = tupPreProcess[2]
			subpopemigration = tupPreProcess[3]
			subpopimmigration = tupPreProcess[4]
			age_size_mean = tupPreProcess[5]
			age_size_std = tupPreProcess[6]			
			xgridpop = tupPreProcess[7]
			ygridpop = tupPreProcess[8]			
			SubpopIN_init = tupPreProcess[9]
			N0 = tupPreProcess[10]
			K_mu = tupPreProcess[11]
			dtype = tupPreProcess[12]
			outsizevals_pass = tupPreProcess[13]
			backsizevals_pass = tupPreProcess[14]
			popmort_out_pass = tupPreProcess[15]
			popmort_back_pass = tupPreProcess[16]
			MgOut_patch_pass = tupPreProcess[17]
			Str_patch_pass = tupPreProcess[18]
			eggmort_pass = tupPreProcess[19]
			setmigrate = tupPreProcess[20]			
			outgrowdays_pass = tupPreProcess[21]
			backgrowdays_pass = tupPreProcess[22]
			Kmu_pass = tupPreProcess[23]			
			Kstd_pass = tupPreProcess[24]
			K_std = tupPreProcess[25]
			popmort_out_sd_pass = tupPreProcess[26]
			popmort_back_sd_pass = tupPreProcess[27]
			eggmort_sd_pass = tupPreProcess[28]
			outsizevals_sd_pass = tupPreProcess[29]
			backsizevals_sd_pass = tupPreProcess[30]
			outgrowdays_sd_pass = tupPreProcess[31]
			backgrowdays_sd_pass = tupPreProcess[32]			
			pop_capture_back_pass = tupPreProcess[33]
			pop_capture_out_pass = tupPreProcess[34]
			pop_capture_back = tupPreProcess[35]
			natal_patches = tupPreProcess[36]
			cor_mat = tupPreProcess[37]
			migrate_patches = tupPreProcess[38]
			N0_pass = tupPreProcess[39]
			allefreqfiles_pass = tupPreProcess[40]
			classvarsfiles_pass = tupPreProcess[41]
			PopTag = tupPreProcess[42] # To pass into AddIndividuals, Emigration, Immigration
			comp_coef_pass = tupPreProcess[43]
			xvars_betas_pass = tupPreProcess[44]
			tempbetas_selection = tupPreProcess[45]
			outhabvals_pass = tupPreProcess[46]
			backhabvals_pass = tupPreProcess[47]
			
			#disperse_patches = tupPreProcess[48] # Pass into DoImmigration / 
			#stray_patches = tupPreProcess[49] # Pass into DoImmigration
			
			MgBack_patch_prob_pass = tupPreProcess[48] # CDClimate then DoImmigrtion
			Disperse_patch_prob_pass = tupPreProcess[49] # Cdclimate thenDoImmigration
			
			# Grab first one only
			K = K_mu # Initialize K with mu		
						
			# Print to log
			stringout = 'DoPreProcess(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			if multiprocessing.current_process().name == "S0" or multiprocessing.current_process().name == "MainProcess":
				print(('DoPreProcess(): ',str(datetime.datetime.now() -start_time1),''))
			
			# ---------------------------------
			# Call GetMetrics()
			# ---------------------------------
			
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			GetMetrics(SubpopIN_init,K,Track_N_Init_pop,Track_K,loci,alleles,0,Track_Ho,Track_Alleles,Track_He,Track_p1,Track_p2,Track_q1,Track_q2,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,Track_N_Init_age,sizeans,age_size_mean,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,packans,RDispersers,IDispersers,xvars_betas_pass,tempbetas_selection,maxfit,minfit,cdevolveans)
			# Print to log
			stringout = 'GetMetrics() Initial: '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
				
			# ---------------------------------
			# Error statements
			# ---------------------------------			
			# Error statement here in case no females or males, then break
			if Track_N_Init_pop[0][0] == 0:
				print('There are no individuals to begin time loop.\n')
				
				#ThisSystem = psutil.Process(os.getpid())
				#ThisSystem.terminate()
				#os.killpg(current_system_pid, signal.SIGTERM)
				
				sys.exit(-1)
				
			# ----------------------------------------------------
			# Call DoUpdate() - output initial file here ind-1.csv
			# ----------------------------------------------------
						
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			DoUpdate(packans,SubpopIN_init,K,xgridpop,ygridpop,-1,nthfile,ithmcrundir,loci,alleles,logfHndl,'Initial')
			# Print to log
			stringout = 'DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			
			# -------------------------------------------
			# Start Generation Looping 
			# -------------------------------------------
			# Begin generation loop
			for gen in range(looptime):
				
				# Timing events: start
				start_timeGen = datetime.datetime.now()
									
				# If initial generation - update with initial populations
				if gen == 0:
					SubpopIN = SubpopIN_init
					del SubpopIN_init
					# Use NatalPop in first generation 
					sourcePop = 'NatalPop'		
				else: # This was for versions previous v1.37
					sourcePop = 'ImmiPop' # For generations greater than 0, this is where the individual is located.
					#sourcePop = 'NatalPop'
				
				# Exit the system if population is 0 or 1
				checkPopN = [len(SubpopIN[x]) for x in range(0,len(SubpopIN))] 
				if sum(checkPopN) == 0:
					stringout = 'Species is extinct.'
					logMsg(logfHndl,stringout)
					# If not extinct from previous generation
					if temp_extinct == 0:					
						print(('Species ' + str(spcNO) + ' went extinct.'))
						temp_extinct = 1
					# Track extinctions
					extinctQ.put(0)
				else:
					# Track extinctions
					extinctQ.put(1)
				# List to track extinctions
				if multiprocessing.current_process().name == "S0" or multiprocessing.current_process().name == "MainProcess":
					ext_list = []	
					for ispecies in range(nspecies):
						ext_list.append(extinctQ.get(block=True))
					# If all species extinct, exit
					if sum(ext_list) == 0:
						# Exit the system if population is 0 or 1
						stringout = 'All species extinct after generation '+str(gen-1)+', program ended.\n'
						logMsg(logfHndl,stringout)
						if multiprocessing.current_process().name == "S0" or multiprocessing.current_process().name == "MainProcess":
							print('All species extinct')
						#print('Population went extinct after generation '+str(gen-1)+'.\n')
						for ispecies in range(nspecies):
							global_extinctQ.put(1)
					else:
						for ispecies in range(nspecies):
							global_extinctQ.put(0)	
				if global_extinctQ.get() == 1:
					break
				#pdb.set_trace()
				# ---------------------------------
				# Call CDClimate()
				# ---------------------------------			
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				# Check gen time equal to cdclimgentime
				for icdtime in range(len(cdclimgentime)): 
					if gen == int(cdclimgentime[icdtime]):
						tupClimate = DoCDClimate(datadir,icdtime,cdclimgentime,matecdmatfile,dispOutcdmatfile,\
						dispBackcdmatfile,straycdmatfile,matemoveno,dispmoveOutno,dispmoveBackno,StrBackno,matemovethreshval,dispmoveOutthreshval,dispmoveBackthreshval,StrBackthreshval,matemoveparA,matemoveparB,matemoveparC,dispmoveOutparA,dispmoveOutparB,dispmoveOutparC,dispmoveBackparA,dispmoveBackparB,dispmoveBackparC,StrBackparA,StrBackparB,StrBackparC,MgOut_patch_pass,Str_patch_pass,Kmu_pass,outsizevals_pass,backsizevals_pass,outgrowdays_pass,backgrowdays_pass,fitvals_pass,popmort_back_pass,popmort_out_pass,eggmort_pass,Kstd_pass,popmort_back_sd_pass,popmort_out_sd_pass,eggmort_sd_pass,outsizevals_sd_pass,backsizevals_sd_pass,outgrowdays_sd_pass,backgrowdays_sd_pass,pop_capture_back_pass,pop_capture_out_pass,cdevolveans,N0_pass,allefreqfiles_pass,classvarsfiles_pass,assortmateModel_pass,assortmateC_pass,subpopmort_pass,PopTag,dispLocalcdmatfile,dispLocalno,dispLocalparA,dispLocalparB,dispLocalparC,dispLocalthreshval,comp_coef_pass,betaFile_selection,xvars_betas_pass,outhabvals_pass,backhabvals_pass,plastic_signalresp_pass,plastic_behaviorresp_pass,plasticans,muterate_pass,sexchromo,MgBack_patch_prob_pass,Disperse_patch_prob_pass)
						
						# Cdmatrix values
						cdmatrix_mate,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,cdmatrix_FXXBack,cdmatrix_MXYBack, cdmatrix_MYYBack, cdmatrix_FYYBack, cdmatrix_FXXStr, cdmatrix_MXYStr, cdmatrix_MYYStr, cdmatrix_FYYStr, cdmatrix_FXXLD, 	cdmatrix_MXYLD, cdmatrix_MYYLD, cdmatrix_FYYLD = tupClimate[0],tupClimate[1],tupClimate[2],tupClimate[3],tupClimate[4], tupClimate[5],tupClimate[6],tupClimate[7],tupClimate[8],tupClimate[9],tupClimate[10],tupClimate[11],tupClimate[12],tupClimate[13],tupClimate[14],tupClimate[15],tupClimate[16]						
						# Threshold values
						thresh_mate, thresh_FXXOut, thresh_MXYOut, thresh_MYYOut, thresh_FYYOut, thresh_FXXBack, thresh_MXYBack, thresh_MYYBack, thresh_FYYBack, thresh_FXXStr, thresh_MXYStr, thresh_MYYStr, thresh_FYYStr, thresh_FXXLD, thresh_MXYLD, thresh_MYYLD, thresh_FYYLD = tupClimate[17], tupClimate[18], tupClimate[19],tupClimate[20],tupClimate[21], tupClimate[22],tupClimate[23],tupClimate[24],tupClimate[25],tupClimate[26],tupClimate[27],tupClimate[28],tupClimate[29],tupClimate[30],tupClimate[31],tupClimate[32],tupClimate[33]						
						# Scale Min
						scalemin_mate,scalemin_FXXOut,scalemin_MXYOut,scalemin_MYYOut,scalemin_FYYOut,scalemin_FXXBack,scalemin_MXYBack,scalemin_MYYBack,scalemin_FYYBack,scalemin_FXXStr,scalemin_MXYStr,scalemin_MYYStr,scalemin_FYYStr,scalemin_FXXLD,scalemin_MXYLD,scalemin_MYYLD,scalemin_FYYLD = tupClimate[34], tupClimate[35], tupClimate[36],tupClimate[37],tupClimate[38], tupClimate[39],tupClimate[40],tupClimate[41],tupClimate[42],tupClimate[43],tupClimate[44],tupClimate[45],tupClimate[46],tupClimate[47],tupClimate[48],tupClimate[49],tupClimate[50]
						# Scale Max
						scalemax_mate, scalemax_FXXOut, scalemax_MXYOut, scalemax_MYYOut,scalemax_FYYOut ,scalemax_FXXBack,scalemax_MXYBack,scalemax_MYYBack,scalemax_FYYBack,scalemax_FXXStr,scalemax_MXYStr,scalemax_MYYStr,scalemax_FYYStr,scalemax_FXXLD,scalemax_MXYLD,scalemax_MYYLD,scalemax_FYYLD = tupClimate[51],tupClimate[52],tupClimate[53],tupClimate[54],tupClimate[55],tupClimate[56],tupClimate[57],tupClimate[58],tupClimate[59],tupClimate[60],tupClimate[61],tupClimate[62],tupClimate[63],tupClimate[64],tupClimate[65],tupClimate[66],tupClimate[67]
						# ParA
						parA_mate,parA_FXXOut,parA_MXYOut,parA_MYYOut,parA_FYYOut,parA_FXXBack,parA_MXYBack,parA_MYYBack,parA_FYYBack,parA_FXXStr,parA_MXYStr,parA_MYYStr,parA_FYYStr,parA_FXXLD,parA_MXYLD,parA_MYYLD,parA_FYYLD = tupClimate[68],tupClimate[69],tupClimate[70],tupClimate[71],tupClimate[72],tupClimate[73],tupClimate[74],tupClimate[75],tupClimate[76],tupClimate[77],tupClimate[78],tupClimate[79],tupClimate[80],tupClimate[81],tupClimate[82],tupClimate[83],tupClimate[84]
						# ParB
						parB_mate,parB_FXXOut,parB_MXYOut,parB_MYYOut,parB_FYYOut,parB_FXXBack,parB_MXYBack,parB_MYYBack,parB_FYYBack,parB_FXXStr,parB_MXYStr,parB_MYYStr,parB_FYYStr,parB_FXXLD,parB_MXYLD,parB_MYYLD,parB_FYYLD = tupClimate[85],tupClimate[86],tupClimate[87],tupClimate[88],tupClimate[89],tupClimate[90],tupClimate[91],tupClimate[92],tupClimate[93],tupClimate[94],tupClimate[95],tupClimate[96],tupClimate[97],tupClimate[98],tupClimate[99],tupClimate[100],tupClimate[101]
						# ParC
						parC_mate,parC_FXXOut,parC_MXYOut,parC_MYYOut,parC_FYYOut,parC_FXXBack,parC_MXYBack,parC_MYYBack,parC_FYYBack,parC_FXXStr,parC_MXYStr,parC_MYYStr,parC_FYYStr,parC_FXXLD,parC_MXYLD,parC_MYYLD,parC_FYYLD = tupClimate[102],tupClimate[103],tupClimate[104],tupClimate[105],tupClimate[106],tupClimate[107],tupClimate[108],tupClimate[109],tupClimate[110],tupClimate[111],tupClimate[112],tupClimate[113],tupClimate[114],tupClimate[115],tupClimate[116],tupClimate[117],tupClimate[118]
						# Movement No
						moveno_mate,moveno_FXXOut,moveno_MXYOut,moveno_MYYOut,moveno_FYYOut,moveno_FXXBack,moveno_MXYBack,moveno_MYYBack,moveno_FYYBack,moveno_FXXStr,moveno_MXYStr,moveno_MYYStr,moveno_FYYStr,moveno_FXXLD,moveno_MXYLD,moveno_MYYLD,moveno_FYYLD = tupClimate[119],tupClimate[120],tupClimate[121],tupClimate[122],tupClimate[123],tupClimate[124],tupClimate[125],tupClimate[126],tupClimate[127],tupClimate[128],tupClimate[129],tupClimate[130],tupClimate[131],tupClimate[132],tupClimate[133],tupClimate[134],tupClimate[135]						
						MgOut_patch_prob = tupClimate[136]
						Str_patch_prob = tupClimate[137]						
						outsizevals_mu = tupClimate[138]
						backsizevals_mu = tupClimate[139]
						outgrowdays_mu = tupClimate[140]
						backgrowdays_mu = tupClimate[141]
						fitvals = tupClimate[142]
						K_mu = tupClimate[143]
						popmort_back_mu = tupClimate[144]
						popmort_out_mu = tupClimate[145]
						eggmort_mu = tupClimate[146]
						K_std = tupClimate[147]
						popmort_back_sd = tupClimate[148]
						popmort_out_sd = tupClimate[149]
						eggmort_sd = tupClimate[150]
						outsizevals_sd = tupClimate[151]
						backsizevals_sd = tupClimate[152]
						outgrowdays_sd = tupClimate[153]
						backgrowdays_sd = tupClimate[154]
						pop_capture_back = tupClimate[155]
						pop_capture_out = tupClimate[156]
						tempN0 = tupClimate[157]
						tempAllelefile = tupClimate[158]
						tempClassVarsfile = tupClimate[159]
						assortmateModel = tupClimate[160]
						assortmateC = tupClimate[161]
						subpopmort_mat = tupClimate[162]
						comp_coef = tupClimate[163]
						betas_selection = tupClimate[164]
						xvars_betas = tupClimate[165]
						outhabvals = tupClimate[166]
						backhabvals = tupClimate[167]
						plastic_signalresp = tupClimate[168] #ts added
						plastic_behaviorresp = tupClimate[169] #ts added, note that plasticans not listed here, but present in tupClimate. in 1.72, present in tupClimate as last piece but not defined down here.
						muterate = tupClimate[170]						
						if gen == 0:
							age_percmort_out_mu,age_percmort_out_sd,age_percmort_back_mu,age_percmort_back_sd,size_percmort_out_mu,size_percmort_out_sd,size_percmort_back_mu,size_percmort_back_sd,age_MgOUT, age_MgBACK,age_S,age_DispProb,age_mature,age_mu,age_sigma,f_leslie_mu,f_leslie_std,age_capture_out,age_capture_back = tupClimate[171],tupClimate[172],tupClimate[173],tupClimate[174],tupClimate[175],tupClimate[176],tupClimate[177],tupClimate[178],tupClimate[179],tupClimate[180],tupClimate[181],tupClimate[182],tupClimate[183],tupClimate[184],tupClimate[185],tupClimate[186],tupClimate[187],tupClimate[188],tupClimate[189]
						MgBack_patch_prob = tupClimate[190]
						Disperse_patch_prob = tupClimate[191]
	
						# ----------------------------------------
						# Introduce new individuals
						# ----------------------------------------
						if (gen != 0 and len(N0_pass[0].split('|')) > 1):							
							SubpopIN = AddIndividuals(SubpopIN,tempN0,tempAllelefile,tempClassVarsfile,datadir,loci,alleles,sizeans,cdinfect,cdevolveans,burningen_cdevolve,fitvals,dtype,N0,natal_patches,gen,PopTag,sexans,logfHndl,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,eggFreq_mu,eggFreq_sd)											
				# -------------------------------------------
				# Update stochastic parameters each year here
				# -------------------------------------------
				tupStoch = DoStochasticUpdate(K_mu,K_std,popmort_back_mu,popmort_back_sd,popmort_out_mu,popmort_out_sd,eggmort_mu,eggmort_sd,outsizevals_mu,outsizevals_sd,backsizevals_mu,backsizevals_sd,outgrowdays_mu,outgrowdays_sd,backgrowdays_mu,backgrowdays_sd,age_percmort_out_mu,age_percmort_out_sd,age_percmort_back_mu,age_percmort_back_sd,size_percmort_out_mu,size_percmort_out_sd,size_percmort_back_mu,size_percmort_back_sd,egg_percmort_mu,egg_percmort_sd,cor_mat,age_mu,age_sigma,f_leslie_mu,f_leslie_std,sexchromo)
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
				f_ind = tupStoch[13]
				f_leslie = tupStoch[14]
				
				# Print to log
				stringout = 'DoCDClimate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				
				# ---------------------------------------
				# Call DoMate() - DoOffspring()
				# ---------------------------------------
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()				
				
				Bearpairs_temp,noOffspring_temp = DoMate(SubpopIN,K,\
				freplace,mreplace,moveno_mate,thresh_mate,\
				cdmatrix_mate,Track_MateDistCD,xgridpop,\
				ygridpop,Track_MateDistCDstd,Track_FAvgMate,Track_MAvgMate,Track_FSDMate,Track_MSDMate,Track_BreedEvents,gen,sourcePop,scalemax_mate,scalemin_mate,parA_mate,parB_mate,parC_mate,Femalepercent_egg,sexans,selfing,assortmateC,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,assortmateModel,subpopmort_mat,Track_BreedFemales,Track_BreedMales,Track_BreedYYMales,Track_BreedYYFemales,Track_MatureCount, Track_ImmatureCount,Track_ToTFemales,Track_ToTMales,Track_ToTYYMales,Track_ToTYYFemales,egg_delay,Bearpairs_temp,natal_patches,offno,transmissionprob,f_ind,age_sigma,sizeans,egg_mean_1,egg_mean_2,egg_mean_ans,equalClutch,dtype,eggmort_patch,Track_EggDeaths,eggmort_pop,noOffspring_temp,Track_Births,Track_BirthsMYY,Track_BirthsFYY,constMortans,outputans)
				
				# Print to log
				stringout = 'DoMate() and DoOffspring: '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				if Track_ToTFemales[gen][0]==0 or (Track_ToTMales[gen][0] + Track_ToTYYMales[gen][0] + Track_ToTYYFemales[gen][0])==0:
					if temp_extinct == 0:
						print(('There are no more females or males left from species ' + str(spcNO) + ' after year '+str(gen)+'.\n'))
					#break
<<<<<<< Updated upstream
						
=======
					
>>>>>>> Stashed changes
				# --------------------------------------------------------------------------------------------------------
				# Call 2nd DoUpdate() - grow, age (selection option),egglay,capture, output ind.csv file;no Age0s; ind.csv
				# --------------------------------------------------------------------------------------------------------
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()
				SubpopIN = DoUpdate(packans,SubpopIN,K,xgridpop,ygridpop,gen,nthfile,ithmcrundir,loci,alleles,logfHndl,'Middle',growans,cdevolveans,fitvals,burningen_cdevolve,age_capture_back,pop_capture_back,Track_CaptureCount_Back,Track_CaptureCount_ClassBack,sizeans,age_size_mean,Track_N_back_age,eggFreq_mu,eggFreq_sd,backsizevals,sizeLoo,sizeR0,size_eqn_1,size_eqn_2,size_eqn_3,backgrowdays,plasticans,burningen_plastic,timeplastic,plastic_signalresp,geneswap,backhabvals,sexchromo)
												
				# Print to log
				stringout = 'Second DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				
				# --------------------------------------------------------
				# Call DoEmigration() - Age0s adding into population here
				# --------------------------------------------------------
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()
<<<<<<< Updated upstream
				SubpopIN = DoEmigration(SubpopIN,K,gen,F_EmiDist,M_EmiDist,cdevolveans,fitvals,F_EmiDist_sd,M_EmiDist_sd,subpopemigration,SelectionDeathsEmi,DisperseDeathsEmi,burningen_cdevolve,MgOut_patch_prob,MgSuccess,AdultNoMg,age_MgOUT,N_Emigration_pop,sourcePop,dtype,setmigrate,sizeans,age_size_mean,PackingDeathsEmi,N_Emigration_age,loci,muterate,mtdna,mutationans,packans,PackingDeathsEmiAge,packpar1,timecdevolve,migrate_patches,outsizevals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,noOffspring_temp,Bearpairs_temp,age_size_std,Femalepercent_egg,transmissionprob,age_mature,alleles,geneswap,allelst,assortmateModel,inheritans_classfiles,eggFreq_mu,eggFreq_sd,sexans,N_beforePack_pop,N_beforePack_age,SelectionDeaths_Age0s,comp_coef,XQs,Track_KadjEmi,Track_KadjImmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,thresh_FXXOut,thresh_MXYOut,thresh_MYYOut,thresh_FYYOut,scalemin_FXXOut,scalemin_MXYOut,scalemin_MYYOut,scalemin_FYYOut,scalemax_FXXOut,scalemax_MXYOut,scalemax_MYYOut,scalemax_FYYOut,parA_FXXOut,parA_MXYOut,parA_MYYOut,parA_FYYOut,parB_FXXOut,parB_MXYOut,parB_MYYOut,parB_FYYOut,parC_FXXOut,parC_MXYOut,parC_MYYOut,parC_FYYOut,moveno_FXXOut,moveno_MXYOut,moveno_MYYOut,moveno_FYYOut,egg_add,outputans)
				
=======

				SubpopIN = DoEmigration(SubpopIN,K,gen,F_EmiDist,M_EmiDist,cdevolveans,fitvals,F_EmiDist_sd,M_EmiDist_sd,subpopemigration,SelectionDeathsEmi,DisperseDeathsEmi,burningen_cdevolve,MgOut_patch_prob,MgSuccess,AdultNoMg,age_MgOUT,N_Emigration_pop,sourcePop,dtype,setmigrate,sizeans,age_size_mean,PackingDeathsEmi,N_Emigration_age,loci,muterate,mtdna,mutationans,packans,PackingDeathsEmiAge,packpar1,timecdevolve,migrate_patches,outsizevals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,noOffspring_temp,Bearpairs_temp,age_size_std,Femalepercent_egg,age_mature,alleles,geneswap,allelst,assortmateModel,inheritans_classfiles,eggFreq_mu,eggFreq_sd,sexans,N_beforePack_pop,N_beforePack_age,SelectionDeaths_Age0s,comp_coef,XQs,Track_KadjEmi,Track_KadjImmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,thresh_FXXOut,thresh_MXYOut,thresh_MYYOut,thresh_FYYOut,scalemin_FXXOut,scalemin_MXYOut,scalemin_MYYOut,scalemin_FYYOut,scalemax_FXXOut,scalemax_MXYOut,scalemax_MYYOut,scalemax_FYYOut,parA_FXXOut,parA_MXYOut,parA_MYYOut,parA_FYYOut,parB_FXXOut,parB_MXYOut,parB_MYYOut,parB_FYYOut,parC_FXXOut,parC_MXYOut,parC_MYYOut,parC_FYYOut,moveno_FXXOut,moveno_MXYOut,moveno_MYYOut,moveno_FYYOut,egg_add,outputans,age_percmort_out, f_leslie,f_leslie_std,disease_vars,Track_DiseaseStates_AddAge0s)
							
>>>>>>> Stashed changes
				# Delete the noOffspring_temp and Bearpairs_temp egg_delay spots used: the first spot in list
				if len(noOffspring_temp) != 0: # But check for extinction
					del(noOffspring_temp[0])
					del(Bearpairs_temp[0])
					# Append a new empty spot for next years cohort
					noOffspring_temp.append(np.asarray([]))
					Bearpairs_temp.append([[-9999,-9999]])				
				
				# Print to log
				stringout = 'DoEmigration(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
						
				# ----------------------------------------
				# Call DoMortality() - when 'Out'
				# ----------------------------------------			
				#pdb.set_trace()
				start_time1 = datetime.datetime.now() # Timing events: start
				SubpopIN = DoMortality(SubpopIN,K,PopDeathsOUT,	popmort_out,age_percmort_out,gen,N_EmiMortality,AgeDeathsOUT,sizeans,age_size_mean,size_percmort_out,SizeDeathsOUT,constMortans,packans,'OUT',sexchromo)
				
				# Print to log
				stringout = 'DoOutMortality(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				
				# ----------------------------------------------------
				# Call DoUpdate() - grow, mature, capture, and optional output indSample.csv
				# ----------------------------------------------------
				#pdb.set_trace()
				start_time1 = datetime.datetime.now() # Timing events: start
				SubpopIN = DoUpdate(packans,SubpopIN,K,xgridpop,ygridpop,gen,nthfile,ithmcrundir,loci,alleles,logfHndl,gridsample,growans,cdevolveans,fitvals,burningen_cdevolve,age_capture_out,pop_capture_out,Track_CaptureCount_Out,Track_CaptureCount_ClassOut,sizeans,age_size_mean,Track_N_out_age,eggFreq_mu,eggFreq_sd,outsizevals,sizeLoo,sizeR0,size_eqn_1,size_eqn_2,size_eqn_3,outgrowdays,plasticans,burningen_plastic,timeplastic,plastic_signalresp,geneswap,outhabvals,sexchromo,age_mature,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope)
			
				# Print to log
				stringout = 'Third DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				
				# ------------------------------------------
				# Call DoImmigration()
				# ------------------------------------------			
				#pdb.set_trace()
				start_time1 = datetime.datetime.now() # Timing events: start
				SubpopIN = DoImmigration(SubpopIN,K,natal_patches,gen,cdevolveans,fitvals,subpopimmigration,SelectionDeathsImm,DisperseDeathsImm,burningen_cdevolve,Str_patch_prob,StrSuccess,age_S,N_Immigration_pop,dtype,sizeans,age_size_mean,PackingDeathsImm,N_Immigration_age,packans,PackingDeathsImmAge,packpar1,homeattempt,timecdevolve,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,backsizevals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,age_percmort_back,comp_coef,XQs,Track_KadjImmi,Track_KadjEmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,f_leslie,f_leslie_std,age_DispProb,cdmatrix_FXXBack,cdmatrix_MXYBack,cdmatrix_MYYBack,cdmatrix_FYYBack,thresh_FXXBack,thresh_MXYBack,thresh_MYYBack,thresh_FYYBack,scalemin_FXXBack,scalemin_MXYBack,scalemin_MYYBack,scalemin_FYYBack,scalemax_FXXBack,scalemax_MXYBack,scalemax_MYYBack,scalemax_FYYBack,parA_FXXBack,parA_MXYBack,parA_MYYBack,parA_FYYBack,parB_FXXBack,parB_MXYBack,parB_MYYBack,parB_FYYBack,parC_FXXBack,parC_MXYBack,parC_MYYBack,parC_FYYBack,moveno_FXXBack,moveno_MXYBack,moveno_MYYBack,moveno_FYYBack,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr,thresh_FXXStr,thresh_MXYStr,thresh_MYYStr,thresh_FYYStr,scalemin_FXXStr,scalemin_MXYStr,scalemin_MYYStr,scalemin_FYYStr,scalemax_FXXStr,scalemax_MXYStr,scalemax_MYYStr,scalemax_FYYStr,parA_FXXStr,parA_MXYStr,parA_MYYStr,parA_FYYStr,parB_FXXStr,parB_MXYStr,parB_MYYStr,parB_FYYStr,parC_FXXStr,parC_MXYStr,parC_MYYStr,parC_FYYStr,moveno_FXXStr,moveno_MXYStr,moveno_MYYStr,moveno_FYYStr,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD,thresh_FXXLD,thresh_MXYLD,thresh_MYYLD,thresh_FYYLD,scalemin_FXXLD,scalemin_MXYLD,scalemin_MYYLD,scalemin_FYYLD,scalemax_FXXLD,scalemax_MXYLD,scalemax_MYYLD,scalemax_FYYLD,parA_FXXLD,parA_MXYLD,parA_MYYLD,parA_FYYLD,parB_FXXLD,parB_MXYLD,parB_MYYLD,parB_FYYLD,parC_FXXLD,parC_MXYLD,parC_MYYLD,parC_FYYLD,moveno_FXXLD,moveno_MXYLD,moveno_MYYLD,moveno_FYYLD,sexchromo,age_MgBACK,MgBack_patch_prob,Disperse_patch_prob,MgOut_patch_prob,age_MgOUT,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,migrate_patches,egg_add,outputans)
								
				# Print to log
				stringout = 'DoImmigration(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
							
				# ------------------------------------------
				# Call DoMortality() - when 'Back'
				# ------------------------------------------
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()
				SubpopIN = DoMortality(SubpopIN,K,PopDeathsIN,popmort_back,age_percmort_back,gen,N_ImmiMortality,AgeDeathsIN,sizeans,age_size_mean,size_percmort_back,SizeDeathsIN,constMortans,packans,'BACK',sexchromo)
				
				# Print to log
				stringout = 'DoInMortality(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				
				# ---------------------------------
				# Call GetMetrics()
				# ---------------------------------
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()
				GetMetrics(SubpopIN,K,Track_N_Init_pop,Track_K,loci,alleles,gen+1,Track_Ho,Track_Alleles,Track_He,Track_p1,Track_p2,Track_q1,Track_q2,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,Track_N_Init_age,sizeans,age_size_mean,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,packans,RDispersers,IDispersers,xvars_betas,betas_selection,maxfit,minfit,cdevolveans)
				
				# Print to log
				stringout = 'GetMetrics(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
										
				# Print to log
				stringout = 'End Generation/Year Loop'+str(gen)+': '+str(datetime.datetime.now() -start_timeGen) + '\n'
				logMsg(logfHndl,stringout)
				if multiprocessing.current_process().name == "S0" or multiprocessing.current_process().name == "MainProcess":
					print(stringout)
			
			# End::generation loop
						
			# ------------------------------------------
			# Call DoPostProcess()
			# ------------------------------------------
			#pdb.set_trace()
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
			Track_EggDeaths,Track_K,Track_N_Init_pop,N_Emigration_pop,N_EmiMortality,N_Immigration_pop,N_ImmiMortality,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,PackingDeathsEmi,PackingDeathsImm,Track_N_Init_age,N_Emigration_age,N_Immigration_age,AgeDeathsOUT,AgeDeathsIN,PackingDeathsEmiAge,PackingDeathsImmAge,Track_MatureCount,Track_ImmatureCount,Track_N_back_age,Track_N_out_age,outputans,gen,Track_CaptureCount_Back,Track_CaptureCount_ClassBack,Track_CaptureCount_Out,Track_CaptureCount_ClassOut,age_size_mean,sizeans,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,SizeDeathsOUT,SizeDeathsIN,N_beforePack_pop,N_beforePack_age,SelectionDeaths_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,F_EmiDist,M_EmiDist,F_EmiDist_sd,M_EmiDist_sd,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,Track_ToTYYMales,Track_BreedYYMales,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,RDispersers,IDispersers,Track_BirthsMYY,Track_KadjEmi,Track_KadjImmi,Track_ToTYYFemales,Track_BirthsFYY,Track_BreedYYFemales)
			
			# Print to log
			stringout = 'DoPostProcess(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			if multiprocessing.current_process().name == "S0" or multiprocessing.current_process().name == "MainProcess":
				print(stringout)
			
			# Print to log
			stringout = 'End Monte Carlo Loop'+str(ithmcrun)+': '+str(datetime.datetime.now() -start_timeMC) + '\n'
			logMsg(logfHndl,stringout)
			if multiprocessing.current_process().name == "S0" or multiprocessing.current_process().name == "MainProcess":
				print(stringout)
			
		# End::Monte Carlo Loop
		
	# End::Batch Loop
