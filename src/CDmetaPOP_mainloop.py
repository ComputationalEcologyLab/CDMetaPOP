# -------------------------------------------------------------------------------------------------
# CDmetaPOP_mainloop.py
# Author: Erin L Landguth, Casey Day
# Created: October 2020
# Description: This is the function/module file for target execution of main loop in parallel runs
# --------------------------------------------------------------------------------------------------

# Python specific functions
import numpy as np 
import sys, datetime
import multiprocessing as mp
from dataclasses import dataclass
from dataclasses import replace
from typing import Any

# CDmetaPOP functions
import CDmetaPOP_Modules as modules
import CDmetaPOP_PostProcess as postprocess
import CDmetaPOP_PreProcess as preprocess
import CDmetaPOP_Mate as mate
import CDmetaPOP_Emigration as emigration
import CDmetaPOP_Immigration as immigration
import CDmetaPOP_Mortality as mortality

@dataclass
class MCArgs:
	ibatch: int
	ithmcrun: int

	spcNO: int
	irun: int

	datadir: str
	outdir: str

	sizeans: str
	constMortans: str
	looptime: int

	gridformat: str
	gridsample: str
	outputans: str

	startcomp: int
	implementcomp: str

	XQs: Any
	nspecies: int
	extinctQ: Any
	global_extinctQ: Any

	passlogfHndl: str

	xyfilename: str

	matecdmatfile: str
	dispOutcdmatfile: str
	dispBackcdmatfile: str
	straycdmatfile: str
	dispLocalcdmatfile: str

	matemoveno: str
	matemoveparA: Any
	matemoveparB: Any
	matemoveparC: Any
	matemovethreshval: Any

	freplace: str
	mreplace: str
	selfing: str

	sexchromo: int
	sexans: str

	assortmateModel_pass: str
	assortmateC_pass: Any

	dispmoveOutno: str
	dispmoveOutparA: Any
	dispmoveOutparB: Any
	dispmoveOutparC: Any
	dispmoveOutthreshval: Any

	dispmoveBackno: str
	dispmoveBackparA: Any
	dispmoveBackparB: Any
	dispmoveBackparC: Any
	dispmoveBackthreshval: Any

	StrBackno: str
	StrBackparA: Any
	StrBackparB: Any
	StrBackparC: Any
	StrBackthreshval: Any

	dispLocalno: str
	dispLocalparA: Any
	dispLocalparB: Any
	dispLocalparC: Any
	dispLocalthreshval: Any

	homeattempt: Any

	offno: Any
	inheritans_classfiles: str
	equalClutch: Any

	muterate_pass: Any
	mutationans: str

	loci: int
	alleles: Any

	mtdna: Any

	geneswap: int

	cdevolveans: str
	burningen_cdevolve: int
	timecdevolve: Any
	betaFile_selection: Any

	plasticans: str
	plastic_signalresp_pass: Any
	plastic_behaviorresp_pass: Any
	burningen_plastic: int
	timeplastic: Any

	growans: str

	sizeLoo: Any
	sizeR0: Any
	size_eqn_1: Any
	size_eqn_2: Any
	size_eqn_3: Any

	eggFreq_mu: float
	eggFreq_sd: float

	egg_mean_ans: Any
	egg_mean_1: float
	egg_mean_2: float

	egg_percmort_mu: Any
	egg_percmort_sd: Any

	Femalepercent_egg: Any

	packans: str
	packpar1: float

	cor_mat_ans: Any
	subpopmort_pass: Any

	egg_delay: int
	egg_add: str

	implementdisease: str

	nthfile: Any
	cdclimgentime: Any

	defaultMature_pass: Any
	mat_slope_pass: Any
	mat_int_pass: Any

	ncores: int
	parallel: str

# --------------------------------------------------------------------------------------------------------------------
def main_loop(inputs, context, XQs, extinctQ, global_extinctQ):

	sizeans = inputs.sizeans
	constMortans = inputs.constMortans
	mcruns = inputs.mcruns
	looptime = inputs.looptime
	nthfile_out = inputs.nthfile_out
	gridformat = inputs.gridformat
	gridsample = inputs.gridsample
	outputans = inputs.outputans
	cdclimgentimelist = inputs.cdclimgentimelist
	startcomp = inputs.startcomp
	implementcomp = inputs.implementcomp
	ncores = inputs.ncores

	appName = context.appName
	appVers = context.appVers
	appRele = context.appRele
	authorNames = context.authorNames
	popvarsfile = context.popvarsfile
	spcNO = context.spcNO
	fileans = context.fileans
	irun = context.irun
	datadir = context.datadir
	outdir = context.outdir
	passlogfHndl = context.passlogfHndl
	nspecies = context.nspecies
	current_system_pid = context.current_system_pid


	with open(passlogfHndl, 'a') as infile:

		modules.logMsg(infile,"\n%s Release %s Version %s\n"%(appName,appRele,appVers), msgVerbose = True)
		modules.logMsg(infile,"Author(s): %s"%(authorNames)+'\n', msgVerbose = True)
		modules.logMsg(infile,"Session runtime inputs from: %s"%(fileans)+'\n\n', msgVerbose = True) 
		modules.logMsg(infile,"Session popvars inputs from: %s"%(popvarsfile[spcNO])+'\n\n', msgVerbose = True)
		modules.logMsg(infile,"On run: %s"%(str(irun))+'\n\n', msgVerbose = True)

	'''Main loop here'''
	# Call function and store PopVar variables
	batchVars,batchVarsIndex,nSimulations = preprocess.loadFile(fileans,1,',',True)

	# ----------------------------------------	
	# Begin Batch Looping 
	# ----------------------------------------
	# This loop is defined by the number of rows in PopVars.csv
	for ibatch in range(nSimulations):
		# Open the logfHndl file
		logfHndl = open(passlogfHndl,'a')
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
		
		growans = batchVars['growth_option'][ibatch]
		sizeLoo = batchVars['growth_Loo'][ibatch] # Check sex ratios
		sizeR0 = batchVars['growth_R0'][ibatch]# Check sex ratios
		size_eqn_1 = batchVars['growth_temp_max'][ibatch]# Check sex ratios
		size_eqn_2 = batchVars['growth_temp_CV'][ibatch]# Check sex ratios
		size_eqn_3 = batchVars['growth_temp_t0'][ibatch]# Check sex ratios
		
		#Maturation pars
		defaultMature_pass = batchVars['mature_default'][ibatch]
		mat_slope_pass = batchVars['mature_eqn_slope'][ibatch]
		mat_int_pass = batchVars['mature_eqn_int'][ibatch]
		
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
		subpopmort_pass = batchVars['subpopmort_file'][ibatch]
		egg_delay = int(batchVars['egg_delay'][ibatch])
		egg_add = batchVars['egg_add'][ibatch]
		
		# PopVars Disease Pars
		implementdisease = batchVars['implement_disease'][ibatch]		
			
		# ------------------------------------------------------------
		# Get nthfile, cdclimgen, alleles in correct format and checks
		# ------------------------------------------------------------
		# Grab the nthfile list range specific to user input, list or sequence
		if not isinstance(nthfile_out, (list,tuple)):
			nthfile_out = int(nthfile_out)
			if nthfile_out > 0:
				nthfile = list(range(0, looptime, nthfile_out))
			else:
				nthfile = [0]		
		else: # If specified years with |
			# Split up list, removing space values, and appending to nthfile
			nthfile = [int(value) for value in nthfile_out if value]
		
		# Error check on nthfile, must be 1 less than looptime for indexing
		modules.validate(max(nthfile) >= looptime, 'nthfile max value must be less than to looptime.')
				
		# Store cdmat file information for tuple or string option passed in from loadFile()
		cdclimgentime = [cdclimgentimelist] if not isinstance(cdclimgentimelist, (list,tuple)) else cdclimgentimelist 
		
		# Create allele array depending on the input format
		alleles = int(batchVars['alleles'][ibatch])*np.ones(loci,int) if len(alleles.split(':')) == 1 else np.asarray(alleles.split(':'),dtype = int)
		
		# ----------------------------
		# For Sex ratios option splits
		# ----------------------------
		modules.validate(sexchromo not in [2, 3, 4], 'Number of sex chromosomes must be 2, 3,or 4.')
				
		# -------------------------------------
		# Error checking 
		# -------------------------------------		
		# DoEmigration() skipped and warning for if selection on
		if cdevolveans != 'N' and timecdevolve.find('Out') != -1:
			if dispOutcdmatfile == 'N':
				modules.logMsg(logfHndl,'Warning: DoEmigration module skipped and spatial selection during this time frame specified, which will also be skipped.')
					 	
		# Constant mortality checks
		modules.validate(constMortans not in ['1','2'], 'Constant mortalities are compounded using option 1 or 2. Enter correct value here.')
			
		# Warning check on multiplicative mortlaity	
		if constMortans == '2':
			modules.logMsg(logfHndl,'Warning: Using multiplicative mortality option, double check for patch and class values when 0%')
				
		# Check on cdevolve answer input
		valid_values = ['1', '2', '1_mat', '2_mat', 'N', 'M', 'G', 'MG_ind', 'MG_link', 'stray', '1_G_ind', '1_G_link', 'Hindex', 'F', 'Plastic', 'Multilocus','runtiming']
		modules.validate(cdevolveans not in valid_values, 'Check CDEvolve answer options.')
		
		# For mature and size ans
		valid_values = ['M', 'MG_ind', 'MG_link', 'G', '1_G_ind', '1_G_link']
		modules.validate(cdevolveans in valid_values and sizeans != 'Y','CDEVOLVE answer is M or G and size answer must be Y.')
				
		# For Hindex answer and each function
		modules.validate('Hindex' in cdevolveans and 'Gauss' in cdevolveans and len(cdevolveans.split('_')[2].split(':')) != 6, 'CDEVOLVE answer is Hindex and 6 parameters for the Gaussian function must be specified, see user manual and example files.')
		modules.validate('Hindex' in cdevolveans and 'Para' in cdevolveans and len(cdevolveans.split('_')[2].split(':')) != 3, 'CDEVOLVE answer is Hindex and 3 parameters for the Parabolic function must be specified, see user manual and example files.')
		modules.validate('Hindex' in cdevolveans and 'Step' in cdevolveans and len(cdevolveans.split('_')[2].split(':')) != 3, 'CDEVOLVE answer is Hindex and 3 parameters for the Step function must be specified, see user manual and example files.')
		modules.validate('Hindex' in cdevolveans and 'Linear' in cdevolveans and len(cdevolveans.split('_')[2].split(':')) != 3, 'CDEVOLVE answer is Hindex and 2 parameters for the Linear function must be specified, see user manual and example files.')
				
		# If cdevolve is turned on must have 2 alleles
		modules.validate((cdevolveans != 'N' or plasticans != 'N') and alleles[0] != 2,'More than 2 alleles per locus specified. CDEVOLVE Or plastic on and 2 alleles should be used needed.')

		# For Plastic answer checks
		if plasticans != 'N':
			modules.validate((plasticans.split('_')[0] != 'Temp') and (plasticans.split('_')[0] != 'Hab'),'Plastic type (Temp/Hab) not entered corectly, check user manual and example files.')
			modules.validate((plasticans.split('_')[1] != 'dom') and (plasticans.split('_')[1] != 'rec') and (plasticans.split('_')[1] != 'codom'), 'Plastic type (dom/codom/rec) not entered corectly, check user manual and example files.')
			modules.validate((timeplastic.find('Out') == -1) and (timeplastic.find('Back') == -1),'Plastic timing must be specified (e.g., Out or Back).')
								
		# Must have more than 1 loci
		modules.validate(loci <= 1, 'CDmetaPOP needs more than 1 locus to run.')
		if cdevolveans != 'N':
			modules.validate((timecdevolve.find('Eggs') == -1) and (timecdevolve.find('Out') == -1) and (timecdevolve.find('Back') == -1) and timecdevolve.find('packing') == -1, 'CDEVOLVE timing must be specified (e.g., Out, Back or Eggs).')
							
		# grid format
		modules.validate((gridformat == 'cdpop' or gridformat == 'general' or gridformat == 'genalex' or gridformat == 'genepop' or gridformat == 'structure') == False, 'Grid format parameter not an option.')
				
		# Check burn in times
		if cdevolveans != 'N' and burningen_cdevolve < geneswap:
			modules.logMsg(logfHndl,'Warning: Selection burnin time < time at which genetic exchange is to initialize, setting burnin time = start genetic exchange time.')			
			burningen_cdevolve = geneswap			
		if plasticans != 'N' and burningen_plastic < geneswap:
			modules.logMsg(logfHndl,'Warning: Selection burnin time < time at which genetic exchange is to initialize, setting burnin time = start genetic exchange time.')
			burningen_plastic = geneswap
			
		# Inherit answer can only be:
		modules.validate(not (inheritans_classfiles == 'random' or inheritans_classfiles == 'Hindex' or inheritans_classfiles == 'mother'),'Inherit answer for multiple class files is not correct: enter either random or Hindex.')
			
		# If inherit answer uses Hindex, mutation can't be on
		#if isinstance(muterate_pass, list):
		if len(muterate_pass) > 1:
			modules.validate(sum(np.asarray(muterate_pass,dtype=float)) != 0.0 and (inheritans_classfiles == 'Hindex' or inheritans_classfiles == 'mother'),'Mutation is not operating with Hindex inheritance options in this version.')
		else:
			modules.validate(sum(np.asarray([muterate_pass],dtype=float)) != 0.0 and (inheritans_classfiles == 'Hindex' or inheritans_classfiles == 'mother'),'Mutation is not operating with Hindex inheritance options in this version.')
			
		# If egg_delay is gretter than 1 v2.68 turned off for further testing
		modules.validate(egg_delay > 0,'Currently, egg delay is not operating beyond 0 year/time unit. Testing sourcePop and gen = 0 issue.')
					
		# Reproduction answers checks 	
		modules.validate((sexans == 'N' or sexans == 'Y' or sexans == 'H') == False,'Reproduction choices either N, Y or H, check user manual.')
		modules.validate(sexans == 'H' and (selfing == 'N' or selfing == 'Y'),'Hermaphroditic mating structure specified - H - then must specify the selfing probability.')
									
		# Spot to add age 0s / eggs / pups /liter, etc
		modules.validate(egg_add not in ['mating', 'nonmating'], 'Egg add choice incorrect.')
		#modules.validate(egg_add not in ['mating', 'nonmating'] and packans.split('_')[0] not in 'logistic', 'Nonmating option should use logistic model.')
		modules.validate(egg_add == 'nonmating' and packans != 'logistic_back', 'Nonmating option must be used with logistic in the back location.')
					
		valid_values = ['N','Both','Back','Out']
		modules.validate(implementdisease not in valid_values, 'Implement disease value incorrect.')
		modules.validate(implementdisease!='N' and alleles[0] != 2,'More than 2 alleles per locus specified with disease defense options.')
		
		valid_values = ['N','Both','Back','Out']
		modules.validate(implementcomp not in valid_values, 'Implement competition value incorrect.')
		
		# Check number of processors requested is available
		if ncores > mp.cpu_count() - 1:
			print("Error: Processors requested (ncores) is greater than number of processors available.")
			sys.exit()
		# Prevent attempts to split cores across MCs when cihld processes have already spawned to handle multiple species and set type of parallel processing.		
		if nspecies > 1:
			parallel = 'species' # Species parallel processing
			if ncores > 1:
				ncores = 1
				if mp.current_process().name == "S0":
					print("Note: Parallel processing of Monte Carlo replicates is unavailable for multispecies runs. Proceeding with replicates in series.")		
		elif mcruns == 1 or ncores == 1:
			parallel = 'N' # No parallel processing
		elif ncores > 1 and mcruns > 1 and nspecies == 1:
			parallel = 'mc' # Monte Carlo parallel processing
		else:
			print("Error: Check inputs for ncores and mcruns in RunVars file.")
			sys.exit()		
		# ---------------------------------------------	
		# Begin Monte-Carlo Looping
		# ---------------------------------------------		
		# Create set of variable arguments to pass into mc_loop(). mp.pool() requires a single argument to be passed in.
		logfHndl.close() # Close log file because you can't pass an open file to mp.pool
		mc_args = MCArgs(
			ibatch = ibatch,
			ithmcrun = 0,
			spcNO = spcNO,
			irun = irun,
			datadir = datadir,
			sizeans = sizeans,
			constMortans = constMortans,
			looptime = looptime,
			gridformat = gridformat,
			gridsample = gridsample,
			outputans = outputans,
			outdir = outdir,
			startcomp = startcomp,
			implementcomp = implementcomp,
			XQs = XQs,
			nspecies = nspecies,
			extinctQ = extinctQ,
			global_extinctQ = global_extinctQ,
			passlogfHndl = passlogfHndl,
			xyfilename = xyfilename,
			matecdmatfile = matecdmatfile,
			dispOutcdmatfile = dispOutcdmatfile,
			dispBackcdmatfile = dispBackcdmatfile,
			straycdmatfile = straycdmatfile,
			dispLocalcdmatfile = dispLocalcdmatfile,
			matemoveno = matemoveno,
			matemoveparA = matemoveparA,
			matemoveparB = matemoveparB,
			matemoveparC = matemoveparC,
			matemovethreshval = matemovethreshval,
			freplace = freplace,
			mreplace = mreplace,
			selfing = selfing,
			sexchromo = sexchromo,
			sexans = sexans,
			assortmateModel_pass = assortmateModel_pass,
			assortmateC_pass = assortmateC_pass,
			dispmoveOutno = dispmoveOutno,
			dispmoveOutparA = dispmoveOutparA,
			dispmoveOutparB = dispmoveOutparB,
			dispmoveOutparC = dispmoveOutparC,
			dispmoveOutthreshval = dispmoveOutthreshval,
			dispmoveBackno = dispmoveBackno,
			dispmoveBackparA = dispmoveBackparA,
			dispmoveBackparB = dispmoveBackparB,
			dispmoveBackparC = dispmoveBackparC,
			dispmoveBackthreshval = dispmoveBackthreshval,
			StrBackno = StrBackno,
			StrBackparA = StrBackparA,
			StrBackparB = StrBackparB,
			StrBackparC = StrBackparC,
			StrBackthreshval = StrBackthreshval,
			dispLocalno = dispLocalno,
			dispLocalparA = dispLocalparA,
			dispLocalparB = dispLocalparB,
			dispLocalparC = dispLocalparC,
			dispLocalthreshval = dispLocalthreshval,
			homeattempt = homeattempt,
			offno = offno,
			inheritans_classfiles = inheritans_classfiles,
			equalClutch = equalClutch,
			muterate_pass = muterate_pass,
			mutationans = mutationans,
			loci = loci,
			alleles = alleles,
			mtdna = mtdna,
			geneswap = geneswap,
			cdevolveans = cdevolveans,
			burningen_cdevolve = burningen_cdevolve,
			timecdevolve = timecdevolve,
			betaFile_selection = betaFile_selection,
			plasticans = plasticans,
			plastic_signalresp_pass = plastic_signalresp_pass,
			plastic_behaviorresp_pass = plastic_behaviorresp_pass,
			burningen_plastic = burningen_plastic,
			timeplastic = timeplastic,
			growans = growans,
			sizeLoo = sizeLoo,
			sizeR0 = sizeR0,
			size_eqn_1 = size_eqn_1,
			size_eqn_2 = size_eqn_2,
			size_eqn_3 = size_eqn_3,
			eggFreq_mu = eggFreq_mu,
			eggFreq_sd = eggFreq_sd,
			egg_mean_ans = egg_mean_ans,
			egg_mean_1 = egg_mean_1,
			egg_mean_2 = egg_mean_2,
			egg_percmort_mu = egg_percmort_mu,
			egg_percmort_sd = egg_percmort_sd,
			Femalepercent_egg = Femalepercent_egg,
			packans = packans,
			packpar1 = packpar1,
			cor_mat_ans = cor_mat_ans,
			subpopmort_pass = subpopmort_pass,
			egg_delay = egg_delay,
			egg_add = egg_add,
			implementdisease = implementdisease,
			nthfile = nthfile,
			cdclimgentime = cdclimgentime,
			defaultMature_pass = defaultMature_pass,
			mat_slope_pass = mat_slope_pass,
			mat_int_pass = mat_int_pass,
			ncores = ncores,
			parallel = parallel
			#'''FXXmat_set':FXXmat_set,
			#'MXYmat_set':MXYmat_set,
			#'MYYmat_set':MYYmat_set,
			#'FYYmat_set':FYYmat_set,
			#'FXXmat_slope':FXXmat_slope,
			#'MXYmat_slope':MXYmat_slope,
			#'MYYmat_slope':MYYmat_slope,
			#'FYYmat_slope':FYYmat_slope,
			#'FXXmat_int':FXXmat_int,
			#'MXYmat_int':MXYmat_int,
			#'MYYmat_int':MYYmat_int,
			#'FYYmat_int':FYYmat_int,'''
		)
		# If using multipsecies, or single species without multiprocessing, run as for loop
		if parallel in ['species','N']:
			# Monte carlo replicate loop
			for ithmcrun in range(mcruns):
				# Need to pass in which mc is running
				run_args = replace(mc_args, ithmcrun=ithmcrun)
				# Run the mc loop
				mc_loop(run_args)
		# If single species with multiprocessing mc replicates, use with Pool()
		elif parallel == 'mc':
			# List comprehension that creates an iterable list of dictionaries to loop through, with an added element of 'ithmcrun' for the mc number
			mc_args_list = [replace(mc_args, ithmcrun=ithmcrun) for ithmcrun in range(mcruns)]
			# Reset counter to name worker id for every new pool creation
			shared_counter = mp.Value('i', 0)
			# Start parallel processing
			with mp.Pool(processes = ncores, initializer = worker_init, initargs = (shared_counter,)) as pool:
				pool.map(mc_loop, mc_args_list)				
				#pool.starmap(mc_loop, [(mc_args,ithmcrun) for ithmcrun in range(mcruns)])				
		else:
			print("Error: Check values for multiprocessing and number of cores.")
			sys.exit()
	# Empty queues before ending batch run to prevent hanging 
	for q in XQs:
		for q2 in q: 
			while not q2.empty():
				q2.get()
			q2.close()
			q2.join_thread()
	logfHndl.close()
	# End::Batch Loop
	
# Function to run monte carlo loop
def mc_loop(args: MCArgs):
	
	# Timing events: start
	start_timeMC = datetime.datetime.now()
	# Keep track so extinction message are printed only once
	temp_extinct = 0	
	# Reopen log file
	logfHndl = open(args.passlogfHndl,'a')
	
	# If multi-species, queues were created before the multiprocess split, so they can't be closed and created again here.
	if args.parallel == 'species':
		pass
	else:
		# Re-create Queues to keep track of extinctions, since Queues cannot be passed to multiprocessing
		args.extinctQ = mp.Queue() # To track extinction. If all species extinct, exit system			
		args.global_extinctQ = mp.Queue() # To track global extinction
	# -----------------------------------------
	# Create storage variables
	# ------------------------------------------	
	# These variables will be stored in output.csv at the end of the simulation						
	
	# GetMetrics()
	Track_p1, Track_p2, Track_q1, Track_q2, Track_Alleles, Track_He, Track_Ho, Track_N_Init_pop, Track_N_Init_age, Track_N_Init_class, Track_K, Track_CaptureCount_Out, Track_CaptureCount_ClassOut, Track_CaptureCount_Back, Track_CaptureCount_ClassBack, maxfit, minfit,Track_DiseaseStates_pop, Track_DiseaseStates_EnvRes  = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [],[],[]
	
	# AddIndividuals()
	Track_DiseaseStates_AddedInds = []
							
	# DoMate()
	Track_FAvgMate, Track_MAvgMate, Track_FSDMate, Track_MSDMate, Track_MateDistCD, Track_MateDistCDstd, Track_BreedEvents, Track_AAaaMates, Track_AAAAMates, Track_aaaaMates, Track_AAAaMates, Track_aaAaMates, Track_AaAaMates, Track_BreedFemales, Track_BreedMales, Track_BreedYYMales, Track_BreedYYFemales, Track_MatureCount, Track_ImmatureCount, Track_ToTFemales, Track_ToTMales, Track_ToTYYMales, Track_ToTYYFemales = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
	
	# DoOffspring
	Track_Births, Track_EggDeaths, Track_BirthsMYY, Track_BirthsFYY, Track_DiseaseStates_AddAge0s = [], [], [], [], []
	
	# DoUpdate
	Track_N_back_age, Track_N_out_age, Track_DiseaseStates_SecondUpdate, Track_DiseaseStates_ThirdUpdate,Track_DiseaseStates_AfterDeaths_SecondUpdate,Track_DiseaseStates_AfterDeaths_ThirdUpdate = [], [],[],[],[],[]
	
	# Emigration()
	N_Emigration_pop, N_Emigration_age, subpopemigration, F_EmiDist, M_EmiDist, F_EmiDist_sd, M_EmiDist_sd, SelectionDeathsEmi, DisperseDeathsEmi, PackingDeathsEmi, PackingDeathsEmiAge, MgSuccess, AdultNoMg, Track_YYSelectionPackDeathsEmi, Track_WildSelectionPackDeathsEmi, SelectionDeaths_Age0s, N_beforePack_pop, N_beforePack_age, Track_KadjEmi = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
	
	# Mortlity after Emigration
	N_EmiMortality, PopDeathsOUT, AgeDeathsOUT, SizeDeathsOUT = [], [], [], []
	
	# Immigration
	N_Immigration_pop, N_Immigration_age, subpopimmigration, F_HomeDist, M_HomeDist, F_HomeDist_sd, M_HomeDist_sd, F_StrayDist, M_StrayDist, F_StrayDist_sd, M_StrayDist_sd, F_ZtrayDist, M_ZtrayDist, F_ZtrayDist_sd, M_ZtrayDist_sd, SelectionDeathsImm, DisperseDeathsImm, PackingDeathsImmAge, PackingDeathsImm, StrSuccess, Track_YYSelectionPackDeathsImmi, Track_WildSelectionPackDeathsImmi, Track_KadjImmi = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

	# Mortality after immigration
	N_ImmiMortality, PopDeathsIN, AgeDeathsIN, SizeDeathsIN = [], [], [], []
	
	# DoOutput()
	Residors, Strayers1, Strayers2, Immigrators, IDispersers, RDispersers, PopSizes_Mean, PopSizes_Std, AgeSizes_Mean, AgeSizes_Std, ClassSizes_Mean, ClassSizes_Std = [], [], [], [], [], [], [], [], [], [], [], []
	
	# Non-tracking variables - create empty 2-D list
	noOffspring_temp = [np.asarray([]),np.asarray([])] 
	Bearpairs_temp = [[[-9999,-9999]],[[-9999,-9999]]]
			
	# ------------------------------------	
	# Call DoPreProcess()
	# ------------------------------------
	
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Call function
	tupPreProcess = preprocess.DoPreProcess(args, 0, logfHndl, subpopemigration, subpopimmigration)


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
	MgBack_patch_prob_pass = tupPreProcess[48] # CDClimate then DoImmigrtion
	Disperse_patch_prob_pass = tupPreProcess[49] # Cdclimate thenDoImmigration
	alldiseaseVars_files = tupPreProcess[50]
	# These are temporary disease vars to pass into the first GetMetrics, will use the cdclimate returned variables for the time loop
	disease_vars_pass = tupPreProcess[51] 
	pathogen_load_pass = tupPreProcess[52]
	disease_fitvals_pass = tupPreProcess[53]
				
	K = K_mu # Initialize K with mu	, Grab first one only	
				
	# Print to log
	stringout = 'DoPreProcess(): '+str(datetime.datetime.now() -start_time1) + ''
	modules.logMsg(logfHndl,stringout)	
	if mp.current_process().name in ("S0", "MainProcess") or mp.current_process().name.endswith("-1"):
		print(('DoPreProcess(): ',str(datetime.datetime.now() -start_time1),''))
	
	# ---------------------------------
	# Call GetMetrics()
	# ---------------------------------
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	modules.GetMetrics(SubpopIN_init,K,Track_N_Init_pop,Track_K,args.loci,args.alleles,0,Track_Ho,Track_Alleles,Track_He,Track_p1,Track_p2,Track_q1,Track_q2,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,Track_N_Init_age,args.sizeans,age_size_mean,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,args.packans,RDispersers,IDispersers,xvars_betas_pass,tempbetas_selection,maxfit,minfit,args.cdevolveans,disease_vars_pass,Track_DiseaseStates_pop,Track_DiseaseStates_EnvRes)
	
	# Print to log
	stringout = 'GetMetrics() Initial: '+str(datetime.datetime.now() -start_time1) + ''
	modules.logMsg(logfHndl,stringout)
	
	# ---------------------------------
	# Error statements
	# ---------------------------------			
	# Error statement here in case no females or males, then break
	modules.validate(Track_N_Init_pop[0][0] == 0, 'There are no individuals to begin time loop for this species ',str(args.spcNO))
				
	# ----------------------------------------------------
	# Call DoUpdate() - output initial file here ind-1.csv
	# ----------------------------------------------------
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	modules.DoUpdate(args.packans,SubpopIN_init,K,xgridpop,ygridpop,-1,args.nthfile,ithmcrundir,args.loci,args.alleles,logfHndl,'Initial')
	
	# Print to log
	stringout = 'DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
	modules.logMsg(logfHndl,stringout)
	
	# -------------------------------------------
	# Start Generation Looping 
	# -------------------------------------------
	# Begin generation loop

	for gen in range(args.looptime):				
		# Timing events: start
		start_timeGen = datetime.datetime.now()
							
		# If initial generation - update with initial populations
		if gen == 0:
			SubpopIN = SubpopIN_init
			del SubpopIN_init
			# Use NatalPop in first generation 
			sourcePop = 'NatalPop'		
		else: # This was for versions previous v1.37
			sourcePop = 'ImmiPop' # For generations greater than 0, this is where the individual is located after DoEmi and DoImmi
			#sourcePop = 'NatalPop'
		
		# Exit the system if population is 0 or 1
		checkPopN = [len(SubpopIN[x]) for x in range(0,len(SubpopIN))] 
		if sum(checkPopN) == 0:
			stringout = 'Species is extinct.'
			modules.logMsg(logfHndl,stringout)
			# If not extinct from previous generation
			if temp_extinct == 0:					
				print(('Species ' + str(args.spcNO) + ' went extinct.'))
				temp_extinct = 1
			# Track extinctions
			args.extinctQ.put(0)
		else:
			# Track extinctions
			args.extinctQ.put(1)
		# List to track extinctions - S0 so only one species handles queues in a multispecies application (i.e., NOT 'S1'), MainProcess for a single species, single core run,
		# and multiprocess for single-species, multi-core runs
		
		if mp.current_process().name == "S0" or mp.current_process().name == "MainProcess" or args.parallel == 'mc':
			ext_list = []	
			for ispecies in range(args.nspecies):
				ext_list.append(args.extinctQ.get(block=True))
			# If all species extinct, exit
			if sum(ext_list) == 0:
				# Exit the system if population is 0 or 1
				stringout = 'All species extinct after generation '+str(gen-1)+', program ended.\n'
				modules.logMsg(logfHndl,stringout)
				if mp.current_process().name in ("S0", "MainProcess") or mp.current_process().name.endswith("-1"):
					print('All species extinct')
				#print('Population went extinct after generation '+str(gen-1)+'.\n')
				for ispecies in range(args.nspecies):
					args.global_extinctQ.put(1)
			else:
				for ispecies in range(args.nspecies):
					args.global_extinctQ.put(0)	
		if args.global_extinctQ.get() == 1:
			break
		
		# ---------------------------------
		# Call CDClimate()
		# ---------------------------------			
		# Timing events: start
		start_time1 = datetime.datetime.now()

		# Check gen time equal to args.cdclimgentime
		for icdtime in range(len(args.cdclimgentime)): 
			if gen == int(args.cdclimgentime[icdtime]):
				climate = preprocess.DoCDClimate(args, icdtime, MgOut_patch_pass, Str_patch_pass, Kmu_pass, outsizevals_pass, backsizevals_pass, outgrowdays_pass, backgrowdays_pass, fitvals_pass, popmort_back_pass, popmort_out_pass, eggmort_pass, Kstd_pass, popmort_back_sd_pass, popmort_out_sd_pass, eggmort_sd_pass, outsizevals_sd_pass, backsizevals_sd_pass, outgrowdays_sd_pass, backgrowdays_sd_pass, pop_capture_back_pass, pop_capture_out_pass, N0_pass, allefreqfiles_pass, classvarsfiles_pass, PopTag, comp_coef_pass, xvars_betas_pass, outhabvals_pass, backhabvals_pass, MgBack_patch_prob_pass, Disperse_patch_prob_pass, alldiseaseVars_files,  pathogen_load_pass, disease_fitvals_pass)
#				# Cdmatrix values
#				cdmatrix_mate, cdmatrix_FXXOut, cdmatrix_MXYOut, cdmatrix_MYYOut, cdmatrix_FYYOut, cdmatrix_FXXBack, cdmatrix_MXYBack, cdmatrix_MYYBack, cdmatrix_FYYBack, cdmatrix_FXXStr, cdmatrix_MXYStr, cdmatrix_MYYStr, cdmatrix_FYYStr, cdmatrix_FXXLD, cdmatrix_MXYLD, cdmatrix_MYYLD, cdmatrix_FYYLD = tupClimate[:17]
#				# Threshold values
#				thresh_mate, thresh_FXXOut, thresh_MXYOut, thresh_MYYOut, thresh_FYYOut, thresh_FXXBack, thresh_MXYBack, thresh_MYYBack, thresh_FYYBack, thresh_FXXStr, thresh_MXYStr, thresh_MYYStr, thresh_FYYStr, thresh_FXXLD, thresh_MXYLD, thresh_MYYLD, thresh_FYYLD = tupClimate[17:34]						
#				# Scale Min
#				scalemin_mate,scalemin_FXXOut,scalemin_MXYOut,scalemin_MYYOut,scalemin_FYYOut,scalemin_FXXBack,scalemin_MXYBack,scalemin_MYYBack,scalemin_FYYBack,scalemin_FXXStr,scalemin_MXYStr,scalemin_MYYStr,scalemin_FYYStr,scalemin_FXXLD,scalemin_MXYLD,scalemin_MYYLD,scalemin_FYYLD = tupClimate[34:51]
#				# Scale Max
#				scalemax_mate, scalemax_FXXOut, scalemax_MXYOut, scalemax_MYYOut,scalemax_FYYOut ,scalemax_FXXBack,scalemax_MXYBack,scalemax_MYYBack,scalemax_FYYBack,scalemax_FXXStr,scalemax_MXYStr,scalemax_MYYStr,scalemax_FYYStr,scalemax_FXXLD,scalemax_MXYLD,scalemax_MYYLD,scalemax_FYYLD = tupClimate[51:68]
#				# ParA
#				parA_mate,parA_FXXOut,parA_MXYOut,parA_MYYOut,parA_FYYOut,parA_FXXBack,parA_MXYBack,parA_MYYBack,parA_FYYBack,parA_FXXStr,parA_MXYStr,parA_MYYStr,parA_FYYStr,parA_FXXLD,parA_MXYLD,parA_MYYLD,parA_FYYLD = tupClimate[68:85]
#				# ParB
#				parB_mate,parB_FXXOut,parB_MXYOut,parB_MYYOut,parB_FYYOut,parB_FXXBack,parB_MXYBack,parB_MYYBack,parB_FYYBack,parB_FXXStr,parB_MXYStr,parB_MYYStr,parB_FYYStr,parB_FXXLD,parB_MXYLD,parB_MYYLD,parB_FYYLD = tupClimate[85:102]
#				# ParC
#				parC_mate,parC_FXXOut,parC_MXYOut,parC_MYYOut,parC_FYYOut,parC_FXXBack,parC_MXYBack,parC_MYYBack,parC_FYYBack,parC_FXXStr,parC_MXYStr,parC_MYYStr,parC_FYYStr,parC_FXXLD,parC_MXYLD,parC_MYYLD,parC_FYYLD = tupClimate[102:119]
#				# Movement No
#				moveno_mate,moveno_FXXOut,moveno_MXYOut,moveno_MYYOut,moveno_FYYOut,moveno_FXXBack,moveno_MXYBack,moveno_MYYBack,moveno_FYYBack,moveno_FXXStr,moveno_MXYStr,moveno_MYYStr,moveno_FYYStr,moveno_FXXLD,moveno_MXYLD,moveno_MYYLD,moveno_FYYLD = tupClimate[119:136]						
#				MgOut_patch_prob = tupClimate[136]
#				Str_patch_prob = tupClimate[137]						
#				# Temp/grow vars (mu)
#				outsizevals_mu,backsizevals_mu,outgrowdays_mu,backgrowdays_mu = tupClimate[138:142]
#				fitvals, K_mu, popmort_back_mu, popmort_out_mu, eggmort_mu, K_std, popmort_back_sd, popmort_out_sd, eggmort_sd = tupClimate[142:151]
#				# Temp/grow vars (sd)
#				outsizevals_sd,backsizevals_sd,outgrowdays_sd,backgrowdays_sd = tupClimate[151:155]
#				pop_capture_back, pop_capture_out, tempN0, tempAllelefile, tempClassVarsfile, assortmateModel, assortmateC, subpopmort_mat, comp_coef, betas_selection, xvars_betas = tupClimate[155:166]
#				# Plastic selection vars
#				outhabvals,backhabvals,plastic_signalresp,plastic_behaviorresp = tupClimate[166:170] #ts added, note that plasticans not listed here, but present in tupClimate. in 1.72, present in tupClimate as last piece but not defined down here.
#				muterate = tupClimate[170]						
#				if gen == 0:
#					age_percmort_out_mu,age_percmort_out_sd,age_percmort_back_mu,age_percmort_back_sd,size_percmort_out_mu,size_percmort_out_sd,size_percmort_back_mu,size_percmort_back_sd,age_MgOUT, age_MgBACK,age_S,age_DispProb,age_mature,age_mu,age_sigma,f_leslie_mu,f_leslie_std,age_capture_out,age_capture_back = tupClimate[171:190]
#				MgBack_patch_prob, Disperse_patch_prob =  tupClimate[190:192]
#				disease_vars = tupClimate[192]
#				# Matruation vars
#				FXXmat_set, MXYmat_set, MYYmat_set, FYYmat_set,FXXmat_slope,MXYmat_slope,MYYmat_slope,FYYmat_slope,FXXmat_int,MXYmat_int,MYYmat_int,FYYmat_int = tupClimate[193:205]

				# =====================================================================
				# 1. Matrices
				# =====================================================================
				cdmatrix_mate		  = climate.matecdmatrix
				cdmatrix_FXXOut		= climate.FXXdispOutcdmatrix
				cdmatrix_MXYOut		= climate.MXYdispOutcdmatrix
				cdmatrix_MYYOut		= climate.MYYdispOutcdmatrix
				cdmatrix_FYYOut		= climate.FYYdispOutcdmatrix
				cdmatrix_FXXBack	   = climate.FXXdispBackcdmatrix
				cdmatrix_MXYBack	   = climate.MXYdispBackcdmatrix
				cdmatrix_MYYBack	   = climate.MYYdispBackcdmatrix
				cdmatrix_FYYBack	   = climate.FYYdispBackcdmatrix
				cdmatrix_FXXStr		= climate.FXXStrcdmatrix
				cdmatrix_MXYStr		= climate.MXYStrcdmatrix
				cdmatrix_MYYStr		= climate.MYYStrcdmatrix
				cdmatrix_FYYStr		= climate.FYYStrcdmatrix
				cdmatrix_FXXLD		 = climate.FXXdispLocalcdmatrix
				cdmatrix_MXYLD		 = climate.MXYdispLocalcdmatrix
				cdmatrix_MYYLD		 = climate.MYYdispLocalcdmatrix
				cdmatrix_FYYLD		 = climate.FYYdispLocalcdmatrix

				# =====================================================================
				# 2. Threshold values
				# =====================================================================
				thresh_mate			= climate.matemovethresh
				thresh_FXXOut		  = climate.FXXdispOutthresh
				thresh_MXYOut		  = climate.MXYdispOutthresh
				thresh_MYYOut		  = climate.MYYdispOutthresh
				thresh_FYYOut		  = climate.FYYdispOutthresh
				thresh_FXXBack = climate.FXXdispOutthresh
				thresh_MXYBack = climate.MXYdispOutthresh
				thresh_MYYBack = climate.MYYdispOutthresh
				thresh_FYYBack = climate.FYYdispOutthresh
				# Maps matching the unique properties sequentially:
				thresh_FXXStr		  = climate.FXXStrthresh
				thresh_MXYStr		  = climate.MXYStrthresh
				thresh_MYYStr		  = climate.MYYStrthresh
				thresh_FYYStr		  = climate.FYYStrthresh
				thresh_FXXLD		   = climate.FXXdispLocalthresh
				thresh_MXYLD		   = climate.MXYdispLocalthresh
				thresh_MYYLD		   = climate.MYYdispLocalthresh
				thresh_FYYLD		   = climate.FYYdispLocalthresh

				# =====================================================================
				# 3. Scale Min
				# =====================================================================
				scalemin_mate		  = climate.mate_ScaleMin
				scalemin_FXXOut		= climate.FXXdispOut_ScaleMin
				scalemin_MXYOut		= climate.MXYdispOut_ScaleMin
				scalemin_MYYOut		= climate.MYYdispOut_ScaleMin
				scalemin_FYYOut		= climate.FYYdispOut_ScaleMin
				scalemin_FXXBack	   = climate.FXXdispBack_ScaleMin
				scalemin_MXYBack	   = climate.MXYdispBack_ScaleMin
				scalemin_MYYBack	   = climate.MYYdispBack_ScaleMin
				scalemin_FYYBack	   = climate.FYYdispBack_ScaleMin
				scalemin_MXYLD = climate.MXYdispBack_ScaleMin
				scalemin_MYYLD = climate.MYYdispBack_ScaleMin
				scalemin_FYYLD = climate.FYYdispBack_ScaleMin
				scalemin_FXXStr		= climate.FXXStr_ScaleMin
				scalemin_MXYStr		= climate.MXYStr_ScaleMin
				scalemin_MYYStr		= climate.MYYStr_ScaleMin
				scalemin_FYYStr		= climate.FYYStr_ScaleMin
				scalemin_FXXLD		 = climate.FXXdispLocal_ScaleMin

				# =====================================================================
				# 4. Scale Max
				# =====================================================================
				scalemax_mate		  = climate.mate_ScaleMax
				scalemax_FXXOut		= climate.FXXdispOut_ScaleMax
				scalemax_MXYOut		= climate.MXYdispOut_ScaleMax
				scalemax_MYYOut		= climate.MYYdispOut_ScaleMax
				scalemax_FYYOut		= climate.FYYdispOut_ScaleMax
				scalemax_FXXBack	   = climate.FXXdispBack_ScaleMax
				scalemax_MXYBack	   = climate.MXYdispBack_ScaleMax
				scalemax_MYYBack	   = climate.MYYdispBack_ScaleMax
				scalemax_FYYBack	   = climate.FYYdispBack_ScaleMax
				scalemax_MXYLD = climate.MXYdispBack_ScaleMax
				scalemax_MYYLD = climate.MYYdispBack_ScaleMax
				scalemax_FYYLD = climate.FYYdispBack_ScaleMax
				scalemax_FXXStr		= climate.FXXStr_ScaleMax
				scalemax_MXYStr		= climate.MXYStr_ScaleMax
				scalemax_MYYStr		= climate.MYYStr_ScaleMax
				scalemax_FYYStr		= climate.FYYStr_ScaleMax
				scalemax_FXXLD		 = climate.FXXdispLocal_ScaleMax

				# =====================================================================
				# 5. Parameter Models (A)
				# =====================================================================
				parA_mate			  = climate.matemoveparA
				parA_FXXOut			= climate.FXXdispmoveOutparA
				parA_MXYOut			= climate.MXYdispmoveOutparA
				parA_MYYOut			= climate.MYYdispmoveOutparA
				parA_FYYOut			= climate.FYYdispmoveOutparA
				parA_FXXBack		   = climate.FXXdispmoveBackparA
				parA_MXYBack		   = climate.MXYdispmoveBackparA
				parA_MYYBack		   = climate.MYYdispmoveBackparA
				parA_FYYBack		   = climate.FYYdispmoveBackparA
				parA_FXXStr			= climate.FXXStrparA
				parA_MXYStr			= climate.MXYStrparA
				parA_MYYStr			= climate.MYYStrparA
				parA_FYYStr			= climate.FYYStrparA
				parA_FXXLD			 = climate.FXXdispLocalparA
				parA_MXYLD			 = climate.MXYdispLocalparA
				parA_MYYLD			 = climate.MYYdispLocalparA
				parA_FYYLD			 = climate.FYYdispLocalparA

				# =====================================================================
				# 6. Parameter Models (B)
				# =====================================================================
				parB_mate			  = climate.matemoveparB
				parB_FXXOut			= climate.FXXdispmoveOutparB
				parB_MXYOut			= climate.MXYdispmoveOutparB
				parB_MYYOut			= climate.MYYdispmoveOutparB
				parB_FYYOut			= climate.FYYdispmoveOutparB
				parB_FXXBack		   = climate.FXXdispmoveBackparB
				parB_MXYBack		   = climate.MXYdispmoveBackparB
				parB_MYYBack		   = climate.MYYdispmoveBackparB
				parB_FYYBack		   = climate.FYYdispmoveBackparB
				parB_FXXStr			= climate.FXXStrparB
				parB_MXYStr			= climate.MXYStrparB
				parB_MYYStr			= climate.MYYStrparB
				parB_FYYStr			= climate.FYYStrparB
				parB_FXXLD			 = climate.FXXdispLocalparB
				parB_MXYLD			 = climate.MXYdispLocalparB
				parB_MYYLD			 = climate.MYYdispLocalparB
				parB_FYYLD			 = climate.FYYdispLocalparB

				# =====================================================================
				# 7. Parameter Models (C)
				# =====================================================================
				parC_mate			  = climate.matemoveparC
				parC_FXXOut			= climate.FXXdispmoveOutparC
				parC_MXYOut			= climate.MXYdispmoveOutparC
				parC_MYYOut			= climate.MYYdispmoveOutparC
				parC_FYYOut			= climate.FYYdispmoveOutparC
				parC_FXXBack		   = climate.FXXdispmoveBackparC
				parC_MXYBack		   = climate.MXYdispmoveBackparC
				parC_MYYBack		   = climate.MYYdispmoveBackparC
				parC_FYYBack		   = climate.FYYdispmoveBackparC
				parC_FXXStr			= climate.FXXStrparC
				parC_MXYStr			= climate.MXYStrparC
				parC_MYYStr			= climate.MYYStrparC
				parC_FYYStr			= climate.FYYStrparC
				parC_FXXLD			 = climate.FXXdispLocalparC
				parC_MXYLD			 = climate.MXYdispLocalparC
				parC_MYYLD			 = climate.MYYdispLocalparC
				parC_FYYLD			 = climate.FYYdispLocalparC

				# =====================================================================
				# 8. Movement Numbers
				# =====================================================================
				moveno_mate			= climate.matemoveno
				moveno_FXXOut		  = climate.FXXdispmoveOutno
				moveno_MXYOut		  = climate.MXYdispmoveOutno
				moveno_MYYOut		  = climate.MYYdispmoveOutno
				moveno_FYYOut		  = climate.FYYdispmoveOutno
				moveno_FXXBack		 = climate.FXXdispmoveBackno
				moveno_MXYBack		 = climate.MXYdispmoveBackno
				moveno_MYYBack		 = climate.MYYdispmoveBackno
				moveno_FYYBack		 = climate.FYYdispmoveBackno
				moveno_FXXStr		  = climate.FXXStrno
				moveno_MXYStr		  = climate.MXYStrno
				moveno_MYYStr		  = climate.MYYStrno
				moveno_FYYStr		  = climate.FYYStrno
				moveno_FXXLD		   = climate.FXXdispLocalno
				moveno_MXYLD		   = climate.MXYdispLocalno
				moveno_MYYLD		   = climate.MXYdispLocalno
				moveno_FYYLD		   = climate.FYYdispLocalno

				# =====================================================================
				# 9. Individual Patch, Size & Mortality Variables
				# =====================================================================
				MgOut_patch_prob	   = climate.tempMgOut_patch_prob
				Str_patch_prob		 = climate.tempStr_patch_prob
				outsizevals_mu		 = climate.tempoutsize
				backsizevals_mu		= climate.tempbacksize
				outgrowdays_mu		 = climate.tempoutgrow
				backgrowdays_mu		= climate.tempbackgrow
				fitvals				= climate.tempfitvals
				K_mu				   = climate.tempK
				popmort_back_mu		= climate.temppopmort_back
				popmort_out_mu		 = climate.temppopmort_out
				eggmort_mu			 = climate.tempeggmort
				K_std				  = climate.tempKstd
				popmort_back_sd		= climate.temppopmort_back_sd
				popmort_out_sd		 = climate.temppopmort_out_sd
				eggmort_sd			 = climate.tempeggmort_sd
				outsizevals_sd		 = climate.tempoutsize_sd
				backsizevals_sd		= climate.tempbacksize_sd
				outgrowdays_sd		 = climate.tempoutgrow_sd
				backgrowdays_sd		= climate.tempbackgrow_sd
				pop_capture_back	   = climate.temppopCapBack
				pop_capture_out		= climate.temppopCapOut
				tempN0				 = climate.tempN0
				tempAllelefile		 = climate.tempAllelefile
				tempClassVarsfile	  = climate.tempClassVarsfile
				assortmateModel		= climate.assortmateModel
				assortmateC			= climate.assortmateC
				subpopmort_mat		 = climate.subpopmort_mat
				comp_coef			  = climate.tempcompcoef
				betas_selection		= climate.tempbetas_selection
				xvars_betas			= climate.tempxvars_betas
				outhabvals			 = climate.tempouthabvals
				backhabvals			= climate.tempbackhabvals
				plastic_signalresp	 = climate.plastic_signalresp
				plastic_behaviorresp   = climate.plastic_behaviorresp
				muterate			   = climate.muterate

				# =====================================================================
				# 10. Demographic Variables
				# =====================================================================
				age_percmort_out_mu	= climate.age_percmort_out
				age_percmort_out_sd	= climate.age_percmort_out_sd
				age_percmort_back_mu   = climate.age_percmort_back
				age_percmort_back_sd   = climate.age_percmort_back_sd
				size_percmort_out_mu   = climate.size_percmort_out
				size_percmort_out_sd   = climate.size_percmort_out_sd
				size_percmort_back_mu  = climate.size_percmort_back
				size_percmort_back_sd  = climate.size_percmort_back_sd
				age_MgOUT			  = climate.age_MgOUT
				age_MgBACK			 = climate.age_MgBACK
				age_S				  = climate.age_S
				age_DispProb		   = climate.age_DispProb
				age_mature			 = climate.age_mature
				age_mu				 = climate.age_mu
				age_sigma			  = climate.age_sigma
				f_leslie_mu			= climate.f_leslie
				f_leslie_std		   = climate.f_leslie_std
				age_capture_out		= climate.age_cap_out
				age_capture_back	   = climate.age_cap_back
				MgBack_patch_prob	  = climate.tempMgBack_patch_prob
				Disperse_patch_prob	= climate.tempDisperse_patch_prob
				disease_vars		   = climate.tupDisease_Vars

				# =====================================================================
				# 11. Maternal & Maturation Settings
				# =====================================================================
				FXXmat_set			 = climate.FXXmat_set
				MXYmat_set			 = climate.MXYmat_set
				MYYmat_set			 = climate.MYYmat_set
				FYYmat_set			 = climate.FYYmat_set
				FXXmat_slope		   = climate.FXXmat_slope
				MXYmat_slope		   = climate.MXYmat_slope
				MYYmat_slope		   = climate.MYYmat_slope
				FYYmat_slope		   = climate.FYYmat_slope
				FXXmat_int			 = climate.FXXmat_int
				MXYmat_int			 = climate.MXYmat_int
				MYYmat_int			 = climate.MYYmat_int
				FYYmat_int			 = climate.FYYmat_int
								
				# ----------------------------------------
				# Introduce new individuals
				# ----------------------------------------
				if (gen != 0 and len(N0_pass[0].split('|')) > 1):							
					SubpopIN = preprocess.AddIndividuals(SubpopIN,tempN0,tempAllelefile,tempClassVarsfile,args.datadir,args.loci,args.alleles,args.sizeans,args.cdevolveans,args.burningen_cdevolve,fitvals,dtype,N0,natal_patches,gen,PopTag,args.sexans,logfHndl,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,args.sexchromo,args.eggFreq_mu,args.eggFreq_sd,disease_vars)

		# -------------------------------------------
		# Update stochastic parameters each year here
		# -------------------------------------------
		tupStoch = preprocess.DoStochasticUpdate(K_mu,K_std,popmort_back_mu,popmort_back_sd,popmort_out_mu,popmort_out_sd,eggmort_mu,eggmort_sd,outsizevals_mu,outsizevals_sd,backsizevals_mu,backsizevals_sd,outgrowdays_mu,outgrowdays_sd,backgrowdays_mu,backgrowdays_sd,age_percmort_out_mu,age_percmort_out_sd,age_percmort_back_mu,age_percmort_back_sd,size_percmort_out_mu,size_percmort_out_sd,size_percmort_back_mu,size_percmort_back_sd,args.egg_percmort_mu,args.egg_percmort_sd,cor_mat,age_mu,age_sigma,f_leslie_mu,f_leslie_std,args.sexchromo,disease_vars)
		K, popmort_back, popmort_out, eggmort_patch, outsizevals, backsizevals, outgrowdays, backgrowdays, age_percmort_out, age_percmort_back, size_percmort_out, size_percmort_back, eggmort_pop, f_ind, f_leslie = tupStoch[:15]
		
		# Print to log
		stringout = 'DoCDClimate(): '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)
	
		# ---------------------------------------
		# Call DoMate() - DoOffspring()
		# ---------------------------------------

		# Timing events: start
		start_time1 = datetime.datetime.now()				
		
		Bearpairs_temp,noOffspring_temp = mate.DoMate(SubpopIN,K,args.freplace,args.mreplace,moveno_mate,thresh_mate,cdmatrix_mate,Track_MateDistCD,xgridpop,ygridpop,Track_MateDistCDstd,Track_FAvgMate,Track_MAvgMate,Track_FSDMate,Track_MSDMate,Track_BreedEvents,gen,sourcePop,scalemax_mate,scalemin_mate,parA_mate,parB_mate,parC_mate,args.Femalepercent_egg,args.sexans,args.selfing,assortmateC,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,assortmateModel,subpopmort_mat,Track_BreedFemales,Track_BreedMales,Track_BreedYYMales,Track_BreedYYFemales,Track_MatureCount, Track_ImmatureCount,Track_ToTFemales,Track_ToTMales,Track_ToTYYMales,Track_ToTYYFemales,args.egg_delay,Bearpairs_temp,natal_patches,args.offno,f_ind,age_sigma,args.sizeans,args.egg_mean_1,args.egg_mean_2,args.egg_mean_ans,args.equalClutch,dtype,eggmort_patch,Track_EggDeaths,eggmort_pop,noOffspring_temp,Track_Births,Track_BirthsMYY,Track_BirthsFYY,args.constMortans,args.outputans,Track_DiseaseStates_AddedInds,disease_vars)
		
		# Print to log
		stringout = 'DoMate() and DoOffspring: '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)
		if Track_ToTFemales[gen][0]==0 or (Track_ToTMales[gen][0] + Track_ToTYYMales[gen][0] + Track_ToTYYFemales[gen][0])==0:
			if temp_extinct == 0:
				print(('There are no more females or males left from species ' + str(args.spcNO) + ' after year '+str(gen)+'.\n'))
			#break
			
		# -------------------------------------------------------------
		# Call 2nd DoUpdate() - grow, age (selection option),egglay,capture, output ind.csv file;no Age0s; ind.csv, disease
		# -------------------------------------------------------------

		# Timing events: start
		start_time1 = datetime.datetime.now()
		SubpopIN = modules.DoUpdate(args.packans,SubpopIN,K,xgridpop,ygridpop,gen,args.nthfile,ithmcrundir,args.loci,args.alleles,logfHndl,'Middle',args.growans,args.cdevolveans,fitvals,args.burningen_cdevolve,age_capture_back,pop_capture_back,Track_CaptureCount_Back,Track_CaptureCount_ClassBack,args.sizeans,age_size_mean,Track_N_back_age,args.eggFreq_mu,args.eggFreq_sd,backsizevals,args.sizeLoo,args.sizeR0,args.size_eqn_1,args.size_eqn_2,args.size_eqn_3,backgrowdays,args.plasticans,args.burningen_plastic,args.timeplastic,plastic_signalresp,args.geneswap,backhabvals,args.sexchromo,Track_DiseaseStates_SecondUpdate,Track_DiseaseStates_AfterDeaths_SecondUpdate,disease_vars)
										
		# Print to log
		stringout = 'Second DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)

		# --------------------------------------------------------------
		# Call DoEmigration() - Age0s adding into population here - OUT vars here
		# --------------------------------------------------------------

		# Timing events: start
		start_time1 = datetime.datetime.now()

		SubpopIN = emigration.DoEmigration(SubpopIN,K,gen,F_EmiDist,M_EmiDist,args.cdevolveans,fitvals,F_EmiDist_sd,M_EmiDist_sd,subpopemigration,SelectionDeathsEmi,DisperseDeathsEmi,args.burningen_cdevolve,MgOut_patch_prob,MgSuccess,AdultNoMg,age_MgOUT,N_Emigration_pop,sourcePop,dtype,setmigrate,args.sizeans,age_size_mean,PackingDeathsEmi,N_Emigration_age,args.loci,muterate,args.mtdna,args.mutationans,args.packans,PackingDeathsEmiAge,args.packpar1,args.timecdevolve,migrate_patches,outsizevals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,args.plasticans,args.burningen_plastic,args.timeplastic,plastic_behaviorresp,noOffspring_temp,Bearpairs_temp,age_size_std,args.Femalepercent_egg,age_mature,args.alleles,args.geneswap,allelst,assortmateModel,args.inheritans_classfiles,args.eggFreq_mu,args.eggFreq_sd,args.sexans,N_beforePack_pop,N_beforePack_age,SelectionDeaths_Age0s,comp_coef,args.XQs,Track_KadjEmi,Track_KadjImmi,args.startcomp,args.spcNO,args.implementcomp,betas_selection,xvars_betas,maxfit,minfit,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,args.sexchromo,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,thresh_FXXOut,thresh_MXYOut,thresh_MYYOut,thresh_FYYOut,scalemin_FXXOut,scalemin_MXYOut,scalemin_MYYOut,scalemin_FYYOut,scalemax_FXXOut,scalemax_MXYOut,scalemax_MYYOut,scalemax_FYYOut,parA_FXXOut,parA_MXYOut,parA_MYYOut,parA_FYYOut,parB_FXXOut,parB_MXYOut,parB_MYYOut,parB_FYYOut,parC_FXXOut,parC_MXYOut,parC_MYYOut,parC_FYYOut,moveno_FXXOut,moveno_MXYOut,moveno_MYYOut,moveno_FYYOut,args.egg_add,args.outputans,age_percmort_out, f_leslie,f_leslie_std,disease_vars,Track_DiseaseStates_AddAge0s)		
					
		# Delete the noOffspring_temp and Bearpairs_temp egg_delay spots used: the first spot in list
		if len(noOffspring_temp) != 0: # But check for extinction
			del(noOffspring_temp[0])
			del(Bearpairs_temp[0])
			# Append a new empty spot for next years cohort
			noOffspring_temp.append(np.asarray([]))
			Bearpairs_temp.append([[-9999,-9999]])				
		
		# Print to log
		stringout = 'DoEmigration(): '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)
				
		# ----------------------------------------
		# Call DoMortality() - when 'Out'
		# ----------------------------------------			

		start_time1 = datetime.datetime.now() # Timing events: start
		SubpopIN = mortality.DoMortality(SubpopIN,K,PopDeathsOUT,	popmort_out,age_percmort_out,gen,N_EmiMortality,AgeDeathsOUT,args.sizeans,age_size_mean,size_percmort_out,SizeDeathsOUT,args.constMortans,args.packans,'OUT',args.sexchromo)
		
		# Print to log
		stringout = 'DoOutMortality(): '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)
		
		# ----------------------------------------------------
		# Call DoUpdate() - grow, mature, capture, and optional output indSample.csv
		# ----------------------------------------------------

		start_time1 = datetime.datetime.now() # Timing events: start
		SubpopIN = modules.DoUpdate(args.packans,SubpopIN,K,xgridpop,ygridpop,gen,args.nthfile,ithmcrundir,args.loci,args.alleles,logfHndl,args.gridsample,args.growans,args.cdevolveans,fitvals,args.burningen_cdevolve,age_capture_out,pop_capture_out,Track_CaptureCount_Out,Track_CaptureCount_ClassOut,args.sizeans,age_size_mean,Track_N_out_age,args.eggFreq_mu,args.eggFreq_sd,outsizevals,args.sizeLoo,args.sizeR0,args.size_eqn_1,args.size_eqn_2,args.size_eqn_3,outgrowdays,args.plasticans,args.burningen_plastic,args.timeplastic,plastic_signalresp,args.geneswap,outhabvals,args.sexchromo,Track_DiseaseStates_ThirdUpdate,Track_DiseaseStates_AfterDeaths_ThirdUpdate,disease_vars,age_mature,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope)
	
		# Print to log
		stringout = 'Third DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)
		
		# ------------------------------------------
		# Call DoImmigration()
		# ------------------------------------------			
		
		start_time1 = datetime.datetime.now() # Timing events: start
		SubpopIN = immigration.DoImmigration(SubpopIN,K,natal_patches,gen,args.cdevolveans,fitvals,subpopimmigration,SelectionDeathsImm,DisperseDeathsImm,args.burningen_cdevolve,Str_patch_prob,StrSuccess,age_S,N_Immigration_pop,dtype,args.sizeans,age_size_mean,PackingDeathsImm,N_Immigration_age,args.packans,PackingDeathsImmAge,args.packpar1,args.homeattempt,args.timecdevolve,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,backsizevals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,args.plasticans,args.burningen_plastic,args.timeplastic,plastic_behaviorresp,age_percmort_back,comp_coef,args.XQs,Track_KadjImmi,Track_KadjEmi,args.startcomp,args.spcNO,args.implementcomp,betas_selection,xvars_betas,maxfit,minfit,f_leslie,f_leslie_std,age_DispProb,cdmatrix_FXXBack,cdmatrix_MXYBack,cdmatrix_MYYBack,cdmatrix_FYYBack,thresh_FXXBack,thresh_MXYBack,thresh_MYYBack,thresh_FYYBack,scalemin_FXXBack,scalemin_MXYBack,scalemin_MYYBack,scalemin_FYYBack,scalemax_FXXBack,scalemax_MXYBack,scalemax_MYYBack,scalemax_FYYBack,parA_FXXBack,parA_MXYBack,parA_MYYBack,parA_FYYBack,parB_FXXBack,parB_MXYBack,parB_MYYBack,parB_FYYBack,parC_FXXBack,parC_MXYBack,parC_MYYBack,parC_FYYBack,moveno_FXXBack,moveno_MXYBack,moveno_MYYBack,moveno_FYYBack,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr,thresh_FXXStr,thresh_MXYStr,thresh_MYYStr,thresh_FYYStr,scalemin_FXXStr,scalemin_MXYStr,scalemin_MYYStr,scalemin_FYYStr,scalemax_FXXStr,scalemax_MXYStr,scalemax_MYYStr,scalemax_FYYStr,parA_FXXStr,parA_MXYStr,parA_MYYStr,parA_FYYStr,parB_FXXStr,parB_MXYStr,parB_MYYStr,parB_FYYStr,parC_FXXStr,parC_MXYStr,parC_MYYStr,parC_FYYStr,moveno_FXXStr,moveno_MXYStr,moveno_MYYStr,moveno_FYYStr,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD,thresh_FXXLD,thresh_MXYLD,thresh_MYYLD,thresh_FYYLD,scalemin_FXXLD,scalemin_MXYLD,scalemin_MYYLD,scalemin_FYYLD,scalemax_FXXLD,scalemax_MXYLD,scalemax_MYYLD,scalemax_FYYLD,parA_FXXLD,parA_MXYLD,parA_MYYLD,parA_FYYLD,parB_FXXLD,parB_MXYLD,parB_MYYLD,parB_FYYLD,parC_FXXLD,parC_MXYLD,parC_MYYLD,parC_FYYLD,moveno_FXXLD,moveno_MXYLD,moveno_MYYLD,moveno_FYYLD,args.sexchromo,age_MgBACK,MgBack_patch_prob,Disperse_patch_prob,MgOut_patch_prob,age_MgOUT,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,migrate_patches,args.egg_add,args.outputans)
						
		# Print to log
		stringout = 'DoImmigration(): '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)
					
		# ------------------------------------------
		# Call DoMortality() - when 'Back'
		# ------------------------------------------
		
		# Timing events: start
		start_time1 = datetime.datetime.now()
		SubpopIN = mortality.DoMortality(SubpopIN,K,PopDeathsIN,popmort_back,age_percmort_back,gen,N_ImmiMortality,AgeDeathsIN,args.sizeans,age_size_mean,size_percmort_back,SizeDeathsIN,args.constMortans,args.packans,'BACK',args.sexchromo)
		
		# Print to log
		stringout = 'DoInMortality(): '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)
		
		# ---------------------------------
		# Call GetMetrics()
		# ---------------------------------
		
		# Timing events: start
		start_time1 = datetime.datetime.now()
		modules.GetMetrics(SubpopIN,K,Track_N_Init_pop,Track_K,args.loci,args.alleles,gen+1,Track_Ho,Track_Alleles,Track_He,Track_p1,Track_p2,Track_q1,Track_q2,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,Track_N_Init_age,args.sizeans,age_size_mean,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,args.packans,RDispersers,IDispersers,xvars_betas,betas_selection,maxfit,minfit,args.cdevolveans,disease_vars,Track_DiseaseStates_pop,Track_DiseaseStates_EnvRes)
		
		# Print to log
		stringout = 'GetMetrics(): '+str(datetime.datetime.now() -start_time1) + ''
		modules.logMsg(logfHndl,stringout)
								
		# Print to log
		stringout = 'End Generation/Year Loop'+str(gen)+': '+str(datetime.datetime.now() -start_timeGen) + '\n'
		modules.logMsg(logfHndl,stringout)
		if mp.current_process().name in ("S0", "MainProcess") or mp.current_process().name.endswith("-1"):
			print(stringout)
	
	# End::generation loop
				
	# ------------------------------------------
	# Call DoPostProcess()
	# ------------------------------------------

	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	postprocess.DoPostProcess(ithmcrundir,args.loci,args.alleles,args.looptime,\
	Track_ToTFemales,Track_ToTMales,Track_BreedFemales,Track_BreedMales,Track_Births,PopDeathsIN,\
	PopDeathsOUT,Track_Alleles,Track_He,Track_Ho,Track_MateDistCD,Track_MateDistCDstd,args.nthfile,logfHndl,\
	Track_p1,Track_p2,Track_q1,Track_q2,subpopemigration,\
	subpopimmigration,Track_FAvgMate,Track_MAvgMate,Track_FSDMate,Track_MSDMate,\
	SelectionDeathsEmi,SelectionDeathsImm,\
	DisperseDeathsEmi,DisperseDeathsImm,\
	Track_BreedEvents,args.gridformat,\
	MgSuccess,AdultNoMg,StrSuccess,\
	Track_EggDeaths,Track_K,Track_N_Init_pop,N_Emigration_pop,N_EmiMortality,N_Immigration_pop,N_ImmiMortality,Track_DiseaseStates_pop,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,PackingDeathsEmi,PackingDeathsImm,Track_N_Init_age,N_Emigration_age,N_Immigration_age,AgeDeathsOUT,AgeDeathsIN,PackingDeathsEmiAge,PackingDeathsImmAge,Track_MatureCount,Track_ImmatureCount,Track_N_back_age,Track_N_out_age,args.outputans,gen,Track_CaptureCount_Back,Track_CaptureCount_ClassBack,Track_CaptureCount_Out,Track_CaptureCount_ClassOut,age_size_mean,args.sizeans,ClassSizes_Mean,ClassSizes_Std,Track_N_Init_class,SizeDeathsOUT,SizeDeathsIN,N_beforePack_pop,N_beforePack_age,SelectionDeaths_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,F_EmiDist,M_EmiDist,F_EmiDist_sd,M_EmiDist_sd,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,Track_ToTYYMales,Track_BreedYYMales,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,RDispersers,IDispersers,Track_BirthsMYY,Track_KadjEmi,Track_KadjImmi,Track_ToTYYFemales,Track_BirthsFYY,Track_BreedYYFemales,disease_vars['ImpDisease'],Track_DiseaseStates_SecondUpdate,Track_DiseaseStates_ThirdUpdate,Track_DiseaseStates_AddAge0s,Track_DiseaseStates_AddedInds,Track_DiseaseStates_AfterDeaths_SecondUpdate,Track_DiseaseStates_AfterDeaths_ThirdUpdate,Track_DiseaseStates_EnvRes,disease_vars)
	# Print to log
	stringout = 'DoPostProcess(): '+str(datetime.datetime.now() -start_time1) + ''
	modules.logMsg(logfHndl,stringout)
	if mp.current_process().name in ("S0", "MainProcess") or mp.current_process().name.endswith("-1"):
		print(stringout)
	
	# Print to log
	stringout = 'End Monte Carlo Loop'+str(args.ithmcrun)+': '+str(datetime.datetime.now() -start_timeMC) + '\n'
	modules.logMsg(logfHndl,stringout)
	if mp.current_process().name in ("S0", "MainProcess") or mp.current_process().name.endswith("-1"):
		print(stringout)
	# End::Monte Carlo Loop

# Function to custom name the multiprocessing processes to be able to identify which process will print messages to terminal.
def worker_init(counter):
	global worker_id
	with counter.get_lock():
		counter.value += 1
		worker_id = counter.value
	# Optional: Rename the process for easier debugging
	mp.current_process().name = f"spawnworker-{worker_id}"
