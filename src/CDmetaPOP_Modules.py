# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Modules.py
# Author: Erin L Landguth
# Created: June 2010
# Description: This is the function/module file for CDmetaPOP vX
# --------------------------------------------------------------------------------------------------

# Import Modules with Except/Try statements

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."

# CDmetaPOP functions
try:
	from CDmetaPOP_PostProcess import DoOutput
except ImportError:
	raise ImportError, "CDmetaPOP PreProcess required."

# Scipy functions
try:
	import scipy.stats
except ImportError:
	raise ImportError, "Scipy required."
	
# Python specific functions
import os, random, copy, pdb, sys
from ast import literal_eval 

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False


# --------------------------------------------------------------------------
def PrepTextFile(textpath):
	'''
	PrepTextFile() - Prepare the input files
	'''
	
	return textpath.strip('\n').strip('\r')
	
	# End::PrepTextFile()

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
		print("%s"%(msg))
		
	# End::logMsg()

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()	
	
# ---------------------------------------------------------------------------------------------------	 
def w_choice_general(lst):
	'''
	w_choice_general()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(x[1] for x in lst)
	n=random.uniform(0,wtotal)
	count = 0
	for item, weight in lst:
		if n < weight:
			break
		n = n-weight
		count = count + 1
	return item,count
	
	#End::w_choice_general()

# ---------------------------------------------------------------------------------------------------	
def Do1LocusSelection(fitvals,genes,location):
	'''
	Do1LocusSelection()
	This function calculates individual differential mortality, ie,
	viability selection, for the 1-locus selection model.
	'''
	
	# If L0A0|L0A0 -- loci under selection:
	if int(genes[0]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][0])/100.
																
	# If L0A0|L0A1 -- loci under selection:
	elif int(genes[0]) == 1 and int(genes[1]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][1])/100.
																															
	# If L0A1|L0A1 -- loci under selection
	elif int(genes[1]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals[location][2])/100.
		
	# Another genotype
	else:
		differentialmortality = 0.0		
	
	return differentialmortality
	
	# End::Do1LocusSelection()
	
# ---------------------------------------------------------------------------------------------------	
def Do2LocusSelection(fitvals,genes,location):
	'''
	Do2LocusSelection()
	This function calculates individual differential mortality, ie,
	viability selection, for the 2-locus selection model.
	'''
	
	# If L0A0|L0A0|L1A0|L1A0 - AABB -- loci under selection:
	if int(genes[0][0]) == 2 and int(genes[1][0]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][0])/100.
									
	# If L0A0|L0A1|L1A0|L1A0 - AaBB -- loci under selection:
	elif int(genes[0][0]) == 1 and int(genes[0][1]) == 1 and int(genes[1][0]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][1])/100.																															
	# If L0A1|L0A1|L1A0|L1A0 - aaBB -- loci under selection
	elif int(genes[0][1]) == 2 and int(genes[1][0]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals[location][2])/100.
									
	# If L0A0|L0A0|L1A0|L1A1 - AABb -- loci under selection:
	elif int(genes[0][0]) == 2 and int(genes[1][0]) == 1 and int(genes[1][1]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][3])/100.
									
	# If L0A0|L0A1|L1A0|L1A1 - AaBb -- loci under selection:
	elif int(genes[0][0]) == 1 and int(genes[0][1]) == 1 and int(genes[1][0]) == 1 and int(genes[1][1]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][4])/100.
																	
	# If L0A1|L0A1|L1A0|L1A1 - aaBb -- loci under selection
	elif int(genes[0][1]) == 2 and int(genes[1][0]) == 1 and int(genes[1][1]) == 1:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals[location][5])/100.
	
	# If L0A0|L0A0|L1A1|L1A1 - AAbb -- loci under selection:
	elif int(genes[0][0]) == 2 and int(genes[1][1]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][6])/100.
									
	# If L0A0|L0A1|L1A1|L1A1 - Aabb -- loci under selection:
	elif int(genes[0][0]) == 1 and int(genes[0][1]) == 1 and int(genes[1][1]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][7])/100.
																															
	# If L0A1|L0A1|L1A1|L1A1 - aabb -- loci under selection
	elif int(genes[0][1]) == 2 and int(genes[1][1]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals[location][8])/100.
	
	# Another genotype
	else:
		differentialmortality = 0.0		
	
		
	return differentialmortality
	
	# End::Do2LocusSelection()
	
# ---------------------------------------------------------------------------------------------------	 
def GetMetrics(SubpopIN,K,Population,K_track,loci,alleles,gen,Ho,Alleles,He,p1,p2,q1,q2,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,ToTMales,ToTFemales,BreedMales,BreedFemales,N_Age,MatureCount,ImmatureCount,sizecall,size_mean,ClassSizes_Mean,ClassSizes_Std,N_Class,sexans,SNPans):
	'''
	GetMetrics()
	This function summarizes the genotypes and
	produces genetic metrics.
	Ho - Observed heterozygosity per generation
	He - Expected heterozygoisty per generation
	Alleles - Total number of unique alleles in genotype*individuals
	'''
		
	# List for total, left, and right
	unique_alleles = Alleles
	
	# Get allele location as seqence from alleles array
	allele_numbers = []
	for i in xrange(loci):
		for j in xrange(alleles[i]):
			allele_numbers.append(j)
	allele_numbers = np.asarray(allele_numbers)
	
	# Get length of classes
	classno = len(size_mean[0])
	
	# Add spots for Age tracking
	N_Age.append([])
	N_Class.append([])
	AgeSizes_Mean.append([])
	AgeSizes_Std.append([])
	ClassSizes_Mean.append([])
	ClassSizes_Std.append([])
	N_Age[gen] = [[] for x in xrange(0,classno)]
	N_Class[gen] = [[] for x in xrange(0,classno)]
	AgeSizes_Mean[gen] = [[] for x in xrange(0,classno)]
	AgeSizes_Std[gen] = [[] for x in xrange(0,classno)]
	ClassSizes_Mean[gen] = [[] for x in xrange(0,classno)]
	ClassSizes_Std[gen] = [[] for x in xrange(0,classno)]
		
	# Extract the genes information from SubpopIN
	tempgenes = []
	tempgenesPop = []
	Population.append([]) # Storage numbers - add spot for generation
	K_track.append([]) # Storage numbers - add spot for generation
	# Create a list to store the subpopulation grid number information
	subgrids = []
	all_freq_sub = []
	ho_count_sub = []
	ho_sub = []
	all_freq_sq_sub = []
	homozygosity_sub = []
	he_sub = []
	sumsubpopsHo = []
	sumsubpopsHe = []
	alleles_sub = []
	tempInf = []
	Residors.append([])
	Strayers1.append([])
	Strayers2.append([])
	Immigrators.append([])
	PopSizes_Mean.append([])
	PopSizes_Std.append([])
	ToTMales.append([]) #Storage add spot for generation
	ToTFemales.append([]) #Storage add spot for generation
	BreedMales.append([]) #Storage add spot for generation
	BreedFemales.append([]) #Storage add spot for generation
	MatureCount.append([])
	ImmatureCount.append([])
	
	# For each supopulation	
	for isub in xrange(len(K)):
		tempgenesPop.append([])
		# For each individual
		for iind in xrange(len(SubpopIN[isub])):
			# Extract genes and convert back to list
			tempgenes.append(literal_eval(SubpopIN[isub][iind]['genes']))
			tempgenesPop[isub].append(literal_eval(SubpopIN[isub][iind]['genes']))
		
		# Add information to Population tracker
		Population[gen].append(len(SubpopIN[isub]))
		K_track[gen].append(K[isub])
		PopSizes_Mean[gen].append(np.mean(SubpopIN[isub]['size']))
		PopSizes_Std[gen].append(np.std(SubpopIN[isub]['size']))
		tempInf.append(list(SubpopIN[isub]['infection']))
		# Extract the disperser type
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'R' in val])
		Residors[gen].append(len(tempname))
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'I' in val])
		Immigrators[gen].append(len(tempname))
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'S' in val])
		Strayers1[gen].append(len(tempname))
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'Z' in val])
		Strayers2[gen].append(len(tempname))
		indexF = np.where(SubpopIN[isub]['sex']==0)[0]
		indexM = np.where(SubpopIN[isub]['sex']==1)[0]
		allfemales = SubpopIN[isub][indexF]
		allmales = SubpopIN[isub][indexM]
		# Get reproduction age individuals
		indexFage = np.where(allfemales['mature'] == 1)[0]
		indexMage = np.where(allmales['mature'] == 1)[0]
		
		if sexans == 'Y':
			# Storage tracking
			ToTMales[gen].append(len(indexM)) 
			ToTFemales[gen].append(len(indexF))
			BreedMales[gen].append(len(indexMage))
			BreedFemales[gen].append(len(indexFage))
		else:
			# Storage tracking
			ToTMales[gen].append(len(indexM)+len(indexF)) 
			ToTFemales[gen].append(len(indexM)+len(indexF))
			BreedMales[gen].append(len(indexMage)+len(indexFage))
			BreedFemales[gen].append(len(indexMage)+len(indexFage))
		MatureCount[gen].append(sum(SubpopIN[isub]['mature']))
		ImmatureCount[gen].append(len(SubpopIN[isub]['mature'])-sum(SubpopIN[isub]['mature']))
		
		# Size class counting		
		# Switch here for size or age control
		if sizecall == 'Y':
			size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
			age_adjusted = np.searchsorted(size_mean_middles, SubpopIN[isub]['size'])
		else:
			# Count up each uniages
			age_adjusted = SubpopIN[isub]['age']
	
		# Tracking age N
		for iage in xrange(len(AgeSizes_Mean[gen])):
			sizeindex = np.where(age_adjusted==iage)[0]
			ageindex = np.where(SubpopIN[isub]['age'] == iage)[0]
			if len(sizeindex) == 0:
				ClassSizes_Mean[gen][iage].append([0])
				N_Class[gen][iage].append(0)
			else:
				ClassSizes_Mean[gen][iage].append(SubpopIN[isub][sizeindex]['size'].tolist())
				N_Class[gen][iage].append(len(SubpopIN[isub][sizeindex]['size'].tolist()))
			if len(ageindex) == 0:
				AgeSizes_Mean[gen][iage].append([0])
				N_Age[gen][iage].append(0)
			else:				
				AgeSizes_Mean[gen][iage].append(SubpopIN[isub][ageindex]['size'].tolist())	
				N_Age[gen][iage].append(len(SubpopIN[isub][ageindex]['size'].tolist()))
		
		# Special case where age class is greater than lastage
		sizeindex = np.where(age_adjusted > iage)[0]
		if len(sizeindex) == 0:
			ClassSizes_Mean[gen][iage].append([0])
			N_Class[gen][iage].append(0)
		else: # Add them to last class
			ClassSizes_Mean[gen][iage].append(SubpopIN[isub][sizeindex]['size'].tolist())
			N_Class[gen][iage].append(len(SubpopIN[isub][sizeindex]['size'].tolist()))
		ageindex = np.where(SubpopIN[isub]['age'] > iage)[0]
		if len(ageindex) == 0:
			AgeSizes_Mean[gen][iage].append([0])
			N_Age[gen][iage].append(0)
		else: # Add them to last class
			AgeSizes_Mean[gen][iage].append(SubpopIN[isub][ageindex]['size'].tolist())	
			N_Age[gen][iage].append(len(SubpopIN[isub][ageindex]['size'].tolist()))
		
		# Temp storage for GetMetrics()
		subgrids.append([])
		all_freq_sub.append([])
		ho_count_sub.append([])
		ho_sub.append([])
		all_freq_sq_sub.append([])
		homozygosity_sub.append([])
		he_sub.append([])
		alleles_sub.append([])
	# Add Population total
	Population[gen].insert(0,sum(Population[gen]))
	K_track[gen].insert(0,sum(K))
	tempInf = sum(tempInf,[])
	Infected.append(sum(tempInf))
	# Add Population totals
	ToTMales[gen].insert(0,sum(ToTMales[gen]))
	ToTFemales[gen].insert(0,sum(ToTFemales[gen]))
	BreedMales[gen].insert(0,sum(BreedMales[gen]))
	BreedFemales[gen].insert(0,sum(BreedFemales[gen]))
	
	# Add Count totals
	MatureCount[gen] = sum(MatureCount[gen])
	ImmatureCount[gen] = sum(ImmatureCount[gen])
	
	# Age tracking	
	for iage in xrange(len(AgeSizes_Mean[gen])):
		tempagesize = np.asarray(sum(AgeSizes_Mean[gen][iage],[]))
		tempagesize = tempagesize[np.where(tempagesize != 0)[0]]
		tempsize = np.asarray(sum(ClassSizes_Mean[gen][iage],[]))
		tempsize = tempsize[np.where(tempsize != 0)[0]]
		N_Age[gen][iage] = sum(N_Age[gen][iage])
		N_Class[gen][iage] = sum(N_Class[gen][iage])		
		AgeSizes_Mean[gen][iage] = np.mean(tempagesize)
		AgeSizes_Std[gen][iage] = np.std(tempagesize)
		ClassSizes_Mean[gen][iage] = np.mean(tempsize)
		ClassSizes_Std[gen][iage] = np.std(tempsize)
	
	# Cast genes as an numpy array as byte type
	genes_array_woNA = np.asarray(tempgenes,dtype='float')
	
	# The total number of alleles
	total_alleles = len(allele_numbers)
	
	# Get unique number of subpops
	nosubpops = len(K)
	
	# And then get the number of filled grids
	filledgrids = Population[gen][0]		
		
	# Get allele frequency for total
	if filledgrids != 0:
		all_freq_tot = np.asarray(np.sum(genes_array_woNA,axis=0),dtype = 'float').reshape(total_alleles)
		if SNPans == 'N':
			all_freq_tot = all_freq_tot/(2*filledgrids)
		else:
			all_freq_tot = all_freq_tot/(filledgrids)
	else:
		all_freq_tot = np.zeros(total_alleles,float)
	
	# Get allele frequency for subpopulations
	for isub in xrange(nosubpops):
		if Population[gen][isub+1] != 0:
			# Cast genes as an numpy array as byte type
			genes_array_subpop = np.asarray(tempgenesPop[isub],dtype='float')
			all_freq_sub[isub].append(np.asarray(np.sum(genes_array_subpop,axis=0),dtype = 'float').reshape(total_alleles))
			if SNPans == 'N':
				all_freq_sub[isub] = all_freq_sub[isub][0]/(2*Population[gen][isub+1])
			else:
				all_freq_sub[isub] = all_freq_sub[isub][0]/(Population[gen][isub+1])
		else:
			all_freq_sub[isub].append(np.zeros(total_alleles,float))
			all_freq_sub[isub] = all_freq_sub[isub][0]
	
	# Create an array to fill up with allele frequencies - only for total
	all_freq_list = np.zeros((total_alleles,2))		
	all_freq_list[:,0] = allele_numbers
	all_freq_list[:,1] = all_freq_tot	
	
	#Calculate the number of homogenous alleles for total
	# Note: this is not calculated for Ho
	ho_count_tot = np.array(genes_array_woNA==2).sum()
	# Calculate the number of homogenous alleles in each subpop
	for isub in xrange(nosubpops):
		# Cast genes as an numpy array as byte type
		genes_array_subpop = np.asarray(tempgenesPop[isub],dtype='float')
		ho_count_sub[isub].append(np.array(genes_array_subpop==2).sum())
	
	# Calculate the observed het for total
	if filledgrids != 0:
		ho_tot = (float(filledgrids*loci - ho_count_tot)/(loci*filledgrids))
	else:
		ho_tot = 0.0
	# Append Ho information (Observed Het)
	Ho.append([ho_tot])		
	# Calculate the observed het in each subpop
	for isub in xrange(nosubpops):
		if Population[gen][isub+1] != 0:
			ho_sub[isub].append((float(Population[gen][isub+1]*loci - ho_count_sub[isub][0])/(loci*Population[gen][isub+1])))
		else:
			ho_sub[isub].append(0.0)
		# Append Ho information (Observed Het)
		Ho[gen].append(ho_sub[isub][0])
		
	# Get the sqare of the allele frequency for total
	all_freq_sq_tot = all_freq_tot**2
	# Calculate the homozygosity for total populations
	homozygosity_tot = sum(all_freq_sq_tot)/loci
	# Get total number of alleles
	alleles_tot = np.array(all_freq_tot>0.).sum()
	# Append allele total information
	unique_alleles.append([alleles_tot])
	# Get the square of the allele frequency for subpops
	# Calculate the homozygosity for subpopulations
	# Get the total number of alleles in each subpop
	for isub in xrange(nosubpops):
		all_freq_sq_sub[isub].append(all_freq_sub[isub]**2)
		homozygosity_sub[isub].append(sum(all_freq_sq_sub[isub][0])/loci)
		alleles_sub[isub].append(np.array(all_freq_sub[isub]>0.).sum())
		# Append allele total information
		unique_alleles[gen].append(alleles_sub[isub][0])
		
	# Store He for [Total]
	if filledgrids != 0:
		he_tot = (1. - homozygosity_tot)
	else:
		he_tot = 0.0
	# Append He information (Expected Het)
	He.append([he_tot])
	# Store He for subpopulations
	for isub in xrange(nosubpops):
		if Population[gen][isub+1] != 0:
			he_sub[isub].append(1. - homozygosity_sub[isub][0])
		else:
			he_sub[isub].append(0.0)
		# Append He information (Expected Het)
		He[gen].append(he_sub[isub][0])
			
	# Get allele frequency totals for selection section
	p1.append(all_freq_tot[0])
	p2.append(all_freq_tot[1])
	q1.append(all_freq_tot[2])
	q2.append(all_freq_tot[3])
	
	# Here we exit function if there are are no Total Females or Males
	if ToTFemales[gen][0]==0:						
		return 
	if ToTMales[gen][0]==0:
		return
			
	#End::GetMetrics()
	
# ---------------------------------------------------------------------------------------------------	 
def InheritGenes(gen,offspring,loci,muterate,mtdna,mutationans,K,dtype,geneswap,allelst,SNPans):
	'''
	InheritGenes()
	Pass along gentic information to survived offspring from parents
	Input: offspring [femalegenes,malegenes,NatalPop,EmiPop,ImmiPop,age,sex,size,infection,name]
	Output: SubpopIN_Age0	[NatalPop,EmiPop,ImmiPop,age,sex,size,infection,name,genes]		
	'''		

	# Create list for appending
	Age0_keep = []
	
	# If there are offspring
	if len(offspring) != 0:
	
		# Begin loop through offspring
		for i in xrange(len(offspring)):
			
			# Check for geneswap time
			if gen >= geneswap and geneswap != 'N':
				
				# Temp storage for i's mother's and father's genes, then offspring genes
				mothergenes=literal_eval(offspring[i]['Mother'])
				fathergenes=literal_eval(offspring[i]['Father'])
				offgenes = []
				
				# Loop through each locus
				for jspot in xrange(loci):
					
					# Temporary index storage
					tempindfather = []
					tempindmother = []
					offgenes.append([])
										
					# Loop through each allele-mother and father have same len of alleles
					for kspot in xrange(len(fathergenes[jspot])):
												
						# Check the homogeneous 2 case - in father genes
						if int(fathergenes[jspot][kspot])==2:
							
							# Get a random number for allele mutation
							mutationrandno = rand()
							
							# Check if random number is less than or equal to muterate
							if mutationrandno <= muterate:
							
								# If backward and forward mutation
								if mutationans == 'random':
									# Randomly choose another allele								
									randallelespot = int(len(fathergenes[jspot])*rand())
									while randallelespot == kspot:
										randallelespot = int(len(fathergenes[jspot])*rand())
									tempindfather.append(randallelespot)
									
								# If just forward mutation
								elif mutationans == 'forward':
									if kspot != len(fathergenes[jspot])-1:
										tempindfather.append(kspot+1)
									else:
										# THen no mutation
										tempindfather.append(kspot)

								# If just forward mutation
								elif mutationans == 'backward':
									if kspot != 0:
										tempindfather.append(kspot-1)
									else:
										# THen no mutation
										tempindfather.append(kspot)
										
								# If forward and backward mutation
								elif mutationans == 'forwardbackward':
									# Then random forward or backward step
									randstep = rand()
									# To go left, but it can't be the first allele
									if randstep < 0.5 and kspot != 0:
										tempindfather.append(kspot-1)
									# To go right, but it can't be the last allele
									elif randstep >= 0.5 and kspot != len(fathergenes[jspot])-1:
										tempindfather.append(kspot+1)
									else:
										# THen no mutation
										tempindfather.append(kspot)
								
								# If forward mutation in A and backward mutation for b (A -> a, b -> B)
								elif mutationans == 'forwardAbackwardBrandomN':
									if jspot == 0 and kspot == 0:
										tempindfather.append(kspot+1)
									elif jspot == 1 and kspot == 1:
										tempindfather.append(kspot-1)
									elif jspot != 0 and jspot != 1:
										# Randomly choose another allele								
										randallelespot = int(len(fathergenes[jspot])*rand())
										while randallelespot == kspot:
											randallelespot = int(len(fathergenes[jspot])*rand())
										tempindfather.append(randallelespot)
									else:
										# THen no mutation
										tempindfather.append(kspot)	
								
								# No other mutation models matched
								else:
									print('The mutation model does not exist.')
									sys.exit(-1)								
							
							# and if random number is not less than or equal to muterate
							else:
								tempindfather.append(kspot)
						
						# Check the homogeneous 2 case - in mother genes	
						if int(mothergenes[jspot][kspot])==2:
							
							# Get a random number for allele mutation
							mutationrandno = rand()
							
							# Check if random number is less than or equal to muterate
							if mutationrandno <= muterate:
							
								# If backward and forward mutation
								if mutationans == 'random':
									# Randomly choose another allele
									randallelespot = int(len(mothergenes[jspot])*rand())
									while randallelespot == kspot:
										randallelespot = int(len(mothergenes[jspot])*rand())
									tempindmother.append(randallelespot)
									
								# If just forward mutation
								elif mutationans == 'forward':
									if kspot != len(mothergenes[jspot])-1:
										tempindmother.append(kspot+1)
									else:
										# THen no mutation
										tempindmother.append(kspot)

								# If just backward mutation
								elif mutationans == 'backward':
									if kspot != 0:
										tempindmother.append(kspot-1)
									else:
										# THen no mutation
										tempindmother.append(kspot)
										
								# If forward or backward step mutation
								elif mutationans == 'forwardbackward':
									# Random draw for forward or backward
									randstep = rand()
									# TO go left, but it can't be the first allele
									if randstep < 0.5 and kspot != 0:
										tempindmother.append(kspot-1)
									# To go right, but it can't be the last allele
									elif randstep >= 0.5 and kspot != len(fathergenes[jspot])-1:
										tempindmother.append(kspot+1)
									else:
										# THen no mutation
										tempindmother.append(kspot)
								
								# If forward mutation in A and backward mutation in b
								elif mutationans == 'forwardAbackwardBrandomN':
									if jspot == 0 and kspot == 0:
										tempindmother.append(kspot+1)
									elif jspot == 1 and kspot == 1:
										tempindmother.append(kspot-1)
									if jspot != 0 and jspot != 1:
										# Randomly choose another allele
										randallelespot = int(len(mothergenes[jspot])*rand())
										while randallelespot == kspot:
											randallelespot = int(len(mothergenes[jspot])*rand())
										tempindmother.append(randallelespot)
									else:
										# THen no mutation
										tempindmother.append(kspot)
								
								# No other mutation models matched
								else:
									print('The mutation model does not exist.')
									sys.exit(-1)					

							# and if random number is not less than or equal to muterate
							else:
								tempindmother.append(kspot)
							
						# Check the hetero 1 case for father genes
						if int(fathergenes[jspot][kspot])==1:
							tempindfather.append(kspot)
							
						# Check the hetero 1 case for mother genes
						if int(mothergenes[jspot][kspot])==1:
							tempindmother.append(kspot)
					
					# Check if the tempindex has a length of 2 (which means it was not homo at
					#	at this locus), then randomly select one of them, and then check for mutation
					# Check from father genes
					if len(tempindfather) == 2:
					
						# Then randomly select one of the homo alleles
						temprandnofather = int(2*rand())
						
						# Delete from list
						del tempindfather[temprandnofather]
						thealleleselected = tempindfather[0]
						
						# Get a random number for allele mutation
						mutationrandno = rand()
						
						# Check if random number is less than or equal to muterate
						if mutationrandno <= muterate:
							
							# If backward and forward mutation
							if mutationans == 'random':
								# Randomly choose another allele
								randallelespot = int(len(fathergenes[jspot])*rand())
								while randallelespot == thealleleselected:
									randallelespot = int(len(fathergenes[jspot])*rand())
								# and then reassign this spot
								tempindfather[0] = randallelespot
								
							# Else if just forward mutation
							elif mutationans == 'forward':
								if thealleleselected != len(fathergenes[jspot])-1:
									# and then reassign this spot
									tempindfather[0] = thealleleselected+1
																
							# Else if just backward mutation
							elif mutationans == 'backward':
								if thealleleselected != 0:
									# and then reassign this spot
									tempindfather[0] = thealleleselected-1
														
							# Else if forward backward mutation
							elif mutationans == 'forwardbackward':
								# Random draw for forward or backward
								randstep = rand()
								# TO go left, but it can't be the first allele
								if randstep < 0.5 and thealleleselected != 0:
									# and then reassign this spot
									tempindfather[0] = thealleleselected-1
								# To go right, but it can't be the last allele
								elif randstep >= 0.5 and thealleleselected != len(fathergenes[jspot])-1:
									# and then reassign this spot
									tempindfather[0] = thealleleselected+1
							
							# Else if forward mutation for A and backward mutation for b
							elif mutationans == 'forwardAbackwardBrandomN':
								if jspot == 0 and thealleleselected == 0:
									# and then reassign this spot
									tempindfather[0] = thealleleselected+1
								elif jspot == 1 and thealleleselected == 1:
									tempindfather.append(kspot-1)
								elif jspot != 0 and jspot != 1:
									# Randomly choose another allele
									randallelespot = int(len(fathergenes[jspot])*rand())
									while randallelespot == thealleleselected:
										randallelespot = int(len(fathergenes[jspot])*rand())
									# and then reassign this spot
									tempindfather[0] = randallelespot
								else:
									# THen no mutation
									tempindfather[0] = thealleleselected	
										
							# No other mutation models matched
							else:
								print('The mutation model does not exist.')
								sys.exit(-1)		
					
					# Check from mother genes
					if len(tempindmother) == 2:
					
						# THen randomly select on of the homo alleles
						temprandnomother = int(2*rand())
						# Delete from list
						del tempindmother[temprandnomother]
						thealleleselected = tempindmother[0]
						
						# Get a random number for allele mutation
						mutationrandno = rand()
						
						# Check if random number is less than or equal to muterate
						if mutationrandno <= muterate:
						
							# If backward and forward mutation
							if mutationans == 'random':
								# Randomly choose another allele
								randallelespot = int(len(mothergenes[jspot])*rand())
								while randallelespot == thealleleselected:
									randallelespot = int(len(mothergenes[jspot])*rand())
								# and then reassign this spot
								tempindmother[0] = randallelespot
								
							# Else if just forward mutation
							elif mutationans == 'forward':
								if thealleleselected != len(mothergenes[jspot])-1:
									# and then reassign this spot
									tempindmother[0] = thealleleselected+1	
																
							# Else if just backward mutation
							elif mutationans == 'backward':
								if thealleleselected != 0:
									# and then reassign this spot
									tempindmother[0] = thealleleselected-1
																
							# Else if forward backward mutation
							elif mutationans == 'forwardbackward':
								# Random draw for forward or backward
								randstep = rand()
								# TO go left, but it can't be the first allele
								if randstep < 0.5 and thealleleselected != 0:
									# and then reassign this spot
									tempindmother[0] = thealleleselected-1
								# To go right, but it can't be the last allele
								elif randstep >= 0.5 and thealleleselected != len(mothergenes[jspot])-1:
									# and then reassign this spot
									tempindmother[0] = thealleleselected+1
							
							# Else if forward mutation for A and backward mutation for b
							elif mutationans == 'forwardAbackwardBrandomN':
								if jspot == 0 and thealleleselected == 0:
									# and then reassign this spot
									tempindmother[0] = thealleleselected+1
								elif jspot == 1 and thealleleselected == 1:
									tempindmother.append(kspot-1)
								elif jspot != 0 and jspot != 1:
									# Randomly choose another allele
									randallelespot = int(len(mothergenes[jspot])*rand())
									while randallelespot == thealleleselected:
										randallelespot = int(len(mothergenes[jspot])*rand())
									# and then reassign this spot
									tempindmother[0] = randallelespot
								else:
									# THen no mutation
									tempindmother[0] = thealleleselected
							
							# No other mutation models matched
							else:
								print('The mutation model does not exist.')
								sys.exit(-1)		
					
					# For SNPs only, randomly choose mother or father
					if SNPans == 'Y':
						snprand = rand()
						if snprand < 0.5:
							tempind = tempindmother
						else:
							tempind = tempindfather					
					
					# Now write to offspring genes array the selected alleles in locus j
					for kspot in xrange(len(fathergenes[jspot])):
						
						if SNPans == 'N':
							# Need a storage temp for genes
							# Hetero case 1 AB
							if tempindfather[0] == kspot and tempindmother[0] != kspot:
								offgenes[jspot].append(1)
							# Homo case AA or BB
							elif tempindfather[0] == kspot and tempindmother[0] == kspot:
								offgenes[jspot].append(2)
							# Hetero case 2 BA
							elif tempindmother[0] == kspot and tempindfather[0] != kspot:
								offgenes[jspot].append(1)
							# Or nothing there at all
							elif tempindmother[0] != kspot and tempindfather[0] != kspot:
								offgenes[jspot].append(0)
						else:
							if tempind[0] == kspot:
								offgenes[jspot].append(1)
							else:
								offgenes[jspot].append(0)
				
				# If mtdna is turned on, then erase the last loci and force it to be mothergenes
				if mtdna == 'Y' and jspot == loci-1:
					print('Check this module')
					sys.exit(-1)
					# Force last locus to be mothergenes
					for imtdna in xrange(len(mothergenes[loci-1])):
						offgenes[jspot] = mothergenes[loci-1][imtdna]
			
			# If not in geneswap time, then initialize with allelst
			else:
				
				# Get genes - For each loci:
				sourcepop = int(offspring[i]['NatalPop'])-1
				offgenes = [] # Storage
				
				# Loop through each locus
				for j in xrange(loci):
									
					# Store genes loci spot
					offgenes.append([])
					
					# Take a random draw from the w_choice function at jth locus
					rand1 = w_choice_general(allelst[sourcepop][j])[0]
					rand2 = w_choice_general(allelst[sourcepop][j])[0]
										
					# 	1s = heterozygous at that locus
					#	2s = homozygous at that locus
					#	0s = absence of allele
					for k in xrange(len(allelst[sourcepop][j])):
											
						# For microsats
						if SNPans == 'N':
							# Assignment of 2, the rest 0
							if rand1 == rand2: 
								if k < rand1 or k > rand1:
									tempindall = 0
								elif k == rand1:
									tempindall = 2
									
							# Assignment of 1s, the rest 0
							if rand1 != rand2:
								if k < min(rand1,rand2) or k > max(rand1,rand2):
									tempindall = 0
								elif k == rand1 or k == rand2:
									tempindall = 1
								else:
									tempindall = 0
						# For SNPs
						else:
							if k == rand1:
								tempindall = 1
							else:
								tempindall = 0
								
						# And to genes list
						offgenes[j].append(tempindall)
			
			# Then record new offspring information to Subpop location [subpop-ofmother,subpop of mother,NASubpop,EmiCD,ImmiCD,age,sex,size,mataure,infection,name,capture,layeggs,genes]
			offpop = offspring[i]['NatalPop']
			name = offspring[i]['name']
			recd = (offpop,offpop,offpop,0.0,-9999,offspring[i]['age'],offspring[i]['sex'],offspring[i]['size'],offspring[i]['mature'],offspring[i]['newmature'],offspring[i]['infection'],name,0,0,0,repr(offgenes))
					
			# Record offspring information to SubpopIN 
			Age0_keep.append(recd)
		
	# Delete offspring variable
	del offspring	
	
	# Return SubpopIN back to array with types
	# Organize type data in SubpopIN - read in from global
	tempN = []
	
	# array and add dtype
	Age0_keep = np.array(Age0_keep,dtype=dtype)
	
	# Return variables from this argument
	return Age0_keep
	
	# End::InheritGenes()

# ---------------------------------------------------------------------------------------------------	
def growInd(Indloc,SubpopIN,sizeLoo,sizeR0,size_1,size_2,size_3,size_4,sizevals,isub,iind,growans,size_mean,gridsample):
	'''
	Growth options
	'''	
	
	# Grow based on von Bertalanffy
	if growans == 'vonB':
		newsize = sizeLoo * (1. - np.exp(-sizeR0*(SubpopIN[isub][iind]['age']+1)))
		if newsize <= 0.:
			print('Warning: von Bertalanffy growth producing negative values.')
			sys.exit(-1)
		SubpopIN[isub][iind]['size'] = newsize	
	# Grow based on temp fit len/size
	elif growans == 'temperature':
		if sizevals[int(Indloc) - 1] != 'N':
			tempval = float(sizevals[int(Indloc) - 1])
			grow = float(size_4[int(Indloc) - 1])
							
			int_R = -sizeR0 * ((scipy.stats.norm(size_1,size_2*size_1).pdf(tempval)) / (scipy.stats.norm(size_1,size_2*size_1).pdf(size_1)))
			
			L_inc = sizeLoo * (1. - np.exp(int_R * (SubpopIN[isub][iind]['age']+1-size_3))) * ((scipy.stats.norm(size_1,size_2*size_1).pdf(tempval)) / (scipy.stats.norm(size_1,size_2*size_1).pdf(size_1)))
			# Get the incremental growth
			L_inc_age = L_inc * np.exp((SubpopIN[isub][iind]['age']+1) * int_R)
			# Update the new size for this individual		
			newsize = SubpopIN[isub][iind]['size'] + (L_inc_age * (grow/365.))
			if newsize <= 0.:
				print('Warning: temperature growth producing negative values.')
				sys.exit(-1)
			SubpopIN[isub][iind]['size'] = newsize				
	# Grow based on bioenergetics
	elif growans == 'bioenergetics':
		print('Bioenergetics equation is not currently implemented.')
		sys.exit(-1)
	# Grow based on known
	elif growans == 'known':
		if gridsample == 'N': # Only apply at 3rd DoUpdate
			if SubpopIN[isub][iind]['age'] > (classno-1): # Age check
				tempage = classno-1 # use last age class
			else:
				tempage = SubpopIN[isub][iind]['age']
			newsize = size_mean[tempage]
			SubpopIN[isub][iind]['size'] = newsize
	# Error check
	else:
		print('Growth options include, vonB, temperature, or bioenergetics. Check that you have entered the correct formate in growth_option in Popvars.csv field.')
		sys.exit(-1)
	
	#End::growInd()
# ---------------------------------------------------------------------------------	
def ageInd(lastage,SubpopIN,isub,iind,sizeans,F_mature,M_mature,Fmat_int,Fmat_slope,Mmat_int,Mmat_slope,eggFreq,Mmat_set,Fmat_set,cdevolveans,fitvals,burningen,gen):
	'''
	Age, mature, and check egg frequency interval here
	'''
	
	# Age here
	# --------
	SubpopIN[isub][iind]['age'] = SubpopIN[isub][iind]['age'] + 1
	
	# Check if the individual is lastage + , then use last index
	if SubpopIN[isub][iind]['age'] >= lastage:
		Indage = lastage-1
	else: # Use the individuals age
		Indage = SubpopIN[isub][iind]['age']
	
	# Check if becomes mature
	# -----------------------
	if SubpopIN[isub][iind]['mature'] == 0:		
		if sizeans == 'N': # Age control
			if SubpopIN[isub][iind]['sex'] == 0: # Female
				if Fmat_set == 'N': # Use prob value
					matval = F_mature[isub][Indage]
				else: # Use set age
					if Indage >= int(Fmat_set): # Age check
						matval = 1.0
					else:
						matval = 0.0				
			else: # Male			
				if Mmat_set == 'N': # Use prob value
					matval = M_mature[isub][Indage]
				else: # Use set age
					if Indage >= int(Mmat_set): # Age check
						matval = 1.0
					else:
						matval = 0.0				
		elif sizeans == 'Y': # Size control
			if (cdevolveans == 'M' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link') and burningen <= gen:
				tempgenes = literal_eval(SubpopIN[isub][iind]['genes'])
				if tempgenes[0][0] == 2: # AA
					tempvals = fitvals[isub][0] # First spot AA
					# Then replace Fmat/Mmat values
					Fmat_int = float(tempvals[3])
					Fmat_slope = float(tempvals[2])
					Mmat_int = float(tempvals[1])
					Mmat_slope = float(tempvals[0])
				elif tempgenes[0][0] == 1 and tempgenes[0][1] == 1: # Aa
					tempvals = fitvals[isub][1] # Second spot Aa
					# Then replace Fmat/Mmat values
					Fmat_int = float(tempvals[3])
					Fmat_slope = float(tempvals[2])
					Mmat_int = float(tempvals[1])
					Mmat_slope = float(tempvals[0])
				elif tempgenes[0][1] == 2: # aa
					tempvals = fitvals[isub][2] # third spot aa
					# Then replace Fmat/Mmat values
					Fmat_int = float(tempvals[3])
					Fmat_slope = float(tempvals[2])
					Mmat_int = float(tempvals[1])
					Mmat_slope = float(tempvals[0])
				else: # Other genotype
					# Then replace Fmat/Mmat values
					Fmat_int = float(Fmat_int)
					Fmat_slope = float(Fmat_slope)
					Mmat_int = float(Mmat_int)
					Mmat_slope = float(Mmat_slope)
				
			if SubpopIN[isub][iind]['sex'] == 0: # Female
				if Fmat_set == 'N': # Use equation - size
					matval = np.exp(Fmat_int + Fmat_slope * SubpopIN[isub][iind]['size']) / (1 + np.exp(Fmat_int + Fmat_slope * SubpopIN[isub][iind]['size']))
				else: # Use set size
					if SubpopIN[isub][iind]['size'] >= int(Fmat_set):
						matval = 1.0
					else:
						matval = 0.0				
			else: # Male			
				if Mmat_set == 'N': # Use equation - size
					matval = np.exp(Mmat_int + Mmat_slope * SubpopIN[isub][iind]['size']) / (1 + np.exp(Mmat_int + Mmat_slope * SubpopIN[isub][iind]['size']))
				else: # Use set size
					if SubpopIN[isub][iind]['size'] >= int(Mmat_set):
						matval = 1.0
					else:
						matval = 0.0						
		else:
			print('Size control option not correct, N or Y.')
			sys.exit(-1)
			
		randmat = rand()
		if randmat < matval:
			SubpopIN[isub][iind]['mature'] = 1# Becomes mature	
			SubpopIN[isub][iind]['newmature'] = 1# Becomes new mature
		else:
			SubpopIN[isub][iind]['mature'] = 0 # Does not mature
			
	# Default check for if age == 6, then make mature
	if SubpopIN[isub][iind]['age'] == 6:
		SubpopIN[isub][iind]['mature'] = 1# Becomes mature	
	
	# Check if mature female, then chance it lays eggs
	if SubpopIN[isub][iind]['mature'] and SubpopIN[isub][iind]['sex'] == 0:
		randegglay = rand()				
		if randegglay < eggFreq:
			SubpopIN[isub][iind]['layeggs'] = 1 # Lays eggs next year
		else:
			SubpopIN[isub][iind]['layeggs'] = 0	# Does not lay eggs next year	
	
	#End::ageInd()
	
# ---------------------------------------------------------------------------------------------------	 
def DoUpdate(SubpopIN,K,xgridpop,ygridpop,gen,nthfile,ithmcrundir,loci,alleles,logfHndl,gridsample,growans,cdevolveans,fitvals = None,burningen = None,ClasscapProb=None,PopcapProb=None,NCap=None,CapClass=None,sizecall=None,size_mean=None,Nclass=None,eggFreq=None,sizevals=None,sizeLoo=None,sizeR0=None,size_1=None,size_2=None,size_3=None,size_4=None,sourcePop=None,sizeans=None,M_mature=None,F_mature=None,Mmat_slope=None,Mmat_int=None,Fmat_slope=None,Fmat_int=None,Mmat_set=None,Fmat_set=None):
	
	'''
	DoUpdate()
	Update Age, Size and some tracker variables.
	Write out information to file.
	'''	
	
	# --------------------------------------------------
	# Tracking numbers for capturing (Middle and Sample and N)
	if gridsample != 'Initial':
		# Get number in class
		classno = len(size_mean[0])
		
		# Patch tracking
		NCap.append([])	
		
		# Add spots for Age tracking
		CapClass.append([]) # If capture 0/1
		CapClass[gen] = [[] for x in xrange(0,classno)]
				
		# Add spots for Age tracking
		Nclass.append([])
		Nclass[gen] = [[] for x in xrange(0,classno)]
		
		# ---------------------------------------------------------------------
		# Begin loop through subpopulations updating tasks at appropriate times
		for isub in xrange(len(K)):
					
			# Get the 'age' adjusted binning - used for capture probability
			# -------------------------------------------------------------
			if sizecall == 'Y':
				size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN[isub]['size'])
			else:
				# Count up each uniages
				age_adjusted = SubpopIN[isub]['age']
			
			# Begin looping through individuals in subpop
			# -------------------------------------------
			for iind in xrange(len(SubpopIN[isub])):
								
				# -----------------------------
				# Grow here - middle and sample
				# -----------------------------
				if growans != 'N':
					Indloc = SubpopIN[isub][iind][sourcePop] # Get location
					# Check for cdevolve growth option - get new growth parameters
					# If MG independent
					if cdevolveans == 'MG_ind':
						Indgenes = literal_eval(SubpopIN[isub][iind]['genes'])
						# BB
						if Indgenes[1][0] == 2:
							genespot = 3
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])							
						# Bb
						elif Indgenes[1][0] == 1 and Indgenes[1][1] == 1:
							genespot = 4
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])
						# bb
						elif Indgenes[1][1] == 2:
							genespot = 5
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])
						else:
							growans = growans
							sizeLoo = float(sizeLoo)
							sizeR0 = float(sizeR0)
							size_1 = float(size_1)
							size_2 = float(size_2)
							size_3 = float(size_3)
					# If MG linked
					if cdevolveans == 'MG_link':
						Indgenes = literal_eval(SubpopIN[isub][iind]['genes'])
						# AA - use BB
						if Indgenes[0][0] == 2:
							genespot = 3
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])							
						# Aa - use Bb
						elif Indgenes[0][0] == 1 and Indgenes[0][1] == 1:
							genespot = 4
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])
						# aa - use bb
						elif Indgenes[0][1] == 2:
							genespot = 5
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])
						else:
							growans = growans
							sizeLoo = float(sizeLoo)
							sizeR0 = float(sizeR0)
							size_1 = float(size_1)
							size_2 = float(size_2)
							size_3 = float(size_3)
					# If just Locus B
					elif cdevolveans == 'G':
						Indgenes = literal_eval(SubpopIN[isub][iind]['genes'])
						# BB
						if Indgenes[1][0] == 2:
							genespot = 0
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])							
						# Bb
						elif Indgenes[1][0] == 1 and Indgenes[1][1] == 1:
							genespot = 1
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])
						# bb
						elif Indgenes[1][1] == 2:
							genespot = 2
							growans = fitvals[int(Indloc)-1][genespot][0]
							sizeLoo = float(fitvals[int(Indloc)-1][genespot][1])
							sizeR0 = float(fitvals[int(Indloc)-1][genespot][2])
							size_1 = float(fitvals[int(Indloc)-1][genespot][3])
							size_2 = float(fitvals[int(Indloc)-1][genespot][4])
							size_3 = float(fitvals[int(Indloc)-1][genespot][5])
						else:
							growans = growans
							sizeLoo = float(sizeLoo)
							sizeR0 = float(sizeR0)
							size_1 = float(size_1)
							size_2 = float(size_2)
							size_3 = float(size_3)
					
					growInd(Indloc,SubpopIN,sizeLoo,sizeR0,size_1,size_2,size_3,size_4,sizevals,isub,iind,growans,size_mean[isub],gridsample)
					
				# --------------------------------------------
				# Age, mature, egg lay frequency here - middle
				# --------------------------------------------
				lastage = classno
				if gridsample == 'Middle':
					
					ageInd(lastage,SubpopIN,isub,iind,sizeans,F_mature,M_mature,Fmat_int,Fmat_slope,Mmat_int,Mmat_slope,eggFreq,Mmat_set,Fmat_set,cdevolveans,fitvals,burningen,gen)
										
				# ---------------------------------
				# Capture here - Middle and Sample
				# ---------------------------------
				# If above last age class
				if age_adjusted[iind] >= lastage:
					Indage = lastage - 1
				else:
					Indage = age_adjusted[iind]
				# Age adjusted capture probability
				capval_age = ClasscapProb[isub][Indage]
				# Patch adjusted capture probability
				capval_pop = PopcapProb[isub]
				
				# Patch level check first
				if capval_pop != 'N':
					capval_pop = float(capval_pop)
					randcapno = rand()
					if randcapno < capval_pop: # Successful patch capture
						# Class level check second
						if capval_age != 'N':					
							capval_age = float(capval_age) # Convert to float
							randcapno = rand()
							#np.random.binomial(Npop,capval)
							if randcapno < capval_age: # Successful capture	
								SubpopIN[isub][iind]['capture'] = 1
								SubpopIN[isub][iind]['recapture'] = SubpopIN[isub][iind]['recapture']+1
			
			# Get the new 'age' adjusted sizes for tracking N
			# -----------------------------------------------
			if sizecall == 'Y':
				size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN[isub]['size'])
			else:
				# Count up each uniages
				age_adjusted = SubpopIN[isub]['age']
			
			# ---------------------------------------------------
			# Track numbers here for capture - Middle and Sample
			# ---------------------------------------------------			
			# Storage tracker for Capture N total
			NCap[gen].append(sum(SubpopIN[isub]['capture']))			
			# Tracking classes: CaptureProb
			for iage in xrange(len(CapClass[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				if len(sizeindex) == 0:
					CapClass[gen][iage].append(0)
				else:
					CapClass[gen][iage].append(sum(SubpopIN[isub][sizeindex]['capture']))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			if len(sizeindex) == 0:
				CapClass[gen][iage].append(0)
			else: # Add them to last class
				CapClass[gen][iage].append(sum(SubpopIN[isub][sizeindex]['capture']))
			
			# ---------------------------------------------------
			# Track numbers here for N - after growth, not Initial
			# ---------------------------------------------------			
			# Tracking classes: CaptureProb
			for iage in xrange(len(Nclass[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				if len(sizeindex) == 0:
					Nclass[gen][iage].append(0)
				else:
					Nclass[gen][iage].append(len(SubpopIN[isub][sizeindex]['size']))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			if len(sizeindex) == 0:
				Nclass[gen][iage].append(0)
			else: # Add them to last class
				Nclass[gen][iage].append(len(SubpopIN[isub][sizeindex]['size']))
		
	# --------------------------------------------------------
	# Track numbers here for N and capture - Middle and Sample
	# --------------------------------------------------------
		# Age tracking
		for iage in xrange(len(CapClass[gen])):		
			CapClass[gen][iage] = sum(CapClass[gen][iage])
			Nclass[gen][iage] = sum(Nclass[gen][iage])
	
	# ------------------------------------------------------------
	# Write out text file for generations specified by nthfile
	# ------------------------------------------------------------	
	# Check if nthfile == generation or "Initial" or "Middle" or "Sample"
	nthfile = np.asarray(nthfile)
	if gridsample == 'Initial':
		DoOutput(SubpopIN,xgridpop,ygridpop,gen,ithmcrundir,loci,alleles,logfHndl,gridsample)
	elif gridsample == 'Middle' or gridsample == 'Sample':
		getyear = np.where(gen == nthfile)[0]
		if len(getyear) != 0:
			DoOutput(SubpopIN,xgridpop,ygridpop,gen,ithmcrundir,loci,alleles,logfHndl,gridsample)
	
	# -----------------------------------------------------
	# Release the captured individuals - Middle and Sample
	# Reset New maturers as well - ALL
	# -----------------------------------------------------
	# Begin loop through subpopulations
	for isub in xrange(len(K)):
		if gridsample != 'Initial':
			SubpopIN[isub]['capture'] = 0
		SubpopIN[isub]['newmature'] = 0

	# Return variables only if updated age
	return SubpopIN
	
	# End::DoUpdate()
	
# ---------------------------------------------------------------------------------------------------
def AddAge0s(SubpopIN_keepAge1plus,K,SubpopIN_Age0,gen,Population,loci,muterate,mtdna,mutationans,dtype,geneswap,allelst,PopulationAge,sizecall,size_mean,cdevolveans,burningen,timecdevolve,fitvals,SelectionDeathsImm_Age0s,SNPans):

	'''
	Add in the Age 0 population.
	'''
	classno = len(size_mean[0])
	# Storage to keep
	SubpopIN_keepK = []
	# Add spot for next generation
	Population.append([]) 
	PopulationAge.append([])
	PopulationAge[gen] = [[] for x in xrange(0,classno)]
	SelectionDeathsImm_Age0s.append([])
	
	for isub in xrange(len(K)):
			
		# Get each SubpopIN pop as array
		SubpopIN_arr = np.array(SubpopIN_keepAge1plus[isub],dtype=dtype)
		
		# Select out the offspring in this pop
		Age0Pop = np.where(SubpopIN_Age0['NatalPop'] == str(isub+1))[0]
		Age0Pop = SubpopIN_Age0[Age0Pop]				
		
		# ----------------
		# InheritGenes()
		# ----------------
		SubpopIN_Age0_temp = InheritGenes(gen,Age0Pop,loci,muterate,mtdna,mutationans,K,dtype,geneswap,allelst,SNPans)	
		
		# --------------------------------
		# Apply spatial selection to Age0s (this might not be the right order)
		# --------------------------------
		if (cdevolveans == '1' or cdevolveans == '1_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
			SubpopIN_Age0_keep = []
			for iind in xrange(len(SubpopIN_Age0_temp)):
				outpool = SubpopIN_Age0_temp[iind]
				# for option 3 in which has to be mature
				if cdevolveans == '1_mat' and outpool['mature'] == 0:
					differentialmortality = 0.0
				else:
					# Call 1-locus selection model
					differentialmortality = Do1LocusSelection(fitvals,literal_eval(outpool['genes'])[0],isub)
				# Then flip the coin to see if outpool survives its location
				randcheck = rand()
				
				# If outpool did not survive: break from loop, move to next outpool
				if randcheck < differentialmortality:					
					continue
				else: # Record if survived
					SubpopIN_Age0_keep.append(outpool)
			# dtype here
			SubpopIN_Age0_keep = np.array(SubpopIN_Age0_keep,dtype=dtype)		
		elif (cdevolveans == '2' or cdevolveans == '2_mat') and gen >= burningen and (timecdevolve == 'Out' or timecdevolve == 'Both'):
			SubpopIN_Age0_keep = []
			for iind in xrange(len(SubpopIN_Age0_temp)):
				outpool = SubpopIN_Age0_temp[iind]
				# for option 3 in which has to be mature
				if cdevolveans == '2_mat' and outpool['mature'] == 0:
					differentialmortality = 0.0
				else:				
					# Call 2-locus selection model
					differentialmortality = Do2LocusSelection(fitvals,literal_eval(outpool['genes'])[0:2],isub)			
				# Then flip the coin to see if outpool survives its location
				randcheck = rand()				
				# If outpool did not survive: break from loop, move to next outpool
				if randcheck < differentialmortality:
					continue
				else: # Record if survived
					SubpopIN_Age0_keep.append(outpool)			
			# dtype here
			SubpopIN_Age0_keep = np.array(SubpopIN_Age0_keep,dtype=dtype)
		else:
			SubpopIN_Age0_keep = SubpopIN_Age0_temp
					
		# Append all information to temp SubpopKeep variable
		SubpopIN_keepK.append(np.concatenate([SubpopIN_arr,SubpopIN_Age0_keep]))				
	
		# Store new N
		Population[gen].append(len(SubpopIN_keepK[isub]))
		SelectionDeathsImm_Age0s[gen].append(len(SubpopIN_Age0_temp)-len(SubpopIN_Age0_keep))
		
		# Age tracking
		# Switch here for size or age control
		if sizecall == 'size':
			size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
			age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])			
		else:
			# Count up each uniages
			age_adjusted = SubpopIN_keepK[isub]['age']
		
		# Tracking age N
		for iage in xrange(len(PopulationAge[gen])):
			sizeindex = np.where(age_adjusted==iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))
			
		# Special case where age class is greater than lastage
		sizeindex = np.where(age_adjusted > iage)[0]
		PopulationAge[gen][iage].append(len(sizeindex))
		
					
	# Add total to N
	Population[gen].insert(0,sum(Population[gen]))
	# Age tracking
	for iage in xrange(len(PopulationAge[gen])):
		PopulationAge[gen][iage] = sum(PopulationAge[gen][iage])	
		
	# add delete here
	return SubpopIN_keepK
	# End::AddAge0s
	