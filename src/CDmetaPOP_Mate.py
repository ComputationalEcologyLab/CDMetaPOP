# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Mate.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for mate processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
	
# Python specific functions
import pdb, random, os, sys, copy
from sets import Set
from ast import literal_eval

# --------------------------------------------------------------------------
def countDuplicatesInList(dupedList):
	'''
	countDuplicatesInList() - Counts dupicates in lists
	'''
	uniqueSet = Set(item for item in dupedList)
	return [dupedList.count(item) for item in uniqueSet]
	
	# End::countDuplicatesInList()

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
#End::count_unique()	
	
# ---------------------------------------------------------------------------------------------------	 
def w_choice_item(lst):
	'''
	w_choice_item()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(lst)
	n=random.uniform(0,wtotal)
	for i in xrange(len(lst)):
		if n < lst[i]:
			break
		n = n-lst[i]
	return i
	
	#End::w_choice_item()

# ---------------------------------------------------------------------------------------------------	
def DoSexual(AAaaMates,AAAAMates,aaaaMates,AAAaMates,aaAaMates,AaAaMates,assortmateC,assortmateModel,xycdmatrix,females,males,matemovethresh,Bearpairs,femalesmated,sourcePop,selfing,subpopmort_mat,count=None):
	'''
	DoSexualYY() and DoSexualNY()
	This function is the mating function for: 
	sexual reproduction
	females	with replacement
	males with replacement.
	Assortative mating checked.
	'''
	
	# For Sexual reproduction NY (count is provided)
	if count != None:
		intfemale = count
	# For Sexual reproduction YY (no count is provided)
	else:
		# Randomly grab a female
		intfemale = int(len(females)*rand())
	
	# Extract the subpopulation this female is in
	femalepop = females[intfemale][sourcePop]

	# Check Assortative Mate model for Hindex or Gene here
	# ----------------------------------------------------
	if len(assortmateModel.split('_')) > 1:
		if assortmateModel.split('_')[1] == 'gene':
			# Get this females genes for assortative mating potential 
			female_genes = females[intfemale]['genes'][0:2]
		elif assortmateModel.split('_')[1] == 'hindex':
			# Get this females genes/hindex for assortive mating potential - round to nearest 10th
			female_hindex = np.around(females[intfemale]['hindex'],1)
		else:
			print('Assortative Mate option entered wrong.')
			sys.exit(-1)
		
	# Extract each male patch probability that female can mate with - careful of indexing
	probarray = xycdmatrix[:,int(femalepop)-1]
	
	# If statement to check if there are patches available in probarray:
	if sum(probarray) != 0.0:
		
		# There are patches in females range, so make sure that a male is in radius
		tempprobarray = copy.deepcopy(probarray) # Copy probarray to change
		
		# continue to check possible patch locations of female to find a male.
		while sum(tempprobarray) != 0.0: 		
			# Select the w_choice item: this is the patch that a male will come from
			itemselect = w_choice_item(tempprobarray)
			
			# Then select all the males in this patch: add one to index to select out subpop
			patchindex = np.where(males[sourcePop]== str(itemselect+1))[0]
			
			# If selfing is not on - only applies to asexual reproduction
			if selfing == 'N':
				# Then remove intfemale from patchindex
				patchindex = patchindex[np.where(patchindex != intfemale)[0]]
		
			# If there are no males in this patch - search other patches
			if len(patchindex) == 0:				
				# Replace probarray with a zero, no males here
				tempprobarray[itemselect] = 0.
				continue
							
			# If subpopulation differential mortality is on
			if not isinstance(subpopmort_mat,str):
								
				# PatchID of female
				female_subpatch = females[intfemale]['popID']
				# PatchID of males patch
				male_subpatch = males[patchindex[0]]['popID']
				
				# if there is a mate selection to another subpatch id
				if female_subpatch != male_subpatch:
					# grab its mortality percentage male moves into female pop (but backwards in matrix (columns are TO)
					differentialmortality = subpopmort_mat[int(female_subpatch)-1][int(male_subpatch)-1]
					
					# check if mating occurs
					continuemate = rand()
					# if randcheck < differentialmortality:
					if continuemate < differentialmortality:
						# Replace probarray with a zero, males can't mate from this patchid
						tempprobarray[itemselect] = 0
						continue
			
			# There are males in this patch, randomly select one while checking for self preference	
			#else:
				
			# Index into the males
			patchmales = males[patchindex]
			
			# Strict self mating option (that is, AA with AA, Aa with Aa, and aa with aa, but using Hindex
			if assortmateModel == '2':
				# Get the males hindex 
				males_hindex = np.around(patchmales['hindex'],1)
				
				# Check for matching males
				males_self = np.where(female_hindex == males_hindex)[0]
				
				# If no males with hindex = females_hindex 
				if len(males_self) == 0:
					# Continue to check next patch					
					tempprobarray[itemselect] = 0. # Replace probarray with a zero, no males here
					continue 
				
				# If males with hindex = females_hindex
				else:
					# Then randomly choose one of the self males
					patchmales = patchmales[males_self]
											
			# Self-preference mating 
			elif assortmateModel == '3a':
				# Get the males Hindex and frequency of each
				males_hindex = np.around(patchmales['hindex'],1)
				males_hindex_count = count_unique(males_hindex)
				males_hindex_fj = males_hindex_count[1]/float(sum(males_hindex_count[1]))
				
				# Calculate probability
				males_hindex_prob = assortmateC ** (female_hindex == males_hindex_count[0]) * males_hindex_fj
				# Normalized
				males_hindex_prob = males_hindex_prob / sum(males_hindex_prob)
				
				# Take a weigthed draw from the 3 genotypes
				selectMaleGenotype = w_choice_item(males_hindex_prob)
					
				# Get selected males of preferential hindex
				patchmales = patchmales[males_hindex == males_hindex_count[0][selectMaleGenotype]]
				
			# Self-preference mating option - multiple species option
			elif assortmateModel == '3b':
				# Get the males Hindex and frequency of each
				males_hindex = np.around(patchmales['hindex'],1)
				males_hindex_count = count_unique(males_hindex)
				males_hindex_fj = males_hindex_count[1]/float(sum(males_hindex_count[1]))
				
				# Calculate probability for females with 1.0 (all males except 0.0)
				if female_hindex == 1.0:
					# Calculate probability
					males_hindex_prob = assortmateC ** (males_hindex_count[0] > 0.0) * males_hindex_fj
				# Calculate probability for females with 0.0 (all males except 1.0)
				elif female_hindex == 0.0:
					# Calculate probability
					males_hindex_prob = assortmateC ** (males_hindex_count[0] < 1.0) * males_hindex_fj
				else: # Hybrid can mate with any
					# Calculate probability
					males_hindex_prob = assortmateC ** (males_hindex_count[0] <= 1.0) * males_hindex_fj	
					
				# Normalized
				males_hindex_prob = males_hindex_prob / sum(males_hindex_prob)
				
				# Take a weigthed draw from the 3 genotypes
				selectMaleGenotype = w_choice_item(males_hindex_prob)
					
				# Get selected males of preferential hindex
				patchmales = patchmales[males_hindex == males_hindex_count[0][selectMaleGenotype]]
				
			# Dominant-preference mating - 
			elif assortmateModel == '4_gene' or assortmateModel == '4_hindex':
				if assortmateModel == '4_hindex':
					# Get the males Hindex and frequency of each
					males_hindex = np.around(patchmales['hindex'],1)
					males_hindex_count = count_unique(males_hindex)
					males_hindex_fj = males_hindex_count[1]/float(sum(males_hindex_count[1]))
				
					# Calculate probability for females greater than 0.0 (all males except 0.0)
					if female_hindex > 0.0:			
						# Calculate probability
						males_hindex_prob = assortmateC ** (males_hindex_count[0] > 0.0) * males_hindex_fj
					# Calculate probability for females with 0.0 (all males except 0.0)
					else: # Hybrid can mate with any
						# Calculate probability
						males_hindex_prob = assortmateC ** (males_hindex_count[0] == 0.0) * males_hindex_fj
					
					# Normalized
					males_hindex_prob = males_hindex_prob / sum(males_hindex_prob)
					
					# Take a weigthed draw from the 3 genotypes
					selectMaleGenotype = w_choice_item(males_hindex_prob)
						
					# Get selected males of preferential hindex
					patchmales = patchmales[males_hindex == males_hindex_count[0][selectMaleGenotype]]
				
				# Special case for sneaker Males - technically not the dominant preference model
				elif assortmateModel == '4_gene':
					# Get the males Genes and frequency of each
					males_genes = patchmales['genes'][:,0:2]
					males_genes_count = count_unique(males_genes[:,0])
					males_genes_fj = males_genes_count[1]/float(sum(males_genes_count[1]))
					# Correct for aa (see Equ in M'Gonigle paper) - kronecker's delta = 0 for sneakers (dgs2020)					
					males_genes_prob = assortmateC ** (males_genes_count[0] > 0) * males_genes_fj
					# Normalized
					males_genes_prob = males_genes_prob / sum(males_genes_prob)
					# Take a weigthed draw from the 3 genotypes
					selectMaleGenotype = w_choice_item(males_genes_prob)
					
					# Get selected males of preferential hindex
					patchmales = patchmales[np.where(males_genes[:,0]==males_genes_count[0][selectMaleGenotype])[0]]

			# Linear hindex preference mating
			elif assortmateModel == '5':
				# Get the males Hindex and frequency of each
				males_hindex = np.around(patchmales['hindex'],1)
				males_hindex_count = count_unique(males_hindex)
				males_hindex_fj = males_hindex_count[1]/float(sum(males_hindex_count[1]))
									
				# Calculate probability
				males_hindex_prob = (1. + (1. - np.abs(males_hindex_count[0] - female_hindex)) * (assortmateC - 1.)) * males_hindex_fj
				# Normalized
				males_hindex_prob = males_hindex_prob / sum(males_hindex_prob)
											
				# Take a weigthed draw from the 3 genotypes
				selectMaleGenotype = w_choice_item(males_hindex_prob)
					
				# Get selected males of preferential hindex
				patchmales = patchmales[males_hindex == males_hindex_count[0][selectMaleGenotype]]
									
			# 'Community' option, select males that are of same 'species' of female
			elif assortmateModel == '6':
				
				# Get female species ID
				female_speciesID = females[intfemale]['species']
				
				# Get the same male species ID
				patchmales = patchmales[patchmales['species'] == female_speciesID]

				# If no males of this species
				if len(patchmales) == 0:
					# Continue to check next patch					
					tempprobarray[itemselect] = 0. # Replace probarray with a zero, no males of this species here
					continue 
								
			# Randomly select a male in patch
			malemate = random.sample(patchmales,1)[0]
			
			# And store the mated pair information.						
			Bearpairs.append([females[intfemale],malemate])
			
			# Tracking
			femalesmated.append(1)
			if len(assortmateModel.split('_')) > 1:
				if assortmateModel.split('_')[1] == 'hindex':
					if (female_hindex == 1.0 and malemate['hindex'] == 0.0) or (female_hindex == 0.0 and malemate['hindex'] == 1.0):
						AAaaMates.append(1)
					elif (female_hindex == 1.0 and malemate['hindex'] == 1.0):
						AAAAMates.append(1)
					elif (female_hindex == 0.0 and malemate['hindex'] == 0.0):
						aaaaMates.append(1)
					elif (female_hindex == 1.0 and (malemate['hindex'] != 0.0 or malemate['hindex'] != 1.0)) or (malemate['hindex'] == 1.0 and (female_hindex != 0.0 or female_hindex != 1.0)):
						AAAaMates.append(1)
					elif (female_hindex == 0.0 and (malemate['hindex'] != 0.0 or malemate['hindex'] != 1.0)) or (malemate['hindex'] == 0.0 and (female_hindex != 0.0 or female_hindex != 1.0)):
						aaAaMates.append(1)
					elif ((female_hindex != 0.0 or female_hindex != 1.0) and (malemate['hindex'] != 0.0 or malemate['hindex'] != 1.0)):
						AaAaMates.append(1)
				else:
					if (female_genes[0] == 2 and malemate['genes'][1] == 2) or (female_genes[1] == 2 and malemate['genes'][0] == 2):
						AAaaMates.append(1)
					elif (female_genes[0] == 2 and malemate['genes'][0] == 2):
						AAAAMates.append(1)
					elif (female_genes[1] == 2 and malemate['genes'][1] == 2):
						aaaaMates.append(1)
					elif (female_genes[0] == 2 and (malemate['genes'][0] == 1 and malemate['genes'][1] == 1)) or (malemate['genes'][0] == 2 and (female_genes[0] == 1 and female_genes[1] == 1)):
						AAAaMates.append(1)
					elif (female_genes[1] == 2 and (malemate['genes'][0] == 1 and malemate['genes'][1] == 1)) or (malemate['genes'][1] == 2 and (female_genes[0] == 1 and female_genes[1] == 1)):
						aaAaMates.append(1)
					elif ((female_genes[0] == 1 and female_genes[1] == 1) and (malemate['genes'][0] == 1 and malemate['genes'][1] == 1)):
						AaAaMates.append(1)
				
			# Then break from patch search loop
			break
					
	# Return Variables from this function
	return Bearpairs,femalesmated
	
	# End::DoSexual()

# ---------------------------------------------------------------------------------------------------	
def DoSexualNN(AAaaMates,AAAAMates,aaaaMates,AAAaMates,aaAaMates,AaAaMates,assortmate,nomales,xycdmatrix,females,\
males,matemovethresh,Bearpairs,femalesmated,subpop,selfing,subpopmort_mat,count=None):
	'''
	DoSexualNN()
	This function is the mating function for
	sexual reproduction
	females	with replacement
	males with replacement
	'''	
	# Assortmate 
	print('Not operating currently')
	sys.exit(-1)
	
	# For Sexual reproduction NY (count is provided)
	if count != None:
		intfemale = count
	# For Sexual reproduction YY (no count is provided)
	else:
		# Randomly grab a female
		intfemale = int(len(females)*rand())
	
	# Extract the subpopulation this female is in: careful of index, subtract 1 for indexing
	femalepop = int(subpop[females[intfemale]]) - 1
	
	# Extract each male patch probability that female can mate in
	probarray = xycdmatrix[:,femalepop]
				
	# If statement to check if there were individuals in probarray:
	if sum(probarray) != 0.0:

		# There are patches in females range, so make sure that a male is in radius
		tempprobarray = copy.deepcopy(probarray) # Copy probarray to change
		
		# continue to check possible patch locations of female to find a male.
		while sum(tempprobarray) != 0.0: 	
		
			# Select the w_choice item: this is the patch that a male will come from
			itemselect = w_choice_item(tempprobarray)
			
			# Then select all the males in this patch: add one to index to select out subpop
			patchindex = np.where(males[sourcePop]== str(itemselect+1))[0]
			
			# If selfing is not on - only applies to asexual reproduction
			if selfing == 'N':
				# Then remove intfemale from patchindex
				patchindex = patchindex[np.where(patchindex != intfemale)[0]]
		
			# Match male index with patchindex
			patchmales = set(males).intersection(patchindex)
		
			# If there are no males in this patch
			if len(patchmales) == 0:				
				# Replace probarray with a zero, no males here
				tempprobarray[itemselect] = 0.
				continue
				
			# If subpopulation differential mortality is on
			if not isinstance(subpopmort_mat,str):
				
				# PatchID of female
				female_subpatch = females[intfemale]['popID']
				# PatchID of males patch
				male_subpatch = males[patchmales[0]]['popID']
				
				# if there is a mate selection to another subpatch id
				if female_subpatch != male_subpatch:
					# grab its mortality percentage male moves into female pop (but backwards in matrix (columns are TO)
					differentialmortality = subpopmort_mat[int(female_subpatch)-1][int(male_subpatch)-1]
					# check if mating occurs
					continuemate = rand()
					# if randcheck < differentialmortality:
					if continuemate < differentialmortality:
						# Replace probarray with a zero, males can't mate from this patchid
						tempprobarray[itemselect] = 0
						continue
				
			#else:			
			# Randomly select a male in patch
			malemate = random.sample(patchmales,1)[0]
		
			# And store the mated pair information.						
			Bearpairs.append([females[intfemale],malemate])
			
			# Then delete that male from the male list
			males = np.delete(males,np.where(males==malemate)[0][0])
			
			# Tracking
			femalesmated.append(1)
			
			# Then break from loop
			break
			
	return Bearpairs,males,femalesmated
	
	# End::DoSexualNN()		

# ---------------------------------------------------------------------------------------------------	 
def DoMate(SubpopIN,K,freplace,mreplace,matemoveno,matemovethresh,xycdmatrix,MateDistCD,xgrid,ygrid,MateDistCDstd,FAvgMate,MAvgMate,FSDMate,MSDMate,Female_BreedEvents,gen,sourcePop,ScaleMax,ScaleMin,A,B,C,Femalepercent,eggFreq,sexans,selfing,assortmateC,AAaaMates,AAAAMates,aaaaMates,AAAaMates,aaAaMates,AaAaMates,assortmateModel,subpopmort_mat,BreedFemales,BreedMales,BreedYYMales,MatureCount,ImmatureCount,ToTFemales,ToTMales,ToTYYMales):

	'''
	DoMate()
	This is the mating function for choosing
	individual mate pairs. 
	Switches for: sexual and asexual mating.	
	'''
	
	# --------------------------------------------------------
	# Preliminary: Needed for both sexual and asexual routines	
	# --------------------------------------------------------
	
	# Get unique number of subpops
	nosubpops = len(K)
	AAaaMates.append([]) # Add for generation split
	AAAAMates.append([])
	aaaaMates.append([])
	AAAaMates.append([])
	aaAaMates.append([])
	AaAaMates.append([])
	BreedMales.append([]) #Storage add spot for generation
	BreedFemales.append([]) #Storage add spot for generation
	BreedYYMales.append([])
	ToTMales.append([]) #Storage add spot for generation
	ToTFemales.append([]) #Storage add spot for generation
	ToTYYMales.append([])
	MatureCount.append([])
	ImmatureCount.append([])
	
	# ---------------------------------------------------
	# Select males and females for mating
	# ---------------------------------------------------

	# Storage variables for breeding age males and females
	#females = []		# These are the zeros or XX
	#males = []			# These are the ones or XY
	
	# Loop through and grab each female and male for probability of mature and probability to lay eggs
	for isub in xrange(len(K)):
		indexF = np.where(SubpopIN[isub]['sex']=='XX')[0]
		indexM = np.where(SubpopIN[isub]['sex']=='XY')[0]
		indexYY = np.where(SubpopIN[isub]['sex']=='YY')[0]
		#indexM = np.concatenate((indexM,indexYY),axis=0) # Here, assume all males the same.
		allfemales = SubpopIN[isub][indexF]
		allmales = SubpopIN[isub][indexM]
		allYYmales = SubpopIN[isub][indexYY]
		# Get reproduction age individuals
		indexFage = np.where(allfemales['layeggs'] == 1)[0]
		if sexans == 'Y':
			indexMage = np.where(allmales['mature'] == 1)[0]
			indexYYage = np.where(allYYmales['mature'] == 1)[0]
		else:
			indexMage = np.where(allmales['layeggs'] == 1)[0]
			indexYYage = np.where(allYYmales['layeggs'] == 1)[0]
		
		# For Tracking
		if sexans == 'Y':
			# Storage tracking
			ToTMales[gen].append(len(indexM)) 
			ToTFemales[gen].append(len(indexF))
			ToTYYMales[gen].append(len(indexYY))
			BreedMales[gen].append(len(indexMage))
			BreedFemales[gen].append(len(indexFage))
			BreedYYMales[gen].append(len(indexYYage))
		else:
			# Storage tracking
			ToTMales[gen].append(len(indexM)+len(indexF)+len(indexYY)) 
			ToTFemales[gen].append(len(indexM)+len(indexF)+len(indexYY))
			ToTYYMales[gen].append(len(indexM)+len(indexF)+len(indexYY))
			BreedMales[gen].append(len(indexMage)+len(indexFage)+len(indexYYage))
			BreedFemales[gen].append(len(indexMage)+len(indexFage)+len(indexYYage))
			BreedYYMales[gen].append(len(indexMage)+len(indexFage)+len(indexYYage))
		MatureCount[gen].append(sum(SubpopIN[isub]['mature']))
		ImmatureCount[gen].append(len(SubpopIN[isub]['mature'])-sum(SubpopIN[isub]['mature']))
		
		# Then get mature females and mature males together for rest of function
		#females.append(list(allfemales[indexFage]))
		# For males, assume YY and XY are the same
		indexM = np.concatenate((indexM,indexYY),axis=0) # Here, assume all males the same.
		allmales = SubpopIN[isub][indexM]		
		if sexans == 'Y':
			indexMage = np.where(allmales['mature'] == 1)[0]
		else:
			indexMage = np.where(allmales['layeggs'] == 1)[0]
		#males.append(list(allmales[indexMage]))
		if isub == 0:	
			females = allfemales[indexFage]
			males = allmales[indexMage]
		else:
			females = np.concatenate((females,allfemales[indexFage]),axis=0)
			males = np.concatenate((males,allmales[indexMage]),axis=0)
	
	# Add Population totals
	ToTMales[gen].insert(0,sum(ToTMales[gen]))
	ToTFemales[gen].insert(0,sum(ToTFemales[gen]))
	ToTYYMales[gen].insert(0,sum(ToTYYMales[gen]))
	BreedMales[gen].insert(0,sum(BreedMales[gen]))
	BreedFemales[gen].insert(0,sum(BreedFemales[gen]))
	BreedYYMales[gen].insert(0,sum(BreedYYMales[gen]))
	
	# Add Count totals
	MatureCount[gen] = sum(MatureCount[gen])
	ImmatureCount[gen] = sum(ImmatureCount[gen])
	
	# Error statement here in case no females or males, then break
	if ToTFemales[gen][0]==0 or (ToTMales[gen][0] + ToTYYMales[gen][0])==0:			
		print('There are no more females or males left in population after year '+str(gen)+'.\n')
		return []
	
	# For sexual reproduction
	if sexans == 'Y':
		# Then get the length of each sex that are reproducing
		nomales = len(males)
		nofemales = len(females)
		# Get the number of times to select mates
		if Femalepercent == 'WrightFisher':
			looptime = nomales + nofemales
		else:
			looptime = nofemales
	# For asexual reproduction
	else:
		allpatches = np.concatenate((females,males),axis=0)
		females = allpatches
		males = allpatches
		del allpatches
		# Then get the length of each sex that are reproducing
		nomales = len(males)
		nofemales = len(females)	
		# Get the number of times to select mates
		looptime = nofemales
	
	# ---------------------------------------------------
	# Choose pairs for mating
	# ---------------------------------------------------
	# Choose mate for each female or individual [female,male]
	Bearpairs = []	
	femalesmated = []
	
	# If there were no reproducing males or females
	if nomales == 0 or nofemales == 0:
		Bearpairs.append([-9999,-9999])
		
	# If there were reproducing males and females
	if nomales != 0 and nofemales != 0:
	
		# For the case of a Female without replacement and a male with replacement
		if freplace == 'N' and mreplace == 'Y':
						
			# Loop through while loop until all females paired up.		
			count = 0		# Initialize the while loop
			while count < looptime:
						
				# Get probability function of user defined input number
				Bearpairs,femalesmated = DoSexual(AAaaMates[gen],AAAAMates[gen],aaaaMates[gen],AAAaMates[gen],aaAaMates[gen],AaAaMates[gen],assortmateC,assortmateModel,xycdmatrix,females,males,matemovethresh,Bearpairs,femalesmated,sourcePop,selfing,subpopmort_mat, count)
												
				# Update count
				count = count + 1
				
		# For the case of a Female with replacement and a male with replacement
		elif freplace == 'Y' and mreplace == 'Y':
				
			# Loop through while loop until all females paired up, but do this total mature individuals times.		
			count = 0		# Initialize the while loop
			while count < looptime:
				
				# Get probability function of user defined input number
				Bearpairs,femalesmated = DoSexual(AAaaMates[gen],AAAAMates[gen],aaaaMates[gen],AAAaMates[gen],aaAaMates[gen],AaAaMates[gen],assortmateC,assortmateModel,xycdmatrix,females,males,matemovethresh,Bearpairs,femalesmated,sourcePop,selfing,subpopmort_mat)
							
				# Update count
				count = count + 1
								
		# For the case of Female with replacement and male without replacement
		elif freplace == 'Y' and mreplace == 'N':
		
			print('Female with replacement and Male without replacement not coded yet.')
			sys.exit(-1)
			
		# For the case of Female without replacement and male without replacement
		elif freplace == 'N' and mreplace == 'N':
			
			# Loop through while loop until all male female pairs occur		
			count = 0		# Initialize the while loop for females
			# Create a temp male to delete from
			tempmales = copy.deepcopy(males)				
			while count < looptime:
							
				# Get probability function of user defined input number
				Bearpairs,tempmales = DoSexualNN(AAaaMates[gen],AAAAMates[gen],aaaaMates[gen],AAAaMates[gen],aaAaMates[gen],AaAaMates[gen],assortmateC,assortmateModel,nomales,xycdmatrix,females,tempmales,matemovethresh,Bearpairs,femalesmated,subpop,selfing,subpopmort_mat,count)
										
				# Update count
				count = count + 1
		
		# Error check
		else:
			print('This Female/Male mating structure does not exist. Must be Y/N combinations.')
			sys.exit(-1)		
	
	# Pairing did not occur
	if len(Bearpairs)==0:
		Bearpairs.append([-9999,-9999])
		
	del females 
	del males
	# ----------------------------------------
	# Summary Stats on Mate functions
	# ----------------------------------------
	# Store the average distance mates were choosen
	tempAvgMateCD = []
	
	# Loop through each CDpair
	for ipair in xrange(len(Bearpairs)):
		
		# else calculate average distances
		if isinstance(Bearpairs[ipair][1],np.void):
		#Bearpairs[ipair][1]!=-9999: # Make sure a mate pair occurred.
			Floc = int(Bearpairs[ipair][0][sourcePop])-1
			Mloc = int(Bearpairs[ipair][1][sourcePop])-1
			probval = xycdmatrix[Floc][Mloc]			
			
			# If panmictic, can't get cost back; also FIDIMO
			if matemoveno == '4' or matemoveno == '6' or matemoveno == '9' or matemoveno == '11':
				cdval = probval
			# IF linear
			elif matemoveno == '1':
				cdval = (probval - 1.) * (-matemovethresh)
				tempAvgMateCD.append(cdval)
			# If inverse square
			elif matemoveno == '2':
				if probval == 1.0:
					cdval = 0.0
				else:	
					cdval = np.sqrt(1. / (probval * (ScaleMax - ScaleMin) + ScaleMin))
			# If neg exp
			elif matemoveno == '5':
				cdval = np.log((probval * (ScaleMax-ScaleMin) + ScaleMin)/float(A)) / (-float(B) * np.log(10))
			
			# If gaussian
			elif matemoveno == '7':				
				cdval = float(B) + np.sqrt(-2*float(C)**2 * np.log((probval*(ScaleMax-ScaleMin)+ScaleMin)/float(A)))
			# If just rescaled
			elif matemoveno == '8':
				cdval = probval*(ScaleMax-ScaleMin)+ScaleMin
			# If pareto
			elif matemoveno == '10':
				if probval == 1.0:
					cdval = 0.0
				else:
					cdval = pow(((float(A)*float(B)**float(A))/probval),(1/(float(A)+1))) - float(B)
			else:
				print('Mate move function does not exist.')
				sys.exit(-1)
			tempAvgMateCD.append(cdval)
			
	# If at least some individuals mated
	if isinstance(Bearpairs[ipair][1],np.void):
	#if Bearpairs[0][0] != -9999:
			
		# And append to MateDistCD
		MateDistCD.append(sum(tempAvgMateCD) / len(Bearpairs))
		MateDistCDstd.append(np.std(tempAvgMateCD))
		
	# No Mates
	else:

		# And append to MateDistCD
		MateDistCD.append(0)
		MateDistCDstd.append(0) 
	
	# Track actual number of breeding events of females and types.
	Female_BreedEvents.append(sum(femalesmated))
	AAaaMates[gen] = sum(AAaaMates[gen])
	AAAAMates[gen] = sum(AAAAMates[gen])
	aaaaMates[gen] = sum(aaaaMates[gen])
	AAAaMates[gen] = sum(AAAaMates[gen])
	aaAaMates[gen] = sum(aaAaMates[gen])
	AaAaMates[gen] = sum(AaAaMates[gen])
	del tempAvgMateCD
	
	# Convert Bearpairs to array with dtype
	Bearpairs = np.asarray(Bearpairs)
	
	# Return variables from this function
	return Bearpairs
	
	#End::DoMate()