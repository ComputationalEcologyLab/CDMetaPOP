# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Mate.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for mate processes.
# --------------------------------------------------------------------------------------------------
	
# Python specific functions
import pdb, os, sys, copy,numbers
from ast import literal_eval
import numpy as np 
from CDmetaPOP_Offspring2 import *
from CDmetaPOP_Mortality import *

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
#End::count_unique()	
	
# ---------------------------------------------------------------------------------------------------
def DoSummaryMate(Bearpairs,sourcePop,xycdmatrix,matemoveno,matemovethresh,ScaleMax, ScaleMin, A, B, C,MateDistCD,MateDistCDstd,Female_BreedEvents, AAaaMates,AAAAMates,aaaaMates,AAAaMates,aaAaMates,AaAaMates,femalesmated,gen,outputans):
	if outputans == 'Y':
		# Store the average distance mates were choosen
		tempAvgMateCD = []	
		# Loop through each CDpair
		for ipair in range(len(Bearpairs)):
			if isinstance(Bearpairs[ipair][1],np.void):
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
			# And append to MateDistCD
			MateDistCD.append(sum(tempAvgMateCD) / len(Bearpairs))
			MateDistCDstd.append(np.std(tempAvgMateCD))		
		# No Mates
		else:
			# And append to MateDistCD
			MateDistCD.append(0)
			MateDistCDstd.append(0)
		del tempAvgMateCD
	else: # Output skipped 
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
	
	#END::DoSummaryMate()
	
# ---------------------------------------------------------------------------------------------------	 
def w_choice_item(lst):
	'''
	w_choice_item()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(lst)
	n=np.random.uniform(0,wtotal)
	for i in range(len(lst)):
		if n < lst[i]:
			break
		n = n-lst[i]
	return i
	
	#End::w_choice_item()

# ---------------------------------------------------------------------------------------------------	
def DoSexual(AAaaMates,AAAAMates,aaaaMates,AAAaMates,aaAaMates,AaAaMates,assortmateC,assortmateModel,xycdmatrix,females,males,Bearpairs,femalesmated,sourcePop,selfing,subpopmort_mat,natal_patches,K,intfemale):
	'''
	DoSexualYY() and DoSexualNY()
	This function is the mating function for: 
	sexual reproduction
	females	with replacement
	males with replacement.
	Assortative mating checked.
	'''
	# Extract the subpopulation this female is in
	femalepop = females[intfemale][sourcePop]
	
	# Extract each male patch probability that female can mate with - careful of indexing - prevents grabbing a male mate in non natal patch
	probarray = xycdmatrix[:,int(femalepop)-1]
	# Where natal grounds = 0, set prob to 0
	probarray[np.where(np.asarray(natal_patches)==0)[0]] = 0.
	# Where K = 0, set prob to 0
	probarray[np.where(np.asarray(K)==0)[0]] = 0.
	
	# And to track if a mate event occurred
	mate_event = 0
	
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
			'''
			# Selfing Options Here
			# --------------------
			if selfing == 'N': # Selfing off
				# Then remove intfemale from patchindex
				patchindex = patchindex[np.where(patchindex != intfemale)[0]]
			elif selfing == 'Y': # Selfing on, but for sexual / asexual reproduction options, only this individual kept.
				# Keep this female in the list for selfing probability
				patchindex = patchindex
			else: # Assume float here a probability value entered and hermaphroditic mating system. 
				selfing = float(selfing)
				# Then check for selfing mate first
				checkselfing = np.random.uniform()
				if checkselfing < selfing:
					# Then selfing occurs - select out this pairing.										
					Bearpairs.append([females[intfemale],females[intfemale]])
					# Tracking
					femalesmated.append(1)
					# Then break from patch search loop
					break
				# else if selfing did not occur continue on for male selection. 
			'''
			# If there are no males in this patch - search other patches
			if len(patchindex) == 0:
				# Replace probarray with a zero, no males here
				tempprobarray[itemselect] = 0.
				continue
			
			# There are males in this patch, randomly select one while checking for self preference	
			# Index into the males
			patchmales = males[patchindex]	
			if assortmateModel == '1': #Random mating
				
				# Randomly select a male in patch
				malemate = np.random.choice(patchmales,1)[0]
				mate_event = 1				
				
			# Strict self mating option (that is, AA with AA, Aa with Aa, and aa with aa, but using Hindex
			if assortmateModel == '2':
				# Get the males hindex 
				female_hindex = np.around(females[intfemale]['hindex'],1)
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
				# Randomly select a male in patch
				malemate = np.random.choice(patchmales,1)[0]
				mate_event = 1
			
			# Self-preference mating 
			elif assortmateModel == '3a':
				# Get the males Hindex and frequency of each
				female_hindex = np.around(females[intfemale]['hindex'],1)
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
				# Randomly select a male in patch
				malemate = np.random.choice(patchmales,1)[0]
				mate_event = 1
				
			# Self-preference mating option - multiple species option
			elif assortmateModel == '3b':
				# Get the males Hindex and frequency of each
				female_hindex = np.around(females[intfemale]['hindex'],1)
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
				# Randomly select a male in patch
				malemate = np.random.choice(patchmales,1)[0]
				mate_event = 1
				
			# Dominant-preference mating - 
			elif assortmateModel == '4_gene' or assortmateModel == '4':
				if assortmateModel == '4': #hindex option
					# Get the males Hindex and frequency of each
					female_hindex = np.around(females[intfemale]['hindex'],1)
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
					# Randomly select a male in patch
					malemate = np.random.choice(patchmales,1)[0]
					mate_event = 1
				
				# Special case for sneaker Males - technically not the dominant preference model
				elif assortmateModel == '4_gene':
					# Get female genes
					#pdb.set_trace()
					female_genes = females[intfemale]['genes'][0:2]
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
					# Randomly select a male in patch
					malemate = np.random.choice(patchmales,1)[0]
					mate_event = 1
			
			# Linear hindex preference mating
			elif assortmateModel == '5':
				female_hindex = np.around(females[intfemale]['hindex'],1)											 
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
				# Randomly select a male in patch
				malemate = np.random.choice(patchmales,1)[0]
				mate_event = 1
			
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
				malemate = np.random.choice(patchmales,1)[0]
				mate_event = 1
			
			# Then break from patch search loop
			break
	
	# Then Store information
	if mate_event == 0:
		Bearpairs.append([-9999,-9999])
		femalesmated.append(0)
		AAaaMates.append(0)
		AAAAMates.append(0)
		aaaaMates.append(0)
		AAAaMates.append(0)
		aaAaMates.append(0)
		AaAaMates.append(0)
	else:
		# And store the mated pair information.						
		Bearpairs.append([females[intfemale],malemate])
		# Get female genes
		female_genes = females[intfemale]['genes'][0:2]
			
		# Tracking
		femalesmated.append(1)
		# Tracking mating to verify if truly random
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
			
	# Return Variables from this function
	return Bearpairs,femalesmated
	
	# End::DoSexual()

# ---------------------------------------------------------------------------------------------------	
def DoSexualNN(AAaaMates,AAAAMates,aaaaMates,AAAaMates,aaAaMates,AaAaMates,assortmate,nomales,xycdmatrix,females,\
males,Bearpairs,femalesmated,subpop,selfing,subpopmort_mat,natal_patches,K,count=None):
	'''
	DoSexualNN()
	This function is the mating function for
	sexual reproduction
	females	without replacement
	males without replacement
	'''	
	# Assortmate 
	print('Sexual NN Not operating currently, add assortmate options and updates v2.69.')
	sys.exit(-1)
	
	# For Sexual reproduction NY (count is provided)
	if count != None:
		intfemale = count
	# For Sexual reproduction YY (no count is provided)
	else:
		# Randomly grab a female
		intfemale = int(len(females)*np.random.uniform())
	
	# Extract the subpopulation this female is in: careful of index, subtract 1 for indexing
	femalepop = int(subpop[females[intfemale]]) - 1
	
	# Extract each male patch probability that female can mate in
	probarray = xycdmatrix[:,femalepop]
	# Where natal grounds = 0, set prob to 0
	probarray[np.where(np.asarray(natal_patches)==0)[0]] = 0.
	# Where K = 0, set prob to 0
	probarray[np.where(np.asarray(K)==0)[0]] = 0.
				
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
			'''
			# Selfing Options Here
			# --------------------
			if selfing == 'N': # Selfing off
				# Then remove intfemale from patchindex
				patchindex = patchindex[np.where(patchindex != intfemale)[0]]
			elif selfing == 'Y': # Selfing on, but for sexual / asexual reproduction options, only this individual kept.
				# Keep this female in the list for selfing probability
				patchindex = patchindex
			else: # Assume float here a probability value entered and hermaphroditic mating system. 
				selfing = float(selfing)
				# Then check for selfing mate first
				checkselfing = np.random.uniform()
				if checkselfing < selfing:
					# Then selfing occurs - select out this pairing.										
					Bearpairs.append([females[intfemale],females[intfemale]])
					# Tracking
					femalesmated.append(1)
					# Then break from patch search loop
					break
				# else if selfing did not occur continue on for male selection. 
			'''
			# Match male index with patchindex
			patchmales = set(males).intersection(patchindex)
		
			# If there are no males in this patch
			if len(patchmales) == 0:				
				# Replace probarray with a zero, no males here
				tempprobarray[itemselect] = 0.
				continue
							
			# Randomly select a male in patch
			malemate = np.random.choice(patchmales,1)[0]
		
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
def DoMate(SubpopIN,K,freplace,mreplace,matemoveno,matemovethresh,xycdmatrix,MateDistCD,xgrid,ygrid,MateDistCDstd,FAvgMate,MAvgMate,FSDMate,MSDMate,Female_BreedEvents,gen,sourcePop,ScaleMax,ScaleMin,A,B,C,Femalepercent,sexans,selfing,assortmateC,AAaaMates,AAAAMates,aaaaMates,AAAaMates,aaAaMates,AaAaMates,assortmateModel,subpopmort_mat,BreedFemales,BreedMales,BreedYYMales,BreedYYFemales,MatureCount,ImmatureCount,ToTFemales,ToTMales,ToTYYMales,ToTYYFemales,egg_delay,Bearpairs_temp,natal_patches,offno,transmissionprob,f_ind,age_sigma,sizeans,egg_mean_1,egg_mean_2,egg_mean_ans,equalClutch,dtype,eggmort_patch,Track_EggDeaths,eggmort_pop,noOffspring_temp,Track_Births,Track_BirthsMYY,Track_BirthsFYY,constMortans,outputans):

	'''
	DoMate()
	This is the mating function for choosing
	individual mate pairs. 
	Switches for: sexual and asexual mating.	
	'''
	# Check Assortative Mate model for Hindex or Gene here
	# ----------------------------------------------------
	if assortmateModel not in ['1','2','3a','3b','4_gene','4','5','6']:
		print('Assortative Mate option entered wrong.')
		sys.exit(-1)
	
	# Selfing needs more checking
	if selfing == 'Y':
		print('Selfing disabled currently.')
		sys.exit(-1)
	
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
	BreedYYFemales.append([])
	ToTMales.append([]) #Storage add spot for generation
	ToTFemales.append([]) #Storage add spot for generation
	ToTYYMales.append([])
	ToTYYFemales.append([])
	MatureCount.append([])
	ImmatureCount.append([])
	# Storage for egg deaths
	Track_EggDeaths.append([])
	Track_Births.append([]) # Spot for generation, 
	Track_BirthsMYY.append([]) # Spot for generation, note this is the offspring number from a MYY after egg deaths.
	Track_BirthsFYY.append([]) # Spot for generation, note this is the offspring number after egg deaths. 	
	
	
	# ---------------------------------------------------
	# Select males and females for mating
	# ---------------------------------------------------
	# Loop through and grab each female and male for probability of mature and probability to lay eggs
	for isub in range(len(K)):
		indexF = np.where(SubpopIN[isub]['sex']=='FXX')[0]
		indexM = np.where(SubpopIN[isub]['sex']=='MXY')[0]
		indexMYY = np.where(SubpopIN[isub]['sex']=='MYY')[0]
		indexFYY = np.where(SubpopIN[isub]['sex']=='FYY')[0]
		allfemales = SubpopIN[isub][indexF]
		allmales = SubpopIN[isub][indexM]
		allYYmales = SubpopIN[isub][indexMYY]
		allYYfemales = SubpopIN[isub][indexFYY]
		# Index for mature - breeding individuals
		if sexans == 'Y' or sexans == 'H':
			indexFage = np.where(allfemales['mature'] == 1)[0]
			indexMage = np.where(allmales['mature'] == 1)[0]
			indexMYYage = np.where(allYYmales['mature'] == 1)[0]
			indexFYYage = np.where(allYYfemales['mature'] == 1)[0]
		else:
			indexFage = np.where(allfemales['layeggs'] >= 1)[0]
			indexMage = np.where(allmales['layeggs'] >= 1)[0]
			indexMYYage = np.where(allYYmales['layeggs'] >= 1)[0]
			indexFYYage = np.where(allYYfemales['layeggs'] >= 1)[0]
		
		# Store For Tracking
		if sexans == 'Y' or sexans == 'H':
			# Storage tracking
			ToTMales[gen].append(len(indexM)) 
			ToTFemales[gen].append(len(indexF))
			ToTYYMales[gen].append(len(indexMYY))
			ToTYYFemales[gen].append(len(indexFYY))
			BreedMales[gen].append(len(indexMage))
			BreedFemales[gen].append(len(indexFage))
			BreedYYMales[gen].append(len(indexMYYage))
			BreedYYFemales[gen].append(len(indexFYYage))
		else:
			# Storage tracking
			ToTMales[gen].append(len(indexM)+len(indexF)+len(indexMYY)+len(indexFYY)) 
			ToTFemales[gen].append(len(indexM)+len(indexF)+len(indexMYY)+len(indexFYY))
			ToTYYMales[gen].append(len(indexM)+len(indexF)+len(indexMYY)+len(indexFYY))
			ToTYYFemales[gen].append(len(indexM)+len(indexF)+len(indexMYY)+len(indexFYY))
			BreedMales[gen].append(len(indexMage)+len(indexFage)+len(indexMYYage)+len(indexFYYage))
			BreedFemales[gen].append(len(indexMage)+len(indexFage)+len(indexMYYage)+len(indexFYYage))
			BreedYYMales[gen].append(len(indexMage)+len(indexFage)+len(indexMYYage)+len(indexFYYage))
			BreedYYFemales[gen].append(len(indexMage)+len(indexFage)+len(indexMYYage)+len(indexFYYage))
		MatureCount[gen].append(sum(SubpopIN[isub]['mature']))
		ImmatureCount[gen].append(len(SubpopIN[isub]['mature'])-sum(SubpopIN[isub]['mature']))
		
		# Then get mature females and mature males together for rest of function
		# For males, assume YY and XY are the same
		indexM = np.concatenate((indexM,indexMYY),axis=0) # Here, assume all males the same.
		allmales = SubpopIN[isub][indexM]		
		if sexans == 'Y' or sexans == 'H':
			indexMage = np.where(allmales['mature'] == 1)[0] 
		else:
			indexMage = np.where(allmales['layeggs'] >= 1)[0]
		
		# For females, assume YY and XX are the same
		indexF = np.concatenate((indexF,indexFYY),axis=0)
		allfemales = SubpopIN[isub][indexF]
		# Overwirte indexFage here with 'layeggs'
		indexFage = np.where(allfemales['layeggs'] >= 1)[0] # Use layeggs for choosing breeding Bearpairs
		# Then repeat this index for when layeggs > 1
		indexFage_rep = np.asarray(allfemales['layeggs'][indexFage],dtype=int)
		indexFage = np.array([val for val, rep in zip(indexFage, indexFage_rep) for i in range(rep)],dtype=int)	
		
		# For all patches, append females/males
		if isub == 0:	
			females = allfemales[indexFage]
			males = allmales[indexMage]
		else:
			females = np.concatenate((females,allfemales[indexFage]),axis=0)
			males = np.concatenate((males,allmales[indexMage]),axis=0)
			
		# Storage Births/Egg Deaths for each patch location
		Track_EggDeaths[gen].append([])
		Track_Births[gen].append([]) 
		Track_BirthsMYY[gen].append([]) 
		Track_BirthsFYY[gen].append([]) 		
		
	# Add Population totals
	ToTMales[gen].insert(0,sum(ToTMales[gen]))
	ToTFemales[gen].insert(0,sum(ToTFemales[gen]))
	ToTYYMales[gen].insert(0,sum(ToTYYMales[gen]))
	ToTYYFemales[gen].insert(0,sum(ToTYYFemales[gen]))
	BreedMales[gen].insert(0,sum(BreedMales[gen]))
	BreedFemales[gen].insert(0,sum(BreedFemales[gen]))
	BreedYYMales[gen].insert(0,sum(BreedYYMales[gen]))
	BreedYYFemales[gen].insert(0,sum(BreedYYFemales[gen]))
	
	# Add Count totals
	MatureCount[gen] = sum(MatureCount[gen])
	ImmatureCount[gen] = sum(ImmatureCount[gen])
	#pdb.set_trace()
	# Error statement here in case no females or males, then break from loop, update any tracking left in this function
	if (ToTFemales[gen][0] + ToTYYFemales[gen][0])==0 or (ToTMales[gen][0] + ToTYYMales[gen][0])==0:			
		MateDistCD.append(0)
		MateDistCDstd.append(0)
		Female_BreedEvents.append(0)
		AAaaMates[gen] = AAAAMates[gen] = aaaaMates[gen] = AAAaMates[gen] = aaAaMates[gen] = AaAaMates[gen] = 0
		Track_Births[gen] = [0 for x in range(0,len(K)+1)] 
		Track_BirthsMYY[gen] = [0 for x in range(0,len(K)+1)] 
		Track_BirthsFYY[gen] = [0 for x in range(0,len(K)+1)] 
		Track_EggDeaths[gen] = [0 for x in range(0,len(K)+1)] 
		#pdb.set_trace()
		Bearpairs_temp[egg_delay] = []	
		noOffspring_temp[egg_delay] = []
		return Bearpairs_temp,noOffspring_temp
	
	# For sexual reproduction
	if sexans == 'Y' or sexans == 'H':
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
	# Choose pairs for mating and assign offpsring number
	# ---------------------------------------------------
	# Choose mate for each female or individual [female,male]
	Bearpairs = []	
	femalesmated = []
	
	# Tracking for offspring
	noOffspring = []
	# Get size or age control here
	if sizeans == 'Y':
		sizecall = 'size'
	elif sizeans == 'N':
		sizecall = 'age'
	else:
		print('Specify Y or N for size control parameters.')
		sys.exit(-1)
	
	# If there were no reproducing males or females
	if nomales == 0 or nofemales == 0:
		Bearpairs.append([-9999,-9999])
		noOffspring = [0]
		Track_Births[gen] = [0 for x in range(0,len(K)+1)] 
		Track_BirthsMYY[gen] = [0 for x in range(0,len(K)+1)] 
		Track_BirthsFYY[gen] = [0 for x in range(0,len(K)+1)] 
		Track_EggDeaths[gen] = [0 for x in range(0,len(K)+1)] 
		femalesmated = [0]
	
	# If there were reproducing males and females
	if nomales != 0 and nofemales != 0:
	
		# For the case of a Female without replacement and a male with replacement
		if (freplace == 'N' and mreplace == 'Y'):
						
			# Loop through while loop until all females paired up.		
			count = 0		# Initialize the while loop
			while count < looptime:
				
				# Extract the subpopulation this female is in: if it is not in natal pop (mating ground), then skip this females
				femalepop = int(females[count][sourcePop] ) - 1
				# Checks this females location and will mate only if female in natal/mating ground patch, else skip this female
				if femalepop in np.where(np.asarray(natal_patches)==1)[0]:					

					# Get probability function of user defined input number
					Bearpairs,femalesmated = DoSexual(AAaaMates[gen],AAAAMates[gen],aaaaMates[gen],AAAaMates[gen],aaAaMates[gen],AaAaMates[gen],assortmateC,assortmateModel,xycdmatrix,females,males,Bearpairs,femalesmated,sourcePop,selfing,subpopmort_mat,natal_patches,K,count)
				
					# Pairing occured	
					if not isinstance(Bearpairs[count][0],int):
									
						# --------------------------------------
						# Call DoOffspringNo() for this Bearpair
						# --------------------------------------
						thisBearpair_noOffspring = DoOffspringNo(offno,Bearpairs[count],f_ind,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans)
						
						# -------------------------------------
						# Call DoEggMortality()
						# -------------------------------------	
						thisBearpair_noOffspring = DoIndividualEggMortality(Bearpairs[count],eggmort_patch,Track_EggDeaths,gen,eggmort_pop,Track_Births,Track_BirthsMYY,Track_BirthsFYY,thisBearpair_noOffspring,constMortans)
						
						# Append to noOffspring list
						noOffspring.append(thisBearpair_noOffspring)
					
					# Pairing did not occur
					else:
						noOffspring.append(0)
				# If this female was skipped
				else:
					Bearpairs.append([-9999,-9999])
					femalesmated.append(0)
					AAaaMates[gen].append(0)
					AAAAMates[gen].append(0)
					aaaaMates[gen].append(0)
					AAAaMates[gen].append(0)
					aaAaMates[gen].append(0)
					AaAaMates[gen].append(0)
					noOffspring.append(0)
												
				# Update count
				count = count + 1
				
		# For the case of a Female with replacement and a male with replacement
		elif freplace == 'Y' and mreplace == 'Y':	
			# Loop through while loop until all females paired up, but do this total mature individuals times.		
			count = 0		# Initialize the while loop
			while count < looptime:
				
				intfemale = int(len(females)*np.random.uniform())
				# Extract the subpopulation this female is in if in natal pop == 0 then skip this females
				femalepop = int(females[intfemale][sourcePop] ) - 1
				# Checks this females location and will mate only if female in natal ground patch, else skip this female
				if femalepop in np.where(np.asarray(natal_patches)==1)[0]:					

					# Get probability function of user defined input number
					Bearpairs,femalesmated = DoSexual(AAaaMates[gen],AAAAMates[gen],aaaaMates[gen],AAAaMates[gen],aaAaMates[gen],AaAaMates[gen],assortmateC,assortmateModel,xycdmatrix,females,males,Bearpairs,femalesmated,sourcePop,selfing,subpopmort_mat,natal_patches,K,intfemale)
										
					# Pairing occured
					if not isinstance(Bearpairs[count][0],int):
					#if Bearpairs[count][0] != -9999:
									
						# --------------------------------------
						# Call DoOffspringNo() for this Bearpair
						# --------------------------------------
						thisBearpair_noOffspring = DoOffspringNo(offno,Bearpairs[count],f_ind,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans)
					
						# -------------------------------------
						# Call DoEggMortality()
						# -------------------------------------	
						thisBearpair_noOffspring = DoIndividualEggMortality(Bearpairs[count],eggmort_patch,Track_EggDeaths,gen,eggmort_pop,Track_Births,Track_BirthsMYY,Track_BirthsFYY,thisBearpair_noOffspring,constMortans)
						
						# Append to noOffspring list
						noOffspring.append(thisBearpair_noOffspring)
					
					# Pairing did not occur
					else:
						noOffspring.append(0)
				# If this female was skipped
				else:
					Bearpairs.append([-9999,-9999])
					femalesmated.append(0)
					AAaaMates[gen].append(0)
					AAAAMates[gen].append(0)
					aaaaMates[gen].append(0)
					AAAaMates[gen].append(0)
					aaAaMates[gen].append(0)
					AaAaMates[gen].append(0)
					noOffspring.append(0)
							
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
							
				# Extract the subpopulation this female is in if in natal pop == 0 then skip this females
				femalepop = int(females[count][sourcePop] ) - 1
				# Checks this females location and will mate only if female in natal ground patch, else skip this female
				if femalepop in np.where(np.asarray(natal_patches)==1)[0]:					

					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualNN(AAaaMates[gen],AAAAMates[gen],aaaaMates[gen],AAAaMates[gen],aaAaMates[gen],AaAaMates[gen],assortmateC,assortmateModel,nomales,xycdmatrix,females,tempmales,Bearpairs,femalesmated,subpop,selfing,subpopmort_mat,natal_patches,K,count)
					
					# Pairing occured	
					if not isinstance(Bearpairs[count][0],int):
									
						# --------------------------------------
						# Call DoOffspringNo() for this Bearpair
						# --------------------------------------
						thisBearpair_noOffspring = DoOffspringNo(offno,Bearpairs[count],f_ind,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans)
					
						# -------------------------------------
						# Call DoEggMortality()
						# -------------------------------------	
						thisBearpair_noOffspring = DoIndividualEggMortality(Bearpairs[count],eggmort_patch,Track_EggDeaths,gen,eggmort_pop,Track_Births,Track_BirthsMYY,Track_BirthsFYY,thisBearpair_noOffspring,constMortans)
						
						# Append to noOffspring list
						noOffspring.append(thisBearpair_noOffspring)
						
					# Pairing did not occur
					else:
						noOffspring.append(0)
										
				# If this female was skipped
				else:
					Bearpairs.append([-9999,-9999])
					femalesmated.append(0)
					AAaaMates[gen].append(0)
					AAAAMates[gen].append(0)
					aaaaMates[gen].append(0)
					AAAaMates[gen].append(0)
					aaAaMates[gen].append(0)
					AaAaMates[gen].append(0)
					noOffspring.append(0)
				# Update count
				count = count + 1
		
		# Error check
		else:
			print('This Female/Male mating structure does not exist. Must be Y/N combinations.')
			sys.exit(-1)		
		
	# For testing EBT old vs new model - adding DoEggMortality() check
	#thisBearpair_noOffspring = DoIndividualEggMortality(Bearpairs[count],eggmort_patch,Track_EggDeaths,gen,eggmort_pop,Track_Births,Track_BirthsMYY,Track_BirthsFYY,thisBearpair_noOffspring,constMortans)
	#DoEggMortality(Bearpairs,eggmort_patch,Age0Deaths,gen,K,eggmort_back,noOffspring,Births,BirthsMYY,BirthsFYY)
	#testEggDeaths = []
	#testBirths = []
	#testBirthsMYY = []
	#testBirthsFYY = []
	#testnoOffspring, testEggDeaths, testBirths, testBirthsMYY, testBirthsFYY = DoEggMortality(Bearpairs, eggmort_patch, testEggDeaths, 0, K, eggmort_pop, noOffspring, testBirths, testBirthsMYY, testBirthsFYY) 
	
	# --------------------------------
	# Clean up after loop - some checks
	# --------------------------------
	if femalesmated == [] or len(Bearpairs) == 0 or len(Bearpairs) != len(noOffspring):
		print('Check Mate pairing and offspring number')
		sys.exit(-1)
	
	# To arrays
	noOffspring = np.asarray(noOffspring)
	Bearpairs = np.asarray(Bearpairs)
	del females 
	del males
	
	# -----------------------------------
	# Tracking Updates Egg Deaths, Births
	# -----------------------------------
	if not isinstance(Track_Births[gen][0],int):
		Track_Births[gen] = [sum(sublist) for sublist in Track_Births[gen]]
		Track_Births[gen].insert(0,sum(Track_Births[gen]))
		Track_EggDeaths[gen] = [sum(sublist) for sublist in Track_EggDeaths[gen]]
		Track_EggDeaths[gen].insert(0,sum(Track_EggDeaths[gen]))
		Track_BirthsMYY[gen] = [sum(sublist) for sublist in Track_BirthsMYY[gen]]
		Track_BirthsMYY[gen].insert(0,sum(Track_BirthsMYY[gen]))
		Track_BirthsFYY[gen] = [sum(sublist) for sublist in Track_BirthsFYY[gen]]
		Track_BirthsFYY[gen].insert(0,sum(Track_BirthsFYY[gen]))
	
	# --------------------------------------------
	# Check for Equal clutch size and reduce pairs
	# --------------------------------------------
	# If equal clutch size is turned on
	if equalClutch == 'Y' and freplace == 'Y':
		noOffspring = DoOffspringClutch(Bearpairs,dtype,noOffspring)
		# A quick update to Tracking: Note only updating Births tracking; Egg deaths not accurate 
		Track_Births[gen][0] = sum(noOffspring)
	
	# --------------------------------------------
	# Check for 0 events and reduce pairs
	# --------------------------------------------	
	if len(np.where(noOffspring == 0)[0]) > 0:
		# Get index of 0 births for mothers
		ind0 = np.where(noOffspring == 0)[0]
		# Delete Bearpairs and offspring no index
		Bearpairs = np.delete(Bearpairs,ind0,0)	
		noOffspring = np.delete(noOffspring,ind0)
		if len(Bearpairs) == 0:
			Bearpairs = np.asarray([[-9999,-9999]])
		else: # Dtype potentially lost through -9999 no mated pairs 
			Bearpairs = Bearpairs.astype(dtype)
	
	# ----------------------------------------
	# Summary Stats on Mate functions
	# ----------------------------------------
	DoSummaryMate(Bearpairs,sourcePop,xycdmatrix,matemoveno,matemovethresh,ScaleMax, ScaleMin, A, B, C,MateDistCD,MateDistCDstd,Female_BreedEvents, AAaaMates,AAAAMates,aaaaMates,AAAaMates,aaAaMates,AaAaMates,femalesmated,gen,outputans)
	
	# ---------------------------------------------------------------
	# Update for egg_delay; create n-D noOffspring array for indexing
	# ---------------------------------------------------------------
	Bearpairs_temp[egg_delay] = Bearpairs	
	noOffspring_temp[egg_delay] = noOffspring		
	
	# Return variables from this function
	return Bearpairs_temp, noOffspring_temp
	
	#End::DoMate()