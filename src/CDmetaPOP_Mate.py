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

# --------------------------------------------------------------------------
def countDuplicatesInList(dupedList):
	'''
	countDuplicatesInList() - Counts dupicates in lists
	'''
	uniqueSet = Set(item for item in dupedList)
	return [dupedList.count(item) for item in uniqueSet]
	
	# End::countDuplicatesInList()

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
def DoSexual(xycdmatrix,females,males,matemovethresh,Bearpairs,femalesmated,sourcePop,selfing,count=None):
	'''
	DoSexualYY() and DoSexualNY()
	This function is the mating function for: 
	sexual reproduction
	females	with replacement
	males with replacement.
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
	
	# Extract each male patch probability that female can mate with - careful of indexing
	probarray = xycdmatrix[int(femalepop)-1]
	
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
			
			# If selfing is on - only applies to asexual reproduction
			if selfing == 'Y':
				# Then remove intfemale from patchindex
				patchindex = patchindex[np.where(patchindex != intfemale)[0]]
		
			# If there are no males in this patch - search other patches
			if len(patchindex) == 0:				
				# Replace probarray with a zero, no males here
				tempprobarray[itemselect] = 0.
				continue
			
			# There are males in this patch, randomly select one	
			else:
				# Index into the males
				patchmales = males[patchindex]
			
				# Randomly select a male in patch
				malemate = random.sample(patchmales,1)[0]
			
				# And store the mated pair information.						
				Bearpairs.append([females[intfemale],malemate])
				
				# Tracking
				femalesmated.append(1)
				
				# Then break from loop
				break
					
	# Return Variables from this function
	return Bearpairs,femalesmated
	
	# End::DoSexual()

# ---------------------------------------------------------------------------------------------------	
def DoSexualNN(nomales,xycdmatrix,females,\
males,matemovethresh,Bearpairs,femalesmated,subpop,selfing,count=None):
	'''
	DoSexualNN()
	This function is the mating function for
	sexual reproduction
	females	with replacement
	males with replacement
	'''	
	
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
	probarray = xycdmatrix[femalepop]
				
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
			
			# If selfing is on - only applies to asexual reproduction
			if selfing == 'Y':
				# Then remove intfemale from patchindex
				patchindex = patchindex[np.where(patchindex != intfemale)[0]]
		
			# Match male index with patchindex
			patchmales = set(males).intersection(patchindex)
		
			# If there are no males in this patch
			if len(patchmales) == 0:				
				# Replace probarray with a zero, no males here
				tempprobarray[itemselect] = 0.
				continue
				
			else:			
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
def DoMate(SubpopIN,K,freplace,mreplace,matemoveno,matemovethresh,xycdmatrix,MateDistCD,xgrid,ygrid,MateDistCDstd,FAvgMate,MAvgMate,FSDMate,MSDMate,Female_BreedEvents,gen,sourcePop,dtype,ScaleMax,ScaleMin,A,B,C,Femalepercent,eggFreq,sexans,selfing):

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
	
	# ---------------------------------------------------
	# Select males and females for mating
	# ---------------------------------------------------

	# Storage variables for breeding age males and females
	females = []		# These are the zeros
	males = []			# These are the ones
	
	# Loop through and grab each female and male for probability of mature and probability to lay eggs
	for isub in xrange(len(K)):
		indexF = np.where(SubpopIN[isub]['sex']==0)[0]
		indexM = np.where(SubpopIN[isub]['sex']==1)[0]
		allfemales = SubpopIN[isub][indexF]
		allmales = SubpopIN[isub][indexM]
		# Get reproduction age individuals
		indexFage = np.where(allfemales['layeggs'] == 1)[0]
		indexMage = np.where(allmales['mature'] == 1)[0]
		females.append(list(allfemales[indexFage]))
		males.append(list(allmales[indexMage]))
	del allfemales
	del allmales	
		
	# For sexual reproduction
	if sexans == 'Y':
		# Flatten the subpop list for random grabs		
		females = sum(females,[])
		males = sum(males,[])
		females = np.asarray(females,dtype=dtype)
		males = np.asarray(males,dtype=dtype)	
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
		# Flatten the subpop list for random grabs		
		females = sum(females,[])
		males = sum(males,[])
		females = np.asarray(females,dtype=dtype)
		males = np.asarray(males,dtype=dtype)
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
				Bearpairs,femalesmated = DoSexual(xycdmatrix,females,males,matemovethresh,Bearpairs,femalesmated,sourcePop,selfing,count)
												
				# Update count
				count = count + 1
				
		# For the case of a Female with replacement and a male with replacement
		elif freplace == 'Y' and mreplace == 'Y':
				
			# Loop through while loop until all females paired up, but do this total mature individuals times.		
			count = 0		# Initialize the while loop
			while count < looptime:
				
				# Get probability function of user defined input number
				Bearpairs,femalesmated = DoSexual(xycdmatrix,females,males,matemovethresh,Bearpairs,femalesmated,sourcePop,selfing)
							
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
				Bearpairs,tempmales = DoSexualNN(nomales,xycdmatrix,females,tempmales,matemovethresh,Bearpairs,femalesmated,subpop,selfing,count)
										
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
		if Bearpairs[ipair][1]!=-9999: # Make sure a mate pair occurred.
			Floc = int(Bearpairs[ipair][0][sourcePop])-1
			Mloc = int(Bearpairs[ipair][1][sourcePop])-1
			probval = xycdmatrix[Floc][Mloc]			
			
			# If panmictic, can't get cost back
			if matemoveno == '4' or matemoveno == '6' or matemoveno == '9':
				cdval = 0.0
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
			else:
				print('Mate move function does not exist.')
				sys.exit(-1)
			tempAvgMateCD.append(cdval)
			
	# If at least some individuals mated
	if Bearpairs[0][0] != -9999:
			
		# And append to MateDistCD
		MateDistCD.append(sum(tempAvgMateCD) / len(Bearpairs))
		MateDistCDstd.append(np.std(tempAvgMateCD))
		
	# No Mates
	else:

		# And append to MateDistCD
		MateDistCD.append(0)
		MateDistCDstd.append(0) 
	
	# Track actual number of breeding events of females.
	Female_BreedEvents.append(sum(femalesmated))
	del tempAvgMateCD
		
	# Convert Bearpairs to array with dtype
	Bearpairs = np.asarray(Bearpairs)
	
	# Return variables from this function
	return Bearpairs
	
	#End::DoMate()