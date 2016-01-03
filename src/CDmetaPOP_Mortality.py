# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Mortality.py
# Author: Erin L Landguth
# Created: December 2013
# Description: This is the function/module file for mortality processes
# --------------------------------------------------------------------------------------------------

# Import Modules with Except/Try statements

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
	
# Python specific functions
import os, random, copy, pdb, sys

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ----------------------------------------------------------------------------------------------	 
def ConstantMortality_Add(SubpopIN,K,PopDeaths,age_percmort,pop_percmort,gen,Population,AgeDeaths,sizecall,size_mean,size_percmort,SizeDeaths):
	'''
	ConstantMortality_Add()
	Additive model P(A) + P(B) 
	Constant mortality applied to each age within each population.
	Returns: Updated SubpopIN and new N for each subpop.	
	'''
	 
	# Get number of population deaths 
	SubpopIN_keep = []
	Population.append([]) # Add spot for generation
	PopDeaths.append([])
	AgeDeaths.append([])
	SizeDeaths.append([])
	
	# Add spots for Age tracking
	for iage in xrange(len(age_percmort[0])):
		AgeDeaths[gen].append([])
		SizeDeaths[gen].append([])
	
	# Loop through each patch
	for isub in xrange(len(K)):		
	
		# Get each SubpopIN pop as array
		SubpopIN_arr = SubpopIN[isub]
		
		# Get N here
		Npop = len(SubpopIN_arr)
		
		# If wished entire patch gone
		if pop_percmort[isub] == 'E' and Npop != 0:
			deleteall = np.asarray(range(0,Npop,1))
			SubpopIN_arr = np.delete(SubpopIN_arr,deleteall)
			SubpopIN_keep.append(SubpopIN_arr)
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(age_percmort[isub])):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)
		
		# Stop if nobody in this patch or wished skipped
		elif pop_percmort[isub] != 'E' and Npop == 0:
			# Get tracking numbers
			SubpopIN_keep.append(SubpopIN_arr)
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(age_percmort[isub])):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)		
		
		# Other case
		elif pop_percmort[isub] == 'E' and Npop == 0:
			# Get tracking numbers
			SubpopIN_keep.append(SubpopIN_arr)
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(age_percmort[isub])):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)			
		
		# Keep going if individuals in this patch
		else:
			# ----------------------------
			# Apply patch level mortality
			# ----------------------------
			# If wished skipped
			if pop_percmort[isub] == 'N':
				SubpopIN_keeppop = SubpopIN_arr	
				Nsurvivors = len(SubpopIN_arr)
			else:			
				# Get number of survivors
				Nsurvivors = int(round((1-pop_percmort[isub])*Npop))
				
				# Select random survivors
				sample_survivors = np.asarray(random.sample(range(Npop),Nsurvivors))
				
				# Case in which there are 0 survivors
				if Nsurvivors == 0:
					deleteall = np.asarray(range(0,Npop,1))
					SubpopIN_keeppop = np.delete(SubpopIN_arr,deleteall)	
				# Grab the survived offspring location
				else:
					SubpopIN_keeppop = SubpopIN_arr[sample_survivors]	
			
			# Only keep going if any survived
			if Nsurvivors != 0:
				# ----------------------------
				# Apply size level mortality
				# ----------------------------
				# Get age adjusted size classes
				size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keeppop['size'])
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)		
				
				# Loop through age getting Nage in pop
				Nsize_samp_ind = []		
				for iage in xrange(len(countages[0])):
								
					AgeClass = countages[0][iage]
					# Check for cases in which age is over last age
					if AgeClass > len(size_mean[isub])-1:
						# Then make a temp age for indexing
						indexforAgeClass = len(size_mean[isub])-1
					else:
						indexforAgeClass = AgeClass
					# Where are these size classes
					Nage_index = np.where(age_adjusted==AgeClass)[0]
					
					# Number in this size class
					Nage = len(Nage_index)
					# Percentage mortality for the age class
					MortAge = size_percmort[isub][indexforAgeClass]
					
					# Make sure want to consider this mort percent
					if MortAge == 'N' or MortAge == 0.0:
						# Keep them all
						Nsize_samp_ind = Nsize_samp_ind + Nage_index.tolist()
						tempmort = 0				
					# Get a random number for the length of Nage and check survival
					elif MortAge > 0.0: # If there is a possible mort event
						randnos = np.random.sample(Nage) # Create list of rand nos
						tempmort = 0 # tracker
						for thisone in xrange(Nage): # Loop through each spot
							if randnos[thisone] > MortAge: # Survives
								# Then store index of survivor
								Nsize_samp_ind = Nsize_samp_ind + [Nage_index[thisone]]
							else: # mortality occured	
								tempmort = tempmort + 1
										
					# Add spots to tracking
					SizeDeaths[gen][indexforAgeClass].append(tempmort)
				
				# Get left over in this patch
				SubpopIN_keeppop = SubpopIN_keeppop[Nsize_samp_ind]
			
				# ----------------------------
				# Apply age level mortality
				# ----------------------------
				# Count up each uniages
				countages = count_unique(SubpopIN_keeppop['age'])
				
				# Loop through age getting Nage in pop
				Nage_samp_ind = []		
				for iage in xrange(len(countages[0])):
								
					AgeClass = countages[0][iage]
					# Check for cases in which age is over last age
					if AgeClass > len(size_mean[isub])-1:
						# Then make a temp age for indexing
						indexforAgeClass = len(size_mean[isub])-1
					else: 
						indexforAgeClass = AgeClass
					
					# Find all with this age
					Nage_index = np.where(SubpopIN_keeppop['age']==AgeClass)[0]		
					
					# Number in this age class
					Nage = len(Nage_index)
					# Percentage mortality for the age class
					MortAge = age_percmort[isub][indexforAgeClass]
					
					# Make sure want to consider this mort percent
					if MortAge == 'N' or MortAge == 0.0:
						# Keep them all
						Nage_samp_ind = Nage_samp_ind + Nage_index.tolist()
						tempmort = 0
					
					# Get a random number for the length of Nage and check survival
					elif MortAge > 0.0: # If there is a possible mort event
						randnos = np.random.sample(Nage) # Create list of rand nos
						tempmort = 0 # tracker
						for thisone in xrange(Nage): # Loop through each spot
							if randnos[thisone] > MortAge: # Survives
								# Then store index of survivor
								Nage_samp_ind = Nage_samp_ind + [Nage_index[thisone]]
							else: # mortality occured	
								tempmort = tempmort + 1
										
					# Add spots to tracking
					AgeDeaths[gen][indexforAgeClass].append(tempmort)
				
				# Grab the survived offspring location
				SubpopIN_keeppop = SubpopIN_keeppop[Nage_samp_ind]		
				
			# If no one survived from patch level mortality
			else:
				for iage in xrange(len(age_percmort[isub])):
					AgeDeaths[gen][iage].append(0)	
					SizeDeaths[gen][iage].append(0)	
			
			# Append all information to temp SubpopKeep variable
			SubpopIN_keep.append(SubpopIN_keeppop)
			
			# Store new N 
			Population[gen].append(len(SubpopIN_keep[isub]))
			# And track age deaths
			PopDeaths[gen].append(Npop - len(SubpopIN_keep[isub]))
		
			# Delete variables
			del(SubpopIN_keeppop)
			del(SubpopIN_arr)
		
	# Add Population total
	Population[gen].insert(0,sum(Population[gen]))	
	PopDeaths[gen].insert(0,sum(PopDeaths[gen]))
	# Sum the AgeDeaths/sizeDeaths
	for iage in xrange(len(age_percmort[0])):
		AgeDeaths[gen][iage] = int(sum(AgeDeaths[gen][iage]))
		SizeDeaths[gen][iage] = int(sum(SizeDeaths[gen][iage]))
	
	# Return variables from this argument
	tupMort = SubpopIN_keep,Population[gen][1:] # Do not return total for N
	return tupMort
	
	# End::DoConstantMortality_Add()

	# ----------------------------------------------------------------------------------------------	 
def ConstantMortality_Multiply(SubpopIN,K,PopDeaths,age_percmort,pop_percmort,gen,Population,AgeDeaths,sizecall,size_mean,size_percmort,SizeDeaths):
	'''
	ConstantMortality()
	Constant mortality applied to each age within each population.
	Returns: Updated SubpopIN and new N for each subpop.	
	'''
	 
	# Get number of population deaths 
	SubpopIN_keep = []
	Population.append([]) # Add spot for generation
	PopDeaths.append([])
	AgeDeaths.append([])
	SizeDeaths.append([])
	# Add spots for Age tracking
	for iage in xrange(len(age_percmort[0])):
		AgeDeaths[gen].append([])
		SizeDeaths[gen].append([])
	
	# Loop through each patch
	for isub in xrange(len(K)):
	
		# Get each SubpopIN pop as array
		SubpopIN_arr = SubpopIN[isub]

		# Get N here
		Npop = len(SubpopIN_arr)
		
		# Stop if nobody in this patch
		if pop_percmort[isub] == 'E' and Npop != 0:
			deleteall = np.asarray(range(0,Npop,1))
			SubpopIN_keeppop = np.delete(SubpopIN_arr,deleteall)
			SubpopIN_keep.append(SubpopIN_keeppop)
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(age_percmort[isub])):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)		

		# If wished entire patch gone
		elif pop_percmort[isub] != 'E' and Npop == 0:
			# Get tracking numbers
			SubpopIN_keep.append(SubpopIN_arr)
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(age_percmort[isub])):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)	
									
		# Other case
		elif pop_percmort[isub] == 'E' and Npop == 0:
			# Get tracking numbers
			SubpopIN_keep.append(SubpopIN_arr)
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(age_percmort[isub])):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)	
		
		# Keep going if individuals in this patch
		else:
			# ----------------------------
			# Apply patch level mortality
			# ----------------------------
			# Get this patches mortality
			MortPatch = pop_percmort[isub]
			# If 'N', then ignore this perc
			if MortPatch == 'N':
				MortPatch_use = 1.0
			else:
				MortPatch_use = MortPatch
							
			# ----------------------------
			# Apply size level mortality
			# ----------------------------
			# Get age adjusted size classes
			size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
			age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_arr['size'])
			# Count up each unique 'sizes'
			countsizes = count_unique(age_adjusted)
			
			# Loop through age getting Nage in pop
			N_samp_ind = []		
			for isize in xrange(len(countsizes[0])):
							
				SizeClass = countsizes[0][isize]
				# Check for cases in which age is over last age
				if SizeClass > len(size_mean[isub])-1:
					# Then make a temp age for indexing
					indexforSizeClass = len(size_mean[isub])-1
				else:
					indexforSizeClass = SizeClass
				# Where are these size classes
				Nsize_index = np.where(age_adjusted==SizeClass)[0]
				
				# Number in this size class
				Nsize = len(Nsize_index)
				# Percentage mortality for the size class
				MortSize = size_percmort[isub][indexforSizeClass]
				# If 'N', then ignore this perc
				if MortSize == 'N':
					MortSize_use = 1.0
				else:
					MortSize_use = MortSize
				
				# --------------------------------------------
				# Apply age level mortality to this size class
				# --------------------------------------------
				
				# Just this size group in this patch
				SubpopIN_arr_sizeclass = SubpopIN_arr[Nsize_index]
				
				# Get ages in this size class
				countages = count_unique(SubpopIN_arr_sizeclass['age'])
				
				# Tracking mortality in this size class
				tempmort_size = 0
				
				# Loop through the ages getting Nage in this size class in this patch...
				for iage in xrange(len(countages[0])):
					AgeClass = countages[0][iage]
					# Check for cases in which age is over last age
					if AgeClass > len(size_mean[isub])-1:
						# Then make a temp age for indexing
						indexforAgeClass = len(size_mean[isub])-1
					else: 
						indexforAgeClass = AgeClass
					
					# Find all with this age
					Nage_index = np.where(SubpopIN_arr_sizeclass['age']==AgeClass)[0]

					# But becareful to link back to index for size group
					Nsize_age_index = Nsize_index[Nage_index] # This index will link back to SubpopIN_arr or patch level
					
					# Number in this age class that is in the size class
					Nage = len(Nage_index)
					# Percentage mortality for the age class
					MortAge = age_percmort[isub][indexforAgeClass]
					# If 'N' then ignore
					if MortAge == 'N':
						MortAge_use = 1.0
					else:
						MortAge_use = MortAge						
					
					# Then apply mortality for all perc - check on N cases
					# ----------------------------------------------------
					if MortAge == 'N' and MortPatch == 'N' and MortSize == 'N':
						MortAll = 0.0
					else:
						MortAll = MortPatch_use * MortSize_use * MortAge_use
					
					# If no mortality here
					if MortAll == 0.0:
						# Keep them all
						N_samp_ind = N_samp_ind + Nsize_age_index.tolist()
						tempmort_age = 0
					# Else apply mortality
					else:
						randnos = np.random.sample(Nage) # Create list of rand nos
						
						tempmort_age = 0 # tracker
						for thisone in xrange(Nage): # Loop through each spot
							if randnos[thisone] > MortAll: # Survives
								# Then store index of survivor
								N_samp_ind = N_samp_ind + [Nsize_age_index[thisone]]
							else: # mortality occured	
								tempmort_age = tempmort_age + 1
						
					# Add spots to tracking
					AgeDeaths[gen][indexforAgeClass].append(tempmort_age)
					tempmort_size = tempmort_size + tempmort_age
				# End::Age loop
				
				# Add spots to tracking size
				SizeDeaths[gen][indexforSizeClass].append(tempmort_size)
			# End::Size Loop
		
			# Append all information to temp SubpopKeep variable
			SubpopIN_keep.append(SubpopIN_arr[N_samp_ind])
			
			# Store new N 
			Population[gen].append(len(SubpopIN_keep[isub]))
			# And track age deaths
			PopDeaths[gen].append(Npop - len(SubpopIN_keep[isub]))
	# End::Patch loop
	
	# Add Population total
	Population[gen].insert(0,sum(Population[gen]))	
	PopDeaths[gen].insert(0,sum(PopDeaths[gen]))
	# Sum the AgeDeaths/sizeDeaths
	for iage in xrange(len(age_percmort[0])):
		AgeDeaths[gen][iage] = int(sum(AgeDeaths[gen][iage]))
		SizeDeaths[gen][iage] = int(sum(SizeDeaths[gen][iage]))
	
	# Return variables from this argument
	tupMort = SubpopIN_keep,Population[gen][1:] # Do not return total for N
	return tupMort
	
	# End::DoConstantMortality_Multiply()
	
# ----------------------------------------------------------------------------------------------	 
def DoMortality(SubpopIN,K,PopDeaths,pop_percmort,age_percmort,gen,Population,AgeDeaths,sizeans,size_mean,size_percmort,SizeDeaths,constMortans,packans):
	'''
	DoMortality()
	Mortality functions for population.
	Patch level mortality first, then age mortality
	'''	
	
	# If not logistic - stage-population model
	if packans != 'logistic':
	
		# Get size or age control here
		if sizeans == 'Y':
			sizecall = 'size'
		elif sizeans == 'N':
			sizecall = 'age'
		else:
			print('Specify Y or N for size control parameters.')
			sys.exit(-1)
		
		# Straight mortality percentages split by age
		if constMortans == '1': # This is the additive model
			tupMort = ConstantMortality_Add(SubpopIN,K,PopDeaths,age_percmort,pop_percmort,gen,Population,AgeDeaths,sizecall,size_mean,size_percmort,SizeDeaths)
		else: # This is the multiplicative model
			tupMort = ConstantMortality_Multiply(SubpopIN,K,PopDeaths,age_percmort,pop_percmort,gen,Population,AgeDeaths,sizecall,size_mean,size_percmort,SizeDeaths)	
		
		SubpopIN = tupMort[0]
		N = tupMort[1]
	
	# For tracking numbers
	else:
		PopDeaths.append([])
		AgeDeaths.append([])
		SizeDeaths.append([])
		Population.append([])
		# Add spots for Age tracking
		for iage in xrange(len(age_percmort[0])):
			AgeDeaths[gen].append([])
			SizeDeaths[gen].append([])
		# Loop through each patch
		for isub in xrange(len(K)):
			PopDeaths[gen].append(0)
			Population[gen].append(len(SubpopIN[isub]))
			for iage in xrange(len(age_percmort[isub])):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)
		# Add Population total	
		PopDeaths[gen].insert(0,sum(PopDeaths[gen]))
		Population[gen].insert(0,sum(Population[gen]))
		# Sum the AgeDeaths/sizeDeaths
		for iage in xrange(len(age_percmort[0])):
			AgeDeaths[gen][iage] = int(sum(AgeDeaths[gen][iage]))
			SizeDeaths[gen][iage] = int(sum(SizeDeaths[gen][iage]))
	
	# Return variables from this argument
	tupAMort = SubpopIN
	return tupAMort
	
	# End::DoMortality()
		
# ----------------------------------------------------------------------------------------------	 
def DoEggMortality(Bearpairs,eggmort_patch,Age0Deaths,gen,K,eggmort_back,noOffspring,Births):
	'''
	DoEggMortality()
	Mortality functions for age 0 class.
	'''
	
	# Storage for egg deaths
	Age0Deaths.append([])
	Births.append([]) # Spot for generation
	tempoffloc = []
	# Loop through each patch
	for i in xrange(len(K)):
	
		# Only if pairing occured
		if Bearpairs[0][0] != -9999:
						
			# Get all the mothers in this patch
			mothers_patch = Bearpairs[:,0]
			mothers_patch_ind = np.where(mothers_patch['NatalPop']==str(i+1))[0]
			
			# Get the offspring in patch
			offspring_patch = noOffspring[mothers_patch_ind]
			Popoffspring = sum(offspring_patch)
			
			# Only keep going if anybody in this patch 
			if sum(offspring_patch) != 0:						
				# ----------------------------
				# Apply patch level mortality
				# ----------------------------		
				# Get number of survivors
				offsurvivors_patch = int(round((1.-eggmort_patch[i])*sum(offspring_patch)))
				# Split based on proportion	
				offspring_patch = np.around(offspring_patch * float(offsurvivors_patch)/sum(offspring_patch))
					
				# Only keep going if still survivors 
				if sum(offspring_patch) != 0:
				# ----------------------------
				# Apply class level mortality
				# ----------------------------				
					# Get number of survivors
					offsurvivors_back = int(round((1.-eggmort_back) * sum(offspring_patch)))
					
					# Split based on proportion only if
					offspring_patch = np.around(offspring_patch * float(offsurvivors_back)/sum(offspring_patch))
					offspring_patch = np.asarray(offspring_patch,dtype='int')
				
				# Linux error, needs to be same type
				offspring_patch = np.array(offspring_patch,dtype=int)
				# -------------------------
				# Update offspring numbers
				# -------------------------
				noOffspring[mothers_patch_ind] = offspring_patch
			
			# -------------------
			# Tracking numbers
			# -------------------
			Age0Deaths[gen].append(Popoffspring - sum(offspring_patch))
			Births[gen].append(Popoffspring)
			
		# if a pairing did not occur, for tracking
		else:
			# -------------------
			# Tracking numbers
			# -------------------
			Age0Deaths[gen].append(0)
			Births[gen].append(0)
			
	Age0Deaths[gen].insert(0,sum(Age0Deaths[gen])) # Total the deaths
	Births[gen].insert(0,sum(Births[gen]))
		
	return noOffspring
	
	# End::DoEggMortality()