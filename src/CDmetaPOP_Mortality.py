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
	
	# For age/size tracking
	'''
	bin_min = min(sum(sum(size_mean,[]),[]))
	bin_max = max(sum(sum(size_mean,[]),[]))
	size_bin = [bin_min]
	for ibin in xrange(len(size_mean[0][0])-1):
		size_bin.append(size_bin[ibin]+(bin_max - bin_min)/(len(size_mean[0][0])-1))
	'''
	# Get the middles for finding closest values
	size_bin = size_mean[0][0]
	size_mean_middles_bin = np.asarray(size_bin)[1:] - np.diff(np.asarray(size_bin).astype('f'))/2
	for iage in xrange(len(size_bin)):
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
			# First get the age/size death tracking numbers			
			age_adjusted_tracking_index = np.searchsorted(size_mean_middles_bin, SubpopIN_arr['size'])
			countages = count_unique(age_adjusted_tracking_index)
			for thiscount in countages[0]:
				SizeDeaths[gen][thiscount].append(countages[1][np.where(thiscount == countages[0])[0][0]])
			# For age...
			age_tracking_index = SubpopIN_arr['age']
			countages = count_unique(age_tracking_index)
			for thiscount in countages[0]:
				# Special case when more age/size groups then bins
				if thiscount >= len(size_bin):
					# Group in the last bin
					AgeDeaths[gen][len(size_bin)-1].append(countages[1][np.where(thiscount == countages[0])[0][0]])
				else:
					AgeDeaths[gen][thiscount].append(countages[1][np.where(thiscount == countages[0])[0][0]])
				
			# Then delete all
			deleteall = np.asarray(range(0,Npop,1))
			SubpopIN_arr = np.delete(SubpopIN_arr,deleteall)
			SubpopIN_keep.append(SubpopIN_arr)
			# Store new N 
			Population[gen].append(0)
			PopDeaths[gen].append(Npop)
		
		# Stop if nobody in this patch
		elif pop_percmort[isub] != 'E' and Npop == 0:
			# Get tracking numbers
			SubpopIN_keep.append(SubpopIN_arr)
			# Store new N 
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(size_bin)):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)
			
		# Other case
		elif pop_percmort[isub] == 'E' and Npop == 0:
			# Get tracking numbers
			SubpopIN_keep.append(SubpopIN_arr)
			# Store new N 
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(size_bin)):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)
					
		# Keep going if individuals in this patch
		else:
			
			# ----------------------------
			# Apply patch level mortality
			# ----------------------------
			# If wished skipped
			if pop_percmort[isub] == 'N' or pop_percmort[isub] == 0.0 or pop_percmort[isub] == '0' or pop_percmort[isub] == '0.0':
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
			# -------------------------------
			if Nsurvivors != 0:
				# ----------------------------
				# Apply size level mortality
				# ----------------------------
				# First check for multiple classfiles, split up 
				countfiles = count_unique(SubpopIN_keeppop['classfile'])
				# Storage for keeping individuals
				N_samp_ind = []
				# Then for each unique file type
				for ifile in xrange(len(countfiles[0])):
					# Get this files unique patch and file pointer
					thistype = countfiles[0][ifile]							
					natalP = int(thistype.split('_')[0].split('P')[1])
					theseclasspars = int(thistype.split('_')[1].split('CV')[1])
					# The size breaks
					size_mean_middles = np.asarray(size_mean[natalP][theseclasspars])[1:] - np.diff(np.asarray(size_mean[natalP][theseclasspars]).astype('f'))/2
					# All individuals information - careful to index into keeppop
					thistype_index = np.where(SubpopIN_keeppop['classfile']==thistype)[0]
					SubpopIN_keeppop_thistype = SubpopIN_keeppop[thistype_index]
					# 'Age_adjusted'
					age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keeppop_thistype['size'])
					
					# Count up each unique 'sizes'
					countages = count_unique(age_adjusted)
					
					# Loop through age getting Nage in pop		
					for iage in xrange(len(countages[0])):
									
						AgeClass = countages[0][iage]
						# Check for cases in which age is over last age
						if AgeClass > len(size_mean[natalP][theseclasspars])-1:
							# Then make a temp age for indexing
							indexforAgeClass = len(size_mean[natalP][theseclasspars])-1
						else:
							indexforAgeClass = AgeClass
						# Where are these size classes
						Nage_index = np.where(age_adjusted==AgeClass)[0]
						
						# Number in this size class
						Nage = len(Nage_index)
						# Percentage mortality for the age class
						MortAge = size_percmort[natalP][theseclasspars][indexforAgeClass]
											
						# Make sure want to consider this mort percent, all N values are 0 values
						# If there are all N values (1,2, or 3) - no mortalities here.
						if len(np.asarray([i for i, val in enumerate(MortAge) if 'N' in val])) == len(MortAge.split('~')): # If All Ns entered, then skip
							# Keep them all - careful here to index back to original array
							N_samp_ind = N_samp_ind + thistype_index[Nage_index].tolist()
							tempmort = 0
						else: 
							# If there are no N values, but all 0 mort values - no mortalities
							if len(np.asarray([i for i, val in enumerate(MortAge) if 'N' in val])) == 0 and sum(np.asarray(MortAge.split('~'),dtype=float)) == 0.0:
								# Keep them all - careful here to index back to original array
								N_samp_ind = N_samp_ind + thistype_index[Nage_index].tolist()
								tempmort = 0							
							else: # Possibly some N values and mort events
								randnos = np.random.sample(Nage) # Create list of rand nos
								for thisone in xrange(Nage): # Loop through each spot
									# Check the sex on individual and if sex options for this mortality were given
									Indsex = SubpopIN_keeppop[thistype_index[Nage_index[thisone]]]['sex']
									if Indsex == 'XX': # Female
										thisoneMortAge = MortAge.split('~')[0]
									elif Indsex == 'XY': # Male
										# Check if more than one value given
										if len(MortAge.split('~')) > 1:
											thisoneMortAge = MortAge.split('~')[1]
										else:
											thisoneMortAge = MortAge.split('~')[0]
									else: # YY male
										# Check if more than one value given
										if len(MortAge.split('~')) == 3:
											thisoneMortAge = MortAge.split('~')[2]
										elif len(MortAge.split('~')) == 2:
											thisoneMortAge = MortAge.split('~')[1]
										else:
											thisoneMortAge = MortAge.split('~')[0]
									# Could be a N value
									if thisoneMortAge == 'N': # Survives
										# Then store index of survivor
										N_samp_ind = N_samp_ind + [thistype_index[Nage_index[thisone]]]
									else: # Not an N value
										if randnos[thisone] > float(thisoneMortAge): # Survives
											# Then store index of survivor
											N_samp_ind = N_samp_ind + [thistype_index[Nage_index[thisone]]]
										else: # mortality occured, store index of mortality for tracking	
											age_adjusted_tracking_index = np.searchsorted(size_mean_middles_bin, SubpopIN_keeppop[thistype_index[Nage_index[thisone]]]['size'])
											SizeDeaths[gen][age_adjusted_tracking_index].append(1)
				
				# Get left over in this patch
				SubpopIN_keeppop = SubpopIN_keeppop[N_samp_ind]
				
				# ----------------------------
				# Apply age level mortality
				# ----------------------------
				# Storage for keeping individuals
				N_samp_ind = []
				
				# Note same countfiles as above size mortality
				
				# The for each unique file type
				for ifile in xrange(len(countfiles[0])):
					# Get this files unique patch and file pointer
					thistype = countfiles[0][ifile]
					natalP = int(thistype.split('_')[0].split('P')[1])
					theseclasspars = int(thistype.split('_')[1].split('CV')[1])
					# All individuals information - careful to index into keeppop
					thistype_index = np.where(SubpopIN_keeppop['classfile']==thistype)[0]
					SubpopIN_keeppop_thistype = SubpopIN_keeppop[thistype_index]
					# 'Age_adjusted'
					age_adjusted = SubpopIN_keeppop_thistype['age']
				
					# Count up each uniages
					countages = count_unique(age_adjusted)
				
					# Loop through age getting Nage in pop						
					for iage in xrange(len(countages[0])):
									
						AgeClass = countages[0][iage]
						# Check for cases in which age is over last age
						if AgeClass > len(size_mean[natalP][theseclasspars])-1:
							# Then make a temp age for indexing
							indexforAgeClass = len(size_mean[natalP][theseclasspars])-1
						else: 
							indexforAgeClass = AgeClass
						
						# Find all with this age
						Nage_index = np.where(age_adjusted==AgeClass)[0]		
						
						# Number in this age class
						Nage = len(Nage_index)
						# Percentage mortality for the age class
						MortAge = age_percmort[natalP][theseclasspars][indexforAgeClass]
																	
						# Make sure want to consider this mort percent, all N values are 0 values
						# If there are all N values (1,2, or 3) - no mortalities here.
						if len(np.asarray([i for i, val in enumerate(MortAge) if 'N' in val])) == len(MortAge.split('~')): # If All Ns entered, then skip
							# Keep them all - careful here to index back to original array
							N_samp_ind = N_samp_ind + thistype_index[Nage_index].tolist()
							tempmort = 0
						else:
							# If there are no N values, but all 0 mort values - no mortalities
							if len(np.asarray([i for i, val in enumerate(MortAge) if 'N' in val])) == 0 and sum(np.asarray(MortAge.split('~'),dtype=float)) == 0.0:
									# Keep them all - careful here to index back to original array
									N_samp_ind = N_samp_ind + thistype_index[Nage_index].tolist()
									tempmort = 0							
							else: # Possibly some N values and mort events
								randnos = np.random.sample(Nage) # Create list of rand nos
								for thisone in xrange(Nage): # Loop through each spot
									# Check the sex on individual and if sex options for this mortality were given
									Indsex = SubpopIN_keeppop[thistype_index[Nage_index[thisone]]]['sex']
									if Indsex == 'XX': # Female
										thisoneMortAge = MortAge.split('~')[0]
									elif Indsex == 'XY': # Male
										# Check if more than one value given
										if len(MortAge.split('~')) > 1:
											thisoneMortAge = MortAge.split('~')[1]
										else:
											thisoneMortAge = MortAge.split('~')[0]
									else: # YY male
										# Check if more than one value given
										if len(MortAge.split('~')) == 3:
											thisoneMortAge = MortAge.split('~')[2]
										elif len(MortAge.split('~')) == 2:
											thisoneMortAge = MortAge.split('~')[1]
										else:
											thisoneMortAge = MortAge.split('~')[0]
									# Could be a N value
									if thisoneMortAge == 'N': # Survives
										# Then store index of survivor
										N_samp_ind = N_samp_ind + [thistype_index[Nage_index[thisone]]]
									else: # Not an N value
										if randnos[thisone] > float(thisoneMortAge): # Survives
											# Then store index of survivor
											N_samp_ind = N_samp_ind + [thistype_index[Nage_index[thisone]]]
										else: # mortality occured, store index of mortality for tracking	
											AgeDeaths[gen][indexforAgeClass].append(1)
													
				# Grab the survived offspring location
				SubpopIN_keeppop = SubpopIN_keeppop[N_samp_ind]		
				
			# If no one survived from patch level mortality
			# ---------------------------------------------
			else:
				# Get the ones that died and age_adjusted for tracking
				age_adjusted_tracking_index = np.searchsorted(size_mean_middles_bin, SubpopIN_arr['size'])
				countages = count_unique(age_adjusted_tracking_index)
				for thiscount in countages[0]:
					SizeDeaths[gen][thiscount].append(countages[1][np.where(thiscount == countages[0])[0][0]])
				# For age...
				age_tracking_index = SubpopIN_arr['age']
				countages = count_unique(age_tracking_index)
				for thiscount in countages[0]:
					# Special case when more age/size groups then bins
					if thiscount >= len(size_bin):
						# Group in the last bin
						AgeDeaths[gen][len(size_bin)-1].append(countages[1][np.where(thiscount == countages[0])[0][0]])
					else:
						AgeDeaths[gen][thiscount].append(countages[1][np.where(thiscount == countages[0])[0][0]])
						
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
	
	# Sum the AgeDeaths/sizeDeaths - in the size_bin classes if multiple files
	for iage in xrange(len(size_bin)):
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
	# For age/size tracking
	'''
	bin_min = min(sum(sum(size_mean,[]),[]))
	bin_max = max(sum(sum(size_mean,[]),[]))
	size_bin = [bin_min]
	for ibin in xrange(len(size_mean[0][0])-1):
		size_bin.append(size_bin[ibin]+(bin_max - bin_min)/(len(size_mean[0][0])-1))
	'''
	size_bin = size_mean[0][0]
	# Get the middles for finding closest values
	size_mean_middles_bin = np.asarray(size_bin)[1:] - np.diff(np.asarray(size_bin).astype('f'))/2
	for iage in xrange(len(size_bin)):
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
			# First get the age/size death tracking numbers			
			age_adjusted_tracking_index = np.searchsorted(size_mean_middles_bin, SubpopIN_arr['size'])
			countages = count_unique(age_adjusted_tracking_index)
			for thiscount in countages[0]:
				SizeDeaths[gen][thiscount].append(countages[1][np.where(thiscount == countages[0])[0][0]])
			# For age...
			age_tracking_index = SubpopIN_arr['age']
			countages = count_unique(age_tracking_index)
			for thiscount in countages[0]:
				# Special case when more age/size groups then bins
				if thiscount >= len(size_bin):
					# Group in the last bin
					AgeDeaths[gen][len(size_bin)-1].append(countages[1][np.where(thiscount == countages[0])[0][0]])
				else:
					AgeDeaths[gen][thiscount].append(countages[1][np.where(thiscount == countages[0])[0][0]])
			# The delete all
			deleteall = np.asarray(range(0,Npop,1))
			SubpopIN_keeppop = np.delete(SubpopIN_arr,deleteall)
			SubpopIN_keep.append(SubpopIN_keeppop)
			# Store new N 
			Population[gen].append(0)
			PopDeaths[gen].append(Npop)
			
		# If wished entire patch gone
		elif pop_percmort[isub] != 'E' and Npop == 0:
			# Get tracking numbers
			SubpopIN_keep.append(SubpopIN_arr)
			# Store new N 
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(size_bin)):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)
							
		# Other case
		elif pop_percmort[isub] == 'E' and Npop == 0:
			# Get tracking numbers
			SubpopIN_keep.append(SubpopIN_arr)
			# Store new N 
			Population[gen].append(0)
			PopDeaths[gen].append(0)
			for iage in xrange(len(size_bin)):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)
			
		# Keep going if individuals in this patch
		else:
			# ----------------------------
			# Apply patch level mortality
			# ----------------------------
			# Get this patches mortality
			MortPatch = pop_percmort[isub]
			# If this is 0, then skip the rests of the checks, no mortality in this patch
			if MortPatch == 0.0 or MortPatch == '0.0' or MortPatch == '0':
				# Append all information to temp SubpopKeep variable
				SubpopIN_keep.append(SubpopIN_arr)					
				# Store new N 
				Population[gen].append(Npop)
				# And track age deaths
				PopDeaths[gen].append(0)
				for iage in xrange(len(size_bin)):
					AgeDeaths[gen][iage].append(0)	
					SizeDeaths[gen][iage].append(0)
			else:			
										
				# ------------------------------
				# Check for multiply class types
				# ------------------------------
				countfiles = count_unique(SubpopIN_arr['classfile'])
				
				# Storage for keeping individuals
				N_samp_ind = []
				# Then for each unique file type in this patch
				for ifile in xrange(len(countfiles[0])):
					# Get this files unique patch and file pointer
					thistype = countfiles[0][ifile]							
					natalP = int(thistype.split('_')[0].split('P')[1])
					theseclasspars = int(thistype.split('_')[1].split('CV')[1])	
					
					# ----------------------------
					# Apply size level mortality
					# ----------------------------
					# The size breaks
					size_mean_middles = np.asarray(size_mean[natalP][theseclasspars])[1:] - np.diff(np.asarray(size_mean[natalP][theseclasspars]).astype('f'))/2
					# All individuals information - careful to index into arr
					thistype_index = np.where(SubpopIN_arr['classfile']==thistype)[0]
					SubpopIN_keep_thistype = SubpopIN_arr[thistype_index]
					# 'Age_adjusted'
					age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keep_thistype['size'])
						
					# Count up each unique 'sizes'
					countsizes = count_unique(age_adjusted)
					
					# Loop through size/age getting Nage in pop
					for isize in xrange(len(countsizes[0])):
								
						SizeClass = countsizes[0][isize]
						# Check for cases in which age is over last age
						if SizeClass > len(size_mean[natalP][theseclasspars])-1:
							# Then make a temp age for indexing
							indexforSizeClass = len(size_mean[natalP][theseclasspars])-1
						else:
							indexforSizeClass = SizeClass
					
						# Where are these size classes
						Nsize_index = np.where(age_adjusted==SizeClass)[0]
						
						# Number in this size class
						Nsize = len(Nsize_index)
						# Percentage mortality for the size class
						MortSize = size_percmort[natalP][theseclasspars][indexforSizeClass]
						
						# --------------------------------------------
						# Apply age level mortality to this size class
						# --------------------------------------------
						
						# Just this size group in this patch
						#SubpopIN_arr_sizeclass = SubpopIN_arr[Nsize_index] # Change made below on v2.39
						SubpopIN_arr_sizeclass = SubpopIN_keep_thistype[Nsize_index]
					
						# Get ages in this size class
						countages = count_unique(SubpopIN_arr_sizeclass['age'])
					
						# Loop through the ages getting Nage in this size class in this patch...
						for iage in xrange(len(countages[0])):
							AgeClass = countages[0][iage]
							# Check for cases in which age is over last age
							if AgeClass > len(size_mean[natalP][theseclasspars])-1:
								# Then make a temp age for indexing
								indexforAgeClass = len(size_mean[natalP][theseclasspars])-1
							else: 
								indexforAgeClass = AgeClass
							
							# Find all with this age
							Nage_index = np.where(SubpopIN_arr_sizeclass['age']==AgeClass)[0]

							# But be careful to link back to index for size group and then class file group
							
							Nsize_age_index = thistype_index[Nsize_index[Nage_index]] # This index will link back to SubpopIN_arr or patch level
							
							# Number in this age class that is in the size class
							Nage = len(Nage_index)
							# Percentage mortality for the age class
							MortAge = age_percmort[natalP][theseclasspars][indexforAgeClass]
							
							# Then apply mortality for all perc - check on N cases
							# ----------------------------------------------------
							# If there are all N
							if len(np.asarray([i for i, val in enumerate(MortAge) if 'N' in val])) == len(MortAge.split('~')) and len(np.asarray([i for i, val in enumerate(MortSize) if 'N' in val])) == len(MortSize.split('~')) and MortPatch == 'N':
								# Keep them all   - no mortality
								N_samp_ind = N_samp_ind + Nsize_age_index.tolist()
							# Else possible values and apply mortality
							else:
								randnos = np.random.sample(Nage) # Create list of rand nos
								for thisone in xrange(Nage): # Loop through each spot
									# Check the sex on individual and if sex options for this mortality were given
									Indsex = SubpopIN_arr[thistype_index[Nage_index[thisone]]]['sex']
									if Indsex == 'XX': # Female
										thisoneMortAge = MortAge.split('~')[0]
										thisoneMortSize = MortSize.split('~')[0]
									elif Indsex == 'XY': # Male
										# Check if more than one value given
										if len(MortAge.split('~')) > 1:
											thisoneMortAge = MortAge.split('~')[1]
										else:
											thisoneMortAge = MortAge.split('~')[0]
										# Check if more than one value given
										if len(MortSize.split('~')) > 1:
											thisoneMortSize = MortSize.split('~')[1]
										else:
											thisoneMortSize = MortSize.split('~')[0]
									else: # YY male
										# Check if more than one value given
										if len(MortAge.split('~')) == 3:
											thisoneMortAge = MortAge.split('~')[2]
										elif len(MortAge.split('~')) == 2:
											thisoneMortAge = MortAge.split('~')[1]
										else:
											thisoneMortAge = MortAge.split('~')[0]
										# Check if more than one value given
										if len(MortSize.split('~')) == 3:
											thisoneMortSize = MortSize.split('~')[2]
										elif len(MortSize.split('~')) == 2:
											thisoneMortSize = MortSize.split('~')[1]
										else:
											thisoneMortSize = MortSize.split('~')[0]
									
									# Then get total mort values
									if thisoneMortSize == 'N' and thisoneMortAge == 'N' and MortPatch == 'N':
										thisoneMortAll = 0.0
									elif thisoneMortSize == 'N' and thisoneMortAge == 'N' and MortPatch != 'N':
										thisoneMortAll = float(MortPatch)
									elif thisoneMortSize == 'N' and thisoneMortAge != 'N' and MortPatch == 'N':
										thisoneMortAll = float(thisoneMortAge)
									elif thisoneMortSize != 'N' and thisoneMortAge == 'N' and MortPatch == 'N':
										thisoneMortAll = float(thisoneMortSize)									
									elif thisoneMortSize == 'N' and thisoneMortAge != 'N' and MortPatch != 'N':
										thisoneMortAll = float(thisoneMortAge) * float(MortPatch)
									elif thisoneMortSize != 'N' and thisoneMortAge == 'N' and MortPatch != 'N':
										thisoneMortAll = float(thisoneMortSize) * float(MortPatch)
									elif thisoneMortSize != 'N' and thisoneMortAge != 'N' and MortPatch == 'N':
										thisoneMortAll = float(thisoneMortSize) * float(thisoneMortAge)
									elif thisoneMortSize != 'N' and thisoneMortAge != 'N' and MortPatch != 'N':
										thisoneMortAll = float(thisoneMortSize) * float(thisoneMortAge) * float(MortPatch)
									else:
										print('Something off in DoConstantMortality_Multiply() with input values')
										sys.exit(-1)
										
									# No mortality
									if thisoneMortAll == 0.0:
										# Then store index of survivor
										N_samp_ind = N_samp_ind + [Nsize_age_index[thisone]]
										#N_samp_ind = N_samp_ind + Nsize_age_index.tolist()
									else:
										if randnos[thisone] > thisoneMortAll: # Survives
											# Then store index of survivor
											N_samp_ind = N_samp_ind + [Nsize_age_index[thisone]]
											
										else: # mortality occured, add tracking spots
											age_adjusted_tracking_index = np.searchsorted(size_mean_middles_bin, SubpopIN_arr[Nsize_age_index[thisone]]['size'])
											SizeDeaths[gen][age_adjusted_tracking_index].append(1)
											AgeDeaths[gen][indexforAgeClass].append(1)
														
						# End::Age loop
					# End::Size Loop
				# End class file loop
				
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
	for iage in xrange(len(size_bin)):
		AgeDeaths[gen][iage] = int(sum(AgeDeaths[gen][iage]))
		SizeDeaths[gen][iage] = int(sum(SizeDeaths[gen][iage]))
	
	# Return variables from this argument
	tupMort = SubpopIN_keep,Population[gen][1:] # Do not return total for N
	return tupMort
	
	# End::DoConstantMortality_Multiply()
	
# ----------------------------------------------------------------------------------------------	 
def DoMortality(SubpopIN,K,PopDeaths,pop_percmort,age_percmort,gen,Population,AgeDeaths,sizeans,size_mean,size_percmort,SizeDeaths,constMortans,packans,timingans):
	'''
	DoMortality()
	Mortality functions for population.
	Patch level mortality first, then age mortality
	'''	
	# If logistic, but only when in the 'back phase', all other options can occur both times.
	if packans != 'logistic' or (packans == 'logistic' and timingans == 'BACK'):
	
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
		for iage in xrange(len(size_mean[0][0])):
			AgeDeaths[gen].append([])
			SizeDeaths[gen].append([])
		# Loop through each patch
		for isub in xrange(len(K)):
			PopDeaths[gen].append(0)
			Population[gen].append(len(SubpopIN[isub]))
			for iage in xrange(len(size_mean[0][0])):
				AgeDeaths[gen][iage].append(0)	
				SizeDeaths[gen][iage].append(0)
		# Add Population total	
		PopDeaths[gen].insert(0,sum(PopDeaths[gen]))
		Population[gen].insert(0,sum(Population[gen]))
		# Sum the AgeDeaths/sizeDeaths
		for iage in xrange(len(size_mean[0][0])):
			AgeDeaths[gen][iage] = int(sum(AgeDeaths[gen][iage]))
			SizeDeaths[gen][iage] = int(sum(SizeDeaths[gen][iage]))
	
	# Return variables from this argument
	tupAMort = SubpopIN
	return tupAMort
	
	# End::DoMortality()
		
# ----------------------------------------------------------------------------------------------	 
def DoEggMortality(Bearpairs,eggmort_patch,Age0Deaths,gen,K,eggmort_back,noOffspring,Births,BirthsYY):
	'''
	DoEggMortality()
	Mortality functions for age 0 class.
	'''
	
	# Storage for egg deaths
	Age0Deaths.append([])
	Births.append([]) # Spot for generation
	BirthsYY.append([]) # Spot for generation, note this is the offspring number after egg deaths. 
	tempoffloc = []
	# Loop through each patch
	for i in xrange(len(K)):
	
		# Only if pairing occured
		#if Bearpairs[0][0] != -9999:
		if not isinstance(Bearpairs[0][0],int):
						
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
			BirthsYY[gen].append( sum(noOffspring[mothers_patch_ind[np.where(Bearpairs[mothers_patch_ind,1]['sex']=='YY')[0]]]) )
			#BirthsYY[gen].append(sum(noOffspring[mothers_patch_ind[np.where(Bearpairs[mothers_patch_ind[np.where(Bearpairs[mothers_patch_ind,1]['sex']=='YY')[0]],1]['age']==0)[0]]]) )
		# if a pairing did not occur, for tracking
		else:
			# -------------------
			# Tracking numbers
			# -------------------
			Age0Deaths[gen].append(0)
			Births[gen].append(0)
			BirthsYY[gen].append(0)
			
	Age0Deaths[gen].insert(0,sum(Age0Deaths[gen])) # Total the deaths
	Births[gen].insert(0,sum(Births[gen]))
	BirthsYY[gen].insert(0,sum(BirthsYY[gen]))
		
	return noOffspring
	
	# End::DoEggMortality()