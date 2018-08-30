# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Emigration.py
# Author: Erin L Landguth
# Created: December 2013
# Description: This is the function/module file for emigration processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
	
# Python specific functions
import pdb, random, copy, os, sys
from ast import literal_eval 
from CDmetaPOP_Modules import *

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
		print("%s"%(msg))
		
	# End::logMsg()

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
	# The case where all of the values in lst are the same
	if len(lst) == count:
		count = count-1
	return item,count
	
	#End::w_choice_general()	

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
def GetProbArray(Fxycdmatrix,Mxycdmatrix,offspring,currentsubpop,K,migrate):
	'''
	GetProbArray()
	This function gets indices for F and M specific cdmatrix values
	In direction cost, this is the column value xycdmatrix[0] grabs column in original cdmatrix; xycdmatrix[:,0] grabs row vals in original cdmatrix
	'''
	
	# Index into offspring array
	currentoff = offspring
	
	# Subpopulation current offspring is in - index - 1
	currentsubpop = int(currentsubpop) - 1

	# Get sex of individual
	indSex = currentoff['sex']
	
	# Append the freegrid probabilities for the offspring choice
	if indSex == 'XX': # Female offspring
		probarray = Fxycdmatrix[currentsubpop]
	elif indSex == 'XY': # Male offspring
		probarray = Mxycdmatrix[currentsubpop]
	elif indSex == 'YY': # Assume same as XY
		probarray = Mxycdmatrix[currentsubpop]
	else:
		print('Invalid offspring list.')
		sys.exit(-1)		
	
	# Where K = 0, turn prob to 0
	probarray[np.where(np.asarray(K)==0)[0]] = 0.
	
	# Where migrate = 0, turn prob to 0
	probarray[np.where(np.asarray(migrate)==0)[0]] = 0.

	return probarray
	
	# End::GetProbArray()
	
# ---------------------------------------------------------------------------------------------------	
def Emigration(SubpopIN,K,Fdispmoveno,Mdispmoveno,\
Fxycdmatrix,Mxycdmatrix,gen,\
cdevolveans,fitvals,SelectionDeaths,DisperseDeaths,\
burningen,ProbPatch,ProbSuccess,AdultNoMg,totalA,ProbAge,Population,sourcePop,dtype,setmigrate,sizecall,size_mean,PackingDeaths,PopulationAge,loci,muterate,mtdna,mutationans,packans,PackingDeathsAge,ithmcrundir,packpar1,timecdevolve,age_percmort,migrate,patchvals,PopTag,subpopmort_mat):

	'''
	DoEmigration()
	This function enforces emigration or movement of individuals
	out of natal populations.
	'''	
	
	# Create variable to store outpool 
	SubpopIN_keep = []	
	
	# Get unique number of subpops
	nosubpops = len(K)
	
	# Add spot to track disperses: dispersing deaths for cdevolve
	NoMg = []
	DisperseDeaths.append([]) # Tried to disperse, but no open spots
	SelectionDeaths.append([]) # Dispersing Deaths cdevolve
	PackingDeaths.append([]) # Deaths due to packing density dependence problem
	ProbSuccess.append([]) # Prob switches
	for x in xrange(0,nosubpops):
		DisperseDeaths[gen].append([])
		SelectionDeaths[gen].append([])
		PackingDeaths[gen].append([])
		SubpopIN_keep.append([])
				
	# Decide where everybody moves - loop through subpop, doesn't need to be shuffled.
	for isub in xrange(len(SubpopIN)):
			
		# Loop through individuals in this subpopulation
		for iind in xrange(len(SubpopIN[isub])):
			
			# Get individual, easier for indexing
			outpool = SubpopIN[isub][iind]
			originalpop = outpool[sourcePop]
			Indsex = outpool['sex']
			
			# Get this individuals original ClassVars file and bins for indexing
			natalP = int(SubpopIN[isub][iind]['classfile'].split('_')[0].split('P')[1])
			theseclasspars = int(SubpopIN[isub][iind]['classfile'].split('_')[1].split('CV')[1])
			
			# If setmigration is turned on
			if setmigrate[isub] == 'Y' and 'I' in outpool['name']:
				indProbans = 'Yes'
							
			# Else no setmigration or outpool is not an Immigrator yet
			else: 
							
				# Check for patch migration probability
				Mg_Patch = ProbPatch[int(originalpop)-1] # index 1 minus
				
				# Check for age/size migration
				indexofProb = outpool[sizecall]
				# If size control then get size nearest to values in age file
				# Index to the correct age file of this original individual
				if sizecall == 'size': 
					closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb)) #Closest size
					Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0] #Index
				else:
					Find = indexofProb
				# Check for ages over last age
				if Find > len(ProbAge[natalP][theseclasspars]) - 1:
					Find = len(ProbAge[natalP][theseclasspars]) - 1 # Make last age					
				Mg_Class = ProbAge[natalP][theseclasspars][Find]
				
				# Get sex class options if given
				if Indsex == 'XX': # Female
					Mg_Class = float(Mg_Class.split('~')[0])
				elif Indsex == 'XY': #Male
					# Check if more than value here
					if len(Mg_Class.split('~')) > 1:
						Mg_Class = float(Mg_Class.split('~')[1])
					else:
						Mg_Class = float(Mg_Class.split('~')[0])
				else: # YY male
					# Check if more than value here
					if len(Mg_Class.split('~')) == 3:
						Mg_Class = float(Mg_Class.split('~')[2])
					elif len(Mg_Class.split('~')) == 2:
						Mg_Class = float(Mg_Class.split('~')[1])
					else:
						Mg_Class = float(Mg_Class.split('~')[0])
				
				# Then multiply these together
				indProb = Mg_Patch * Mg_Class
			
				randProb = rand()	# Get a random number				
				# Flip the coin for migration
				if randProb < indProb:					
					# Then migrate out....
					indProbans = 'Yes'				
				# Patch migration not a success
				else:
					indProbans = 'No'					
				
			# Then check migration success
			if indProbans == 'Yes':
				
				# Then Migrate Out....			
				# Create a function here that gets indices for male and female
				probarray = GetProbArray(Fxycdmatrix,Mxycdmatrix,outpool,originalpop,K,migrate)
										
				# If statement to check if there are spots for offpsring to disperse to
				if sum(probarray) != 0.0:
					
					# CDEVOLVE
					if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen) and (timecdevolve.find('Out') != -1):		
												
						# Select the w_choice item
						iteminlist = w_choice_item(probarray)
						
						# for option 3 in which has to be mature
						if cdevolveans == '1_mat' and outpool['mature'] == 0:
							differentialmortality = 0.0
						else:						
							# Call 1-locus selection model
							differentialmortality = Do1LocusSelection(fitvals,literal_eval(outpool['genes'])[0],iteminlist)
												
						# Then flip the coin to see if outpool survives its location
						randcheck = rand()
						dispersingto = iteminlist
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(1)
							continue
													
					# CDEVOLVE - 2 loci
					elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen) and (timecdevolve.find('Out') != -1):
						
						# Select the w_choice item
						iteminlist = w_choice_item(probarray)

						# for option 3 in which has to be mature
						if cdevolveans == '2_mat' and outpool['mature'] == 0:
							differentialmortality = 0.0
						else:
							# Call 2-locus selection model
							differentialmortality = Do2LocusSelection(fitvals,literal_eval(outpool['genes'])[0:2],iteminlist)
												
						# Then flip the coin to see if outpool survives its location
						randcheck = rand()
						dispersingto = iteminlist
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(1)
							continue
							
					# CDEVOLVE - Hindex
					elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen) and (timecdevolve.find('Out') != -1):
						
						# Select the w_choice item
						iteminlist = w_choice_item(probarray)

						# Call Hindex selection model
						differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
												
						# Then flip the coin to see if outpool survives its location
						randcheck = rand()
						dispersingto = iteminlist
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(1)
							continue
					
					# If subpopulation differential mortality is on
					elif not isinstance(subpopmort_mat,str):
						
						# Select the w_choice item
						iteminlist = w_choice_item(probarray)
						
						# What subpatchID is individual coming from
						from_subpatch = PopTag[isub]
												
						# What subpatchID is individual proposing to go to
						to_subpatch = PopTag[iteminlist]
						
						# If it is dispersing to another subpatchID
						if from_subpatch != to_subpatch:
							
							# grab the differential mortality associated with moving into this new subpatchID - from subpatch TO subpatch - cols are TO, read row, then col for correct spot
							differentialmortality = subpopmort_mat[int(to_subpatch)-1][int(from_subpatch)-1]
							# Check if survives
							randcheck = rand()
							dispersingto = iteminlist
							# If it does not survive
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue				
					
					# If not cdevolve or if cdevolve but it is in burn in gen
					else:
						
						# Select the w_choice item
						iteminlist = w_choice_item(probarray)	
						
					# Append information to variable [outpool, pop it dispersed to Str() + 1, and name]
					# Name if previous individual or not an offspring - {Initial,Residor,Immigrant,Emigrant,Stayor}_{Year born}_{Natal Pop}_{Numeric ID}
					tosubpop = str(iteminlist+1)
					outpool_name = outpool['name']
					outpool_name = outpool_name.split('_')
					name = 'E'+str(tosubpop)+'_F'+str(originalpop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]	
					
					# Record string name of OriginalSubpop,ToSubpop,NAsubpop,EmiCD,ImmiCD,age,sex,size,infection,name,capture,layeggs,genes				
					recd = (originalpop,tosubpop,'NA',-9999,-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(tosubpop)-1],outpool['genes'])
								
					# Record outpool disperse information
					SubpopIN_keep[int(tosubpop)-1].append(recd)

					SelectionDeaths[gen][int(tosubpop)-1].append(0)
					DisperseDeaths[gen][int(tosubpop)-1].append(0)
					ProbSuccess[gen].append(1)
					
				# If statement to check if there were not spots to disperse to
				elif sum(probarray) == 0.0:
					
					# Then Break from the loop and move to next outpool
					SelectionDeaths[gen][int(originalpop)-1].append(0)
					DisperseDeaths[gen][int(originalpop)-1].append(1)
					ProbSuccess[gen].append(1)				
					continue
									
			#End::Prob migrate out if statement
			# -------------------------------
			
			# If migration probability not a success
			else:
								
				# CDEVOLVE - 1 locus
				if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen) and (timecdevolve.find('Out') != -1):
											
					# for option 3 in which has to be mature
					if cdevolveans == '1_mat' and outpool['mature'] == 0:
						differentialmortality = 0.0
					else:					
						# Call 1-locus selection model
						differentialmortality = Do1LocusSelection(fitvals,literal_eval(outpool['genes'])[0],int(originalpop)-1)
											
					# Then flip the coin to see if outpool survives its location
					randcheck = rand()
					dispersingto = int(originalpop)-1
					# If outpool did not survive: break from loop, move to next outpool
					if randcheck < differentialmortality:
						SelectionDeaths[gen][dispersingto].append(1)
						DisperseDeaths[gen][dispersingto].append(0)
						ProbSuccess[gen].append(0)
						continue
												
				# CDEVOLVE - 2 loci
				elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen) and (timecdevolve.find('Out') != -1):
					
					# for option 3 in which has to be mature
					if cdevolveans == '2_mat' and outpool['mature'] == 0:
						differentialmortality = 0.0
					else:
						# Call 2-locus selection model
						differentialmortality = Do2LocusSelection(fitvals,literal_eval(outpool['genes'])[0:2],int(originalpop)-1)
										
					# Then flip the coin to see if outpool survives its location
					randcheck = rand()
					dispersingto = int(originalpop)-1
					# If outpool did not survive: break from loop, move to next outpool
					if randcheck < differentialmortality:
						SelectionDeaths[gen][dispersingto].append(1)
						DisperseDeaths[gen][dispersingto].append(0)
						ProbSuccess[gen].append(0)
						continue				
				
				# CDEVOLVE - Hindex
				elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen) and (timecdevolve.find('Out') != -1):

					# Call Hindex selection model
					differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[int(originalpop)-1])
											
					# Then flip the coin to see if outpool survives its location
					randcheck = rand()
					dispersingto = int(originalpop)-1
					# If outpool did not survive: break from loop, move to next outpool
					if randcheck < differentialmortality:
						SelectionDeaths[gen][dispersingto].append(1)
						DisperseDeaths[gen][dispersingto].append(0)
						ProbSuccess[gen].append(0)
						continue				
				
				# If it is an adult: 2 checks, ID is not new and index spots are the same
				outpool_name = outpool['name']
				outpool_name = outpool_name.split('_')
				if len(outpool_name) != 6:
					print('Error in ID field.')
					sys.exit(-1)
				
				name = 'R'+str(originalpop)+'_F'+str(originalpop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]	
					
				# Record string name of OriginalSubpop,ToSubpop,NA,EmiCD,ImmiCD,age,sex,size,infection,name,capture,genes 
				recd = (originalpop,originalpop,'NA',-9999,-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(originalpop)-1],outpool['genes'])
							
				# Record outpool disperse information
				SubpopIN_keep[int(originalpop)-1].append(recd)

				SelectionDeaths[gen][int(originalpop)-1].append(0)
				DisperseDeaths[gen][int(originalpop)-1].append(0)
				ProbSuccess[gen].append(1)
				NoMg.append(1)
				
			#End::Prob migrate stay in else statement
			# -----------------------------------
			
		#End::For loop individual
		#------------------------
		
	#End::For loop Subpop
	#---------------------
	
	# ----------------------------------------------
	# SubpopIN - sort/rank/pack by Kage options
	# ----------------------------------------------	
	# If multiple ClassVars are given then bin min to max
	if sizecall == 'size':
		'''
		bin_min = min(sum(sum(size_mean,[]),[]))
		bin_max = max(sum(sum(size_mean,[]),[]))
		size_bin = [bin_min]
		for ibin in xrange(len(size_mean[0][0])-1):
			size_bin.append(size_bin[ibin]+(bin_max - bin_min)/(len(size_mean[0][0])-1))
		'''
		# Get the middles for finding closest values
		size_bin = size_mean[0][0]
		size_mean_middles = np.asarray(size_bin)[1:] - np.diff(np.asarray(size_bin).astype('f'))/2	
	
	# Temp lists and tracking updates
	SubpopIN_keepK = []
	Population.append([]) #Add spot for next generation
	PackingDeathsAge.append([])
	PopulationAge.append([])
	
	# Age tracking numbers will be for the first size bin given
	PackingDeathsAge[gen] = [[] for x in xrange(0,len(size_mean[0][0]))]
	PopulationAge[gen] = [[] for x in xrange(0,len(size_mean[0][0]))]	
	
	# -------------------
	# Packing is selected
	# -------------------
	if packans == 'packing':

		# Loop through each subpopulation
		for isub in xrange(len(K)):
					
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# Switch here for size or age control
			if sizecall == 'size': # Careful here and use first ClassVars
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_arr[sizecall])
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)
			else:
				# Count up each uniages
				countages = count_unique(SubpopIN_arr['age'])
						
			# K,N for this population
			Kpop = K[isub]
			Npop = sum(countages[1])
			if Npop == 0 or Kpop == 0:
				# Append all information to temp SubpopKeep variable
				Nage_samp_ind = np.arange(0)					
			else:		
				# Check here on numbers
				if Kpop == Npop:
					Kscale = 1.0
				else:
					# Get K_scaling for population
					Kscale = np.log(float(Kpop)/Npop) / (1. - (Npop/float(Kpop)))
				
				# -------------------------
				# Age loop
				# -------------------------
				# Loop through each age class - here recursive, start with largest class first
				Nage_samp_ind = [] 
				Kage_hab_adj = 0.
				Kage_hab_adj_inc = 0.
				for iage in -np.sort(-countages[0]):			
					
					# Age class
					Ageclass = iage
					
					# Special case when age is greater than last age class only used for indexing now
					if Ageclass > len(size_mean[0][0])-1:
						indexforAgeclass = len(size_mean[0][0]) - 1
					else:
						indexforAgeclass = Ageclass
					
					# N for this age coming in
					Nage = countages[1][np.where(countages[0]==Ageclass)[0][0]]
					if sizecall == 'size': # Use the adjusted age classes
						Nage_index = np.where(age_adjusted==Ageclass)[0]
					else:
						Nage_index = np.where(SubpopIN_arr['age']==Ageclass)[0]
					
					# Get Age_scaling for this pop's age class
					Agescale = np.exp(Kscale * (1. - (Nage / float(Kpop))))
					
					# Kage proportion of habitat available (add one to age class) - equation from AB.
					Kage_hab = np.exp(packpar1*(Ageclass+1))
					
					# Class count for redsitributing Khab
					classcount = len(np.where(countages[0]<=Ageclass)[0])
					#classcount = Ageclass+1
					Kage_hab_adj_inc = Kage_hab_adj_inc + (Kage_hab_adj/classcount)
					
					# Adjust the Kage_hab
					Kage_hab = Kage_hab + Kage_hab_adj_inc
					
					# Kage for this population - this is a proportion
					Kage = Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))
										
					# Get Kused
					if Nage <= Kage:
						# Grab all the ages
						Kused = len(Nage_index)
						Nage_samp_ind.append(list(Nage_index))
						PackingDeathsAge[gen][indexforAgeclass].append(0)
					else:
						# Grab a random draw of Kage from numbers - or sort by size and grab the top size classes?
						Kused = int(round(Kage))
						Nage_samp_ind.append(random.sample(Nage_index,Kused))
						PackingDeathsAge[gen][indexforAgeclass].append(Nage-Kused)
					
					# The adjust the habitat or reallocation for next class
					if Kage == 0:
						Kage_hab_adj = Kage_hab
					else:
						Kage_hab_adj = Kage_hab - (Kused * Kage_hab / Kage)
					
			# Clean up index
			Nage_samp_ind = sum(Nage_samp_ind,[])
				
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepK.append(SubpopIN_arr[Nage_samp_ind])
						
			# Store new N - it can be possible to be less than K - its OK - rounding error
			Population[gen].append(len(SubpopIN_keepK[isub]))			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
			PackingDeaths[gen][isub] = Npop - len(Nage_samp_ind)
			
			# Track class size - get age adjusted again. 		
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']
			
			# Tracking age N
			for iage in xrange(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))	
	
	# ---------------------
	# Packing is turned off
	# ---------------------
	elif packans == 'N':
		
		for isub in xrange(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# K,N for this population
			Kpop = K[isub]
			Npop = len(SubpopIN_arr)
			if Npop == 0 or Kpop == 0:
				# Append all information to temp SubpopKeep variable
				Nage_samp_ind = np.arange(len(SubpopIN_arr))
				PackingDeaths[gen][isub] = 0
			
			else:
				# else do nothing if Npop is below Kpop
				if Npop <= Kpop:	
					# Append all information to temp SubpopKeep variable
					Nage_samp_ind = np.arange(len(SubpopIN_arr))
					PackingDeaths[gen][isub] = 0
			
				else:
					# Grab a random draw of Kage from numbers
					Nage_samp_ind = random.sample(np.arange(len(SubpopIN_arr)),Kpop)
					PackingDeaths[gen][isub] = len(SubpopIN_arr) -Kpop					
			
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepK.append(SubpopIN_arr[Nage_samp_ind])
			
			# Get size adjusted age for tracking		
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']
			
			# Tracking age N
			for iage in xrange(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				if len(sizeindex) == 0:
					PopulationAge[gen][iage].append(0)
				else:				
					PopulationAge[gen][iage].append(len(SubpopIN_keepK[isub][sizeindex]['size'].tolist()))
				# Just store 0 for packing deaths age
				PackingDeathsAge[gen][iage].append(0)
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			if len(sizeindex) == 0:
				PopulationAge[gen][iage].append(0)
			else: # Add them to last class
				PopulationAge[gen][iage].append(len(SubpopIN_keepK[isub][sizeindex]['size'].tolist()))
				
			# Store new N - it can be possible to be less than K - its OK - rounding error
			Population[gen].append(len(SubpopIN_keepK[isub]))			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
	
	# ---------------------
	# Packing is turned off - this is a special case in which all inds are kept.
	# ---------------------
	elif packans == 'N_keepeggs':
		
		for isub in xrange(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# K,N for this population
			Kpop = K[isub]
			Npop = len(SubpopIN_arr)
			if Npop == 0 or Kpop == 0:
				# Append all information to temp SubpopKeep variable
				Nage_samp_ind = np.arange(len(SubpopIN_arr))
				PackingDeaths[gen][isub] = 0
			
			else:
				# else do nothing if Npop is below Kpop
				if Npop <= Kpop:	
					# Append all information to temp SubpopKeep variable
					Nage_samp_ind = np.arange(len(SubpopIN_arr))
					PackingDeaths[gen][isub] = 0
			
				else:
					# Grab a random draw of Kage from numbers CHANGE MADE HERE Kpop to Npop
					Nage_samp_ind = random.sample(np.arange(len(SubpopIN_arr)),Npop)
					PackingDeaths[gen][isub] = len(SubpopIN_arr) -Npop					
			
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepK.append(SubpopIN_arr[Nage_samp_ind])
			
			# Get size adjusted age for tracking		
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']
			
			# Tracking age N
			for iage in xrange(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				if len(sizeindex) == 0:
					PopulationAge[gen][iage].append(0)
				else:				
					PopulationAge[gen][iage].append(len(SubpopIN_keepK[isub][sizeindex]['size'].tolist()))
				# Just store 0 for packing deaths age
				PackingDeathsAge[gen][iage].append(0)
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			if len(sizeindex) == 0:
				PopulationAge[gen][iage].append(0)
			else: # Add them to last class
				PopulationAge[gen][iage].append(len(SubpopIN_keepK[isub][sizeindex]['size'].tolist()))
				
			# Store new N - it can be possible to be less than K - its OK - rounding error
			Population[gen].append(len(SubpopIN_keepK[isub]))			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
	
	
	# --------------------------------------------
	# Logistic selected - 
	# -------------------------------------------- 
	elif packans == 'logistic':	
				
		# Loop through each patch
		for isub in xrange(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# Count up each uniages - logistic function uses actual ages even if size control is specified.
			countages = count_unique(SubpopIN_arr['age']) # ages will be 1+
			
			# K,N for this population
			Kpop = K[isub]
			Npop = len(SubpopIN_arr)
			
			# If none in this patch
			if Npop == 0 or Kpop == 0:
				# Append all information to temp SubpopKeep variable
				Nage_samp_ind = np.arange(0)
				
			# If populated - continue with DDMortality
			else:			
				
				# -------------------------
				# Age loop
				# -------------------------
				Nage_samp_ind = []
				
				# Get total number of deaths for each age class
				for iage in countages[0]: # Note only works for ages 1+
					
					# Age class
					Ageclass = iage									
					
					# Current Ni,t or Nage
					Nage = countages[1][np.where(countages[0]==Ageclass)[0][0]]
					Nage_index = np.where(SubpopIN_arr['age']==Ageclass)[0]
											
					# Special case: if last age class, no survivros
					if Ageclass > len(size_mean[0][0])-1:
						mortreturn = Nage # Then apply mortality to all age class.
						indexforAgeclass = len(size_mean[0][0]) - 1
					else: 
						indexforAgeclass = Ageclass
						# Number in next age class
						if len(np.where(countages[0]==(Ageclass+1))[0]) != 0:
							Niplus1 = countages[1][np.where(countages[0]==(Ageclass+1))[0][0]]
						else:
							Niplus1 = 0.
												
						# Ni+1,t * (1 - Nt/K) * (si-1 * Ni,t - Ni+1,t)
						# Careful age_percmort include age 0, so index into agemort is 1 minus iage
						# Careful age_percmort might have multiple files: index to first file age_percmort[isub][0]
						Ntplus1 = Niplus1 + (1. - Npop/float(Kpop))*((1.-float(age_percmort[isub][0][Ageclass-1])) * Nage - float(Niplus1))
						if Ntplus1 <= 0: # If above K, remove to Kpop
							mortreturn = Nage - Kpop
						else:
							mortreturn = Nage - int(Ntplus1)
						
						# If no mortality
						if mortreturn <= 0:
							mortreturn = 0
							# Grab all the ages
							Nage_samp_ind.append(list(Nage_index))
						# If mortality to this age group occurred. 
						else:
							# Grab a random draw of these ages to keep
							Nage_samp_ind.append(random.sample(Nage_index,Nage-mortreturn))
					# Tracker
					PackingDeathsAge[gen][indexforAgeclass].append(mortreturn)
					
			# Clean up index
			Nage_samp_ind = sum(Nage_samp_ind,[])	
			
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepK.append(SubpopIN_arr[Nage_samp_ind])
			
			# Store new N - it can be possible to be less than K - its OK - rounding error
			Population[gen].append(len(SubpopIN_keepK[isub]))			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
			PackingDeaths[gen][isub] = Npop - len(Nage_samp_ind)
			
			# Count up each uniages
			age_adjusted = SubpopIN_keepK[isub]['age']			
			# Tracking age N
			for iage in xrange(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))	
	
	else:
		print('See user manual for population model options.')
		sys.exit(-1)
		
	# Summary numbers
	SelectionDeaths[gen].insert(0,sum(SelectionDeaths[gen]))
	DisperseDeaths[gen].insert(0,sum(DisperseDeaths[gen]))
	PackingDeaths[gen].insert(0,sum(PackingDeaths[gen]))
	ProbSuccess[gen] = sum(ProbSuccess[gen])
	AdultNoMg.append(sum(NoMg))
	# Add Population total
	Population[gen].insert(0,sum(Population[gen]))
	# Age tracking
	for iage in xrange(len(PopulationAge[gen])): 
		PopulationAge[gen][iage] = sum(PopulationAge[gen][iage])
		PackingDeathsAge[gen][iage] = sum(PackingDeathsAge[gen][iage])		
	
	# Return variables from this function
	return SubpopIN_keepK
	
	# End::DoEmigration()
	
# ---------------------------------------------------------------------------------------------------	
def CalculateDispersalMetrics(OffDisperseIN,xgrid,ygrid,\
Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,FDispDistCD,MDispDistCD,\
FDispDistCDstd,MDispDistCDstd,subpopmigration,gen,Fthreshold,Mthreshold,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC):
	'''
	CalculateDispersalMetrics()
	This function calculates how far disperses are moving.
	'''		
	# Store the average dispersal distance offspring went
	# temp variable to store offspring dispersal distance
	FtempAvgDispDistED = []
	MtempAvgDispDistED = []
	FtempAvgDispDistCD = []
	MtempAvgDispDistCD = []
	Fcount = 0
	Mcount = 0
	subpopmigration.append([]) # This adds a spot for next generation
	
	# Loop through each subpop
	for isub in xrange(len(OffDisperseIN)):
		
		# Extract the disperser type
		miIndex = np.asarray([i for i, val in enumerate(OffDisperseIN[isub]['name']) if 'E' in val])
		if len(miIndex) != 0:
			Ind = OffDisperseIN[isub][miIndex]
		else:
			Ind = []
		
		# Loop through each OffDisperseIN
		for ioffspring in xrange(len(Ind)):
			# Grab information from this individual
			indFrom = Ind['NatalPop'][ioffspring]
			indTo = Ind['EmiPop'][ioffspring]
			indSex = Ind['sex'][ioffspring]
			
			# Store migration numbers
			if indTo != indFrom:
				subpopmigration[gen][int(indTo)-1].append(1)			
			
			# If female - CD distances			
			if indSex == 'XX':
				Fcount = Fcount + 1			
				probval = Fxycdmatrix[int(indFrom)-1][int(indTo)-1]
				
				# If panmictic
				if Fdispmoveno == '4' or Fdispmoveno == '6': 
					cdval = 0.
					
				elif Fdispmoveno == '9':
					cdval = probval
				
				# If linear
				elif Fdispmoveno == '1':
					cdval = (probval - 1.) * (-Fthreshold)				
				
				# If inverse square
				elif Fdispmoveno == '2':
					if probval == 1.0:
						cdval = 0.0
					else:	
						cdval = np.sqrt(1. / (probval * (FScaleMax - FScaleMin) + FScaleMin))
					
				# If neg exponetial
				elif Fdispmoveno == '5':
					cdval = np.log((probval * (FScaleMax-FScaleMin) + FScaleMin)/float(FA)) / (-float(FB) * np.log(10))
				# If Gaussian	
				elif Fdispmoveno == '7':
					cdval = float(FB) + np.sqrt(-2*float(FC)**2 * np.log((probval*(FScaleMax-FScaleMin)+FScaleMin)/float(FA)))
				# If matrix
				elif Fdispmoveno == '8':
					cdval = (1. - probval)*(FScaleMax-FScaleMin)+FScaleMin
				# If pareto
				elif Fdispmoveno == '10':
					if probval == 1.0:
						cdval = 0.0
					else:
						cdval = pow(((float(FA)*float(FB)**float(FA))/probval),(1/(float(FA)+1))) - float(FB)
				# Write to temp	
				FtempAvgDispDistCD.append(cdval)				
					
			# Else if Male
			elif indSex == 'XY' or indSex == 'YY': # assume males same movement matrix 			
				Mcount = Mcount + 1
				probval = Mxycdmatrix[int(indFrom)-1][int(indTo)-1]
				
				# If panmictic
				if Mdispmoveno == '4' or Mdispmoveno == '6': 
					cdval = 0.0
				elif Mdispmoveno == '9':
					cdval = probval
				# If linear
				elif Mdispmoveno == '1':					
					cdval = (probval - 1.) * (-Mthreshold)
					
				# If inverse square
				elif Mdispmoveno == '2':
					if probval == 1.0:
						cdval = 0.0
					else:	
						cdval = np.sqrt(1. / (probval * (MScaleMax - MScaleMin) + MScaleMin))
					
				# If neg exponetial
				elif Mdispmoveno == '5':
					cdval = np.log((probval * (MScaleMax-MScaleMin) + MScaleMin)/float(MA)) / (-float(MB) * np.log(10))
				elif Mdispmoveno == '7':
					cdval = float(MB) + np.sqrt(-2*float(MC)**2 * np.log((probval*(MScaleMax-MScaleMin)+MScaleMin)/float(MA)))
				elif Mdispmoveno == '8':
					cdval = (1. - probval)*(MScaleMax-MScaleMin)+MScaleMin
				# If pareto
				elif Mdispmoveno == '10':
					if probval == 1.0:
						cdval = 0.0
					else:
						cdval = pow(((float(MA)*float(MB)**float(MA))/probval),(1/(float(MA)+1))) - float(MB)
				MtempAvgDispDistCD.append(cdval)
			
			# Store the traveled distance - carefully index here
			OffDisperseIN[isub][miIndex[ioffspring]]['EmiCD'] = cdval
			
		# The store the number of disperses to separate subpops
		subpopmigration[gen][isub]=sum(subpopmigration[gen][isub])
		subpopmigration[gen+1].append([0]) # This adds spots for subpops in	
	
	# If at least 1 Female offspring dispersed
	if Fcount > 0:		
		# And append to DispDistCD
		FDispDistCD.append(sum(FtempAvgDispDistCD)/Fcount)
		FDispDistCDstd.append(np.std(FtempAvgDispDistCD))
	else:
		# And append to DispDistCD
		FDispDistCD.append(0)
		FDispDistCDstd.append(0)
	
	# If at least 1 Male offspring dispersed
	if Mcount > 0:		
		# And append to DispDistCD
		MDispDistCD.append(sum(MtempAvgDispDistCD)/Mcount)
		MDispDistCDstd.append(np.std(MtempAvgDispDistCD))
	else:
		# And append to DispDistCD
		MDispDistCD.append(0)
		MDispDistCDstd.append(0)
	
	# End::CalculateDispersalMetrics()
	
# ---------------------------------------------------------------------------------------------------	 
def DoEmigration(SubpopIN,K,Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,gen,xgridcopy,ygridcopy,FDispDistCD,MDispDistCD,cdevolveans,fitvals,FDispDistCDstd,MDispDistCDstd,subpopmigration,SelectionDeaths,DisperseDeaths,burningen,Prob,ProbSuccess,AdultNoMg,totalA,ProbAge,Fthreshold,Mthreshold,Population,sourcePop,dtype,setmigrate,sizeans,size_mean,PackingDeaths,PopulationAge,loci,muterate,mtdna,mutationans,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC,packans,PackingDeathsAge,ithmcrundir,packpar1,timecdevolve,age_percmort,migrate,patchvals,PopTag,subpopmort_mat):
	'''
	DoEmigration()
	Disperse the individuals to patch locations
	Input: Units of dipsersal, movement function,
	SubpopIN, cdmatrix 
	Output: SubpopIN = [subpop,age,sex,infection,name,genes]
	'''		
	
	# Get size or age control here - for Mg Prob value
	if sizeans == 'Y':
		sizecall = 'size'
	elif sizeans == 'N':
		sizecall = 'age'
	else:
		print('Specify Y or N for size control parameters.')
		sys.exit(-1)	
	
	# Check if skip module
	# --------------------
	if not (isinstance(Fxycdmatrix,str) and isinstance(Mxycdmatrix,str)):
		SubpopIN = Emigration(SubpopIN,K,Fdispmoveno,\
		Mdispmoveno,\
		Fxycdmatrix,Mxycdmatrix,gen,\
		cdevolveans,fitvals,SelectionDeaths,DisperseDeaths,burningen,Prob,ProbSuccess,AdultNoMg,totalA,ProbAge,Population,sourcePop,dtype,setmigrate,sizecall,size_mean,PackingDeaths,PopulationAge,loci,muterate,mtdna,mutationans,packans,PackingDeathsAge,ithmcrundir,packpar1,timecdevolve,age_percmort,migrate,patchvals,PopTag,subpopmort_mat)
		
		# Calculate Dispersal Metrics for movers out
		CalculateDispersalMetrics(SubpopIN,xgridcopy,ygridcopy,\
		Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,\
		FDispDistCD,MDispDistCD,FDispDistCDstd,MDispDistCDstd,subpopmigration,\
		gen,Fthreshold,Mthreshold,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC)
	# Skip module, but store tracking numbers from this module and update 'EmiPop' key
	# --------------------------------------------------------------------------------
	else:
		# Patch-based with total
		SelectionDeaths.append(['NA']*(len(SubpopIN)+1))
		DisperseDeaths.append(['NA']*(len(SubpopIN)+1))
		Population.append(['NA']*(len(SubpopIN)+1))
		PackingDeaths.append(['NA']*(len(SubpopIN)+1))
		
		# Population value
		ProbSuccess.append('NA')
		AdultNoMg.append('NA')
		FDispDistCD.append('NA')
		MDispDistCD.append('NA')
		FDispDistCDstd.append('NA')
		MDispDistCDstd.append('NA')
		
		# Age values
		PopulationAge.append(['NA']*len(size_mean[0][0]))
		PackingDeathsAge.append(['NA']*len(size_mean[0][0]))
		
		# Update 'EmiPop'
		for isub in xrange(len(SubpopIN)):
			SubpopIN[isub]['EmiPop'] = SubpopIN[isub]['NatalPop']		
		
	# Return variables from this argument
	return SubpopIN
	
	# End::DoDisperse()