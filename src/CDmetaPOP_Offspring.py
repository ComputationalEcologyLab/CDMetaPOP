# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Offspring.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for offspring processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	from numpy.random import *
	import numpy as np
except ImportError:
	raise ImportError, "Numpy required."
import pdb,sys
from sets import Set
import random,copy
from ast import literal_eval
from scipy.stats import truncnorm
from CDmetaPOP_Mortality import DoEggMortality

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringVars(Bearpairs,Femalepercent,sourcePop,size_mean,transmissionprob,gen,sizecall,M_mature,F_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,noOffspring,size_std,inheritans_classfiles):
	'''
	DoOffspringVars()
	This function assigns the age (0), sex, and size of each offspring.
	'''	
	
	# Create empty variable for storing individual offspring information
	offspring=[]
		
	# Only if pairing occured
	if Bearpairs[0][0] != -9999:
		
		# Error check
		if len(noOffspring) != len(Bearpairs):
			print('Offspring mismatch with Bearpairs.')
			sys.exit(-1)
		count = 0	
		# Loop through each mate pair
		# ----------------------------
		for i in xrange(len(Bearpairs)):
			# -------------------------------------------------------
			# Get parent's information for multiple classfiles
			# if inheritans_classfiles is random use this information
			# -------------------------------------------------------
			mothers_file = Bearpairs[i][0]['classfile']
			mothers_natalP = int(mothers_file.split('_')[0].split('P')[1])
			mothers_theseclasspars = int(mothers_file.split('_')[1].split('CV')[1])
			
			fathers_file = Bearpairs[i][1]['classfile']
			fathers_natalP = int(fathers_file.split('_')[0].split('P')[1])
			fathers_theseclasspars = int(fathers_file.split('_')[1].split('CV')[1])
			
			# ---------------------
			# Get Hindex
			# ---------------------
			offspring_hindex = Bearpairs[i][0]['hindex']/2. + Bearpairs[i][1]['hindex']/2.
			
			# And then loop through each offspring from that mate pair
			# --------------------------------------------------------
			for j in xrange(noOffspring[i]):
				
				# ------------------------
				# Mother's patch location
				# ------------------------
				patchindex = int(Bearpairs[i][0][sourcePop])-1
				
				# ------------------------
				# Get classfile assignment
				# ------------------------
				randno = rand() # Random number
				if inheritans_classfiles == 'random':
					if randno < 0.5:
						natalP = fathers_natalP
						theseclasspars = fathers_theseclasspars
					else:
						natalP = mothers_natalP
						theseclasspars = mothers_theseclasspars	
				else: #Hindex draw
					if randno <= offspring_hindex: # Use 1.0 files
						natalP = 0
						theseclasspars = 0
					else: # Use 0.0 files
						natalP = 0
						theseclasspars = 1					
				
				# --------------------------
				# Assign sex here
				# --------------------------			
				# Case for Wright Fisher or not
				if Femalepercent != 'WrightFisher':
					# Select sex of the jth offspring - select a random number
					randsex = int(100*rand())				
					# If that random number is less the Femalepercent, assign it to be a female
					if randsex < int(Femalepercent):
						offsex = 0				
					# If the random number is greater than the Femalepercent, assign it to be a male
					else:
						offsex = 1
				# Special case for WrightFisher
				else:
					offsex = int(2*rand())
				
				# --------------------------
				# Assign infection here
				# --------------------------			
				# If parent has infection
				if Bearpairs[i][0]['infection'] == 1 or\
				Bearpairs[i][1]['infection'] == 1:			
					# Get a random number
					randinfection = rand()				
					# If random flip is less than transmissionprob
					if randinfection < transmissionprob:				
						# Then append infection status to offspring 
						infect = 1					
					# If offspring does not get infected
					else:
						infect = 0				
				# If offspring does not get infected.
				else:			
					# Then append infection status to offspring
					infect = 0
				
				# --------------------------
				# Assign ID here
				# --------------------------
				mother_name = Bearpairs[i][0]['name']
				mother_name = mother_name.split('_')
				name = 'Age0_'+'F'+Bearpairs[i][0][sourcePop]+'_P'+Bearpairs[i][0][sourcePop]+'_Y'+str(gen)+'_UO'+str(count)
				
				# Unique ID is needed for sorting later, however, this unique name can get large, check length and reset.
				check = name.split('_')[-1]
				if len(check) > 80:
					print('Too many offspring, recheck fecundity values.')
					sys.exit(-1)
				# Assign ID: 'Age0_Y{}_P{}_M{}_O{}
				id = name
				
				# --------------------------
				# Assign size here
				# --------------------------
				mu,sigma = size_mean[natalP][theseclasspars][0],size_std[natalP][theseclasspars][0]			
				# Case here for sigma == 0
				if sigma != 0:
					lower, upper = 0,np.inf
					sizesamp  = truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)			
				else:
					sizesamp = mu
					
				# --------------------------
				# Assign maturity here
				# --------------------------
				# Assign maturity: age or size switch, then male or female mature switch
				if sizecall == 'age':			
					if offsex == 0: # Female check 
						if Fmat_set == 'N': # Use prob value
							matval = F_mature[natalP][theseclasspars][0]
						else: # Use set age
							if int(Fmat_set) == 0: # Age of offspring is 0
								matval = 1.0
							else:
								matval = 0.0
					else: # Male check
						if Mmat_set == 'N': # Use prob value in agevars
							matval = M_mature[natalP][theseclasspars][0]
						else: # Use set age
							if int(Mmat_set) == 0: # Age of offspring is 0
								matval = 1.0
							else:
								matval = 0.0
				else: # If size control				
					if offsex == 0:	# Female check
						if Fmat_set == 'N': # Use equation - size
							matval = np.exp(Fmat_int + Fmat_slope * sizesamp) / (1 + np.exp(Fmat_int + Fmat_slope * sizesamp))
						else: # Use set size
							if sizesamp >= int(Fmat_set):
								matval = 1.0
							else:
								matval = 0.0
					else: # Male check
						if Mmat_set == 'N': # Use equation - size
							matval = np.exp(Mmat_int + Mmat_slope * sizesamp) / (1 + np.exp(Mmat_int + Mmat_slope * sizesamp))
						else: # Use set size
							if sizesamp >= int(Mmat_set):
								matval = 1.0
							else:
								matval = 0.0
							
				randmat = rand()
				if randmat < matval:
					mature = 1
				else:
					mature = 0
				
				# --------------------------
				# REcord information
				# --------------------------			
				# And then recd new information of offspring [Mothergenes,Fathergenes,natalpop,emipop,immipop,emicd,immicd,age0,sex,size,mature,newmature,infection,id,capture,recapture,layeggs,Mothers Hindex, Fathers Hindex, ClassVars File,PopID]
				recd = (Bearpairs[i][0]['genes'],Bearpairs[i][1]['genes'],Bearpairs[i][0][sourcePop],'NA','NA',-9999,-9999,0,offsex,sizesamp,mature,mature,infect,id,0,0,0,Bearpairs[i][0]['hindex'],Bearpairs[i][1]['hindex'],'P'+str(natalP)+'_CV'+str(theseclasspars),Bearpairs[i][0]['popID'])
				offspring.append(recd)
				count = count + 1 # For unique naming tracking				
	# If there was not a pairing
	else:
		offspring.append([])
	
	# Variables returned
	return offspring
	
	# End::DoOffspringVars()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringRandom(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring):
	'''
	DoOffspringRandom()
	This function chooses a random number of 
	offspring for a mated pair.
	'''	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):			
			
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			littersamp = 0
		
		# If females did mate up, then assign random drawn number
		else:
		
			# If size control then use parameters for length age_mu and CV
			if sizecall == 'size':
				if egg_mean_ans == 'linear':
					litter_mu = egg_mean_1 + egg_mean_2 * Bearpairs[i][0]['size']
				elif egg_mean_ans == 'exp':
					litter_mu = egg_mean_1 * np.exp(egg_mean_2*Bearpairs[i][0]['size'])
				elif egg_mean_ans == 'pow':
					litter_mu = egg_mean_1 * Bearpairs[i][0]['size']**egg_mean_2				
				else:
					print('Egg mean answer not an option, enter exp or linear.')
					sys.exit(-1)
			else: # Use the AgeVars given mu and sigma
				# Grab the age or size of the female
				ageF = Bearpairs[i][0]['age']
				# Get the original ClassVars file location.
				natalP = int(Bearpairs[i][0]['classfile'].split('_')[0].split('P')[1])
				theseclasspars = int(Bearpairs[i][0]['classfile'].split('_')[1].split('CV')[1])
				# If age is greater than last age
				if ageF > len(age_mu[natalP][theseclasspars]) - 1:
					ageF = len(age_mu[natalP][theseclasspars]) - 1
				litter_mu = age_mu[natalP][theseclasspars][ageF]
				
			if litter_mu <= 0.:				
				littersamp = 0
			else:
				# Set the litter size
				littersamp = int(int((litter_mu))*rand())
		
		# Append Offspring number to end of Pairs [F,M,#offspring]
		noOffspring.append(littersamp)	
	
	# Variables returned
	return noOffspring
	# End::DoOffspringRandom()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringPoisson(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring):
	'''
	DoOffspringPoisson()
	This function chooses a number of offspring 
	from a Poisson distribution for a mated pair.
	'''		
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			littersamp = 0
		
		# If females did mate up, then assign random drawn number
		else:
		
			
			# If size control then use parameters for length age_mu and CV
			if sizecall == 'size':
				if egg_mean_ans == 'linear':
					litter_mu = egg_mean_1 + egg_mean_2 * Bearpairs[i][0]['size']
				elif egg_mean_ans == 'exp':
					litter_mu = egg_mean_1 * np.exp(egg_mean_2*Bearpairs[i][0]['size'])
				elif egg_mean_ans == 'pow':
					litter_mu = egg_mean_1 * Bearpairs[i][0]['size']**egg_mean_2	
				else:
					print('Egg mean answer not an option, enter exp or linear.')
					sys.exit(-1)				
			else: # Use the AgeVars given mu and sigma
				# Grab the age of the female
				ageF = Bearpairs[i][0]['age']
				# Get the original ClassVars file location.
				natalP = int(Bearpairs[i][0]['classfile'].split('_')[0].split('P')[1])
				theseclasspars = int(Bearpairs[i][0]['classfile'].split('_')[1].split('CV')[1])
				# If age is greater than last age
				if ageF > len(age_mu[natalP][theseclasspars]) - 1:
					ageF = len(age_mu[natalP][theseclasspars]) - 1
				litter_mu = age_mu[natalP][theseclasspars][ageF]
			
			if litter_mu <= 0.:				
				littersamp = 0
			else:
				# Set the litter size
				littersamp = int(poisson(litter_mu))
			
		# Append Offspring number to end of Pairs [F,M,#offspring]	
		noOffspring.append(littersamp)			
		
	# Variables returned
	return noOffspring
	# End::DoOffspringPoisson()
	
# ---------------------------------------------------------------------------------------------------	
def DoOffspringNormal(Bearpairs,age_mu,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring):
	'''
	DoOffspringNormal()
	This function chooses a number of offspring 
	from a Normal distribution for a mated pair.
	'''		
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
		
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			littersamp = 0
		
		# If females did mate up, then assign random drawn number
		else:
		
			# If size control then use parameters for length age_mu and CV
			if sizecall == 'size':
				if i == 0:
					print('Warning: size control specified with offspring number that does not have standard deviation, using sigma from Agevars file.')
				if egg_mean_ans == 'linear':
					litter_mu = egg_mean_1 + egg_mean_2 * Bearpairs[i][0]['size']
				elif egg_mean_ans == 'exp':
					litter_mu = egg_mean_1 * np.exp(egg_mean_2*Bearpairs[i][0]['size'])
				elif egg_mean_ans == 'pow':
					litter_mu = egg_mean_1 * Bearpairs[i][0]['size']**egg_mean_2	
				else:
					print('Egg mean answer not an option, enter exp or linear.')
					sys.exit(-1)
				
				# Get the original ClassVars file location.
				natalP = int(Bearpairs[i][0]['classfile'].split('_')[0].split('P')[1])
				theseclasspars = int(Bearpairs[i][0]['classfile'].split('_')[1].split('CV')[1])
				litter_sigma = np.mean(age_sigma[natalP][theseclasspars])
			else: # Use the AgeVars given mu and sigma
				# Grab the age or size of the female
				ageF = Bearpairs[i][0]['age']
				# Get the original ClassVars file location.
				natalP = int(Bearpairs[i][0]['classfile'].split('_')[0].split('P')[1])
				theseclasspars = int(Bearpairs[i][0]['classfile'].split('_')[1].split('CV')[1])
				# If age is greater than last age
				if ageF > len(age_mu[natalP][theseclasspars]) - 1:
					ageF = len(age_mu[natalP][theseclasspars]) - 1
				litter_mu = age_mu[natalP][theseclasspars][ageF]
				litter_sigma = age_sigma[natalP][theseclasspars][ageF]
				
			if litter_mu <= 0.:				
				littersamp = 0
			else:
				# Set the litter size
				littersamp = int(np.random.normal(litter_mu,litter_sigma))
			if littersamp < 0:
				littersamp = 0
	
		# Append Offspring number to end of Pairs [F,M,#offspring]
		noOffspring.append(littersamp)	
		
	# Variables returned
	return noOffspring
	# End::DoOffspringNormal()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringConstant(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring):
	'''
	DoOffspringConstant()
	This function chooses a constant number of 
	offspring for each mated pair.
	'''	
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			littersamp = 0
		
		# If females did mate up, then assign random drawn number
		else:
			
			# If size control then use parameters for length age_mu and CV
			if sizecall == 'size':
				if egg_mean_ans == 'linear':
					litter_mu = egg_mean_1 + egg_mean_2 * Bearpairs[i][0]['size']
				elif egg_mean_ans == 'exp':
					litter_mu = egg_mean_1 * np.exp(egg_mean_2*Bearpairs[i][0]['size'])
				elif egg_mean_ans == 'pow':
					litter_mu = egg_mean_1 * Bearpairs[i][0]['size']**egg_mean_2	
				else:
					print('Egg mean answer not an option, enter exp or linear.')
					sys.exit(-1)
			else: # Use the AgeVars given mu and sigma
				# Grab the age or size of the female
				ageF = Bearpairs[i][0]['age']
				# Get the original ClassVars file location.
				natalP = int(Bearpairs[i][0]['classfile'].split('_')[0].split('P')[1])
				theseclasspars = int(Bearpairs[i][0]['classfile'].split('_')[1].split('CV')[1])
				# If age is greater than last age
				if ageF > len(age_mu[natalP][theseclasspars]) - 1:
					ageF = len(age_mu[natalP][theseclasspars]) - 1
				litter_mu = age_mu[natalP][theseclasspars][ageF]
				
			if litter_mu <= 0.:				
				littersamp = 0
			else:
				# Set the litter size
				littersamp = int(litter_mu)
	
		# Append Offspring number to end of Pairs [F,M,#offspring]
		noOffspring.append(littersamp)	
			
	# Variables returned
	return noOffspring
	# End::DoOffspringConstant()
	
# ---------------------------------------------------------------------------------------------------	
def DoClutch(Bearpairs,dtype,noOffspring):
	'''
	DoClutch()
	The assigns an equal clutch to each female that potentially mated more than once.
	'''	
	
	# Only if pairing occured
	if len(Bearpairs)!=0:
		# Get just mothers
		mothers = Bearpairs[:,0]
		mothers = np.asarray(mothers,dtype=dtype)
		unimo = count_unique(mothers['name'])
		
		noOffspring = np.asarray(noOffspring)
		
		# Loop over unique mothers
		for imo in xrange(len(unimo[0])):
			# if there is more than one mate pair
			if unimo[1][imo] != 1:
				duplicateLoc = np.where(mothers['name']==unimo[0][imo])[0]
				
				# Divide the births
				for ipairs in xrange(unimo[1][imo]):
					# Then divide the clutch size
					noOffspring[duplicateLoc[ipairs]] = noOffspring[duplicateLoc[ipairs]] / unimo[1][imo]
					# Case for when less births than pairings
		
	return noOffspring
	# End::DoClutch()
	
# ---------------------------------------------------------------------------------------------------	 
def DoOffspring(offno,Bearpairs,Births,transmissionprob,gen,K,sourcePop,\
age_mu,age_sigma,sizeans,egg_mean_1,egg_mean_2,egg_mean_ans,equalClutch,dtype,Mmat_set,Fmat_set,eggmort_patch,EggDeaths,eggmort_back):
	'''
	DoOffspring()
	Choose number of Offspring for each mated pair 
	Input: selection choice for offspring number distribution draw
	offno, Bearpairs, lmbda.
	Age specific recruitment option here (age_mu)
	Output: Bear Pairs + # of offspring [Female,Male,#Offspring]
	offspring = [Femalegenes,Malegenes,age,sex,infection,name{AGE0}]
	'''	
	
	# Get unique number of subpops
	nosubpops = len(K)
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
	
	# Only if pairings occurred
	if Bearpairs[0][0] != -9999:
			
		# Function 1 is a uniform random draw between 0 and lmdba number	
		if (offno=='1'):
			
			noOffspring = DoOffspringRandom(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring)
		
		# Function 2 is a Poisson draw
		elif (offno=='2'):
		
			noOffspring = DoOffspringPoisson(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring)
			
		# Function 3 is a constant of lmbda offspring per each pairing
		elif (offno=='3'):
		
			noOffspring = DoOffspringConstant(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring)
		
		# Function 4 is a normal draw
		elif (offno=='4'):
			noOffspring = DoOffspringNormal(Bearpairs,age_mu,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring)
		
		# Other functions
		else:
			print('This offspring birth rate option (offno) does not exist.')
			sys.exit(-1)
		
		# If equal clutch size is turned on
		if equalClutch == 'Y':		
			noOffspring = DoClutch(Bearpairs,dtype,noOffspring)
	
	# Make sure as array
	noOffspring = np.asarray(noOffspring)
	
	# Check if there were 0 litter size events, delete those Bearpairs
	if len(np.where(noOffspring == 0)[0]) > 0:
		# Get index of 0 births for mothers
		ind0 = np.where(noOffspring == 0)[0]
		# Delete bearpairs and offspring no
		Bearpairs = np.delete(Bearpairs,ind0,0)	
		noOffspring = np.delete(noOffspring,ind0)
		if len(Bearpairs) == 0:
			Bearpairs = [[-9999,-9999]]
	
	# -------------------------------------
	# Call DoEggMortality()
	# -------------------------------------	
	
	noOffspring = DoEggMortality(Bearpairs,eggmort_patch,EggDeaths,gen,K,eggmort_back,noOffspring,Births)
	
	# Check if there were 0 litter size events, delete those Bearpairs
	if len(np.where(noOffspring == 0)[0]) > 0:
		# Get index of 0 births for mothers
		ind0 = np.where(noOffspring == 0)[0]
		# Delete bearpairs and offspring no
		Bearpairs = np.delete(Bearpairs,ind0,0)	
		noOffspring = np.delete(noOffspring,ind0)
		if len(Bearpairs) == 0:
			Bearpairs = [[-9999,-9999]]
	
	# Return variables from this argument
	return noOffspring,Bearpairs
		
	# End::DoOffspring()