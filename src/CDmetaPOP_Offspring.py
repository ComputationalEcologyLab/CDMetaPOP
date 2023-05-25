# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Offspring.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for offspring processes.
# --------------------------------------------------------------------------------------------------

import pdb,sys, copy, numbers
from ast import literal_eval
from scipy.stats import truncnorm
from CDmetaPOP_Mortality import DoEggMortality
from CDmetaPOP_Modules import logMsg
import numpy as np

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringVars(Bearpairs,Femalepercent,sourcePop,size_mean,transmissionprob,gen,sizeans,age_mature,noOffspring,size_std,inheritans_classfiles,eggFreq,sexans,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo):
	'''
	DoOffspringVars()
	This function assigns the age (0), sex, and size of each offspring.
	'''	
	
	# Create empty variable for storing individual offspring information
	offspring=[]
	
	# Only if pairing occured
	#if len(Bearpairs) != 1 and Bearpairs[0][0] != -9999:
	#if len(Bearpairs[0][0]) != 1:
	if not (isinstance(Bearpairs[0][0],numbers.Integral)):
	#if Bearpairs[0][0] != -9999:
		
		# Error check
		if len(noOffspring) != len(Bearpairs):
			print('Offspring mismatch with Bearpairs.')
			sys.exit(-1)
		count = 0	
		# Loop through each mate pair
		# ----------------------------
		for i in range(len(Bearpairs)):
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
			for j in range(noOffspring[i]):
				
				# ------------------------
				# Mother's patch location
				# ------------------------
				patchindex = int(Bearpairs[i][0][sourcePop])-1
				
				# ------------------------
				# Get classfile assignment
				# ------------------------
				randno = np.random.uniform() # Random number
				if inheritans_classfiles == 'random':
					if randno < 0.5:
						natalP = fathers_natalP
						theseclasspars = fathers_theseclasspars
					else:
						natalP = mothers_natalP
						theseclasspars = mothers_theseclasspars	
				elif inheritans_classfiles == 'Hindex': #Hindex draw
					if randno <= offspring_hindex: # Use 1.0 files
						natalP = 0
						theseclasspars = 0
					else: # Use 0.0 files
						natalP = 0
						theseclasspars = 1
				elif inheritans_classfiles == 'mother': # Inherits mother's class pars
					natalP = mothers_natalP
					theseclasspars = mothers_theseclasspars	
				else:
					print('Error in inherit class vars file answer')
					sys.exit(-1)
				
				# --------------------------
				# Assign sex here
				# --------------------------								
				# Special case for WrightFisher
				if Femalepercent == 'WrightFisher':
					# Error with YYs here:
					if Bearpairs[i][1]['sex'] == 'MYY' or Bearpairs[i][0]['sex'] == 'FYY':
						print('Wright Fisher option specified for sex ratios. YY individuals should not be considered; use probability value for Femaleprob_Egg.')
						sys.exit(-1)
					offsex = int(2*np.random.uniform())
					if offsex == 0:
						offsex = 'FXX'
					else:
						offsex = 'MXY'
				# Else Prob given
				elif isinstance(float(Femalepercent),float):
					mothers_sex = Bearpairs[i][0]['sex']
					fathers_sex = Bearpairs[i][1]['sex'] 
					
					# FXX and MXY
					if mothers_sex == 'FXX' and fathers_sex == 'MXY':
						randsex = np.random.uniform()				
						# If that random number is less the Femalepercent/prob, assign it to be a female
						if randsex < float(Femalepercent):
							offsex = 'FXX'
						# If the random number is greater than the Femalepercent/prob, assign it to be a male
						else:
							offsex = 'MXY'
					# FXX and MYY
					elif mothers_sex == 'FXX' and fathers_sex == 'MYY':
						offsex = 'MXY'
					# FYY and MXY
					elif mothers_sex == 'FYY' and fathers_sex == 'MXY':
						randsex = np.random.uniform()				
						# If that random number is less the Femalepercent/prob, assign it to be a female
						if randsex < float(Femalepercent):
							offsex = 'MXY'
						# If the random number is greater than the Femalepercent/prob, assign it to be a male
						else:
							offsex = 'MYY'
					# FYY and MYY
					elif mothers_sex == 'FYY' and fathers_sex == 'MYY':
						offsex = 'MYY'
					elif sexans == 'H':
						if mothers_sex == 'FXX' and fathers_sex == 'FXX':
							offsex = 'FXX'
						else:
							randsex = np.random.uniform()				
							# If that random number is less the Femalepercent/prob, assign it to be a female
							if randsex < float(Femalepercent):
								offsex = 'FXX'
							# If the random number is greater than the Femalepercent/prob, assign it to be a male
							else:
								offsex = 'MXY'
						
					else:
						pdb.set_trace()
						print('Error in sex assignment.')
						sys.exit(-1)
					#print(fathers_sex)	
				# Error check
				else:
					print('Egg_Femaleprob is not correct data type.')
					sys.exit(-1)
				# Indexing used later for sex splits
				if offsex == 'FXX':
					sxspot = 0
				elif offsex == 'MXY':
					sxspot = 1
				elif offsex == 'MYY':
					sxspot = 2
				else:
					sxspot = 3
				
				# --------------------------
				# Assign infection here
				# --------------------------			
				# If parent has infection
				if Bearpairs[i][0]['infection'] == 1 or\
				Bearpairs[i][1]['infection'] == 1:			
					# Get a random number
					randinfection = np.random.uniform()				
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
				father_name = Bearpairs[i][1]['name']
				father_name = father_name.split('_')
				name = 'Age0'+'_F'+Bearpairs[i][0][sourcePop]+'_m'+Bearpairs[i][0][sourcePop]+'f'+Bearpairs[i][1][sourcePop]+'_P'+Bearpairs[i][0][sourcePop]+'_Y'+str(gen)+'_UO'+str(count)
				
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
					#lower, upper = 0,np.inf
					#sizesamp  = truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)	
					sizesamp = np.random.normal(mu,sigma)
					if sizesamp < 0:
						sizesamp = 0
				else:
					sizesamp = mu
				
				# --------------------------
				# Assign maturity here
				# --------------------------
				matval = 0.0 # Initialize
				agetemp = 0
				# Check default age/size for maturity
				if offsex == 'FXX':
					if FXXmat_set != 'N':
						if len(FXXmat_set.split('age')) == 2: # Maturation value set for age
							AgeMature = int(FXXmat_set.split('age')[1])
							if agetemp >= AgeMature: # If the age is > than default mature value, then becomes mature.
								matval = 1.0
							else: 
								matval = 0.0
						elif len(FXXmat_set.split('size')) == 2: # Maturation value set for size
							SizeMature = int(FXXmat_set.split('size')[1])
							if sizesamp >= SizeMature:
								matval = 1.0
							else: 
								matval = 0.0
						else:
							print('Female XX Maturation default set values age or size not specified in PopVars.')
							sys.exit(-1)
				elif offsex == 'MXY':
					if MXYmat_set != 'N':
						if len(MXYmat_set.split('age')) == 2: # Maturation value set for age
							AgeMature = int(MXYmat_set.split('age')[1])
							if agetemp >= AgeMature: # If the age is > than default mature value, then becomes mature.
								matval = 1.0
							else: 
								matval = 0.0
						elif len(MXYmat_set.split('size')) == 2: # Maturation value set for size
							SizeMature = int(MXYmat_set.split('size')[1])
							if sizesamp >= SizeMature:
								matval = 1.0
							else: 
								matval = 0.0
						else:
							print('Male XY Maturation default set values age or size not specified in PopVars.')
							sys.exit(-1)
				elif offsex == 'MYY':
					if MYYmat_set != 'N':
						if len(MYYmat_set.split('age')) == 2: # Maturation value set for age
							AgeMature = int(MYYmat_set.split('age')[1])
							if agetemp >= AgeMature: # If the age is > than default mature value, then becomes mature.
								matval = 1.0
							else: 
								matval = 0.0
						elif len(MYYmat_set.split('size')) == 2: # Maturation value set for size
							SizeMature = int(MYYmat_set.split('size')[1])
							if sizesamp >= SizeMature:
								matval = 1.0
							else: 
								matval = 0.0
						else:
							print('Male YY Maturation default set values age or size not specified in PopVars.')
							sys.exit(-1)
				elif offsex == 'FYY':
					if FYYmat_set != 'N':
						if len(FYYmat_set.split('age')) == 2: # Maturation value set for age
							AgeMature = int(FYYmat_set.split('age')[1])
							if agetemp >= AgeMature: # If the age is > than default mature value, then becomes mature.
								matval = 1.0
							else: 
								matval = 0.0
						elif len(FYYmat_set.split('size')) == 2: # Maturation value set for size
							SizeMature = int(FYYmat_set.split('size')[1])
							if sizesamp >= SizeMature:
								matval = 1.0
							else: 
								matval = 0.0
						else:
							print('Female YY Maturation default set values age or size not specified in PopVars.')
							sys.exit(-1)
					
				# If mat val is not 1, then run size/age probs
				if matval != 1.0:		
					# Check age values for mature
					if sizeans == 'age': # Age control
						if len(age_mature[natalP][theseclasspars][agetemp].split('~')) == 1:
							matval = float(age_mature[natalP][theseclasspars][agetemp].split('~')[0])
						elif len(age_mature[natalP][theseclasspars][agetemp].split('~')) != sexchromo:
							print('ClassVars age maturation probabilities must be length 1 or length of number of sex_chromo specified.')
							sys.exit(-1)							
						else:
							matval = float(age_mature[natalP][theseclasspars][agetemp].split('~')[sxspot])
						
					# If size control specified, grab slope/int values from PopVars	
					elif sizeans == 'size': # Size control							
						if offsex == 'FXX': # Female		
							matval = np.exp(float(FXXmat_int) + float(FXXmat_slope) * sizesamp) / (1 + np.exp(float(FXXmat_int) + float(FXXmat_slope) * sizesamp))
						elif offsex == 'MXY': # Male			
							matval = np.exp(float(MXYmat_int) + float(MXYmat_slope) * sizesamp) / (1 + np.exp(float(MXYmat_int) + float(MXYmat_slope) * sizesamp))
						elif offsex == 'MYY': # Male			
							matval = np.exp(float(MYYmat_int) + float(MYYmat_slope) * sizesamp) / (1 + np.exp(float(MYYmat_int) + float(MYYmat_slope) * sizesamp))
						elif offsex == 'FYY': # Male			
							matval = np.exp(float(FYYmat_int) + float(FYYmat_slope) * sizesamp) / (1 + np.exp(float(FYYmat_int) + float(FYYmat_slope) * sizesamp))	
					# Error check 	
					else:
						print('Size control option not correct, N or Y.')
						sys.exit(-1)
				
				# Check probability and egg laying
				randmat = np.random.uniform()
				if randmat < matval:
					mature = 1
					randegglay = np.random.uniform()
					# If sexans 'Y' and female, check layEggs
					if sexans == 'Y' or sexans == 'H':
						if offsex == 'FXX' or offsex == 'FYY':
							if randegglay < eggFreq:
								offlayeggs = 1
							else:
								offlayeggs = 0
						else:
							offlayeggs = 0
					else:				
						if randegglay < eggFreq:
							offlayeggs = 1
						else:
							offlayeggs = 0
				else:
					mature = 0
					offlayeggs = 0
								
				# --------------------------
				# REcord information
				# --------------------------			
				# And then recd new information of offspring [Mothergenes,Fathergenes,natalpop,emipop,immipop,emicd,immicd,age0,sex,size,mature,newmature,infection,id,motherid,fatherid,capture,recapture,layeggs,Mothers Hindex, Fathers Hindex, ClassVars File,PopID,speciesID]
				recd = (Bearpairs[i][0]['genes'],Bearpairs[i][1]['genes'],Bearpairs[i][0][sourcePop],'NA','NA',-9999,-9999,agetemp,offsex,sizesamp,mature,mature,infect,id,Bearpairs[i][0]['name'],Bearpairs[i][1]['name'],0,0,offlayeggs,Bearpairs[i][0]['hindex'],Bearpairs[i][1]['hindex'],'P'+str(natalP)+'_CV'+str(theseclasspars),Bearpairs[i][0]['popID'],Bearpairs[i][0]['species'])
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
	for i in range(len(Bearpairs)):			
			
		# If female did not mate up, then assign 0 offspring
		#if Bearpairs[i][1] == -9999:
		if isinstance(Bearpairs[i][1],numbers.Integral):
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
				litter_mu = float(age_mu[natalP][theseclasspars][ageF])
				
			if litter_mu <= 0.:				
				littersamp = 0
			else:
				# Set the litter size
				littersamp = int(round((litter_mu)*np.random.uniform()))
		
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
	for i in range(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		#if Bearpairs[i][1] == -9999:
		if isinstance(Bearpairs[i][1],numbers.Integral):
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
				litter_mu = float(age_mu[natalP][theseclasspars][ageF])
			
			if litter_mu <= 0.:				
				littersamp = 0
			else:
				# Set the litter size
				littersamp = int(round(np.random.poisson(litter_mu)))
				if littersamp < 0:
					pdb.set_trace()
	
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
	for i in range(len(Bearpairs)):
		
		# If female did not mate up, then assign 0 offspring
		#if Bearpairs[i][1] == -9999:
		if isinstance(Bearpairs[i][1],numbers.Integral):
			littersamp = 0
		
		# If females did mate up, then assign random drawn number
		else:
		
			# If size control then use parameters for length age_mu and CV
			if sizecall == 'size':
				if i == 0:
					stringout = 'Warning: size control specified with offspring number that does not have standard deviation, using sigma from Agevars file.'
					#logMsg(logfHndl,stringout)
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
				ageF = Bearpairs[i][0]['age']
				#litter_sigma = np.mean(np.asarray(age_sigma[natalP][theseclasspars],dtype=float))
				#pdb.set_trace()
				# Casey edit to account for when age is larger than number of age classes
				if ageF > len(age_sigma[natalP][theseclasspars])-1:
					litter_sigma = float(age_sigma[natalP][theseclasspars][len(age_sigma[natalP][theseclasspars])-1])
				else:
					litter_sigma = float(age_sigma[natalP][theseclasspars][ageF])
				if litter_mu <= 0.:				
					littersamp = 0
				else:
					# Set the litter size
					littersamp = round(np.random.normal(litter_mu,litter_sigma))
			else: # Use the AgeVars given mu and sigma
				# Grab the age or size of the female
				ageF = Bearpairs[i][0]['age']
				# Get the original ClassVars file location.
				natalP = int(Bearpairs[i][0]['classfile'].split('_')[0].split('P')[1])
				theseclasspars = int(Bearpairs[i][0]['classfile'].split('_')[1].split('CV')[1])
				# If age is greater than last age
				if ageF > len(age_mu[natalP][theseclasspars]) - 1:
					ageF = len(age_mu[natalP][theseclasspars]) - 1
				litter_mu = float(age_mu[natalP][theseclasspars][ageF])
				#litter_sigma = float(age_sigma[natalP][theseclasspars][ageF]) # Remove this as updated in DOStochastic
				littersamp = litter_mu
			
			if littersamp < 0:
				littersamp = 0
	
		# Append Offspring number to end of Pairs [F,M,#offspring]
		noOffspring.append(int(round(littersamp)))	
		
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
	for i in range(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		#if Bearpairs[i][1] == -9999:
		if isinstance(Bearpairs[i][1],numbers.Integral):
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
				litter_mu = float(age_mu[natalP][theseclasspars][ageF])
				
			if litter_mu <= 0.:				
				littersamp = 0
			else:
				# Set the litter size
				littersamp = int(round(litter_mu))
		
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
		for imo in range(len(unimo[0])):
			# if there is more than one mate pair
			if unimo[1][imo] != 1:
				duplicateLoc = np.where(mothers['name']==unimo[0][imo])[0]
				
				# Divide the births
				for ipairs in range(unimo[1][imo]):
					# Then divide the clutch size
					noOffspring[duplicateLoc[ipairs]] = noOffspring[duplicateLoc[ipairs]] / unimo[1][imo]
					# Case for when less births than pairings
		
	return noOffspring
	# End::DoClutch()
	
# ---------------------------------------------------------------------------------------------------	 
def DoOffspring(offno,Bearpairs,transmissionprob,gen,K,sourcePop,\
age_mu,age_sigma,sizeans,egg_mean_1,egg_mean_2,egg_mean_ans,equalClutch,dtype,eggmort_patch,EggDeaths,eggmort_back,egg_delay,noOffspring_temp,Births,BirthsMYY,BirthsFYY):
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
	
	# Check for Population extinct here and skip
	# ------------------------------------------
	if len(Bearpairs) != 0:
		
		# Only if pairings occurred
		#if Bearpairs[egg_delay][0][0] != -9999:
		if not isinstance(Bearpairs[egg_delay][0][0],numbers.Integral):
				
			# Function 1 is a uniform random draw between 0 and lmdba number	
			if (offno=='1'):
				
				noOffspring = DoOffspringRandom(Bearpairs[egg_delay],age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring)
			
			# Function 2 is a Poisson draw
			elif (offno=='2'):
			
				noOffspring = DoOffspringPoisson(Bearpairs[egg_delay],age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring)				
				
			# Function 3 is a constant of lmbda offspring per each pairing
			elif (offno=='3'):
			
				noOffspring = DoOffspringConstant(Bearpairs[egg_delay],age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring)
			
			# Function 4 is a normal draw
			elif (offno=='4'):
				noOffspring = DoOffspringNormal(Bearpairs[egg_delay],age_mu,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans,noOffspring)
			
			# Other functions
			else:
				print('This offspring birth rate option (offno) does not exist.')
				sys.exit(-1)
			
			# If equal clutch size is turned on
			if equalClutch == 'Y':		
				noOffspring = DoClutch(Bearpairs[egg_delay],dtype,noOffspring)
		
		# Make sure as array
		noOffspring = np.asarray(noOffspring)
		
		# Check if there were 0 litter size events, delete those Bearpairs[egg_delay]
		if len(np.where(noOffspring == 0)[0]) > 0:
			# Get index of 0 births for mothers
			ind0 = np.where(noOffspring == 0)[0]
			# Delete Bearpairs[egg_delay] and offspring no
			Bearpairs[egg_delay] = np.delete(Bearpairs[egg_delay],ind0,0)	
			noOffspring = np.delete(noOffspring,ind0)
			if len(Bearpairs[egg_delay]) == 0:
				Bearpairs[egg_delay] = [[-9999,-9999]]
		
		# -------------------------------------
		# Call DoEggMortality()
		# -------------------------------------	
		noOffspring = DoEggMortality(Bearpairs[egg_delay],eggmort_patch,EggDeaths,gen,K,eggmort_back,noOffspring,Births,BirthsMYY,BirthsFYY)
		
		# Check if there were 0 litter size events, delete those Bearpairs[egg_delay]
		if len(np.where(noOffspring == 0)[0]) > 0:
			# Get index of 0 births for mothers
			ind0 = np.where(noOffspring == 0)[0]
			# Delete Bearpairs[egg_delay] and offspring no
			Bearpairs[egg_delay] = np.delete(Bearpairs[egg_delay],ind0,0)	
			noOffspring = np.delete(noOffspring,ind0)
			if len(Bearpairs[egg_delay]) == 0:
				Bearpairs[egg_delay] = [[-9999,-9999]]
		
		# ---------------------------------------------------------------
		# Update for egg_delay; create n-D noOffspring array for indexing
		# ---------------------------------------------------------------
		noOffspring_temp[egg_delay] = noOffspring	
		
	# Population extinct
	else:
		noOffspring_temp = []
		Births.append([0 for x in range(0,len(K)+1)] )
		BirthsMYY.append([ 0 for x in range(0,len(K)+1)] )
		BirthsFYY.append([ 0 for x in range(0,len(K)+1)] )
		EggDeaths.append( [0 for x in range(0,len(K)+1)] )
	
	return noOffspring_temp, Bearpairs
	# End::DoOffspring()