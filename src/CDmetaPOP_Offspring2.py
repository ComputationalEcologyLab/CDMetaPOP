# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Offspring.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for offspring processes.
# --------------------------------------------------------------------------------------------------

import pdb,sys, copy, numbers
from ast import literal_eval
from scipy.stats import truncnorm
from CDmetaPOP_Mortality import DoEggMortality, DoIndividualEggMortality
from CDmetaPOP_Modules import *
import numpy as np

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringVars(Bearpairs,Femalepercent,sourcePop,size_mean,gen,sizeans,age_mature,noOffspring,size_std,inheritans_classfiles,eggFreq_mu,eggFreq_sd,sexans,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,egg_add,SubpopIN_keepAge1plus,PopTag,disease_vars,K):
	'''
	DoOffspringVars()
	This function assigns the age (0), sex, and size of each offspring.
	'''	
	#pdb.set_trace()	
	# Create empty variable for storing individual offspring information
	offspring=[]
	#pdb.set_trace() 
	# Need to check gen 1+ for 'NA; in EmiPop location - think only in gen 0
	# if nonmating option specified, then need to find where each Bearpair mother went or survived after Emigration
	# Also a pairing has to occur
	if egg_add == 'nonmating' and isinstance(Bearpairs[0][0],np.void):
		#Bearpairs[0][0]['name'] = 'this_Ptesting' # THis is for testing situation where mother did not make it through to SubpopIN
		Bearpairs[:,0]['EmiPop'] = 'NA' # For generation 1+ We can overwrite the EmiPop for the mothers
		# For each ID in Bearpair mothers
		mothers_Bearpairs_names = [s.split('_P')[1] for s in Bearpairs[:,0]['name'].tolist()]		
		#pdb.set_trace()		
		# Find where they are in SubpopIN_keepAge1plus
		#tempspot = []
		for isub in range(len(SubpopIN_keepAge1plus)):	
			FemalesInThisSubPOP = np.where(SubpopIN_keepAge1plus[isub]['sex'] == 'FXX')[0] # Indexes into SubpopIN
			thissubpop = [s.split('_P')[1] for s in SubpopIN_keepAge1plus[isub][FemalesInThisSubPOP]['name'].tolist()] # Female only names
			# Loop through the females names in this patch (may or maynot be a mother)
			for imove,itemmove in enumerate(thissubpop):
				# Else loop through the mothers names where mated
				for imate,itemmate in enumerate(mothers_Bearpairs_names):										
					if itemmate == itemmove: 
						#tempspot.append([isub,imate,itemmate,FemalesInThisSubPOP[imove],itemmove])
						# Update the 'EmiPop' location in Bearpairs
						Bearpairs[imate][0]['EmiPop'] = str(isub + 1)
						continue
	
	# -----------------------------------------------------------------------------
	# Only if pairing occured - Continue through initialization of individual loop
	if sum(noOffspring) > 0:
		
		# Error check
		if len(noOffspring) != len(Bearpairs):
			print('Offspring mismatch with Bearpairs.')
			sys.exit(-1)
		count = 0	
		# Loop through each mate pair
		# ----------------------------
		for i in range(len(Bearpairs)):			
			# Error check for eggdelay and sourcePop
			if Bearpairs[i][0][sourcePop] == 'NA':
				print('Issue with Egg Delay = 1.')
				sys.exit(-1)
			# The mother didn't make it so assume offspring didn't either - skip this pair and offspring
			if egg_add == 'nonmating' and Bearpairs[i][0]['EmiPop'] == 'NA':
				continue
			
			# ---------------------
			# Get Hindex
			# ---------------------
			offspring_hindex = Bearpairs[i][0]['hindex']/2. + Bearpairs[i][1]['hindex']/2.
			
			# --------------------------------------------------------
			# And then loop through each offspring from that mate pair
			# --------------------------------------------------------
			for j in range(noOffspring[i]):
				
				# ------------------------
				# Mother's patch location - note sourcePop updated if nonmating option
				# ------------------------
				if egg_add == 'mating':
					patchindex = int(Bearpairs[i][0][sourcePop])-1 # Assume gen = 0 is natalpop and gen>0 immipop
				else:
					patchindex = int(Bearpairs[i][0]['EmiPop'])-1					
				
				# ------------------------
				# Get classfile assignment
				# ------------------------
				# Extract and parse classfile details for mother and father
				mother_info = Bearpairs[i][0]['classfile'].split('_')
				father_info = Bearpairs[i][1]['classfile'].split('_')

				mothers_natalP = int(mother_info[0].split('P')[1])
				mothers_theseclasspars = int(mother_info[1].split('CV')[1])

				fathers_natalP = int(father_info[0].split('P')[1])
				fathers_theseclasspars = int(father_info[1].split('CV')[1])
				
				# Generate random number for inheritance logic
				randno = np.random.uniform()
				# Determine inheritance based on specified method
				if inheritans_classfiles == 'random':
					natalP, theseclasspars = (
						(fathers_natalP, fathers_theseclasspars)
						if randno < 0.5
						else (mothers_natalP, mothers_theseclasspars)
					)
				elif inheritans_classfiles == 'Hindex':  # Hindex logic
					natalP, theseclasspars = (0, 0) if randno <= offspring_hindex else (0, 1)
				elif inheritans_classfiles == 'mother':  # Inherits mother's class pars
					natalP, theseclasspars = mothers_natalP, mothers_theseclasspars
				else:
					print('Error in inherit class vars file answer')
					sys.exit(-1)

				# --------------------------
				# Assign sex here
				# --------------------------								
				# Special case for WrightFisher
				if Femalepercent == 'WrightFisher':
					# Error check for YY individuals
					if Bearpairs[i][1]['sex'] == 'MYY' or Bearpairs[i][0]['sex'] == 'FYY':
						print('Wright Fisher option specified for sex ratios. YY individuals should not be considered; use probability value for Femaleprob_Egg.')
						sys.exit(-1)
					# Assign sex based on random uniform distribution
					offsex = 'FXX' if int(2 * np.random.uniform()) == 0 else 'MXY'

				# Probabilistic sex assignment
				elif isinstance(float(Femalepercent), float):
					mothers_sex = Bearpairs[i][0]['sex']
					fathers_sex = Bearpairs[i][1]['sex']
					randsex = np.random.uniform()

					# Sex determination logic
					if mothers_sex == 'FXX' and fathers_sex == 'MXY':
						offsex = 'FXX' if randsex < float(Femalepercent) else 'MXY'
					elif mothers_sex == 'FXX' and fathers_sex == 'MYY':
						offsex = 'MXY'
					elif mothers_sex == 'FYY' and fathers_sex == 'MXY':
						offsex = 'MXY' if randsex < float(Femalepercent) else 'MYY'
					elif mothers_sex == 'FYY' and fathers_sex == 'MYY':
						offsex = 'MYY'
					elif sexans == 'H':
						offsex = 'FXX' if mothers_sex == 'FXX' and fathers_sex == 'FXX' else (
							'FXX' if randsex < float(Femalepercent) else 'MXY'
						)
					else:
						print('Error in sex assignment. Possible selfing on with sexual reproduction.')
						sys.exit(-1)

				# Error check for invalid data type
				else:
					print('Egg_Femaleprob is not correct data type.')
					sys.exit(-1)
				'''
				# ChatGPT generated section
				# Indexing used later for sex splits
				sex_map = {'FXX': 0, 'MXY': 1, 'MYY': 2}
				sxspot = sex_map.get(offsex, 3)
				'''
								
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
				if disease_vars['ImpDisease'] != 'N':
					if disease_vars['OffAns'][patchindex].lower() == 'susceptible':
						# Make all S state = 0
						initialOffState = 0
				else: # Else no disease wished.
					# Make all S state = 0
					initialOffState = 0
				
				#pdb.set_trace()
				
				# --------------------------
				# Assign ID here
				# --------------------------				
				mother_name = Bearpairs[i][0]['name'].split('_')
				father_name = Bearpairs[i][1]['name'].split('_')
				# Note: Age0 + FROM(location mother is laying pup/eggs) + mother(location of mother in mating location)father(location of father in mating location)+ PATCH(patch location born in) + YEAR(born) + Unique ID number
				
				id = 'Age0'+'_F'+str(patchindex+1)+'_m'+Bearpairs[i][0][sourcePop]+'f'+Bearpairs[i][1][sourcePop]+'_P'+str(patchindex+1)+'_Y'+str(gen)+'_UO'+str(count)
				
				# Unique ID is needed for sorting later, however, this unique name can get large, check length and reset.
				if len(id.split('_')[-1]) > 80:
					print('Too many offspring, recheck fecundity values.')
					sys.exit(-1)
								
				# --------------------------
				# Assign size here
				# --------------------------
				mu,sigma = size_mean[natalP][theseclasspars][0],size_std[natalP][theseclasspars][0]						
				if sigma != 0: # Case here for sigma == 0
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
					# Get Egg laying rate 
					tempEggFreq=[] # Temp list value to store egg lay events
					stochastic_update(eggFreq_mu,eggFreq_sd,tempEggFreq)
					tempEggFreq = tempEggFreq[0] # Note indexing into first spot since list created above
					# If sexans 'Y' and female, check layEggs
					if sexans == 'Y' or sexans == 'H':
						if offsex == 'FXX' or offsex == 'FYY':
							if tempEggFreq < 1: # If egg laying is less than 1 event per year
								randegglay = np.random.uniform()
								if randegglay < tempEggFreq:
									offlayeggs = 1
								else:
									offlayeggs = 0
							else: # egg laying is greater than 1 event per year
								offlayeggs = np.round(tempEggFreq,0)
						else:
							offlayeggs = 0
					else:
						if tempEggFreq < 1: # If egg laying is less than 1 event per year
							randegglay = np.random.uniform()
							if randegglay < eggFreq:
								offlayeggs = 1
							else:
								offlayeggs = 0
						else: # egg laying is greater than 1 event per year
							offlayeggs = np.round(tempEggFreq,0)
				else:
					mature = 0
					offlayeggs = 0					
					
				# --------------------------
				# REcord information
				# --------------------------			
				# And then recd new information of offspring [Mothergenes,Fathergenes,natalpop or mating pop,emipop or nonmating pop (where born in nonmating),immipop,emicd,immicd,age0,sex,size,mature,newmature,id,motherid,fatherid,capture,recapture,layeggs,Mothers Hindex, Fathers Hindex, ClassVars File,PopID,speciesID]
				if egg_add == 'mating':	
					recd = (Bearpairs[i][0]['genes'],Bearpairs[i][1]['genes'],str(patchindex+1),'NA','NA',-9999,-9999,agetemp,offsex,sizesamp,mature,mature,initialOffState,id,Bearpairs[i][0]['name'],Bearpairs[i][1]['name'],0,0,offlayeggs,Bearpairs[i][0]['hindex'],Bearpairs[i][1]['hindex'],'P'+str(natalP)+'_CV'+str(theseclasspars),PopTag[patchindex],Bearpairs[i][0]['species'])
				else:
					recd = (Bearpairs[i][0]['genes'],Bearpairs[i][1]['genes'],'NA',str(patchindex+1),'NA',-9999,-9999,agetemp,offsex,sizesamp,mature,mature,initialOffState,id,Bearpairs[i][0]['name'],Bearpairs[i][1]['name'],0,0,offlayeggs,Bearpairs[i][0]['hindex'],Bearpairs[i][1]['hindex'],'P'+str(natalP)+'_CV'+str(theseclasspars),PopTag[patchindex],Bearpairs[i][0]['species'])
				offspring.append(recd)
				count = count + 1 # For unique naming tracking				
	# If there was not a pairing
	else:
		offspring = [[]]*len(K)
	
	# Variables returned
	return offspring	
	# End::DoOffspringVars()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringRandom(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans):
	'''
	DoOffspringRandom()
	'''
		
	# If size control then use parameters for length age_mu and CV
	if sizecall == 'size':
		if egg_mean_ans == 'linear':
			litter_mu = egg_mean_1 + egg_mean_2 * Bearpairs[0]['size']
		elif egg_mean_ans == 'exp':
			litter_mu = egg_mean_1 * np.exp(egg_mean_2*Bearpairs[0]['size'])
		elif egg_mean_ans == 'pow':
			litter_mu = egg_mean_1 * Bearpairs[0]['size']**egg_mean_2				
		else:
			print('Egg mean answer not an option, enter exp or linear.')
			sys.exit(-1)
	else: # Use the AgeVars given mu and sigma
		# Grab the age or size of the female
		ageF = Bearpairs[0]['age']
		# Get the original ClassVars file location.
		natalP = int(Bearpairs[0]['classfile'].split('_')[0].split('P')[1])
		theseclasspars = int(Bearpairs[0]['classfile'].split('_')[1].split('CV')[1])
		# If age is greater than last age
		if ageF > len(age_mu[natalP][theseclasspars]) - 1:
			ageF = len(age_mu[natalP][theseclasspars]) - 1
		litter_mu = float(age_mu[natalP][theseclasspars][ageF])
		
	if litter_mu <= 0.:				
		littersamp = 0
	else:
		# Set the litter size
		littersamp = int(round((litter_mu)*np.random.uniform()))	
	
	# Variables returned
	return littersamp
	# End::DoOffspringRandom()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringPoisson(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans):
	'''
	DoOffspringPoisson()
	'''		
			
	# If size control then use parameters for length age_mu and CV
	if sizecall == 'size':
		if egg_mean_ans == 'linear':
			litter_mu = egg_mean_1 + egg_mean_2 * Bearpairs[0]['size']
		elif egg_mean_ans == 'exp':
			litter_mu = egg_mean_1 * np.exp(egg_mean_2*Bearpairs[0]['size'])
		elif egg_mean_ans == 'pow':
			litter_mu = egg_mean_1 * Bearpairs[0]['size']**egg_mean_2	
		else:
			print('Egg mean answer not an option, enter exp or linear.')
			sys.exit(-1)				
	else: # Use the AgeVars given mu and sigma
		# Grab the age of the female
		ageF = Bearpairs[0]['age']
		# Get the original ClassVars file location.
		natalP = int(Bearpairs[0]['classfile'].split('_')[0].split('P')[1])
		theseclasspars = int(Bearpairs[0]['classfile'].split('_')[1].split('CV')[1])
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
		
	# Variables returned
	return littersamp
	# End::DoOffspringPoisson()
	
# ---------------------------------------------------------------------------------------------------	
def DoOffspringNormal(Bearpairs,age_mu,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans):
	'''
	DoOffspringNormal()
	'''		
	
	# If size control then use parameters for length age_mu and CV
	if sizecall == 'size':
		if i == 0:
			stringout = 'Warning: size control specified with offspring number that does not have standard deviation, using sigma from Agevars file.'
			#logMsg(logfHndl,stringout)
		if egg_mean_ans == 'linear':
			litter_mu = egg_mean_1 + egg_mean_2 * Bearpairs[0]['size']
		elif egg_mean_ans == 'exp':
			litter_mu = egg_mean_1 * np.exp(egg_mean_2*Bearpairs[0]['size'])
		elif egg_mean_ans == 'pow':
			litter_mu = egg_mean_1 * Bearpairs[0]['size']**egg_mean_2	
		else:
			print('Egg mean answer not an option, enter exp or linear.')
			sys.exit(-1)
		
		# Get the original ClassVars file location.
		natalP = int(Bearpairs[0]['classfile'].split('_')[0].split('P')[1])
		theseclasspars = int(Bearpairs[0]['classfile'].split('_')[1].split('CV')[1])
		ageF = Bearpairs[0]['age']
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
		ageF = Bearpairs[0]['age']
		# Get the original ClassVars file location.
		natalP = int(Bearpairs[0]['classfile'].split('_')[0].split('P')[1])
		theseclasspars = int(Bearpairs[0]['classfile'].split('_')[1].split('CV')[1])
		# If age is greater than last age
		if ageF > len(age_mu[natalP][theseclasspars]) - 1:
			ageF = len(age_mu[natalP][theseclasspars]) - 1
		litter_mu = float(age_mu[natalP][theseclasspars][ageF])
		#litter_sigma = float(age_sigma[natalP][theseclasspars][ageF]) # Remove this as updated in DOStochastic
		littersamp = litter_mu
	
	if littersamp < 0:
		littersamp = 0
			
	# Variables returned
	return int(round(littersamp))
	# End::DoOffspringNormal()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringConstant(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans):
	'''
	DoOffspringConstant()
	'''	

	# If size control then use parameters for length age_mu and CV
	if sizecall == 'size':
		if egg_mean_ans == 'linear':
			litter_mu = egg_mean_1 + egg_mean_2 * Bearpairs[0]['size']
		elif egg_mean_ans == 'exp':
			litter_mu = egg_mean_1 * np.exp(egg_mean_2*Bearpairs[0]['size'])
		elif egg_mean_ans == 'pow':
			litter_mu = egg_mean_1 * Bearpairs[0]['size']**egg_mean_2	
		else:
			print('Egg mean answer not an option, enter exp or linear.')
			sys.exit(-1)
	else: # Use the AgeVars given mu and sigma
		# Grab the age or size of the female
		ageF = Bearpairs[0]['age']
		# Get the original ClassVars file location.
		natalP = int(Bearpairs[0]['classfile'].split('_')[0].split('P')[1])
		theseclasspars = int(Bearpairs[0]['classfile'].split('_')[1].split('CV')[1])
		# If age is greater than last age
		if ageF > len(age_mu[natalP][theseclasspars]) - 1:
			ageF = len(age_mu[natalP][theseclasspars]) - 1
		litter_mu = float(age_mu[natalP][theseclasspars][ageF])
		
	if litter_mu <= 0.:				
		littersamp = 0
	else:
		# Set the litter size
		littersamp = int(round(litter_mu))
				
	# Variables returned
	return littersamp
	# End::DoOffspringConstant()
	
# ---------------------------------------------------------------------------------------------------	
def DoOffspringClutch(Bearpairs,dtype,noOffspring):
	'''
	DoClutch()
	The assigns an equal clutch to each female that potentially mated more than once.
	'''	
	
	# Only if pairing occured
	if len(Bearpairs)!=0:
		# Get just mothers
		mothers = Bearpairs[:,0]
		mothers = np.asarray(mothers,dtype=dtype)
		unimo = count_unique(mothers['name']) # Grabs unique mothers name unimo[0] and the count of that name unimo[1]
		
		# For each unique mother, get sum of each instance of that unique mother / count of that unique mother
		for imo in range(len(unimo[0])):
			
			# if there is more than one mate pair
			if unimo[1][imo] != 1:
				# The index locations (maps to both noOffspring and mothers/Bearpairs) of this unique spot
				duplicateIndexLocations = np.where(mothers['name']==unimo[0][imo])[0]
				
				# Get the mean eggs this female will lay
				thisfemale_meaneggs = int(np.round(sum(noOffspring[duplicateIndexLocations]) / unimo[1][imo]))
								
				# Zero the noOffspring locations for this females
				noOffspring[duplicateIndexLocations] = 0
				
				# Draw from the mated pairs array, thisfemale_meaneggs times
				matedpair_IndexLocations = np.random.choice(duplicateIndexLocations, thisfemale_meaneggs,replace=True)
				
				# Then loop through the mated pair locations and add 1 to each noOffspring index spot
				for ipairs in range(len(matedpair_IndexLocations)):
					noOffspring[matedpair_IndexLocations[ipairs]] = noOffspring[matedpair_IndexLocations[ipairs]] + 1
						
	return noOffspring
	# End::DoClutch()
	
# ---------------------------------------------------------------------------------------------------	 
def DoOffspringNo(offno,Bearpairs,age_mu,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans):
	'''
	DoOffspring()
	Get number of Offspring for mated pair
	'''	
	
	# First, call offno, and assign number of eggs
	# --------------------------------------------
	if (offno=='1'): # Function 1 is a uniform random draw between 0 and lmdba number		
		noOffspring = DoOffspringRandom(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans)
	elif (offno=='2'): # Function 2 is a Poisson draw	
		noOffspring = DoOffspringPoisson(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans)				
	elif (offno=='3'): # Function 3 is a constant of lmbda offspring per each pairing	
		noOffspring = DoOffspringConstant(Bearpairs,age_mu,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans)
	elif (offno=='4'): # Function 4 is a normal draw
		noOffspring = DoOffspringNormal(Bearpairs,age_mu,age_sigma,sizecall,egg_mean_1,egg_mean_2,egg_mean_ans)
	else: # Other functions
		print('This offspring birth rate option (offno) does not exist.')
		sys.exit(-1)
	
	return noOffspring		
	
	# End::DoOffspring()
	
# ---------------------------------------------------------------------------------------------------
def AddAge0s(SubpopIN_keepAge1plus,K,SubpopIN_Age0,gen,Population,loci,muterate,mtdna,mutationans,dtype,geneswap,allelst,PopulationAge,sizecall,size_mean,cdevolveans,burningen_cdevolve,timecdevolve,fitvals,SelectionDeaths_Age0s,assortmateModel,patchvals,packans,noalleles,plasticans,sexans,eggFreq_mu,eggFreq_sd,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,egg_add,disease_vars,Track_DiseaseStates_AddAge0s,EHom=None):

	'''
	Add in the Age 0 population.
	'''
	#pdb.set_trace()
	classno = len(size_mean[0][0])
	# Storage to keep
	SubpopIN_keepK = []
	# Add spot for next generation
	Population.append([]) 
	PopulationAge.append([])
	PopulationAge[gen] = [[] for x in range(0,classno)]
	SelectionDeaths_Age0s.append([])
	Track_DiseaseStates_AddAge0s.append([])
	if egg_add == 'mating': # Born in mating grounds
		egg_add_call = 'NatalPop'
	else: # Born in nonmating grounds
		egg_add_call = 'EmiPop'	
				
	# Loop through each patch
	for isub in range(len(K)):
		
		# Get each SubpopIN pop as array and Age0s array
		SubpopIN_arr = np.array(SubpopIN_keepAge1plus[isub],dtype=dtype)
		# Make sure there are offspring to add here
		if len(SubpopIN_Age0[isub]) != 0:
			Age0Pop = SubpopIN_Age0[np.where(SubpopIN_Age0[egg_add_call] == str(isub+1))[0]]
		else:
			Age0Pop = SubpopIN_Age0[isub]
				
		# ----------------
		# InheritGenes()
		# ----------------		
		SubpopIN_Age0_temp = InheritGenes(gen,Age0Pop,loci,muterate,mtdna,mutationans,K,dtype,geneswap,allelst,assortmateModel,noalleles,plasticans,cdevolveans,egg_add_call,egg_add)
		
		# --------------------------------
		# Apply spatial selection to Age0s (this might not be the right order)
		# --------------------------------
		# 1-locus selection model
		if (cdevolveans in {'1', '1_mat', '1_G_ind', '1_G_link'}) and (gen >= burningen_cdevolve) and ('Eggs' in timecdevolve):			
			SubpopIN_Age0_keep = []
			for iind in range(len(SubpopIN_Age0_temp)):
				outpool = SubpopIN_Age0_temp[iind]
				# for option 3 in which has to be mature
				if cdevolveans == '1_mat' and outpool['mature'] == 0:
					differentialmortality = 0.0
				else:
					# Call 1-locus selection model
					differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],isub)
				# Then flip the coin to see if outpool survives its location
				randcheck = np.random.uniform()
				
				# If outpool did not survive: break from loop, move to next outpool
				if randcheck < differentialmortality:					
					continue
				else: # Record if survived
					SubpopIN_Age0_keep.append(outpool)
			# dtype here
			SubpopIN_Age0_keep = np.array(SubpopIN_Age0_keep,dtype=dtype)		
		
		# 2-locus model
		elif (cdevolveans in {'2','2_mat'}) and (gen >= burningen_cdevolve) and ('Eggs' in timecdevolve):
			SubpopIN_Age0_keep = []
			for iind in range(len(SubpopIN_Age0_temp)):
				outpool = SubpopIN_Age0_temp[iind]
				# for option 3 in which has to be mature
				if cdevolveans == '2_mat' and outpool['mature'] == 0:
					differentialmortality = 0.0
				else:				
					# Call 2-locus selection model
					differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],isub)			
				# Then flip the coin to see if outpool survives its location
				randcheck = np.random.uniform()				
				# If outpool did not survive: break from loop, move to next outpool
				if randcheck < differentialmortality:
					continue
				else: # Record if survived
					SubpopIN_Age0_keep.append(outpool)			
			# dtype here
			SubpopIN_Age0_keep = np.array(SubpopIN_Age0_keep,dtype=dtype)
		
		# Hindex cdevolveans
		elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and ('Eggs' in timecdevolve):
			SubpopIN_Age0_keep = []
			for iind in range(len(SubpopIN_Age0_temp)):
				outpool = SubpopIN_Age0_temp[iind]				
				
				# Call 2-locus selection model
				differentialmortality =	DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[isub])
							
				# Then flip the coin to see if outpool survives its location
				randcheck = np.random.uniform()				
				# If outpool did not survive: break from loop, move to next outpool
				if randcheck < differentialmortality:
					continue
				else: # Record if survived
					SubpopIN_Age0_keep.append(outpool)			
			# dtype here
			SubpopIN_Age0_keep = np.array(SubpopIN_Age0_keep,dtype=dtype)
					
		# CDEVOLVE - Inbreeding F
		elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and ('Eggs' in timecdevolve):			
			SubpopIN_Age0_keep = []
			for iind in range(len(SubpopIN_Age0_temp)):
				outpool = SubpopIN_Age0_temp[iind]
				
				# Call F selection model
				differentialmortality = DoFSelection(fitvals,outpool['genes'],isub,EHom,cdevolveans)
									
				# Then flip the coin to see if outpool survives its location
				randcheck = np.random.uniform()				
				# If outpool did not survive: break from loop, move to next outpool
				if randcheck < differentialmortality:
					continue
				else: # Record if survived
					SubpopIN_Age0_keep.append(outpool)			
			# dtype here
			SubpopIN_Age0_keep = np.array(SubpopIN_Age0_keep,dtype=dtype)
		
		# CDEVOLVE - Inbreeding * Outbreeding
		elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and ('Eggs' in timecdevolve):
			SubpopIN_Age0_keep = []
			for iind in range(len(SubpopIN_Age0_temp)):
				outpool = SubpopIN_Age0_temp[iind]			
				
				# Call Hindex selection model
				differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],isub,EHom,cdevolveans,outpool['hindex'],patchvals[isub])
										
				# Then flip the coin to see if outpool survives its location
				randcheck = np.random.uniform()
				# If outpool did not survive: break from loop, move to next outpool
				if randcheck < differentialmortality:
					continue
				else: # Record if survived
					SubpopIN_Age0_keep.append(outpool)			
			# dtype here
			SubpopIN_Age0_keep = np.array(SubpopIN_Age0_keep,dtype=dtype)
				
		# MLocus Selection 
		elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and ('Eggs' in timecdevolve):
			SubpopIN_Age0_keep = []
			for iind in range(len(SubpopIN_Age0_temp)):
				outpool = SubpopIN_Age0_temp[iind]				
				
				# Call MLocus selection model
				differentialmortality =	DoMLocusSelection(outpool['genes'],isub,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
							
				# Then flip the coin to see if outpool survives its location
				randcheck = np.random.uniform()				
				# If outpool did not survive: break from loop, move to next outpool
				if randcheck < differentialmortality:
					continue
				else: # Record if survived
					SubpopIN_Age0_keep.append(outpool)			
			# dtype here
			SubpopIN_Age0_keep = np.array(SubpopIN_Age0_keep,dtype=dtype)
		
		# Maturation values need to be updated here for cdevolveans M
		elif (cdevolveans in {'M', 'MG_ind','MG_link'}) and (gen >= burningen_cdevolve): # cdevolve answer mature			
			print('These cdevolveans debugging and checking still in progress: M, MG, etc.')
			sys.exit(-1)
						
		else:
			SubpopIN_Age0_keep = SubpopIN_Age0_temp
		
		# Append all information to temp SubpopKeep variable
		SubpopIN_keepK.append(np.concatenate([SubpopIN_arr,SubpopIN_Age0_keep]))
		
		# Tracking for disease
		indstates_inthispatch = SubpopIN_Age0_keep['states'] 
		updated_countstates = np.bincount(indstates_inthispatch, minlength=disease_vars['noStates'][isub])		
		Track_DiseaseStates_AddAge0s[gen].append([])
		Track_DiseaseStates_AddAge0s[gen][isub].extend(updated_countstates)
		
		# Store new Ns 
		Population[gen].append(len(SubpopIN_keepK[isub]))
		SelectionDeaths_Age0s[gen].append(len(SubpopIN_Age0_temp)-len(SubpopIN_Age0_keep))
				
		# Age tracking
		# Switch here for size or age control
		if sizecall == 'size' and packans.split('_')[0] != 'logistic': # Use min and max for tracking numbers.
			# Get the middles for finding closest values
			size_bin = size_mean[0][0]
			size_mean_middles = np.asarray(size_bin)[1:] - np.diff(np.asarray(size_bin).astype('f'))/2
			age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])			
		else:
			# Count up each uniages
			age_adjusted = SubpopIN_keepK[isub]['age']
		
		# Tracking age N
		for iage in range(len(PopulationAge[gen])):
			sizeindex = np.where(age_adjusted==iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))
			
		# Special case where age class is greater than lastage
		sizeindex = np.where(age_adjusted > iage)[0]
		PopulationAge[gen][iage].append(len(sizeindex))		
					
	# Add total to Ns
	Population[gen].insert(0,sum(Population[gen]))
	Track_DiseaseStates_AddAge0s[gen].insert(0,np.sum(np.asarray(Track_DiseaseStates_AddAge0s[gen]),axis=0).tolist())
	# Age tracking
	for iage in range(len(PopulationAge[gen])):
		PopulationAge[gen][iage] = sum(PopulationAge[gen][iage])	
		
	# add delete here
	return SubpopIN_keepK
	# End::AddAge0s
		