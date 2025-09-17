# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Immigration.py
# Author: Erin L Landguth
# Created: January 2013
# Description: This is the function/module file for immigration processes.
# --------------------------------------------------------------------------------------------------

# Python specific functions
import pdb, copy, os, sys, multiprocessing
from ast import literal_eval 
from CDmetaPOP_Modules import *
import numpy as np 

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
		print(("%s"%(msg)))
		
	# End::logMsg()

# ---------------------------------------------------------------------------------------------------	 
def w_choice_general(lst):
	'''
	w_choice_general()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(x[1] for x in lst)
	n=np.random.uniform(0,wtotal)
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
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()
	
# ---------------------------------------------------------------------------------------------------	
def GetProbArray(offspring,answer,K,natal_patches,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXX,cdmatrix_MXY,cdmatrix_MYY,cdmatrix_FYY,fitvals=None):
	'''
	GetProbArray()
	This function gets indices for F and M specific cdmatrix values
	In direction cost, this is the row value; xycdmatrix[0] grabs column in original cdmatrix; xycdmatrix[:,0] grabs row vals in original cdmatrix or use xycdmatrix[currentpop][natalpop]
	'''
	#if answer == 'immigrator':
	#	pdb.set_trace()
	# Index into offspring array
	currentoff = offspring
	
	# Subpopulation current offspring is in - index - 1
	currentsubpop = int(currentoff['EmiPop']) - 1
	natalsubpop = int(currentoff['NatalPop']) - 1
	indSex = currentoff['sex']	

	# If straying process - use full cdmatrix
	if answer == 'strayer_emiPop':
		# Append the freegrid probabilities for the offspring choice
		if indSex == 'FXX': # Female offspring
			probarray = copy.deepcopy(cdmatrix_FXX[currentsubpop])
		elif indSex == 'MXY': # Male offspring
			probarray = copy.deepcopy(cdmatrix_MXY[currentsubpop])
		elif indSex == 'MYY': # Male YY
			probarray = copy.deepcopy(cdmatrix_MYY[currentsubpop])
		elif indSex == 'FYY': # FeMale YY
			probarray = copy.deepcopy(cdmatrix_FYY[currentsubpop])
		else:
			print('Invalid offspring list.')
			sys.exit(-1)
		
		# Where natal grounds = 0, set prob to 0
		probarray[np.where(np.asarray(natal_patches)==0)[0]] = 0.
		# Where K = 0, set prob to 0
		probarray[np.where(np.asarray(K)==0)[0]] = 0.
	
	# If straying process - use full cdmatrix
	elif answer == 'strayer_natalPop':
		
		# Append the freegrid probabilities for the offspring choice
		if indSex == 'FXX': # Female offspring
			probarray = copy.deepcopy(cdmatrix_FXX[natalsubpop])
		elif indSex == 'MXY': # Male offspring
			probarray = copy.deepcopy(cdmatrix_MXY[natalsubpop])
		elif indSex == 'MYY': # Male YY
			probarray = copy.deepcopy(cdmatrix_MYY[natalsubpop])
		elif indSex == 'FYY': # FeMale YY
			probarray = copy.deepcopy(cdmatrix_FYY[natalsubpop])
		else:
			print('Invalid offspring list.')
			sys.exit(-1)		
				
		# Where natal grounds = 0, set prob to 0
		probarray[np.where(np.asarray(natal_patches)==0)[0]] = 0.
		# Where K = 0, set prob to 0
		probarray[np.where(np.asarray(K)==0)[0]] = 0.
		
		# Where currentpop to 0 probabilities set to zero
		if indSex == 'FXX': # Female 
			nothere = copy.deepcopy(cdmatrix_FXX[currentsubpop])
		elif indSex == 'MXY': # Male 
			nothere = copy.deepcopy(cdmatrix_MXY[currentsubpop])
		elif indSex == 'MYY': # Male 
			nothere = copy.deepcopy(cdmatrix_MYY[currentsubpop])
		elif indSex == 'FYY': # FeMale YY
			nothere = copy.deepcopy(cdmatrix_FYY[currentsubpop])
		probarray[np.where(nothere == 0.)[0]] = 0.
		
	# If local dispersal process
	elif answer == 'localDispersal':
		
		# Append the freegrid probabilities for the offspring choice
		if indSex == 'FXX': # Female offspring
			probarray = copy.deepcopy(cdmatrix_FXX[natalsubpop])
		elif indSex == 'MXY': # Male offspring
			probarray = copy.deepcopy(cdmatrix_MXY[natalsubpop])
		elif indSex == 'MYY': # Male YY
			probarray = copy.deepcopy(cdmatrix_MYY[natalsubpop])
		elif indSex == 'FYY': # FeMale YY
			probarray = copy.deepcopy(cdmatrix_FYY[natalsubpop])
		else:
			print('Invalid offspring list.')
			sys.exit(-1)
		
		# Where natal grounds = 0, set prob to 0
		probarray[np.where(np.asarray(natal_patches)==0)[0]] = 0.
		# Where K = 0, set prob to 0
		probarray[np.where(np.asarray(K)==0)[0]] = 0.
		
	# Return home attempt
	elif answer == 'immigrator':
		checkback = True
		# Check for runtiming option
		if cdevolveans == 'runtiming':
			# Get this individuals genes
			indGenes = currentoff['genes']
			if indGenes[0] == 2:
				fitvalspot = 0				
			elif indGenes[0] == 1:
				fitvalspot = 1			
			elif indGenes[0] == 0:
				fitvalspot = 2
			
			# Get the probability this ind will move to based on genotype
			natalGeneProb = float(fitvals[natalsubpop][fitvalspot].split(';')[0])
			randback = np.random.uniform()
			if randback < natalGeneProb: # Gene matches natal patch location to go back to
				checkback = True
			else: # Gene does not match natal patch location, check for another spot to go to
				checkback = False
				# Grab the patch-based values that match this inds Gene (given in 'fitvals')
				temptheseGenePatch =  np.asarray(fitvals)[:,fitvalspot]
				theseGenePatch = np.asarray([float(entry.split(';')[0]) for entry in temptheseGenePatch])
				
				# Append the freegrid probabilities for the offspring choice
				if indSex == 'FXX': # Female offspring
					probarray = copy.deepcopy(cdmatrix_FXX[currentsubpop])
				elif indSex == 'MXY': # Male offspring
					probarray = copy.deepcopy(cdmatrix_MXY[currentsubpop])
				elif indSex == 'MYY': # Male YY
					probarray = copy.deepcopy(cdmatrix_MYY[currentsubpop])
				elif indSex == 'FYY': # FeMale YY
					probarray = copy.deepcopy(cdmatrix_FYY[currentsubpop])
				else:
					print('Invalid offspring list.')
					sys.exit(-1)
				
				# Update probability matrix based on the patch values gene match
				probarray = probarray * theseGenePatch				
				# Where natal grounds = 0, set prob to 0
				probarray[np.where(np.asarray(natal_patches)==0)[0]] = 0.
				# Where K = 0, set prob to 0
				probarray[np.where(np.asarray(K)==0)[0]] = 0.
				
		if checkback == True:	# If not runtiming or special case for runtiming	
			# If K is 0, then can not return home # Added natal_patches[natalsubpop] check on 2 Jan 2024
			if K[natalsubpop] == 0 or natal_patches[natalsubpop] == 0:
				probarray = [0.0]
			else:
				# Append the freegrid probabilities for the offspring choice
				if indSex == 'FXX': # Female 
					# Get probarray from current to natal - only one number!
					probarray_BA = [copy.deepcopy(cdmatrix_FXX[currentsubpop][natalsubpop])]
					probarray_AB = [copy.deepcopy(cdmatrix_FXX[natalsubpop][currentsubpop])]
				elif indSex == 'MXY': # Male 
					# Get probarray from current to natal - only one number!
					probarray_BA = [copy.deepcopy(cdmatrix_MXY[currentsubpop][natalsubpop])]
					probarray_AB = [copy.deepcopy(cdmatrix_MXY[natalsubpop][currentsubpop])]
				elif indSex == 'MYY': # Male 
					# Get probarray from current to natal - only one number!
					probarray_BA = [copy.deepcopy(cdmatrix_MYY[currentsubpop][natalsubpop])]
					probarray_AB = [copy.deepcopy(cdmatrix_MYY[natalsubpop][currentsubpop])]
				elif indSex == 'FYY': # FeMale YY 
					# Get probarray from current to natal - only one number!
					probarray_BA = [copy.deepcopy(cdmatrix_FYY[currentsubpop][natalsubpop])]
					probarray_AB = [copy.deepcopy(cdmatrix_FYY[natalsubpop][currentsubpop])]			
				else:
					print('Invalid offspring list.')
					sys.exit(-1)
				# Check asymmetrical cost - get relative difference for probability
				if probarray_BA[0] == 0.0: # Absolutely can not make it back, e.g., complete barrier
					probarray = [0.0]
				if probarray_BA[0] < probarray_AB[0]: # Harder to make it back - higher cost lower prob, e.g., partial barrier
					probarray = [1. - ((probarray_AB[0] - probarray_BA[0]) / probarray_AB[0])]
					# Here check if makes it back
					randback = np.random.uniform()
					if randback >= probarray[0]: # Does not make it back 
						probarray[0] = 0.0
					else: # Does make it back
						probarray[0] = 1.0 
				else: # Same probability/cost to make it back, then assume it migrates back
					probarray = [1.0]
	else:
		print('wrong answer in GetProbArray.')
		sys.exit(-1)
	
# Check plastic response here for temperature response and dominant "allele"
	if (plasticans != 'N') and (gen >= burningen_plastic) and (timeplastic.find('Out') != -1) and (plasticans.split('_')[0] == 'Temp') and (plasticans.split('_')[1] == "dom"):

		# Get location in genes array for plastic region
		# ----------------------------------------------
		Indgenes = currentoff['genes']
		# If cdevolve is on
		if cdevolveans != 'N':
			# Then the first l loci are for selection, next for plastic region
			if cdevolveans.split('_')[0] == 'P': # This is for multilocus selection, not currently implemented, to be moved over from cdpop
				selloci = int(cdevolveans.split('_')[2].split('L')[1])
			elif cdevolveans == '1' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans == 'stray' or cdevolveans == 'Hindex':
				selloci = 1
			elif cdevolveans == '2' or cdevolveans == 'MG' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == '2_mat':
				selloci = 2
			else:
				print('CDEVOLVEANS not entered correctly; DoUpdate() error.')
				sys.exit(-1)
		# If selection is not on
		else:
			selloci = 0 # zero loci in selection
		# Get number of plastic loci
		plaloci = 1
		# Get index for plastic region
		plaloci_index = range(selloci*2,selloci*2+plaloci*2)
		
		# Get the plastic behaviorial response threshold
		# ----------------------------------------------
		
		# Does this individual have the allele for response (2,2) or (2,0) or (0,2)
		if Indgenes[plaloci_index][0] == 2 or Indgenes[plaloci_index][1] == 2:
					
			if answer == 'immigrator': # probarray is length 1
				# Only check patchvals in natalsubpop, which is where this individual is going if probarray is 1
				if patchvals[natalsubpop] >= plastic_behaviorresp:
					probarray = [0.0]
				
			else: # probarray is length patch dim
				# Whereever temperature threshold, turn prob to 0
				probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2]))
	
    	# Check plastic response here for temperature response and recessive "allele"
	if (plasticans != 'N') and (gen >= burningen_plastic) and (timeplastic.find('Out') != -1) and (plasticans.split('_')[0] == 'Temp') and (plasticans.split('_')[1] == "rec"):

		# Get location in genes array for plastic region
		# ----------------------------------------------
		Indgenes = currentoff['genes']
		# If cdevolve is on
		if cdevolveans != 'N':
			# Then the first l loci are for selection, next for plastic region
			if cdevolveans.split('_')[0] == 'P': # This is for multilocus selection, not currently implemented, to be moved over from cdpop
				selloci = int(cdevolveans.split('_')[2].split('L')[1])
			elif cdevolveans == '1' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans == 'stray' or cdevolveans == 'Hindex':
				selloci = 1
			elif cdevolveans == '2' or cdevolveans == 'MG' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == '2_mat':
				selloci = 2
			else:
				print('CDEVOLVEANS not entered correctly; DoUpdate() error.')
				sys.exit(-1)
		# If selection is not on
		else:
			selloci = 0 # zero loci in selection
		# Get number of plastic loci
		plaloci = 1
		# Get index for plastic region
		plaloci_index = range(selloci*2,selloci*2+plaloci*2)
		
		# Get the plastic behaviorial response threshold
		# ----------------------------------------------
		
		# Does this individual have the allele for response (2,2) or (2,0) or (0,2)
		if Indgenes[plaloci_index][0] == 2 and Indgenes[plaloci_index][1] == 2:
					
			if answer == 'immigrator': # probarray is length 1
				# Only check patchvals in natalsubpop, which is where this individual is going if probarray is 1
				if patchvals[natalsubpop] >= plastic_behaviorresp:
					probarray = [0.0]
				
			else: # probarray is length patch dim
				# Whereever temperature threshold, turn prob to 0
				probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2]))
    
        	# Check plastic response here for temperature response and intermediate response
	if (plasticans != 'N') and (gen >= burningen_plastic) and (timeplastic.find('Out') != -1) and (plasticans.split('_')[0] == 'Temp') and (plasticans.split('_')[1] == "codom"):

		# Get location in genes array for plastic region
		# ----------------------------------------------
		Indgenes = currentoff['genes']
		# If cdevolve is on
		if cdevolveans != 'N':
			# Then the first l loci are for selection, next for plastic region
			if cdevolveans.split('_')[0] == 'P': # This is for multilocus selection, not currently implemented, to be moved over from cdpop
				selloci = int(cdevolveans.split('_')[2].split('L')[1])
			elif cdevolveans == '1' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans == 'stray' or cdevolveans == 'Hindex':
				selloci = 1
			elif cdevolveans == '2' or cdevolveans == 'MG' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == '2_mat':
				selloci = 2
			else:
				print('CDEVOLVEANS not entered correctly; DoUpdate() error.')
				sys.exit(-1)
		# If selection is not on
		else:
			selloci = 0 # zero loci in selection
		# Get number of plastic loci
		plaloci = 1
		# Get index for plastic region
		plaloci_index = range(selloci*2,selloci*2+plaloci*2)
		
		# Get the plastic behaviorial response threshold
		# ----------------------------------------------
		
		# Does this individual have the allele to intiate response
		# 2 in both = 0, one 2 = 50% reduction
		if Indgenes[plaloci_index][0] == 2 and Indgenes[plaloci_index][1] == 2:
		
			# Whereever temperature threshold, turn prob to value input as third part of plasticans
			probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2]))
            
		if Indgenes[plaloci_index][0] == 1 and Indgenes[plaloci_index][1] == 2:
		
			# Whereever temperature threshold, reduce probably y 50% of the value input as the third part of plasticans
			probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = ( (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2])) + probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)] ) /2
            
		if Indgenes[plaloci_index][0] == 2 and Indgenes[plaloci_index][1] == 1:
		
			# Whereever temperature threshold, reduce probably y 50% of the value input as the third part of plasticans
			probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = ( (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2])) + probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)] ) /2
    
    
	if (plasticans != 'N') and (gen >= burningen_plastic) and (timeplastic.find('Back') != -1) and (plasticans.split('_')[0] == 'Hab') and (plasticans.split('_')[1] == "dom"):
		
		# Get location in genes array for plastic region
		# ----------------------------------------------
		Indgenes = currentoff['genes']
		# If cdevolve is on
		if cdevolveans != 'N':
			# Then the first l loci are for selection, next for plastic region
			if cdevolveans.split('_')[0] == 'P': # This is for multilocus selection, not currently implemented, to be moved over from cdpop
				selloci = int(cdevolveans.split('_')[2].split('L')[1])
			elif cdevolveans == '1' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans == 'stray' or cdevolveans == 'Hindex':
				selloci = 1
			elif cdevolveans == '2' or cdevolveans == 'MG' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == '2_mat':
				selloci = 2
			else:
				print('CDEVOLVEANS not entered correctly; DoUpdate() error.')
				sys.exit(-1)
		# If selection is not on
		else:
			selloci = 0 # zero loci in selection
		# Get number of plastic loci
		plaloci = 1
		# Get index for plastic region
		plaloci_index = range(selloci*2,selloci*2+plaloci*2)
		
		# Get the plastic behaviorial response threshold
		# ----------------------------------------------
		
		# Does this individual have the allele for response (2,2) or (2,0) or (0,2)
		if Indgenes[plaloci_index][0] == 2 or Indgenes[plaloci_index][1] == 2:
					
			if answer == 'immigrator': # probarray is length 1
				# Only check patchvals in natalsubpop, which is where this individual is going if probarray is 1
				if patchvals[natalsubpop] >= plastic_behaviorresp:
					probarray = [0.0]
				
			else: # probarray is length patch dim
				# Whereever temperature threshold, turn prob to 0
				probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2]))
                
	if (plasticans != 'N') and (gen >= burningen_plastic) and (timeplastic.find('Back') != -1) and (plasticans.split('_')[0] == 'Hab') and (plasticans.split('_')[1] == "rec"):
		
		# Get location in genes array for plastic region
		# ----------------------------------------------
		Indgenes = currentoff['genes']
		# If cdevolve is on
		if cdevolveans != 'N':
			# Then the first l loci are for selection, next for plastic region
			if cdevolveans.split('_')[0] == 'P': # This is for multilocus selection, not currently implemented, to be moved over from cdpop
				selloci = int(cdevolveans.split('_')[2].split('L')[1])
			elif cdevolveans == '1' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans == 'stray' or cdevolveans == 'Hindex':
				selloci = 1
			elif cdevolveans == '2' or cdevolveans == 'MG' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == '2_mat':
				selloci = 2
			else:
				print('CDEVOLVEANS not entered correctly; DoUpdate() error.')
				sys.exit(-1)
		# If selection is not on
		else:
			selloci = 0 # zero loci in selection
		# Get number of plastic loci
		plaloci = 1
		# Get index for plastic region
		plaloci_index = range(selloci*2,selloci*2+plaloci*2)
		
		# Get the plastic behaviorial response threshold
		# ----------------------------------------------
		
		# Does this individual have the allele for response (2,2)
		if Indgenes[plaloci_index][0] == 2 and Indgenes[plaloci_index][1] == 2:
					
			if answer == 'immigrator': # probarray is length 1
				# Only check patchvals in natalsubpop, which is where this individual is going if probarray is 1
				if patchvals[natalsubpop] >= plastic_behaviorresp:
					probarray = [0.0]
				
			else: # probarray is length patch dim
				# Whereever temperature threshold, turn prob to 0
				probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2]))
                
	if (plasticans != 'N') and (gen >= burningen_plastic) and (timeplastic.find('Back') != -1) and (plasticans.split('_')[0] == 'Hab') and (plasticans.split('_')[1] == "codom"):
		
		# Get location in genes array for plastic region
		# ----------------------------------------------
		Indgenes = currentoff['genes']
		# If cdevolve is on
		if cdevolveans != 'N':
			# Then the first l loci are for selection, next for plastic region
			if cdevolveans.split('_')[0] == 'P': # This is for multilocus selection, not currently implemented, to be moved over from cdpop
				selloci = int(cdevolveans.split('_')[2].split('L')[1])
			elif cdevolveans == '1' or cdevolveans == 'M' or cdevolveans == 'G' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link' or cdevolveans == 'stray' or cdevolveans == 'Hindex':
				selloci = 1
			elif cdevolveans == '2' or cdevolveans == 'MG' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == '2_mat':
				selloci = 2
			else:
				print('CDEVOLVEANS not entered correctly; DoUpdate() error.')
				sys.exit(-1)
		# If selection is not on
		else:
			selloci = 0 # zero loci in selection
		# Get number of plastic loci
		plaloci = 1
		# Get index for plastic region
		plaloci_index = range(selloci*2,selloci*2+plaloci*2)
		
		# Get the plastic behaviorial response threshold
		# ----------------------------------------------
		
		# Does this individual have the allele to intiate response
		# 2 in both = 0, one 2 = 50% reduction
		if Indgenes[plaloci_index][0] == 2 and Indgenes[plaloci_index][1] == 2:
		
			# Whereever temperature threshold, turn prob to value input as third part of plasticans
			probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2]))
            
		if Indgenes[plaloci_index][0] == 1 and Indgenes[plaloci_index][1] == 2:
		
			# Whereever temperature threshold, reduce probably y 50% (calc by mean) of the value input as the third part of plasticans
			probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = ( (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2])) + probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)] ) /2
            
		if Indgenes[plaloci_index][0] == 2 and Indgenes[plaloci_index][1] == 1:
		
			# Whereever temperature threshold, reduce probably y 50% (calc by mean) of the value input as the third part of plasticans
			probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)[0]] = ( (probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)]) * (float(plasticans.split('_')[2])) + probarray[np.where(np.asarray(patchvals) >= plastic_behaviorresp)] ) /2
                

	return probarray
	
	# End::GetProbArray()
	
# ---------------------------------------------------------------------------------------------------	
def Immigration(SubpopIN,K,natal_patches,gen,cdevolveans,fitvals,SelectionDeaths,DisperseDeaths,burningen_cdevolve,Str_patch_prob,StrSuccess,age_S,Population,dtype,sizecall,size_mean,PackingDeaths,PopulationAge,packans,PackingDeathsAge,packpar1,homeattempt,timecdevolve,patchvals,PopTag,subpopmort_mat,Track_YYSelectionPackDeaths,Track_WildSelectionPackDeaths,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,age_percmort,comp_coef,XQs,Kadj_track,Track_KadjEmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,f_leslie,f_leslie_std,ageDispProb,cdmatrix_FXXBack,cdmatrix_MXYBack,cdmatrix_MYYBack,cdmatrix_FYYBack,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD,sexchromo,age_MgBACK,MgBack_patch_prob,Disperse_patch_prob,MgOut_patch_prob,age_MgOUT,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,migrate,egg_add):
	
	'''
	Immigration()
	This function settles individuals to a location through residency, straying, immigration.
	SubpopIN - (NatalPop,EmiPop,ImmiPop,age,sex,infection,name,genes)
	original dic name 'NatalPop' = CurrentPop or where they will spawn, could be reassigned from Local dispersal movement
	natalP var is pulled from the classfile to grab actual Natal Pop where this individual was born - natalP does NOT change, but 'NatalPop' could change
	EmiPop - where it moved to
	ImmiPop - where it moved from
	
	'''			
		
	# Store the keepers
	SubpopIN_keep = []
			
	# Get unique number of subpops
	nosubpops = len(K)
	
	# Add spot to track numbers of individuals
	DisperseDeaths.append([])
	SelectionDeaths.append([])
	StrSuccess.append([]) # Prob switches
	PackingDeaths.append([])
	Track_YYSelectionPackDeaths.append([]) # YY selection deaths in packing algorithm
	Track_WildSelectionPackDeaths.append([]) # wild type selection deaths in packing algorithm
	for x in range(0,nosubpops):
		DisperseDeaths[gen].append([]) 
		SelectionDeaths[gen].append([])
		PackingDeaths[gen].append([])
		Track_YYSelectionPackDeaths[gen].append([])
		Track_WildSelectionPackDeaths[gen].append([])
		SubpopIN_keep.append([])	
	
	# If F Selection is On, then calculation of EHom here for entire population
	if (cdevolveans.split('_')[0] == 'F' or cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
		EHom = calc_EHom(SubpopIN)
	else:
		EHom = [] # Make empty to pass null into callDiffMortality funcitons
	#pdb.set_trace()
	# Decide where everybody moves - loop through subpop, doesn't need to be shuffled.
	for isub in range(len(SubpopIN)):

		# Loop through individuals in this subpopulation
		for iind in range(len(SubpopIN[isub])):
			# Get individual, easier for indexing
			outpool = SubpopIN[isub][iind]
			originalpop = outpool['NatalPop'] # OR patch trying to get back to
			emipop = outpool['EmiPop']
			outpool_genes = outpool['genes']
			indSex = outpool['sex']
			# For indexing into 
			if indSex == 'FXX':
				sxspot = 0
			elif indSex == 'MXY':
				sxspot = 1
			elif indSex == 'MYY':
				sxspot = 2
			else:
				sxspot = 3

			# Get this individuals original ClassVars file and bins for indexing
			natalP = int(SubpopIN[isub][iind]['classfile'].split('_')[0].split('P')[1])
			theseclasspars = int(SubpopIN[isub][iind]['classfile'].split('_')[1].split('CV')[1])
			
			# -----------------------------------------------------------------------------------------
			# Check if Migrant - 'E' assumes they have moved OUT or 'EO' they have been out for a while
			# -----------------------------------------------------------------------------------------
			#if isub + 1 != int(emipop): # If current patch is NOT the same as the individuals 'emipop' == 'Resident'
			if outpool['name'][0] == 'E' or (outpool['name'][0:4] == 'Age0' and not packans == 'anadromy' and egg_add == 'nonmating'):
			#if outpool['name'][0] == 'E' or outpool['name'][0:4] == 'Age0':
			#if outpool['name'][0] == 'E':
				
				#if outpool['name'][0:4] == 'Age0' and not packans == 'anadromy' and egg_add == 'nonmating':
				#	print('What patch do return age0 migrants go?')
				#	sys.exit(-1)
				
				#----------------------------------------------------
				# Check if migrate back (age_MgBACK) * (patch_MgBack)
				#---------------------------------------------------- 			
				# Add here the check for patch migration back probability
				MgBACK_Patch = MgBack_patch_prob[isub]
				
				# Check for age or size based migration back
				indexofProb = outpool[sizecall]
				# If size control then get size nearest to values in age file
				if sizecall == 'size':
					closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb))
					Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0]
				else:
					Find = indexofProb
				# Check for ages over last age (-1 for length counting)
				if Find > len(age_MgBACK[natalP][theseclasspars]) - 1:
					Find = len(age_MgBACK[natalP][theseclasspars]) - 1 # Make last age
				MgBACK_Class = age_MgBACK[natalP][theseclasspars][Find]
				
				# Get sex class options if given
				if len(MgBACK_Class.split('~')) == 1:
					MgBACK_Class = float(MgBACK_Class.split('~')[0])
				elif len(MgBACK_Class.split('~')) != sexchromo:
					print('Number of age-specific probability parameters must match sex_chromo.')
					sys.exit(-1)
				else:
					MgBACK_Class = float(MgBACK_Class.split('~')[sxspot])
				
				# Then multiply these together or update when patch-based Mg Back is added
				indProb = MgBACK_Patch * MgBACK_Class
				
				# Decide if immigrant - or individual that migrates back
				randProb = np.random.uniform()	# Get a random number			
				# Flip the coin for moving back
				if randProb < indProb:										
					# Then True
					indProbans_MgBACK = 'Yes'							
				# Not true
				else:
					indProbans_MgBACK = 'No'				
				
				# ------------------------------------------------------------------------------
				# Individual is migrating back (either through straying or back to natal grounds)
				# -------------------------------------------------------------------------------
				if indProbans_MgBACK == 'Yes':	
												
					# -----------------------------------
					# Check for strayer first
					# -----------------------------------
					# No genotype tied to stray rate
					if cdevolveans != 'stray':
											
						# Check for patch stray probability - here we select the stray patch prob of it's origin pop
						Str_Patch = Str_patch_prob[int(originalpop)-1] # index 1 minus
						#Str_Patch = Str_patch_prob[isub] # changed 12/22/2023; stray from patch in or natal?
						
						# Check for age/size stray
						indexofProb = outpool[sizecall]
						# If size control then get size nearest to values in age file
						if sizecall == 'size':
							closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb))
							Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0]
						else:
							Find = indexofProb
						# Check for ages over last age
						if Find > len(age_S[natalP][theseclasspars]) - 1:
							Find = len(age_S[natalP][theseclasspars]) - 1 # Make last age
						Str_Class = age_S[natalP][theseclasspars][Find]
													
						# Get sex class options if given
						if len(Str_Class.split('~')) == 1:
							Str_Class = float(Str_Class.split('~')[0])
						elif len(Str_Class.split('~')) != sexchromo:
							print('Number of age-specific probability parameters must match sex_chromo.')
							sys.exit(-1)
						else:
							Str_Class = float(Str_Class.split('~')[sxspot])

						# Then multiply these together
						indProb = Str_Patch * Str_Class
											
					# If genotype tied to stray rate
					else:
						#pdb.set_trace()
						if outpool_genes[0] == 2: #AA
							indProb = float(fitvals[isub][0])
						elif outpool_genes[0] == 1 and outpool_genes[1] == 1: #Aa
							indProb = float(fitvals[isub][1])
						elif outpool_genes[1] == 2: #aa
							indProb = float(fitvals[isub][2])
						else: #Something else
							print('Straying genotype in use - only single locus model used.')
							sys.exit(-1)							

					# Decide if strayer?
					randProb = np.random.uniform()	# Get a random number			
					# Flip the coin for patch stray
					if randProb < indProb:			
										
						# Then stray
						indProbans_stray = 'Yes'
								
					# Patch stray not a success
					else:
						indProbans_stray = 'No'
					
					# ---------------------------------
					# Straying occurred for the Migrant
					# ---------------------------------
					if indProbans_stray == 'Yes':				
						# Stray success - go anywhere from current location and use full cdmatrix, minus where K = 0 and natal grounds are 0
						# ---------------------------------------------------------
						
						# Get the probability array
						probarray = GetProbArray(outpool,'strayer_emiPop',K,natal_patches,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr)
												
						# If spot available to stray to
						if sum(probarray) != 0.0:
						
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)
							namethis = 'S'
						
						# If spot was not available to stray to ('Move Deaths')
						elif sum(probarray) == 0.0:
							# Then Break from the loop and move to next outpool	
							print("Nowhere to stray to. This individual dies.")
							SelectionDeaths[gen][int(emipop)-1].append(0)
							DisperseDeaths[gen][int(emipop)-1].append(1)
							StrSuccess[gen].append(1)
							continue
							
						# ----------------------------------------------------------
						# Check CDEVOLVE Selection Death and spatial mortality Death
						# ----------------------------------------------------------
						differentialmortality = callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,'Back',outpool,fitvals,iteminlist,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom)
													
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][iteminlist].append(1)
							DisperseDeaths[gen][iteminlist].append(0)
							StrSuccess[gen].append(1)
							continue						
												
						# For recording purposes
						# -----------------------
						disppop = str(iteminlist+1)
						immipop = str(iteminlist+1)
						outpool_name = outpool['name']
						outpool_name = outpool_name.split('_')
						
						name = namethis+str(disppop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
						
						# Note 'NatalPop' gets updated here to the new straypop now
						recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
									
						# Record outpool disperse information	
						SubpopIN_keep[int(disppop)-1].append(recd)		
						
						# Record outpool disperse information	
						SelectionDeaths[gen][int(disppop)-1].append(0)
						DisperseDeaths[gen][int(disppop)-1].append(0)
						StrSuccess[gen].append(1)
						continue													
														
					#End::Prob success of stray
					# -------------------------------
					
					# -----------------------
					# No straying for Migrant
					# -----------------------
					else:
						# -------------------------------------
						# Attempt to move back to natal grounds
						# -------------------------------------
						probarray = GetProbArray(outpool,'immigrator',K,natal_patches,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXXBack,cdmatrix_MXYBack,cdmatrix_MYYBack,cdmatrix_FYYBack,fitvals)
												
						# If spots available to move to:
						if sum(probarray) != 0.0:													
							
							# Special case for runtiming - could be an array here, then random grab new natal patch to move back to
							if len(probarray) > 1:
								
								# Select the w_choice item to move to
								iteminlist = w_choice_item(probarray)
								namethis = 'I'
								
								# ----------------------------------------------------------
								# Check CDEVOLVE Selection Death and spatial mortality Death
								# ----------------------------------------------------------
								differentialmortality = callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,'Back',outpool,fitvals,iteminlist,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom)
																								
								# Then flip the coin to see if outpool survives its location
								randcheck = np.random.uniform()
								# If outpool did not survive: break from loop, move to next outpool
								if randcheck < differentialmortality:
									SelectionDeaths[gen][iteminlist].append(1)
									DisperseDeaths[gen][iteminlist].append(0)
									StrSuccess[gen].append(0)
									continue
																	
								# For recording purposes 
								# -----------------------
								disppop = str(iteminlist+1)
								immipop = str(iteminlist+1)							
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = namethis+str(disppop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]	
																		
								recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
													
								# Record outpool disperse information	
								SubpopIN_keep[int(disppop)-1].append(recd)				
							
								SelectionDeaths[gen][int(disppop)-1].append(0)
								DisperseDeaths[gen][int(disppop)-1].append(0)
								StrSuccess[gen].append(0)
								continue	

							# Else probarray returned one value either [1.0] and it will move back to its natal patch and also attempt to local disperse 
							else:
								iteminlist = int(originalpop)-1
							
								# ------------------------------------------------------------------------------
								# Attempt Local Dispersal
								# ------------------------------------------------------------------------------
								Disp_Patch = Disperse_patch_prob[iteminlist]						
								
								# Next check if it decides to local dispersal given population and age
								indexofProb = outpool[sizecall] # Check for age/size 
								# If size control then get size nearest to values in age file
								if sizecall == 'size':
									closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb))
									Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0]
								else:
									Find = indexofProb
								# Check for ages over last age
								if Find > len(ageDispProb[natalP][theseclasspars]) - 1:
									Find = len(ageDispProb[natalP][theseclasspars]) - 1 # Make last age					
								Disp_Class = ageDispProb[natalP][theseclasspars][Find]
								# Get sex class options if given
								if len(Disp_Class.split('~')) == 1:
									Disp_Class = float(Disp_Class.split('~')[0])
								elif len(Disp_Class.split('~')) != sexchromo:
									print('Number of age-specific capture probability parameters must match sex_chromo.')
									sys.exit(-1)
								else:
									Disp_Class = float(Disp_Class.split('~')[sxspot])
								
								# Then multiply these together
								indProb = Disp_Patch * Disp_Class
								
								# Decide if disperse?
								randProb = np.random.uniform()	# Get a random number			
								
								# Flip the coin for patch stray
								if randProb < indProb:									
									# Then Disperse
									indProbans_localD = True
									
								# Patch stray not a success
								else:
									indProbans_localD = False
								
								# ------------------------------------------------------------------------
								# Local disperse this individual from it's natal patch once it moved back
								# ------------------------------------------------------------------------
								if indProbans_localD:
									# Then locally move this individual - go anywhere from natal patch and use full cdmatrix, minus where K = 0 and natal grounds are 0
									# ---------------------------------------------------------
									
									# Get the probability array - be careful not to overwrite previous if statement!
									probarray_LD = GetProbArray(outpool,'localDispersal',K,natal_patches,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD)
									
									# If spots available to move to:
									if sum(probarray_LD) != 0.0:
									
										# Select the w_choice item
										iteminlist = w_choice_item(probarray_LD)
										namethis = 'ID'									
										
									# No spots available to move to, then stays in natal patch
									elif sum(probarray_LD) == 0.0:
										
										# Then Break from the loop and move to next outpool	
										print("Nowhere for this immigrant to disperse to. This individual dies.")
										SelectionDeaths[gen][int(emipop)-1].append(0)
										DisperseDeaths[gen][int(originalpop)-1].append(1)
										StrSuccess[gen].append(0)
										continue
								# ---------------------------------------------------------------------------------
								# Local dispersal is turned off - stays in it's natal patch area int(originalpop)-1
								# ---------------------------------------------------------------------------------	
								else:
									# Then it makes it back to original pop - or pop it came from
									iteminlist = int(originalpop)-1
									namethis = 'I'
								#pdb.set_trace()	
								# ----------------------------------------------------------
								# Check CDEVOLVE Selection Death and spatial mortality Death
								# ----------------------------------------------------------
								differentialmortality = callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,'Back',outpool,fitvals,iteminlist,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom)
																								
								# Then flip the coin to see if outpool survives its location
								randcheck = np.random.uniform()
								# If outpool did not survive: break from loop, move to next outpool
								if randcheck < differentialmortality:
									SelectionDeaths[gen][iteminlist].append(1)
									DisperseDeaths[gen][iteminlist].append(0)
									StrSuccess[gen].append(0)
									continue
																	
								# For recording purposes 
								# -----------------------
								disppop = str(iteminlist+1)
								immipop = str(iteminlist+1)							
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = namethis+str(disppop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]	
																		
								recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
													
								# Record outpool disperse information	
								SubpopIN_keep[int(disppop)-1].append(recd)				
							
								SelectionDeaths[gen][int(disppop)-1].append(0)
								DisperseDeaths[gen][int(disppop)-1].append(0)
								StrSuccess[gen].append(0)
								continue							
										
						# It did not make it back to natal grounds
						elif sum(probarray) == 0.0:
							# Could not make it home (either K = 0 or exceded threshold (complete or partial barrier)
							if homeattempt == 'mortality':
								# Store information
								SelectionDeaths[gen][int(emipop)-1].append(0)
								DisperseDeaths[gen][int(emipop)-1].append(1)
								StrSuccess[gen].append(0)
								# Then Break from the loop and move to next outpool
								continue
							
							# Then attempt to stray the individual one more time
							elif homeattempt == 'stray_emiPop':

								# Get the probability array
								probarray = GetProbArray(outpool,'strayer_emiPop',K,natal_patches,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr)
			
							# Then attempt to stray the individual one more time
							elif homeattempt == 'stray_natalPop':
								
								# Get the probability array
								probarray = GetProbArray(outpool,'strayer_natalPop',K,natal_patches,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr)	
								
							else: # Error check
								print('Home attempt must be either mortality or stray options. See user manual.')
								sys.exit(-1)								
								
							# If statement to check if there are spots for offpsring to stray to
							if sum(probarray) != 0.0:
								
								# Select the w_choice item
								iteminlist = w_choice_item(probarray)
								namethis = 'Z'
								
							# If statement to check if there were not spots to disperse to in straying
							elif sum(probarray) == 0.0:
								# Then Break from the loop and move to next outpool	
								print("Return migrant cannot make it back, cannot stray, and dies.")
								SelectionDeaths[gen][int(emipop)-1].append(0)
								DisperseDeaths[gen][int(emipop)-1].append(1)
								StrSuccess[gen].append(1)
								continue
									
							# ----------------------------------------------------------
							# Check CDEVOLVE Selection Death and spatial mortality Death
							# ----------------------------------------------------------
							differentialmortality = callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,'Back',outpool,fitvals,iteminlist,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom)
																							
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][iteminlist].append(1)
								DisperseDeaths[gen][iteminlist].append(0)
								StrSuccess[gen].append(1)
								continue
																	
							# For recording purposes 
							# -----------------------
							disppop = str(iteminlist+1)
							immipop = disppop
							outpool_name = outpool['name']
							outpool_name = outpool_name.split('_')
							
							name = namethis+str(disppop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
								
							recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
											
							# Record outpool disperse information	
							SubpopIN_keep[int(disppop)-1].append(recd)		
							
							# Record outpool disperse information	
							SelectionDeaths[gen][int(disppop)-1].append(0)
							DisperseDeaths[gen][int(disppop)-1].append(0)
							StrSuccess[gen].append(1)
							continue										
								
						
				#End::Migrant Checks
				# ------------------------
				
				# ------------------------------------------------------------
				# Individual is not migrating back - consider local dispersal 
				# ------------------------------------------------------------
				else:
					# Check local dispersal, otherwise resident, not currently implemented
					#pdb.set_trace()
					Disp_Patch = Disperse_patch_prob[int(isub)]
					# Next check if it decides to local dispersal given population and age
					indexofProb = outpool[sizecall] # Check for age/size 
					# If size control then get size nearest to values in age file
					if sizecall == 'size':
						closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb))
						Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0]
					else:
						Find = indexofProb
					# Check for ages over last age
					if Find > len(ageDispProb[natalP][theseclasspars]) - 1:
						Find = len(ageDispProb[natalP][theseclasspars]) - 1 # Make last age					
					Disp_Class = ageDispProb[natalP][theseclasspars][Find]
					# Get sex class options if given
					if len(Disp_Class.split('~')) == 1:
						Disp_Class = float(Disp_Class.split('~')[0])
					elif len(Disp_Class.split('~')) != sexchromo:
						print('Number of age-specific capture probability parameters must match sex_chromo.')
						sys.exit(-1)
					else:
						Disp_Class = float(Disp_Class.split('~')[sxspot])
					
					# Then multiply these together
					indProb = Disp_Patch * Disp_Class					
					
					# Decide if disperse?
					randProb = np.random.uniform()	# Get a random number			
					
					# Flip the coin for patch stray
					if randProb < indProb:									
						# Then Disperse
						indProbans_localD = True
						
					# Patch stray not a success
					else:
						indProbans_localD = False
					
					# ------------------------------------------------------------------------
					# Local disperse this individual from it's migration patch
					# ------------------------------------------------------------------------
					if indProbans_localD:
						
						# Then locally move this individual - go anywhere from natal patch and use full cdmatrix, minus where K = 0 and natal grounds are 0
						# ---------------------------------------------------------
						
						# Get the probability array - be careful not to overwrite previous if statement!
						probarray_LD = GetProbArray(outpool,'localDispersal',K,natal_patches,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD)
						
						#-----------------------
						# Open spots to move to
						#----------------------
						if sum(probarray_LD) != 0.0:
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray_LD)
							
							# Age0s should not turn Es here - RDs
							if outpool['name'][0:4] == 'Age0':
								namethis = 'RD'
							# Not Age0s
							else:
								namethis = 'ED' 
							
						#-----------------------
						# NO Open spots to move to
						#----------------------
						elif sum(probarray_LD) == 0.0:
							
							print('Nowhere for this migrant to locally disperse to. This individual will die.')
							# Then Break from the loop and move to next outpool	
							SelectionDeaths[gen][int(emipop)-1].append(0)
							DisperseDeaths[gen][int(emipop)-1].append(1)
							StrSuccess[gen].append(0)
							continue
					# ---------------------------------------------------------------------------------
					# Local dispersal is turned off - stays in it's migration patch
					# ---------------------------------------------------------------------------------	
					else:
						# Set the iteminlist to not move
						iteminlist = isub
						
						# Age0s should not turn Es here - Rs
						if outpool['name'][0:4] == 'Age0':
							namethis = 'R'
						# Not Age0s
						else:
							namethis = 'EO' 
						 
						
					# ----------------------------------------------------------
					# Check CDEVOLVE Selection Death and spatial mortality Death
					# ----------------------------------------------------------
					differentialmortality = callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,'Back',outpool,fitvals,iteminlist,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom)
																					
					# Then flip the coin to see if outpool survives its location
					randcheck = np.random.uniform()
					# If outpool did not survive: break from loop, move to next outpool
					if randcheck < differentialmortality:
						SelectionDeaths[gen][iteminlist].append(1)
						DisperseDeaths[gen][iteminlist].append(0)
						StrSuccess[gen].append(0)
						continue
														
					# For recording purposes 
					# -----------------------
					disppop = str(iteminlist+1)
					immipop = disppop							# CHECK. Is this the patch emigrated from?
					outpool_name = outpool['name']
					outpool_name = outpool_name.split('_')
					
					#name = namethis+str(disppop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]	
					# Not sure if this is right to get the individual to migrate back to its original patch after potentially dispersing to new emigration patch
					name = namethis+str(originalpop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]	
															
					#recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
					# Because this individual is in the ocean, we keep original pop and emipop becomes disppop
					recd = (originalpop,disppop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
										
					# Record outpool disperse information	
					SubpopIN_keep[int(disppop)-1].append(recd)				
				
					SelectionDeaths[gen][int(disppop)-1].append(0)
					DisperseDeaths[gen][int(disppop)-1].append(0)
					StrSuccess[gen].append(0)
					continue							
					
			# ------------------------------------------------------------------------------
			# Else not a migrant (resident or Age0) - 'emipop' == 'Resident' or isub + 1 == 'emipop'
			# ------------------------------------------------------------------------------
			else:
				
				# -----------------------------
				# Check for a Migration OUT
				#------------------------------
				MgOut_Patch = MgOut_patch_prob[isub]
				# Next check on age level
				indexofProb = outpool[sizecall] # Check for age/size 
				# If size control then get size nearest to values in age file
				if sizecall == 'size':
					closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb))
					Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0]
				else:
					Find = indexofProb
				# Check for ages over last age
				if Find > len(age_MgOUT[natalP][theseclasspars]) - 1:
					Find = len(age_MgOUT[natalP][theseclasspars]) - 1 # Make last age					
				MgOut_Class = age_MgOUT[natalP][theseclasspars][Find]
				# Get sex class options if given
				if len(MgOut_Class.split('~')) == 1:
					MgOut_Class = float(MgOut_Class.split('~')[0])
				elif len(MgOut_Class.split('~')) != sexchromo:
					print('Number of age-specific capture probability parameters must match sex_chromo.')
					sys.exit(-1)
				else:
					MgOut_Class = float(MgOut_Class.split('~')[sxspot])
				
				# Then multiply these together
				indProb = MgOut_Patch * MgOut_Class
								
				# Decide if disperse?
				randProb = np.random.uniform()	# Get a random number			
				# Flip the coin for patch stray
				if randProb < indProb:									
					# Then Disperse
					indProbans_MgOut = True
				# Patch stray not a success
				else:
					indProbans_MgOut = False
				
				# -----------------------------------------------------------
				# Age0 individual that migrates Out for anadromy option only
				# -----------------------------------------------------------
				if outpool['name'][0:4] == 'Age0' and packans == 'anadromy' and indProbans_MgOut:
					#if outpool['name'][0:4] == 'Age0' and indProbans_MgOut and packans == 'anadromy':
					#pdb.set_trace()
					# -------------------------------------------------------------	--------------Then migrate Out				
					# Get the probability array - call cdmatrix from current pop isub or natal pop to move to migration grounds only
					probarray = GetProbArray(outpool,'localDispersal',K,migrate,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut)				
				
					#-----------------------
					# Open spots to move to
					#----------------------
					if sum(probarray) != 0.0:
						
						# Select the w_choice item
						iteminlist = w_choice_item(probarray)	
						
						# ----------------------------------------------------------
						# Check CDEVOLVE Selection Death and spatial mortality Death - Note selection option changed to OUT
						# ----------------------------------------------------------
						differentialmortality = callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,'Out',outpool,fitvals,iteminlist,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom)													
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][iteminlist].append(1)
							DisperseDeaths[gen][iteminlist].append(0)
							StrSuccess[gen].append(1)
							continue					
						
						# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-get in DoCal,age,sex,size,states,capture,name
						# ---------------------------------------
						disppop = str(iteminlist+1)
						immipop = disppop
						outpool_name = outpool['name']
						outpool_name = outpool_name.split('_')
						name = 'E'+str(disppop)+'_F'+str(isub+1)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
						
						recd = (str(isub+1),disppop,'NA',outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
									
						# Record outpool disperse information	
						SubpopIN_keep[int(disppop)-1].append(recd)				
					
						SelectionDeaths[gen][int(disppop)-1].append(0)
						DisperseDeaths[gen][int(disppop)-1].append(0)
						StrSuccess[gen].append(0)
						continue
							
					#-------------------------
					# NO Open spots to move to
					#------------------------- 
					elif sum(probarray) == 0.0:
						print('Nowhere to locally disperse to. All individuals will die.')
						# Then Break from the loop and move to next outpool	
						SelectionDeaths[gen][int(emipop)-1].append(0)
						DisperseDeaths[gen][int(emipop)-1].append(1)
						StrSuccess[gen].append(1)
						continue
										
				# -------------------------------------------------------
				# All IDs (including Age0s)
				# -------------------------------------------------------
				else:
					# Local dispersal check
					Disp_Patch = Disperse_patch_prob[isub]							
								
					# Next check if it decides to local dispersal given population and age
					indexofProb = outpool[sizecall] # Check for age/size 
					# If size control then get size nearest to values in age file
					if sizecall == 'size':
						closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb))
						Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0]
					else:
						Find = indexofProb
					# Check for ages over last age
					if Find > len(ageDispProb[natalP][theseclasspars]) - 1:
						Find = len(ageDispProb[natalP][theseclasspars]) - 1 # Make last age					
					Disp_Class = ageDispProb[natalP][theseclasspars][Find]
					# Get sex class options if given
					if len(Disp_Class.split('~')) == 1:
						Disp_Class = float(Disp_Class.split('~')[0])
					elif len(Disp_Class.split('~')) != sexchromo:
						print('Number of age-specific capture probability parameters must match sex_chromo.')
						sys.exit(-1)
					else:
						Disp_Class = float(Disp_Class.split('~')[sxspot])
					
					# Then multiply these together
					indProb = Disp_Patch * Disp_Class
									
					# Decide if disperse?
					randProb = np.random.uniform()	# Get a random number			
					# Flip the coin for patch stray
					if randProb < indProb:									
						# Then Disperse
						indProbans_RLD = True
					# Patch stray not a success
					else:
						indProbans_RLD = False
					
					# ------------------------------
					# Local disperse this individual
					# ------------------------------
					if indProbans_RLD:
						# Then locally move this individual - go anywhere from current location and use full cdmatrix, minus where K = 0 and natal grounds are 0
						# ---------------------------------------------------------					
						# Get the probability array
						probarray = GetProbArray(outpool,'localDispersal',K,natal_patches,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD)			
					
						#-----------------------
						# Open spots to move to
						#----------------------
						if sum(probarray) != 0.0:
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)	
							
							# ----------------------------------------------------------
							# Check CDEVOLVE Selection Death and spatial mortality Death
							# ----------------------------------------------------------
							differentialmortality = callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,'Back',outpool,fitvals,iteminlist,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom)													
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][iteminlist].append(1)
								DisperseDeaths[gen][iteminlist].append(0)
								StrSuccess[gen].append(1)
								continue					
							
							# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-get in DoCal,age,sex,size,states,capture,name
							# ---------------------------------------
							disppop = str(iteminlist+1)
							immipop = disppop
							outpool_name = outpool['name']
							outpool_name = outpool_name.split('_')
							name = 'RD'+str(disppop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
							
							recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
										
							# Record outpool disperse information	
							SubpopIN_keep[int(disppop)-1].append(recd)				
						
							SelectionDeaths[gen][int(disppop)-1].append(0)
							DisperseDeaths[gen][int(disppop)-1].append(0)
							StrSuccess[gen].append(0)
							continue
								
						#-----------------------
						# NO Open spots to move to
						#---------------------- 
						elif sum(probarray) == 0.0:
							print('Nowhere to locally disperse to. All individuals will die.')
							# Then Break from the loop and move to next outpool	
							SelectionDeaths[gen][int(emipop)-1].append(0)
							DisperseDeaths[gen][int(emipop)-1].append(1)
							StrSuccess[gen].append(1)
							continue
							
					# -------------------------------------------------------
					# Local dispersal is turned off - Stays in Resident patch
					# -------------------------------------------------------
					else:					
						iteminlist = isub
						# ------------------------------
						# Check CDEVOLVE Selection Death
						# ------------------------------
						differentialmortality = callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,'Back',outpool,fitvals,iteminlist,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom)												
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][iteminlist].append(1)
							DisperseDeaths[gen][iteminlist].append(0)
							StrSuccess[gen].append(0)
							continue
							
						# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-get in DoCal,age,sex,size,states,capture,name
						# ---------------------------------------
						disppop = str(iteminlist+1)
						immipop = disppop # or Currentpop
						outpool_name = outpool['name']
						outpool_name = outpool_name.split('_')
						name = 'R'+str(immipop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
						
						recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['states']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(immipop)-1],outpool['species'],outpool['genes'])
									
						# Record outpool disperse information	
						SubpopIN_keep[int(immipop)-1].append(recd)				
					
						SelectionDeaths[gen][int(immipop)-1].append(0)
						DisperseDeaths[gen][int(immipop)-1].append(0)
						StrSuccess[gen].append(0)
						continue				
						
			#End::Resident Checks 
			# -------------------
				
		#End::For loop individual
		#------------------------
			
	#End::For loop Subpop
	#---------------------
	#pdb.set_trace()
	# ----------------------------------------------
	# SubpopIN - sort/rank/pack by Kage options
	# ----------------------------------------------	
	# If multiple ClassVars are given then bin min to max
	if sizecall == 'size':
		'''
		bin_min = min(sum(sum(size_mean,[]),[]))
		bin_max = max(sum(sum(size_mean,[]),[]))
		size_bin = [bin_min]
		for ibin in range(len(size_mean[0][0])-1):
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
	PackingDeathsAge[gen] = [[] for x in range(0,len(size_mean[0][0]))]
	PopulationAge[gen] = [[] for x in range(0,len(size_mean[0][0]))]	
	Kadj_track.append([])
	#pdb.set_trace()
	# -------------------
	# Packing is selected
	# -------------------
	if packans == 'packing':

		# ------------------------------------------
		# Get other species Ns from all Patches here
		# ------------------------------------------
		if implementcomp == 'Back':
			Nself_pop = [len(SubpopIN_keep[x]) for x in range(0,len(SubpopIN_keep))]
			# Ignore queus if only one species
			if len(XQs) > 0:
				# Loop through queue spots, Putting Nself_pop for grabbing for other species.
				for ispecies in range(len(XQs[spcNO])):
					XQs[spcNO][ispecies].put(Nself_pop) 
			Nother_pop = []
			popcount = 0
			# Ignore queues if only one species
			if len(XQs) > 0:
				for ispecies in range(len(XQs[spcNO])):			
					if ispecies != spcNO: 
						Nother_pop.append(XQs[ispecies][spcNO].get(block = True))
						# Error check here for equal patches
						if len(Nother_pop[popcount]) != len(Nself_pop):
							print("Species systems must have the same number of patches.")
							sys.exit(-1)
						popcount = popcount + 1	
		
		# Loop through each subpopulation
		for isub in range(len(K)):
					
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# -----------------------------------
			# Get countages - for 'age' adjusted
			# -----------------------------------
			# Switch here for size or age control
			if sizecall == 'size': # Careful here and use first ClassVars
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_arr[sizecall])
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)
			else:
				# Count up each uniages
				countages = count_unique(SubpopIN_arr['age'])
			
			# ------------------------
			# Get packing parameters
			# ------------------------			
			# K,N for this population
			Kpop = K[isub]
			Npop = sum(countages[1])
			
			# Overwrite Kpop with K adjusted from DoEmigration that was updated from competition - careful of indexing Tracking variable is patch + 1, first value is total
			if implementcomp == 'Out':
				Kpop = Track_KadjEmi[gen][isub+1]
						
			if implementcomp == 'Back':
				#if multiprocessing.current_process().name == "S1":
				#	pdb.set_trace()	
				# -----------------
				# Adjustment for K if more than one species
				# -----------------
				#Ignore K adjustment if only one species
				if len(XQs) > 0:
					alphas = comp_coef[isub].split(';') # Extracting alphas
					if gen >= startcomp: # Check timing 
						tempspeciesSum = [] # Create sum list for  >= 2 species
						for ispecies in range(len(Nother_pop)):
							tempspeciesSum.append(Nother_pop[ispecies][isub]*float(alphas[ispecies]))
						Kpop = int(Kpop - sum(tempspeciesSum)) # Adjustment or K		
			Kadj_track[gen].append(Kpop) # For Tracking
			
			# Continue packing
			if Npop == 0 or Kpop <= 0:
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
				#Kij_proportion = np.exp(packpar1*(countages[0]+1)) # Equ 8
				Kij_proportion = np.exp(packpar1*(np.asarray(list(range(len(size_mean[0][0]))))+1))
				Kij_proportion = Kij_proportion/sum(Kij_proportion) # Equ 8 rescaled
				Kij_proportion = np.flipud(Kij_proportion) # Flip to match next age loop (old to young)
				AgeClass_reverse = np.flipud(np.asarray(list(range(len(size_mean[0][0])))))
				#for iage in -np.sort(-countages[0]):
				for iage in AgeClass_reverse:
					
					# Age class
					Ageclass = iage
					
					# Special case when age is greater than last age class only used for indexing now
					if Ageclass > len(size_mean[0][0])-1:
						indexforAgeclass = len(size_mean[0][0]) - 1
					else:
						indexforAgeclass = Ageclass
					
					# N for this age coming in
					if len(np.where(countages[0]==Ageclass)[0]) == 0:
						Nage = 0
					else:
						Nage = countages[1][np.where(countages[0]==Ageclass)[0][0]]
					if sizecall == 'size': # Use the adjusted age classes
						Nage_index = np.where(age_adjusted==Ageclass)[0]
					else:
						Nage_index = np.where(SubpopIN_arr['age']==Ageclass)[0]
					
					# Get Age_scaling for this pop's age class (*Stage R)
					Agescale = np.exp(Kscale * (1. - (Nage / float(Kpop))))
					
					# Kage proportion of habitat available (add one to age class) - equation from AB.
					#Kage_hab = np.exp(packpar1*(Ageclass+1))
					Kage_hab = Kij_proportion[np.where(iage==AgeClass_reverse)][0]
					
					# Class count for redsitributing Khab
					classcount = len(np.where(AgeClass_reverse<=Ageclass)[0])
					#classcount = Ageclass+1
					Kage_hab_adj_inc = Kage_hab_adj_inc + (Kage_hab_adj/classcount)
					#Kage_hab_adj_inc = Kage_hab_adj_inc[0]
					
					# Adjust the Kage_hab
					Kage_hab = Kage_hab + Kage_hab_adj_inc
					
					# Kage for this population 
					#Kage = Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))
					#Added this check because this expression can result in large 'inf' values
					if np.isnan(Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))):
						Kage = 1000000000 #This should safely be above any Nage value
					else:
						Kage = Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))
					# Get Kused
					if Nage <= Kage:
						# Grab all the ages - they all survive
						Kused = len(Nage_index)
						Nage_samp_ind.append(list(Nage_index))
						# Tracking numbers ------------------------------
						PackingDeathsAge[gen][indexforAgeclass].append(0)
						Track_YYSelectionPackDeaths[gen][isub].append(0)
						Track_WildSelectionPackDeaths[gen][isub].append(0)
					else: # Some die
						# Get Kused - this is the number that is kept
						Kused = int(round(Kage))
												
						# CDEVOLVE PACKING OPTION
						if timecdevolve == 'packing' and cdevolveans.split('_')[0] == 'Hindex':
							
							Fitness = float(cdevolveans.split('_')[2].split(':')[1])
							X = Nage-Kused # Total number of deaths
							hindex_arr = SubpopIN_arr[Nage_index]['hindex'] # Hindex of all the ind here
							N_w = len(np.where(hindex_arr > 0)[0]) # Number of wild type
							N_yy = len(np.where(hindex_arr == 0.)[0]) # Number of YY type
							Deaths_yy = int(round((X * N_yy) / ((Fitness * N_w) + N_yy)))
							Deaths_w = X - Deaths_yy
							
							# Here, recalculate deaths for the cases in which more deaths calculated than N showed up
							if Deaths_yy > N_yy:
								# Adjust the actuald deaths to wild type
								Deaths_yy = N_yy
								Deaths_w = X - Deaths_yy								
							
							# Then keep the propotional amount of each type
							if Kused == 0: # No one survives given packing parameters
								# Grab a random number of Kused from Nage sample. 
								Nage_samp_ind.append(np.random.choice(Nage_index,Kused,replace=False).tolist())
								# Tracking numbers -------------------------------------
								Track_YYSelectionPackDeaths[gen][isub].append(N_yy)
								Track_WildSelectionPackDeaths[gen][isub].append(N_w)
							
							else: # Some individuals survive given packing parameters
								# As with Nage_index - index for where the YYs are
								Nage_index_yy = Nage_index[np.where(hindex_arr == 0.)[0]]
								# Then get the Kused for just the yy class
								Kused_yy = N_yy - Deaths_yy
																
								# Same for wild type Nage_index - all hindex greater than 0 at this point
								Nage_index_wild = Nage_index[np.where(hindex_arr > 0)[0]]
								#The get Kused for just the wild class
								Kused_wild = N_w - Deaths_w								
								
								# Then sample those that survive and stay - YYs
								Nage_samp_ind.append(np.random.choice(Nage_index_yy,Kused_yy,replace=False).tolist())								
								# Then sample those that survive and stay - Wild type
								Nage_samp_ind.append(np.random.choice(Nage_index_wild,Kused_wild,replace=False).tolist())								
																
								# Tracking numbers ------------------------------------------
								Track_YYSelectionPackDeaths[gen][isub].append(Deaths_yy)
								Track_WildSelectionPackDeaths[gen][isub].append(Deaths_w)					
							
						# CDEVOLVE Packing option off
						else:
							# Grab a random number of Kused from Nage sample. 
							Nage_samp_ind.append(np.random.choice(Nage_index,Kused,replace=False).tolist())
							# Tracking numbers ----------------------------------
							Track_YYSelectionPackDeaths[gen][isub].append(0)
							Track_WildSelectionPackDeaths[gen][isub].append(0)
						
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
			Track_YYSelectionPackDeaths[gen][isub] = sum(Track_YYSelectionPackDeaths[gen][isub])
			Track_WildSelectionPackDeaths[gen][isub] = sum(Track_WildSelectionPackDeaths[gen][isub])
			
			# Track class size - get age adjusted again. 		
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']
			
			# Tracking age N
			for iage in range(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))	
	# -------------------------------------------------
	# Packing without reallocation of space is selected
	# -------------------------------------------------
	elif packans == 'packing_N':

		# ------------------------------------------
		# Get other species Ns from all Patches here
		# ------------------------------------------
		if implementcomp == 'Back':
			Nself_pop = [len(SubpopIN_keep[x]) for x in range(0,len(SubpopIN_keep))]
			# Ignore queus if only one species
			if len(XQs) > 0:
				# Loop through queue spots, Putting Nself_pop for grabbing for other species.
				for ispecies in range(len(XQs[spcNO])):
					XQs[spcNO][ispecies].put(Nself_pop) 
			Nother_pop = []
			popcount = 0
			# Ignore queues if only one species
			if len(XQs) > 0:
				for ispecies in range(len(XQs[spcNO])):			
					if ispecies != spcNO: 
						Nother_pop.append(XQs[ispecies][spcNO].get(block = True))
						# Error check here for equal patches
						if len(Nother_pop[popcount]) != len(Nself_pop):
							print("Species systems must have the same number of patches.")
							sys.exit(-1)
						popcount = popcount + 1	
		
		# Loop through each subpopulation
		for isub in range(len(K)):
					
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# -----------------------------------
			# Get countages - for 'age' adjusted
			# -----------------------------------
			# Switch here for size or age control
			if sizecall == 'size': # Careful here and use first ClassVars
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_arr[sizecall])
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)
			else:
				# Count up each uniages
				countages = count_unique(SubpopIN_arr['age'])
			
			# ------------------------
			# Get packing parameters
			# ------------------------			
			# K,N for this population
			Kpop = K[isub]
			Npop = sum(countages[1])
			
			# Overwrite Kpop with K adjusted from DoEmigration that was updated from competition - careful of indexing Tracking variable is patch + 1, first value is total
			if implementcomp == 'Out':
				Kpop = Track_KadjEmi[gen][isub+1]
						
			if implementcomp == 'Back':
				#if multiprocessing.current_process().name == "S1":
				#	pdb.set_trace()	
				# -----------------
				# Adjustment for K if more than one species
				# -----------------
				#Ignore K adjustment if only one species
				if len(XQs) > 0:
					alphas = comp_coef[isub].split(';') # Extracting alphas
					if gen >= startcomp: # Check timing 
						tempspeciesSum = [] # Create sum list for  >= 2 species
						for ispecies in range(len(Nother_pop)):
							tempspeciesSum.append(Nother_pop[ispecies][isub]*float(alphas[ispecies]))
						Kpop = int(Kpop - sum(tempspeciesSum)) # Adjustment or K		
			Kadj_track[gen].append(Kpop) # For Tracking
			
			# Continue packing
			if Npop == 0 or Kpop <= 0:
				# Append all information to temp SubpopKeep variable
				Nage_samp_ind = np.arange(0)					
			else:		
				# Check here on numbers
				if Kpop == Npop:
					Kscale = 1.0
				else:
					# Get K_scaling for population - Equation 6 - R0
					Kscale = np.log(float(Kpop)/Npop) / (1. - (Npop/float(Kpop)))
				
				# -------------------------
				# Age loop
				# -------------------------
				# Loop through each age class - here recursive, start with largest class first
				Nage_samp_ind = [] 
				Kage_hab_adj = 0.
				Kage_hab_adj_inc = 0.
				#Kij_proportion = np.exp(packpar1*(countages[0]+1)) # Equ 8
				Kij_proportion = np.exp(packpar1*(np.asarray(list(range(len(size_mean[0][0]))))+1))
				Kij_proportion = Kij_proportion/sum(Kij_proportion) # Equ 8 rescaled
				Kij_proportion = np.flipud(Kij_proportion) # Flip to match next age loop (old to young)
				AgeClass_reverse = np.flipud(np.asarray(list(range(len(size_mean[0][0])))))
				#for iage in -np.sort(-countages[0]):
				for iage in AgeClass_reverse:
					
					# Age class
					Ageclass = iage
					
					# Special case when age is greater than last age class only used for indexing now
					#if Ageclass > len(size_mean[0][0])-1:
					if iage == AgeClass_reverse[0] and sizecall == 'age':
						indexforAgeclass = len(size_mean[0][0]) - 1
					else:
						indexforAgeclass = Ageclass
					
					
					# N for this age(s) coming in
					if len(np.where(countages[0]==Ageclass)[0]) == 0:
						Nage = 0
					else:
						if iage == AgeClass_reverse[0] and sizecall == 'age': # Special case when age is greater than last age
							Nage = sum(countages[1][np.where(countages[0]>=Ageclass)[0]])
						else:
							Nage = countages[1][np.where(countages[0]==Ageclass)[0][0]]
					
					# Get index of these ages
					if sizecall == 'size': # Use the adjusted age classes
						Nage_index = np.where(age_adjusted==Ageclass)[0]
					else:
						if iage == AgeClass_reverse[0] and sizecall == 'age': # Special case when age is greater than last age
							Nage_index = np.where(SubpopIN_arr['age']>=Ageclass)[0]
						else:
							Nage_index = np.where(SubpopIN_arr['age']==Ageclass)[0]
					
					# Get Age_scaling for this pop's age class (*Stage R) - Equation 7
					Agescale = np.exp(Kscale * (1. - (Nage / float(Kpop))))
					
					# Kage proportion of habitat available (add one to age class) - equation from AB.
					#Kage_hab = np.exp(packpar1*(Ageclass+1))
					Kage_hab = Kij_proportion[np.where(iage==AgeClass_reverse)][0]
					
					# Class count for redsitributing Khab
					classcount = len(np.where(AgeClass_reverse<=Ageclass)[0])
					#classcount = Ageclass+1
					Kage_hab_adj_inc = Kage_hab_adj_inc + (Kage_hab_adj/classcount)
					#Kage_hab_adj_inc = Kage_hab_adj_inc[0]
					
					# Adjust the Kage_hab
					Kage_hab = Kage_hab + Kage_hab_adj_inc
					
					# Kage for this population 
					#Kage = Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))
					#Added this check because this expression can result in large 'inf' values
					if np.isnan(Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))):
						Kage = 1000000000 #This should safely be above any Nage value
					else:
						# Equation 9
						Kage = Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))
					# Get Kused
					if Nage <= Kage:
						# Grab all the ages - they all survive
						Kused = len(Nage_index)
						Nage_samp_ind.append(list(Nage_index))
						# Tracking numbers ------------------------------
						PackingDeathsAge[gen][indexforAgeclass].append(0)
						Track_YYSelectionPackDeaths[gen][isub].append(0)
						Track_WildSelectionPackDeaths[gen][isub].append(0)
					else: # Some die
						# Get Kused - this is the number that is kept
						Kused = int(round(Kage))
												
						# CDEVOLVE PACKING OPTION
						if timecdevolve == 'packing' and cdevolveans.split('_')[0] == 'Hindex':
							
							Fitness = float(cdevolveans.split('_')[2].split(':')[1])
							X = Nage-Kused # Total number of deaths
							hindex_arr = SubpopIN_arr[Nage_index]['hindex'] # Hindex of all the ind here
							N_w = len(np.where(hindex_arr > 0)[0]) # Number of wild type
							N_yy = len(np.where(hindex_arr == 0.)[0]) # Number of YY type
							Deaths_yy = int(round((X * N_yy) / ((Fitness * N_w) + N_yy)))
							Deaths_w = X - Deaths_yy
							
							# Here, recalculate deaths for the cases in which more deaths calculated than N showed up
							if Deaths_yy > N_yy:
								# Adjust the actuald deaths to wild type
								Deaths_yy = N_yy
								Deaths_w = X - Deaths_yy								
							
							# Then keep the propotional amount of each type
							if Kused == 0: # No one survives given packing parameters
								# Grab a random number of Kused from Nage sample. 
								Nage_samp_ind.append(np.random.choice(Nage_index,Kused,replace=False).tolist())
								# Tracking numbers -------------------------------------
								Track_YYSelectionPackDeaths[gen][isub].append(N_yy)
								Track_WildSelectionPackDeaths[gen][isub].append(N_w)
							
							else: # Some individuals survive given packing parameters
								# As with Nage_index - index for where the YYs are
								Nage_index_yy = Nage_index[np.where(hindex_arr == 0.)[0]]
								# Then get the Kused for just the yy class
								Kused_yy = N_yy - Deaths_yy
																
								# Same for wild type Nage_index - all hindex greater than 0 at this point
								Nage_index_wild = Nage_index[np.where(hindex_arr > 0)[0]]
								#The get Kused for just the wild class
								Kused_wild = N_w - Deaths_w								
								
								# Then sample those that survive and stay - YYs
								Nage_samp_ind.append(np.random.choice(Nage_index_yy,Kused_yy,replace=False).tolist())								
								# Then sample those that survive and stay - Wild type
								Nage_samp_ind.append(np.random.choice(Nage_index_wild,Kused_wild,replace=False).tolist())								
																
								# Tracking numbers ------------------------------------------
								Track_YYSelectionPackDeaths[gen][isub].append(Deaths_yy)
								Track_WildSelectionPackDeaths[gen][isub].append(Deaths_w)					
							
						# CDEVOLVE Packing option off
						else:
							# Grab a random number of Kused from Nage sample. 
							Nage_samp_ind.append(np.random.choice(Nage_index,Kused,replace=False).tolist())
							# Tracking numbers ----------------------------------
							Track_YYSelectionPackDeaths[gen][isub].append(0)
							Track_WildSelectionPackDeaths[gen][isub].append(0)
						
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
			Track_YYSelectionPackDeaths[gen][isub] = sum(Track_YYSelectionPackDeaths[gen][isub])
			Track_WildSelectionPackDeaths[gen][isub] = sum(Track_WildSelectionPackDeaths[gen][isub])
			
			# Track class size - get age adjusted again. 		
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']
			
			# Tracking age N
			for iage in range(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))	
	# -------------------------------------------------------------
	# Packing option 1 - Extra space allocated to only 1 size class
	# -------------------------------------------------------------
	elif packans == 'packing_1':
		
		# ------------------------------------------
		# Get other species Ns from all Patches here
		# ------------------------------------------
		if implementcomp == 'Back':
			Nself_pop = [len(SubpopIN_keep[x]) for x in range(0,len(SubpopIN_keep))]
			# Ignore queus if only one species
			if len(XQs) > 0:
				# Loop through queue spots, Putting Nself_pop for grabbing for other species.
				for ispecies in range(len(XQs[spcNO])):
					XQs[spcNO][ispecies].put(Nself_pop) 
			Nother_pop = []
			popcount = 0
			# Ignore queues if only one species
			if len(XQs) > 0:
				for ispecies in range(len(XQs[spcNO])):			
					if ispecies != spcNO: 
						Nother_pop.append(XQs[ispecies][spcNO].get(block = True))
						# Error check here for equal patches
						if len(Nother_pop[popcount]) != len(Nself_pop):
							print("Species systems much have the same number of patches.")
							sys.exit(-1)
					popcount = popcount + 1	
		# Loop through each subpopulation
		for isub in range(len(K)):
					
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# -----------------------------------
			# Get countages - for 'age' adjusted
			# -----------------------------------
			# Switch here for size or age control
			if sizecall == 'size': # Careful here and use first ClassVars
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_arr[sizecall])
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)
			else:
				# Count up each uniages
				countages = count_unique(SubpopIN_arr['age'])
			
			# ------------------------
			# Get packing parameters
			# ------------------------			
			# K,N for this population
			Kpop = K[isub]
			Npop = sum(countages[1])
			
			# Overwrite Kpop with K adjusted from DoEmigration that was updated from competition - careful of indexing Tracking variable is patch + 1, first value is total
			if implementcomp == 'Out':
				Kpop = Track_KadjEmi[gen][isub+1]
			
			if implementcomp == 'Back':
				#if multiprocessing.current_process().name == "S1":
				#	pdb.set_trace()	
				# -----------------
				# Adjustment for K if more than one species
				# -----------------
				#Ignore K adjustment if only one species
				if len(XQs) > 0:
					alphas = comp_coef[isub].split(';') # Extracting alphas
					if gen >= startcomp: # Check timing 
						tempspeciesSum = [] # Create sum list for  >= 2 species
						for ispecies in range(len(Nother_pop)):
							tempspeciesSum.append(Nother_pop[ispecies][isub]*float(alphas[ispecies]))
						Kpop = int(Kpop - sum(tempspeciesSum)) # Adjustment or K		
			Kadj_track[gen].append(Kpop) # For Tracking
			
			# Continue packing
			if Npop == 0 or Kpop == 0:
				# Append all information to temp SubpopKeep variable
				Nage_samp_ind = np.arange(0)					
			else:		
				# Check here on numbers
				if Kpop == Npop:
					Kscale = 1.0
				else:
					# Get K_scaling for population
					Kscale = np.log(float(Kpop)/Npop) / (1. - (Npop/float(Kpop))) # Equation 6
				
				# -------------------------
				# Age loop
				# -------------------------
				# Loop through each age class - here recursive, start with largest class first
				Nage_samp_ind = [] 
				Kage_hab_adj = 0.
				Kage_hab_adj_inc = 0.				
				Kij_proportion = np.exp(packpar1*(np.asarray(list(range(len(size_mean[0][0]))))+1))
				Kij_proportion = Kij_proportion/sum(Kij_proportion) # Equ 8 rescaled
				Kij_proportion = np.flipud(Kij_proportion) # Flip to match next age loop (old to young)
				AgeClass_reverse = np.flipud(np.asarray(list(range(len(size_mean[0][0])))))
				######################## Start edits for packing 1 #####################
				#Goal here is to only distribute available space to one class below the given size class
				AgeClass = np.asarray(list(range(len(size_mean[0][0])))) #Age classes not reversed
				Agescale = []
				Kage = []
				Kage_hab_adj = []
				Kage_hab = []
				Kused = []
				Kused_new = [] # This is for the second calculation of Kused after the adjusted Kage_hab is calculated
				Kage_new = [] # This is the second calculation of Kage after the adjusted Kage_hab is calculated
				
				for iage in AgeClass:					
					#Nage = countages[1][iage]
					# Special case when age is greater than last age class only used for indexing now
					if iage > len(size_mean[0][0])-1:
						indexforAgeclass = len(size_mean[0][0]) - 1 
					else:
						indexforAgeclass = iage
					
					# N for this age coming in
					if len(np.where(countages[0]==iage)[0]) == 0:
						Nage = 0
					else:
						Nage = countages[1][np.where(countages[0]==iage)[0][0]]
					if sizecall == 'size': # Use the adjusted age classes
						Nage_index = np.where(age_adjusted==iage)[0]
					else:
						Nage_index = np.where(SubpopIN_arr['age']==iage)[0]
					# Need to calculate Stage R (agescale) for each age class (equation 7)
					Agescale.append(np.exp(Kscale * (1.-(Nage/float(Kpop))))) 
					
					# Here calculate Kage for a given size class based on equation 9
					# Check for inf value in case of small N and large AgeScale value
					if np.isnan(Nage * np.exp(Agescale[iage] * (1.-(Nage/(Kpop * Kij_proportion[(iage+1)*-1]))))):
						Kage.append(1000000000)
					else:
						Kage.append(Nage * np.exp(Agescale[iage] * (1.-(Nage/(Kpop * Kij_proportion[(iage+1)*-1]))))) 
					
					# Now calculate Kused for each age class so we can calculate Kage_hab_adj
					if Nage <= Kage[iage]: 
						Kused.append(len(Nage_index))
						# Don't append to Nage_samp_ind until second Kused?
					else:
						Kused.append(Kage[iage])						
						
					# Now calculate Kage_hab_adj - here just use Kij_proportion instead of Kage_hab_adj, equation 10
					if Kage[iage] <= 0:
						Kage_hab_adj.append(0)
					else:
						Kage_hab_adj.append(Kij_proportion[(iage+1)*-1] - ((Kused[iage]/Kage[iage])*Kij_proportion[(iage+1)*-1])) # Represents available space in each age class
					#Kage_hab_adj_inc not necessary here because the available space is not divided among classes							
				
				#Reset Kused and Kage now that Kage_hab_adj is already calculated, will calculate new Kused based on new Kage_hab
				for iage in AgeClass:					
					# Special case when age is greater than last age class only used for indexing now
					if iage > len(size_mean[0][0])-1:
						indexforAgeclass = len(size_mean[0][0]) - 1
					else:
						indexforAgeclass = iage
					
					# N for this age coming in
					if len(np.where(countages[0]==iage)[0]) == 0:
						Nage = 0
					else:
						Nage = countages[1][np.where(countages[0]==iage)[0][0]]
					if sizecall == 'size': # Use the adjusted age classes
						Nage_index = np.where(age_adjusted==iage)[0]
					else:
						Nage_index = np.where(SubpopIN_arr['age']==iage)[0]
					# Recalculate Kij_proportion for each age class based on available space from Kage_hab_adj
					if iage == AgeClass[-1]: #If largest class, needs to be fixed to account for old fish.
						Kage_hab.append(Kij_proportion[(iage+1)*-1]) # Just your own Kij_proportion value
					else:
						Kage_hab.append(Kij_proportion[(iage+1)*-1] + Kage_hab_adj[iage+1]) #your Kij_proportion + available space from class above you
					# Now recalculate Kage with the new Kage_hab, and calculate deaths from Kused, equation 9
					# Check for inf value in case of small N and large AgeScale value
					if np.isnan(Nage * np.exp(Agescale[iage] * (1.-(Nage/(Kpop*Kage_hab[iage]))))):
						Kage_new.append(1000000000)
					else:
						Kage_new.append(Nage * np.exp(Agescale[iage] * (1.-(Nage/(Kpop*Kage_hab[iage]))))) 
					
					
										
					#Now calculate Kused again to determine # of deaths
					if Nage <= Kage_new[iage]:
						Kused_new.append(Nage)
						# Add tracking stuff here
						Nage_samp_ind.append(list(Nage_index))
						# Tracking numbers ------------------------------
						PackingDeathsAge[gen][indexforAgeclass].append(0)
						Track_YYSelectionPackDeaths[gen][isub].append(0)
						Track_WildSelectionPackDeaths[gen][isub].append(0)
					else: # some die
						Kused_new.append(int(round(Kage_new[iage])))
						if timecdevolve == 'packing' and cdevolveans.split('_')[0] == 'Hindex':
							Fitness = float(cdevolveans.split('_')[2].split(':')[1])
							X = Nage-Kused_new[iage] # Total number of deaths
							hindex_arr = SubpopIN_arr[Nage_index]['hindex'] # Hindex of all the ind here
							N_w = len(np.where(hindex_arr > 0)[0]) # Number of wild type
							N_yy = len(np.where(hindex_arr == 0.)[0]) # Number of YY type
							Deaths_yy = int(round((X * N_yy) / ((Fitness * N_w) + N_yy)))
							Deaths_w = X - Deaths_yy
							
							# Here, recalculate deaths for the cases in which more deaths calculated than N showed up
							if Deaths_yy > N_yy:
								# Adjust the actuald deaths to wild type
								Deaths_yy = N_yy
								Deaths_w = X - Deaths_yy								
							
							# Then keep the propotional amount of each type
							if Kused_new[iage] == 0: # No one survives given packing parameters
								# Grab a random number of Kused from Nage sample. 
								Nage_samp_ind.append(np.random.choice(Nage_index,Kused_new[iage],replace=False).tolist())
								# Tracking numbers -------------------------------------
								Track_YYSelectionPackDeaths[gen][isub].append(N_yy)
								Track_WildSelectionPackDeaths[gen][isub].append(N_w)
							
							else: # Some individuals survive given packing parameters
								# As with Nage_index - index for where the YYs are
								Nage_index_yy = Nage_index[np.where(hindex_arr == 0.)[0]]
								# Then get the Kused for just the yy class
								Kused_yy = N_yy - Deaths_yy
																
								# Same for wild type Nage_index - all hindex greater than 0 at this point
								Nage_index_wild = Nage_index[np.where(hindex_arr > 0)[0]]
								#The get Kused for just the wild class
								Kused_wild = N_w - Deaths_w								
								
								# Then sample those that survive and stay - YYs
								Nage_samp_ind.append(np.random.choice(Nage_index_yy,Kused_yy,replace=False).tolist())								
								# Then sample those that survive and stay - Wild type
								Nage_samp_ind.append(np.random.choice(Nage_index_wild,Kused_wild,replace=False).tolist())								
																
								# Tracking numbers ------------------------------------------
								Track_YYSelectionPackDeaths[gen][isub].append(Deaths_yy)
								Track_WildSelectionPackDeaths[gen][isub].append(Deaths_w)				
						# CDEVOLVE Packing option off
						else:
							# Grab a random number of Kused from Nage sample. 
							#pdb.set_trace()
							Nage_samp_ind.append(np.random.choice(Nage_index,Kused_new[iage],replace=False).tolist())
							# Tracking numbers ----------------------------------
							Track_YYSelectionPackDeaths[gen][isub].append(0)
							Track_WildSelectionPackDeaths[gen][isub].append(0)
						
						PackingDeathsAge[gen][indexforAgeclass].append(Nage-Kused_new[iage])
						#Then all the tracking stuff can go here						
								
			# Clean up index
			#pdb.set_trace()
			Nage_samp_ind = sum(Nage_samp_ind,[])
				
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepK.append(SubpopIN_arr[Nage_samp_ind])
						
			# Store new N - it can be possible to be less than K - its OK - rounding error
			Population[gen].append(len(SubpopIN_keepK[isub]))			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
			PackingDeaths[gen][isub] = Npop - len(Nage_samp_ind)
			Track_YYSelectionPackDeaths[gen][isub] = sum(Track_YYSelectionPackDeaths[gen][isub])
			Track_WildSelectionPackDeaths[gen][isub] = sum(Track_WildSelectionPackDeaths[gen][isub])
			
			# Track class size - get age adjusted again. 		
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']
			
			# Tracking age N
			for iage in range(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))		
	
	# -------------------
	# Anadromy is selected
	# -------------------
	elif packans == 'anadromy':
		if sizecall == "size": # For now, only allow for age-based option
			print("Currently, only age-based simulations are allowed for the 'anadromy' popmodel")
			sys.exit(-1)
		# ------------------------------------------
		# Get other species Ns from all Patches here
		# ------------------------------------------
		if implementcomp == 'Back':
			Nself_pop = [len(SubpopIN_keep[x]) for x in range(0,len(SubpopIN_keep))]
			# Ignore queus if only one species
			if len(XQs) > 0:
				print("Single species models only for anadromy popmodel option for now")
				sys.exit(-1)
				# Loop through queue spots, Putting Nself_pop for grabbing for other species.
				for ispecies in range(len(XQs[spcNO])):
					XQs[spcNO][ispecies].put(Nself_pop) 
			Nother_pop = []
			popcount = 0
			# Ignore queues if only one species
			if len(XQs) > 0:
				for ispecies in range(len(XQs[spcNO])):			
					if ispecies != spcNO: 
						Nother_pop.append(XQs[ispecies][spcNO].get(block = True))
						# Error check here for equal patches
						if len(Nother_pop[popcount]) != len(Nself_pop):
							print("Species systems must have the same number of patches.")
							sys.exit(-1)
						popcount = popcount + 1	
		
		# Loop through each subpopulation
		for isub in range(len(K)):
					
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# -----------------------------------
			# Get countages - for 'age' adjusted
			# -----------------------------------
			# Switch here for size or age control
			if sizecall == 'size': # Careful here and use first ClassVars
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_arr[sizecall])
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)
			else:
				# Count up each uniages
				countages = count_unique(SubpopIN_arr['age'])
			
			# ------------------------
			# Get packing parameters
			# ------------------------			
			# K,N for this population
			Kpop = K[isub]
			# Add in to popvars as variable eventually e.g., anadromy_1
			packAge = 1 # Age of oldest individuals to be packed
			Npop = sum(countages[1])
			# Only ages 0 and 1 contriubte to Npop for packing
			Npop = sum(countages[1][np.where(countages[0]<=packAge)])
			
			# Overwrite Kpop with K adjusted from DoEmigration that was updated from competition - careful of indexing Tracking variable is patch + 1, first value is total
			if implementcomp == 'Out':
				Kpop = Track_KadjEmi[gen][isub+1]
						
			if implementcomp == 'Back':
				#if multiprocessing.current_process().name == "S1":
				#	pdb.set_trace()	
				# -----------------
				# Adjustment for K if more than one species
				# -----------------
				#Ignore K adjustment if only one species
				if len(XQs) > 0:
					alphas = comp_coef[isub].split(';') # Extracting alphas
					if gen >= startcomp: # Check timing 
						tempspeciesSum = [] # Create sum list for  >= 2 species
						for ispecies in range(len(Nother_pop)):
							tempspeciesSum.append(Nother_pop[ispecies][isub]*float(alphas[ispecies]))
						Kpop = int(Kpop - sum(tempspeciesSum)) # Adjustment or K		
			Kadj_track[gen].append(Kpop) # For Tracking
			
			# Continue packing
			if Npop == 0 or Kpop <= 0:
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
				# Age loop - anadromy option only packs age classes 0 and 1
				# -------------------------
				# Loop through each age class - here recursive, start with largest class first (Age class 1)
				Nage_samp_ind = [] 
				Kage_hab_adj = 0.
				Kage_hab_adj_inc = 0.
				#Kij_proportion = np.exp(packpar1*(countages[0]+1)) # Equ 8
				Kij_proportion = np.exp(packpar1*(np.asarray(list(range(len(size_mean[0][0]))))+1)[0:packAge+1])
				Kij_proportion = Kij_proportion/sum(Kij_proportion) # Equ 8 rescaled
				Kij_proportion = np.flipud(Kij_proportion) # Flip to match next age loop (old to young)
				AgeClass_reverse = np.flipud(np.asarray(list(range(len(size_mean[0][0])))))[-(packAge+1):]
				# Loop through ages 1 and 0
				for iage in AgeClass_reverse:
					
					# Age class
					Ageclass = iage
					
					# Special case when age is greater than last age class only used for indexing now
					# Not needed when only working with first two age classes
					'''if Ageclass > len(size_mean[0][0])-1:
						indexforAgeclass = len(size_mean[0][0]) - 1
					else:
						indexforAgeclass = Ageclass'''
					indexforAgeclass = Ageclass
					# N for this age coming in
					if len(np.where(countages[0]==Ageclass)[0]) == 0:
						Nage = 0
					else: # Does this need a check for if age is higher than highest class?
						Nage = countages[1][np.where(countages[0]==Ageclass)[0][0]]
					if sizecall == 'size': # Use the adjusted age classes
						Nage_index = np.where(age_adjusted==Ageclass)[0]
					else: # Same here, need a check for age higher than last age class?
						Nage_index = np.where(SubpopIN_arr['age']==Ageclass)[0]
					
					# Get Age_scaling for this pop's age class (*Stage R)
					Agescale = np.exp(Kscale * (1. - (Nage / float(Kpop))))
					
					# Kage proportion of habitat available (add one to age class) - equation from AB.
					#Kage_hab = np.exp(packpar1*(Ageclass+1))
					Kage_hab = Kij_proportion[np.where(iage==AgeClass_reverse)][0]
					
					# Class count for redsitributing Khab
					classcount = len(np.where(AgeClass_reverse<=Ageclass)[0])
					#classcount = Ageclass+1
					Kage_hab_adj_inc = Kage_hab_adj_inc + (Kage_hab_adj/classcount)
					#Kage_hab_adj_inc = Kage_hab_adj_inc[0]
					
					# Adjust the Kage_hab
					Kage_hab = Kage_hab + Kage_hab_adj_inc
					
					# Kage for this population 
					#Kage = Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))
					#Added this check because this expression can result in large 'inf' values
					if np.isnan(Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))):
						Kage = 1000000000 #This should safely be above any Nage value
					else:
						Kage = Nage * np.exp(Agescale * (1. - (Nage / (Kpop * Kage_hab))))
					# Get Kused
					if Nage <= Kage:
						# Grab all the ages - they all survive
						Kused = len(Nage_index)
						Nage_samp_ind.append(list(Nage_index))
						# Tracking numbers ------------------------------
						PackingDeathsAge[gen][indexforAgeclass].append(0)
						Track_YYSelectionPackDeaths[gen][isub].append(0)
						Track_WildSelectionPackDeaths[gen][isub].append(0)
					else: # Some die
						# Get Kused - this is the number that is kept
						Kused = int(round(Kage))
												
						# CDEVOLVE PACKING OPTION
						if timecdevolve == 'packing' and cdevolveans.split('_')[0] == 'Hindex':
							print("cdevolve options not ready for anadromy popmodel. Change cdevolve answer to 'N'.")
							sys.exit(-1)
							Fitness = float(cdevolveans.split('_')[2].split(':')[1])
							X = Nage-Kused # Total number of deaths
							hindex_arr = SubpopIN_arr[Nage_index]['hindex'] # Hindex of all the ind here
							N_w = len(np.where(hindex_arr > 0)[0]) # Number of wild type
							N_yy = len(np.where(hindex_arr == 0.)[0]) # Number of YY type
							Deaths_yy = int(round((X * N_yy) / ((Fitness * N_w) + N_yy)))
							Deaths_w = X - Deaths_yy
							
							# Here, recalculate deaths for the cases in which more deaths calculated than N showed up
							if Deaths_yy > N_yy:
								# Adjust the actuald deaths to wild type
								Deaths_yy = N_yy
								Deaths_w = X - Deaths_yy								
							
							# Then keep the propotional amount of each type
							if Kused == 0: # No one survives given packing parameters
								# Grab a random number of Kused from Nage sample. 
								Nage_samp_ind.append(np.random.choice(Nage_index,Kused,replace=False).tolist())
								# Tracking numbers -------------------------------------
								Track_YYSelectionPackDeaths[gen][isub].append(N_yy)
								Track_WildSelectionPackDeaths[gen][isub].append(N_w)
							
							else: # Some individuals survive given packing parameters
								# As with Nage_index - index for where the YYs are
								Nage_index_yy = Nage_index[np.where(hindex_arr == 0.)[0]]
								# Then get the Kused for just the yy class
								Kused_yy = N_yy - Deaths_yy
																
								# Same for wild type Nage_index - all hindex greater than 0 at this point
								Nage_index_wild = Nage_index[np.where(hindex_arr > 0)[0]]
								#The get Kused for just the wild class
								Kused_wild = N_w - Deaths_w								
								
								# Then sample those that survive and stay - YYs
								Nage_samp_ind.append(np.random.choice(Nage_index_yy,Kused_yy,replace=False).tolist())								
								# Then sample those that survive and stay - Wild type
								Nage_samp_ind.append(np.random.choice(Nage_index_wild,Kused_wild,replace=False).tolist())								
																
								# Tracking numbers ------------------------------------------
								Track_YYSelectionPackDeaths[gen][isub].append(Deaths_yy)
								Track_WildSelectionPackDeaths[gen][isub].append(Deaths_w)					
							
						# CDEVOLVE Packing option off
						else:
							# Grab a random number of Kused from Nage sample. 
							Nage_samp_ind.append(np.random.choice(Nage_index,Kused,replace=False).tolist())
							# Tracking numbers ----------------------------------
							Track_YYSelectionPackDeaths[gen][isub].append(0)
							Track_WildSelectionPackDeaths[gen][isub].append(0)
						
						PackingDeathsAge[gen][indexforAgeclass].append(Nage-Kused)
					
					# The adjust the habitat or reallocation for next class
					if Kage == 0:
						Kage_hab_adj = Kage_hab
					else:
						Kage_hab_adj = Kage_hab - (Kused * Kage_hab / Kage)
					
			# Clean up index, double check from here with Erin
			Nage_samp_ind = sum(Nage_samp_ind,[])
			index2plus = np.where(SubpopIN_arr['age'] > 1)[0].tolist()
			# Combine age 1s and 0s that survived with age 2+
			Nage_samp_ind = Nage_samp_ind + index2plus
				
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepK.append(SubpopIN_arr[Nage_samp_ind])
						
			# Store new N - it can be possible to be less than K - its OK - rounding error
			Population[gen].append(len(SubpopIN_keepK[isub]))			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
			PackingDeaths[gen][isub] = sum(countages[1]) - len(Nage_samp_ind)
			Track_YYSelectionPackDeaths[gen][isub] = sum(Track_YYSelectionPackDeaths[gen][isub])
			Track_WildSelectionPackDeaths[gen][isub] = sum(Track_WildSelectionPackDeaths[gen][isub])
			
			# Track class size - get age adjusted again. 		
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']
			
			# Tracking age N
			for iage in range(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))	
	# ---------------------
	# Packing is turned off
	# ---------------------
	elif packans == 'N':
	
		# ------------------------------------------
		# Get other species Ns from all Patches here
		# ------------------------------------------
		if implementcomp == 'Back':
			Nself_pop = [len(SubpopIN_keep[x]) for x in range(0,len(SubpopIN_keep))]
			# Ignore queus if only one species
			if len(XQs) > 0:
				# Loop through queue spots, Putting Nself_pop for grabbing for other species.
				for ispecies in range(len(XQs[spcNO])):
					XQs[spcNO][ispecies].put(Nself_pop) 
			Nother_pop = []
			popcount = 0
			# Ignore queues if only one species
			if len(XQs) > 0:
				for ispecies in range(len(XQs[spcNO])):			
					if ispecies != spcNO: 
						Nother_pop.append(XQs[ispecies][spcNO].get(block = True))
						# Error check here for equal patches
						if len(Nother_pop[popcount]) != len(Nself_pop):
							print("Species systems must have the same number of patches.")
							sys.exit(-1)
						popcount = popcount + 1	
		
		for isub in range(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# K,N for this population
			Kpop = K[isub]
			Npop = len(SubpopIN_arr)
			
			# Overwrite Kpop with K adjusted from DoEmigration that was updated from competition - careful of indexing Tracking variable is patch + 1, first value is total
			if implementcomp == 'Out':
				Kpop = Track_KadjEmi[gen][isub+1]
			
			if implementcomp == 'Back':
				# -----------------
				# Adjustment for K if more than one species
				# -----------------
				#Ignore K adjustment if only one species
				'''if len(XQs) > 0:
					alphas = comp_coef[isub].split(';') # Extracting alphas
					if gen >= startcomp: # Check timing 
						tempspeciesSum = [] # Create sum list for  >= 2 species
						for ispecies in range(len(Nother_pop)):
							tempspeciesSum.append(Nother_pop[ispecies][isub]*float(alphas[ispecies]))
						Kpop = int(Kpop - sum(tempspeciesSum)) # Adjustment or K'''
			Kadj_track[gen].append(Kpop) # For Tracking
			
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
					Nage_samp_ind = np.random.choice(np.arange(len(SubpopIN_arr)),int(Kpop),replace=False).tolist()
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
			for iage in range(len(PopulationAge[gen])):
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
			Track_YYSelectionPackDeaths[gen][isub] = 0 # Zero for now
			Track_WildSelectionPackDeaths[gen][isub] = 0 # Zero fro now
	
	# ---------------------
	# Packing is turned off - this is a special case in which all inds are kept.
	# ---------------------
	elif packans == 'N_keepeggs':
		
		# ------------------------------------------
		# Get other species Ns from all Patches here
		# ------------------------------------------
		if implementcomp == 'Back':
			Nself_pop = [len(SubpopIN_keep[x]) for x in range(0,len(SubpopIN_keep))]
			# Ignore queus if only one species
			if len(XQs) > 0:
				# Loop through queue spots, Putting Nself_pop for grabbing for other species.
				for ispecies in range(len(XQs[spcNO])):
					XQs[spcNO][ispecies].put(Nself_pop) 
			Nother_pop = []
			popcount = 0
			# Ignore queues if only one species
			if len(XQs) > 0:
				for ispecies in range(len(XQs[spcNO])):			
					if ispecies != spcNO: 
						Nother_pop.append(XQs[ispecies][spcNO].get(block = True))
						# Error check here for equal patches
						if len(Nother_pop[popcount]) != len(Nself_pop):
							print("Species systems must have the same number of patches.")
							sys.exit(-1)
						popcount = popcount + 1	
						
		for isub in range(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# K,N for this population
			Kpop = K[isub]
			Npop = len(SubpopIN_arr)
			
			# Overwrite Kpop with K adjusted from DoEmigration that was updated from competition - careful of indexing Tracking variable is patch + 1, first value is total
			if implementcomp == 'Out':
				Kpop = Track_KadjEmi[gen][isub+1]
				
			if implementcomp == 'Back':
				#if multiprocessing.current_process().name == "S1":
				#	pdb.set_trace()	
				# -----------------
				# Adjustment for K if more than one species
				# -----------------
				#Ignore K adjustment if only one species
				if len(XQs) > 0:
					alphas = comp_coef[isub].split(';') # Extracting alphas
					if gen >= startcomp: # Check timing 
						tempspeciesSum = [] # Create sum list for  >= 2 species
						for ispecies in range(len(Nother_pop)):
							tempspeciesSum.append(Nother_pop[ispecies][isub]*float(alphas[ispecies]))
						Kpop = int(Kpop - sum(tempspeciesSum)) # Adjustment or K	
			Kadj_track[gen].append(Kpop) # For Tracking
			
			# Get  indexing numbers
			#Nage_ind_all = np.arange(len(SubpopIN_arr))
			Nage_ind_adults = np.where(SubpopIN_arr['age'] != 0)[0]
			
			if Npop == 0 or Kpop == 0:
				# Append all information to temp SubpopKeep variable
				Nage_samp_ind = np.arange(len(SubpopIN_arr))
				PackingDeaths[gen][isub] = 0
			
			else:
				# Else check only adults, do nothing is adult Npop below Kpop
				if len(Nage_ind_adults) <= Kpop:
					Nage_samp_ind = np.arange(len(SubpopIN_arr))
					PackingDeaths[gen][isub] = 0
				# Else check only adults, apply a mortality to adult if this pop over Kpop
				else:
					Nage_samp_ind_adults = np.random.choice(Nage_ind_adults,Kpop, replace=False).tolist()
					PackingDeaths[gen][isub] = len(Nage_ind_adults) - Kpop
					# Then get the Age0 index
					Nage_samp_ind_age0s = np.where(SubpopIN_arr['age'] == 0)[0]
					# Concate these together
					Nage_samp_ind = np.concatenate((Nage_samp_ind_age0s,np.asarray(Nage_samp_ind_adults)))
				'''# else do nothing if Npop is below Kpop
				if Npop <= Kpop:	
					# Append all information to temp SubpopKeep variable
					Nage_samp_ind = np.arange(len(SubpopIN_arr))
					PackingDeaths[gen][isub] = 0
			
				else:
					# Grab a random draw of Kage from numbers CHANGE MADE HERE Kpop to Npop
					Nage_samp_ind = np.random.choice(np.arange(len(SubpopIN_arr)),Npop, replace=False).tolist()
					PackingDeaths[gen][isub] = len(SubpopIN_arr) -Npop'''					
			
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepK.append(SubpopIN_arr[Nage_samp_ind])
			
			# Get size adjusted age for tracking		
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_keepK[isub]['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']
			
			# Tracking age N
			for iage in range(len(PopulationAge[gen])):
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
			Track_YYSelectionPackDeaths[gen][isub] = 0 # Zero for now
			Track_WildSelectionPackDeaths[gen][isub] = 0 # Zero fro now	
	
	# --------------------------------------------
	# Logistic selected with implementcomp = 'back' 
	# -------------------------------------------- 
	elif packans.split('_')[0] == 'logistic':	
		# Check that popmodel answer is either logistic_out or logistic_back
		try:
			logisticAns = packans.split('_')[1]
		except:
			print("Logistic popmodel option must specify out or back.")
			sys.exit(-1)
		# ------------------------------------------
		# Get other species Ns from all Patches here
		# ------------------------------------------		
		if logisticAns.lower() == 'back':
			if len(XQs) > 0:
				if implementcomp.lower() != 'back':
					print("Competition must be implemented during the same time step (out or back) as the logistic growth model")
					sys.exit(-1)
			Ntself_pop = []
			# Length of each patch without age 0s
			#if multiprocessing.current_process().name == "S0":
				#ForkablePdb().set_trace()	
				# Adjusting N instead of adjusting K
			for isub in range(0,len(SubpopIN_keep)):
				SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
				thisone = SubpopIN_arr['age']
				thisone_ages = count_unique(thisone[np.where(thisone!=0)])
				# Create Nt from countages
				Nt = []
				# Exclude age 0s that were just added
				for iage in range(len(size_mean[0][0])): #Use len(size_mean) for # of ages
					if iage+1 in thisone_ages[0]:
						Nt.append(thisone_ages[1][thisone_ages[0].tolist().index(iage+1)])
					else:
						Nt.append(0)
				# Create list of Nt's, one for each patch
				Ntself_pop.append(Nt)	
			
			# Add to queues if more than one species
			if len(XQs) > 0:
				# Loop through queue spots, Putting Nself_pop for grabbing for other species.
				for ispecies in range(len(XQs[spcNO])): 
					XQs[spcNO][ispecies].put(Ntself_pop) 
			Ntother_pop = []
			popcount = 0
			# Grab from queues only if more than one species
			if len(XQs) > 0:
				for ispecies in range(len(XQs[spcNO])):			
					if ispecies != spcNO: 
						#Nt other pop
						Ntother_pop.append(XQs[ispecies][spcNO].get(block = True))
						# Error check here for equal patches
						if len(Ntother_pop[popcount]) != len(Ntself_pop):
							print("Species systems must have the same number of patches.")
							sys.exit(-1)
						popcount = popcount + 1
				
			# -----------------------
			# Loop through each patch
			for isub in range(len(K)):
				
				# Get each SubpopIN pop as array
				SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
				
				# Count up each uniages - logistic function uses actual ages even if size control is specified.
				countages = count_unique(SubpopIN_arr['age']) # ages will be 1+
				
				# K,N for this population
				Kpop = K[isub]
				# Includes age 0s
				Npop = len(SubpopIN_arr)
				# Excludes age 0s
				Nt = Ntself_pop[isub]
				# Keep track of P for logistic P/K 
				Pisub = sum(Nt)
				#if multiprocessing.current_process().name == "S1":
					#ForkablePdb().set_trace()
				# ----------------------------
				# Parameters for Second species (comp_coef = alpha2_1)
				# ----------------------------
				if len(XQs) > 0 and gen >= startcomp:
					alphas = comp_coef[isub].split(';') # Relabling here, so clear for us right now				
					tempspeciesSum = []
					for ispecies in range(len(Ntother_pop)):
						# if species/patch is extinct, add 0
						#if Ntother_pop[ispecies][isub] == 0:
							#tempspeciesSum.append(0)						
						# Multiply Nt of other species times alpha
						#else:						
						tempspeciesSum.append(sum(Ntother_pop[ispecies][isub])*float(alphas[ispecies]))
							#tempspeciesSum.append(Ntother_pop[ispecies][isub]*float(alphas[ispecies]))
							
					# This is the adjusted N for the N/K ratio in the logistic (P/K, Miller and Ankley 2004)
					Pisub = sum(Nt) + sum(tempspeciesSum)
				
				# Track P/N ratio when using competition
				if sum(Nt) != 0 and len(XQs) > 0:
					Kadj_track[gen].append(Pisub/sum(Nt))
				else:
					Kadj_track[gen].append(np.nan)
				
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
					# Define Leslie matrix M; Jensen 1995 Eco Modelling
					
					# Error check here for more than 1 class vars - move to preprocessing eventually
					if len(f_leslie[isub]) > 1:
						print('Leslie matrix option does not allow for multiple classvars.')
						sys.exit(-1)
					fecund = np.asarray(f_leslie[isub][0], dtype=float) #Import from classvars
					#if multiprocessing.current_process().name == "S0":
					#	ForkablePdb().set_trace()	
									
					Mb = [] # Create N x N-1 matrix				
					for idx in range(len(fecund)-1): 
						for idx2 in range(len(fecund)):
							if idx==idx2:
								if age_percmort[isub][0][idx].split('~')[0] != age_percmort[isub][0][idx].split('~')[1]:
									print('Leslie matrix option does not allow for sex specific mortality.')
									sys.exit(-1)
								Mb.append(1-float(age_percmort[isub][0][idx].split('~')[0]))
							else:
								Mb.append(0.)
					# Add the final row for complete N x N
					M = np.array(fecund.tolist() + Mb).reshape(len(fecund),len(fecund))		
					# dominant eigenvector of M - Stable age distribution or proportion of K allocated to each age class
					vals, vects = np.linalg.eig(M)
					maxcol = list(vals).index(max(vals))
					# Dominant (right) eigenvector
					eigenvect = vects[:,maxcol]
					# Stable Age Distribution
					SAD = eigenvect/sum(eigenvect)
					# Calculate intrinsic growth rate r
					valsabs = np.abs(vals)
					max_index = np.argmax(valsabs)
					# Dominant eigenvalue (greatest magnitude from 0)
					valdom = vals[max_index]
					if valdom == 0:
						print("Invalid Leslie matrix")
						print(M)
						sys.exit(-1)
					# intrinsic growth rate
					# Updated to use absolte value of dominant eigenvalue. I'm not sure if this is ok, but a negative value is impossible to calculate a growth rate. So 
					# either abs(valdom) or max(vals), not sure which as the dominant eigenvalue is meant to be the value with the greatest magnitude
					igr = np.log(abs(valdom))
					# Add a check here for if length of Nt is not len of fecund
					if len(Nt) != len(fecund): 
						print("Check code")
						pdb.set_trace()
					
					# Get total number of deaths for each age class
					Ntplus1 = np.exp((igr*-1)*(Pisub/Kpop))*(np.dot(M,Nt))
					#Ntplus1 = np.rint(Ntplus1) # Round
					# CDay fix for low population sizes compared to age classes
					Ntplus1 = np.divmod(Ntplus1,1)[0] + (np.divmod(Ntplus1,1)[1] >= np.random.rand(len(Ntplus1)))
					Ntplus1[Ntplus1<0] = 0 # Ensure not negative

					#for iage in countages[0]: # Age 0+
					for iage in range(len(Ntplus1)):
						
						# Age class
						Ageclass = iage									
						if Ageclass in countages[0]:
						# Get index for where these ages are
						
							Nage = countages[1][np.where(countages[0]==Ageclass)[0][0]]
						
							Nage_index = np.where(SubpopIN_arr['age']==Ageclass)[0]
							
							# Special case: if last age class, no survivros
							if Ageclass > len(fecund)-1:
								mortreturn = Nage # Then apply mortality to all age class.
								indexforAgeclass = len(fecund) - 1
							else:
								mortreturn = Nage - round(Ntplus1[Ageclass])
								if mortreturn <= 0: # No mortality
									mortreturn = 0
									# Grab all the ages
									Nage_samp_ind.append(list(Nage_index))
								else: # Mortality occurs
									if Nage-mortreturn <= 0: # Kill off all
										Nage_samp_ind.append([])
									else:
										# Grab a random draw of these ages to keep
										Nage_samp_ind.append(np.random.choice(Nage_index,int(Nage-mortreturn),replace=False).tolist())
							# Tracker
							PackingDeathsAge[gen][Ageclass].append(Nage - len(Nage_samp_ind[-1]))
								
						else:
							Nage_samp_ind.append([])
							# Tracker
							PackingDeathsAge[gen][Ageclass].append(0)
					
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
				Track_YYSelectionPackDeaths[gen][isub] = 0 # # No selection in logistic for now
				Track_WildSelectionPackDeaths[gen][isub] = 0 # No selection in logistic for now
				
				# Count up each uniages
				age_adjusted = SubpopIN_keepK[isub]['age']			
				# Tracking age N
				for iage in range(len(PopulationAge[gen])):
					sizeindex = np.where(age_adjusted==iage)[0]
					PopulationAge[gen][iage].append(len(sizeindex))
				# Special case where age class is greater than last age
				sizeindex = np.where(age_adjusted > iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))	
		
		elif logisticAns.lower() == 'out': # What happens if implementcomp = Out
			# Loop through patches for tracking
			for isub in range(len(K)):
				# Get each SubpopIN pop as array
				SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
				
				# Append all information to temp SubpopKeep variable
				SubpopIN_keepK.append(SubpopIN_arr)
				# Store new N - it can be possible to be less than K - its OK - rounding error
				Population[gen].append(len(SubpopIN_keepK[isub]))			
				# K,N for this population
				Kpop = K[isub]
				Npop = len(SubpopIN_arr) # Don't think this is needed
				
				# Track P/N ratio calculated in Emigration, logistic=out
				Kadj_track[gen].append(Track_KadjEmi[gen][isub+1]) # For Tracking
				
				# Store some numbers in this loop too.
				SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
				DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
				PackingDeaths[gen][isub] = 0 # No packing deaths here			
				Track_YYSelectionPackDeaths[gen][isub] = 0 # Zero for now
				Track_WildSelectionPackDeaths[gen][isub] = 0 # Zero fro now
				
				# Count up each uniages
				age_adjusted = SubpopIN_arr['age']
				
				# Tracking numbers - N before packing, packing deaths (0), 
				# -----------------------------------
				
				for iage in range(len(PopulationAge[gen])):
					sizeindex = np.where(age_adjusted==iage)[0]
					PopulationAge[gen][iage].append(len(sizeindex))
					PackingDeathsAge[gen][iage].append(0) # Store 0 for packing deaths age
				# Special case where age class is greater than lastage
				sizeindex = np.where(age_adjusted > iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))
		else:
			print("Logistic popmodel option must specify either out or back.")
			sys.exit(-1)
	else:
		print('See user manual for population model options.')
		sys.exit(-1)
	# Summary numbers
	SelectionDeaths[gen].insert(0,sum(SelectionDeaths[gen]))
	DisperseDeaths[gen].insert(0,sum(DisperseDeaths[gen]))
	PackingDeaths[gen].insert(0,sum(PackingDeaths[gen]))
	Kadj_track[gen].insert(0,sum(Kadj_track[gen]))
	Track_YYSelectionPackDeaths[gen].insert(0,sum(Track_YYSelectionPackDeaths[gen]))
	Track_WildSelectionPackDeaths[gen].insert(0,sum(Track_WildSelectionPackDeaths[gen]))
	StrSuccess[gen] = sum(StrSuccess[gen])
	# Add Population total
	Population[gen].insert(0,sum(Population[gen]))
	# Age tracking
	for iage in range(len(PopulationAge[gen])): 
		PopulationAge[gen][iage] = sum(PopulationAge[gen][iage])
		PackingDeathsAge[gen][iage] = sum(PackingDeathsAge[gen][iage])		
	#pdb.set_trace()	
	return SubpopIN_keepK
	
	# End::DoImmigration()
	
# ---------------------------------------------------------------------------------------------------	
def CalculateDispersalMetrics(OffDisperseIN,FDispDistCD,MDispDistCD,FDispDistCDstd,MDispDistCDstd,subpopmigration,name,gen,cdmatrix_FXX,cdmatrix_MXY,cdmatrix_MYY,cdmatrix_FYY,thresh_FXX,thresh_MXY,thresh_MYY,thresh_FYY,scalemin_FXX,scalemin_MXY,scalemin_MYY,scalemin_FYY,scalemax_FXX,scalemax_MXY,scalemax_MYY,scalemax_FYY,parA_FXX,parA_MXY,parA_MYY,parA_FYY,parB_FXX,parB_MXY,parB_MYY,parB_FYY,parC_FXX,parC_MXY,parC_MYY,parC_FYY,moveno_FXX,moveno_MXY,moveno_MYY,moveno_FYY):
	'''
	CalculateDispersalMetrics()
	This function calculates how far disperses are moving.
	'''		
	# Store the average dispersal distance offspring went
	# temp variable to store offspring dispersal distance
	FtempAvgDispDistCD = []
	MtempAvgDispDistCD = []
	Fcount = 0
	Mcount = 0
	tempN = 0
	
	# Loop through each subpop
	for isub in range(len(OffDisperseIN)):
		
		# Extract the disperser type
		# Add in 'Z' case
		if name == 'S':
			miIndex = np.asarray([i for i, val in enumerate(OffDisperseIN[isub]['name']) if 'S' in val])
			frompop = 'EmiPop'
			topop = 'ImmiPop'
			
		elif name == 'Z':
			miIndex = np.asarray([i for i, val in enumerate(OffDisperseIN[isub]['name']) if 'Z' in val])
			frompop = 'EmiPop'
			topop = 'ImmiPop'
					
		elif name == 'D':
			miIndex = np.asarray([i for i, val in enumerate(OffDisperseIN[isub]['name']) if 'D' in val])
			frompop = 'ImmiPop'
			topop = 'NatalPop'
		else:
			# Get all movers, but not Age0s
			miIndex_temp = np.asarray([i for i, val in enumerate(OffDisperseIN[isub]['name']) if not 'Age' in val])
			if len(miIndex_temp) != 0:
				Ind_temp = OffDisperseIN[isub][miIndex_temp]
				#miIndex = np.asarray([i for i, val in enumerate(Ind_temp['name']) if not 'R' in val])
				miIndex = miIndex_temp
			else:
				miIndex = miIndex_temp
			frompop = 'NatalPop'
			topop = 'ImmiPop'
	
		if len(miIndex) != 0:
			Ind = OffDisperseIN[isub][miIndex]
		else:
			Ind = []
		tempN=tempN+len(Ind) # Count up this disperer type
		
		# Loop through each OffDisperseIN
		for ioff in range(len(Ind)):
			# Grab information from this individual
			indFrom = Ind[frompop][ioff]
			indTo = Ind[topop][ioff]
			indSex = Ind['sex'][ioff]
						
			# Ignore FYY and MYY - only calculate for FXX and MXY
			if indSex == 'FYY' or indSex == 'MYY':
				continue
			
			# Store migration numbers
			if name != 'All':
				if indTo != indFrom:
					subpopmigration[gen][int(indTo)-1].append(1)
			
			# If female - CD distances
			if indSex == 'FXX':
				Fcount = Fcount + 1
				probval = cdmatrix_FXX[int(indFrom)-1][int(indTo)-1]
				
				# If panmictic
				if moveno_FXX == '4' or moveno_FXX == '6': 
					cdval = 0.0
				# If prob matrix or FIDIMO
				elif moveno_FXX == '9' or moveno_FXX == '11':
					cdval = probval
					
				# If linear
				elif moveno_FXX == '1':
					cdval = (probval - 1.) * (-thresh_FXX)				
				
				# If inverse square
				elif moveno_FXX == '2':
					if probval == 1.0:
						cdval = 0.0
					else:	
						cdval = np.sqrt(1. / (probval * (scalemax_FXX - scalemin_FXX) + scalemin_FXX))
					
				# If neg exponetial
				elif moveno_FXX == '5':
					cdval = np.log((probval * (scalemax_FXX-scalemin_FXX) + scalemin_FXX)/float(parA_FXX)) / (-float(parB_FXX) * np.log(10))
					
				elif moveno_FXX == '7':
					cdval = float(parB_FXX) + np.sqrt(-2*float(parC_FXX)**2 * np.log((probval*(scalemax_FXX-scalemin_FXX)+scalemin_FXX)/float(parA_FXX)))
					
				elif moveno_FXX == '8':
					cdval = (1.-probval)*(scalemax_FXX-scalemin_FXX)+scalemin_FXX
				# If pareto
				elif moveno_FXX == '10':
					if probval == max(cdmatrix_FXX[isub]):
						cdval = 0.0
					else:
						cdval = pow(((float(parA_FXX)*float(parB_FXX)**float(parA_FXX))/probval),(1/(float(parA_FXX)+1))) - float(parB_FXX)
				else:
					print('Movement function does not exist')
					sys.exit(-1)
					
				FtempAvgDispDistCD.append(cdval)				
					
			# Else if Male
			elif indSex == 'MXY': # assume males same movement matrix 	
				Mcount = Mcount + 1
				probval = cdmatrix_MXY[int(indFrom)-1][int(indTo)-1]
				
				# If panmictic
				if moveno_MXY == '4' or moveno_MXY == '6': 
					cdval = 0.0
				# If prob matrix or FIDIMO
				elif moveno_MXY == '9' or moveno_MXY == '11':
					cdval = probval
					
				# If linear
				elif moveno_MXY == '1':					
					cdval = (probval - 1.) * (-thresh_MXY)
					
				# If inverse square
				elif moveno_MXY == '2':
					if probval == 1.0:
						cdval = 0.0
					else:	
						cdval = np.sqrt(1. / (probval * (scalemax_MXY - scalemin_MXY) + scalemin_MXY))
					
				# If neg exponetial
				elif moveno_MXY == '5':
					cdval = np.log((probval * (scalemax_MXY-scalemin_MXY) + scalemin_MXY)/float(parA_MXY)) / (-float(parB_MXY) * np.log(10))
				elif moveno_MXY == '7':
					cdval = float(parB_MXY) + np.sqrt(-2*float(parC_MXY)**2 * np.log((probval*(scalemax_MXY-scalemin_MXY)+scalemin_MXY)/float(parA_MXY)))
				elif moveno_MXY == '8':
					cdval = (1. - probval)*(scalemax_MXY-scalemin_MXY)+scalemin_MXY
				# If pareto
				elif moveno_MXY == '10':
					if probval == max(cdmatrix_MXY[isub]):
						cdval = 0.0
					else:
						cdval = pow(((float(parA_MXY)*float(parB_MXY)**float(parA_MXY))/probval),(1/(float(parA_MXY)+1))) - float(parB_MXY)
				else:
					print('Movement function does not exist')
					sys.exit(-1)
				MtempAvgDispDistCD.append(cdval)
			
			else:
				print('Error in sex assignment.')
				sys.exit(-1)
			
			# Store the traveled distance - carefully index here
			if name != 'All':
				OffDisperseIN[isub][miIndex[ioff]]['ImmiCD'] = cdval
		
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
	
	return tempN
	# End::CalculateDispersalMetrics()

# ---------------------------------------------------------------------------------------------------	 
def DoImmigration(SubpopIN,K,natal_patches,gen,cdevolveans,fitvals,subpopmigration,SelectionDeaths,DisperseDeaths,burningen_cdevolve,Str_patch_prob,StrSuccess,age_S,Population,dtype,sizeans,size_mean,PackingDeaths,N_Immigration_age,packans,PackingDeathsAge,packpar1,homeattempt,timecdevolve,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,patchvals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,age_percmort,comp_coef,XQs,Kadj_track,Track_KadjEmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,f_leslie,f_leslie_std,ageDispProb,cdmatrix_FXXBack,cdmatrix_MXYBack,cdmatrix_MYYBack,cdmatrix_FYYBack,thresh_FXXBack,thresh_MXYBack,thresh_MYYBack,thresh_FYYBack,scalemin_FXXBack,scalemin_MXYBack,scalemin_MYYBack,scalemin_FYYBack,scalemax_FXXBack,scalemax_MXYBack,scalemax_MYYBack,scalemax_FYYBack,parA_FXXBack,parA_MXYBack,parA_MYYBack,parA_FYYBack,parB_FXXBack,parB_MXYBack,parB_MYYBack,parB_FYYBack,parC_FXXBack,parC_MXYBack,parC_MYYBack,parC_FYYBack,moveno_FXXBack,moveno_MXYBack,moveno_MYYBack,moveno_FYYBack,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr,thresh_FXXStr,thresh_MXYStr,thresh_MYYStr,thresh_FYYStr,scalemin_FXXStr,scalemin_MXYStr,scalemin_MYYStr,scalemin_FYYStr,scalemax_FXXStr,scalemax_MXYStr,scalemax_MYYStr,scalemax_FYYStr,parA_FXXStr,parA_MXYStr,parA_MYYStr,parA_FYYStr,parB_FXXStr,parB_MXYStr,parB_MYYStr,parB_FYYStr,parC_FXXStr,parC_MXYStr,parC_MYYStr,parC_FYYStr,moveno_FXXStr,moveno_MXYStr,moveno_MYYStr,moveno_FYYStr,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD,thresh_FXXLD,thresh_MXYLD,thresh_MYYLD,thresh_FYYLD,scalemin_FXXLD,scalemin_MXYLD,scalemin_MYYLD,scalemin_FYYLD,scalemax_FXXLD,scalemax_MXYLD,scalemax_MYYLD,scalemax_FYYLD,parA_FXXLD,parA_MXYLD,parA_MYYLD,parA_FYYLD,parB_FXXLD,parB_MXYLD,parB_MYYLD,parB_FYYLD,parC_FXXLD,parC_MXYLD,parC_MYYLD,parC_FYYLD,moveno_FXXLD,moveno_MXYLD,moveno_MYYLD,moveno_FYYLD,sexchromo,age_MgBACK,MgBack_patch_prob,Disperse_patch_prob,MgOut_patch_prob,age_MgOUT,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,migrate,egg_add,outans):

	'''
	DoImmigration()
	Disperse the individual back to patch
	Input: Units of dipsersal, movement function,
	SubpopIN, cdmatrix 
	Output: SubpopIN = [subpop,age,sex,infection,name,genes] 
	'''		
	
	# Population extinct check here
	checkPopN = [len(SubpopIN[x]) for x in range(0,len(SubpopIN))] 
	if sum(checkPopN) != 0: 
	
		# Get size or age control here
		if sizeans == 'Y':
			sizecall = 'size'
		elif sizeans == 'N':
			sizecall = 'age'
		else:
			print('Specify Y or N for size control parameters.')
			sys.exit(-1)
		
		SubpopIN = Immigration(SubpopIN,K,natal_patches,gen,cdevolveans,fitvals,SelectionDeaths,DisperseDeaths,\
		burningen_cdevolve,Str_patch_prob,StrSuccess,age_S,Population,dtype,sizecall,size_mean,PackingDeaths,N_Immigration_age,packans,PackingDeathsAge,packpar1,homeattempt,timecdevolve,patchvals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,age_percmort,comp_coef,XQs,Kadj_track,Track_KadjEmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,f_leslie,f_leslie_std,ageDispProb,cdmatrix_FXXBack,cdmatrix_MXYBack,cdmatrix_MYYBack,cdmatrix_FYYBack,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD,sexchromo,age_MgBACK,MgBack_patch_prob,Disperse_patch_prob,MgOut_patch_prob,age_MgOUT,cdmatrix_FXXOut,cdmatrix_MXYOut,cdmatrix_MYYOut,cdmatrix_FYYOut,migrate,egg_add)
		
		if outans == 'Y':
			# Calculate Dispersal Metrics for strayers 'S' 
			tempStrayS = CalculateDispersalMetrics(SubpopIN,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,subpopmigration,'S',gen,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr,thresh_FXXStr,thresh_MXYStr,thresh_MYYStr,thresh_FYYStr,scalemin_FXXStr,scalemin_MXYStr,scalemin_MYYStr,scalemin_FYYStr,scalemax_FXXStr,scalemax_MXYStr,scalemax_MYYStr,scalemax_FYYStr,parA_FXXStr,parA_MXYStr,parA_MYYStr,parA_FYYStr,parB_FXXStr,parB_MXYStr,parB_MYYStr,parB_FYYStr,parC_FXXStr,parC_MXYStr,parC_MYYStr,parC_FYYStr,moveno_FXXStr,moveno_MXYStr,moveno_MYYStr,moveno_FYYStr)
			
			# Calculate Dispersal Metrics for strayers 'Z' 
			tempStrayZ = CalculateDispersalMetrics(SubpopIN,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,subpopmigration,'Z',gen,cdmatrix_FXXStr,cdmatrix_MXYStr,cdmatrix_MYYStr,cdmatrix_FYYStr,thresh_FXXStr,thresh_MXYStr,thresh_MYYStr,thresh_FYYStr,scalemin_FXXStr,scalemin_MXYStr,scalemin_MYYStr,scalemin_FYYStr,scalemax_FXXStr,scalemax_MXYStr,scalemax_MYYStr,scalemax_FYYStr,parA_FXXStr,parA_MXYStr,parA_MYYStr,parA_FYYStr,parB_FXXStr,parB_MXYStr,parB_MYYStr,parB_FYYStr,parC_FXXStr,parC_MXYStr,parC_MYYStr,parC_FYYStr,moveno_FXXStr,moveno_MXYStr,moveno_MYYStr,moveno_FYYStr)
			
			# Calculate Dispersal Metrics for local dispersers only 'ID' and 'RD'
			tempImmiD = CalculateDispersalMetrics(SubpopIN,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,subpopmigration,'D',gen,cdmatrix_FXXLD,cdmatrix_MXYLD,cdmatrix_MYYLD,cdmatrix_FYYLD,thresh_FXXLD,thresh_MXYLD,thresh_MYYLD,thresh_FYYLD,scalemin_FXXLD,scalemin_MXYLD,scalemin_MYYLD,scalemin_FYYLD,scalemax_FXXLD,scalemax_MXYLD,scalemax_MYYLD,scalemax_FYYLD,parA_FXXLD,parA_MXYLD,parA_MYYLD,parA_FYYLD,parB_FXXLD,parB_MXYLD,parB_MYYLD,parB_FYYLD,parC_FXXLD,parC_MXYLD,parC_MYYLD,parC_FYYLD,moveno_FXXLD,moveno_MXYLD,moveno_MYYLD,moveno_FYYLD)
					
			# Track Subpopulation migration numbers here
			subpopmigration.append([]) # This adds a spot for next generation
			for isub in range(len(K)):
				subpopmigration[gen][isub]=sum(subpopmigration[gen][isub])
				subpopmigration[gen+1].append([0]) # This adds spots for subpops in
		else:
			#ProbAge.append( [0 for x in range(0,len(size_mean[0][0]))] ) # This was
			# Dispersal Metrics
			F_HomeDist.append('NA')
			M_HomeDist.append('NA')
			F_HomeDist_sd.append('NA')
			M_HomeDist_sd.append('NA')
			F_ZtrayDist.append('NA')
			M_ZtrayDist.append('NA')
			F_ZtrayDist_sd.append('NA')
			M_ZtrayDist_sd.append('NA')
			F_StrayDist.append('NA')
			M_StrayDist.append('NA')
			F_StrayDist_sd.append('NA')
			M_StrayDist_sd.append('NA')
			
			# Track Subpopulation migration numbers here
			subpopmigration.append([]) # This adds a spot for next generation
			for isub in range(len(K)):
				subpopmigration[gen][isub]='NA'
				subpopmigration[gen+1].append(['NA']) # This adds spots for subpops in
			
	else: # Population Extinct, return tracking variables 0 only
		# Population variables here
		Population.append( [0 for x in range(0,len(SubpopIN)+1)] )
		#Kadj_track.append( [0 for x in range(0,len(SubpopIN))] )
		Kadj_track.append( [0 for x in range(0,len(SubpopIN)+1)] )
		SelectionDeaths.append( [0 for x in range(0,len(SubpopIN)+1)] )
		DisperseDeaths.append( [0 for x in range(0,len(SubpopIN)+1)] )
		PackingDeaths.append( [0 for x in range(0,len(SubpopIN)+1)] )
		Track_YYSelectionPackDeathsImmi.append( [0 for x in range(0,len(SubpopIN)+1)] )
		Track_WildSelectionPackDeathsImmi.append( [0 for x in range(0,len(SubpopIN)+1)] )
		StrSuccess.append( [0 for x in range(0,len(SubpopIN)+1)] )
		# Age tracking variables here
		PackingDeathsAge.append( [0 for x in range(0,len(size_mean[0][0]))] )
		N_Immigration_age.append( [0 for x in range(0,len(size_mean[0][0]))] )
				
					
		#ProbAge.append( [0 for x in range(0,len(size_mean[0][0]))] ) # This was
		# Dispersal Metrics
		F_HomeDist.append(0)
		M_HomeDist.append(0)
		F_HomeDist_sd.append(0)
		M_HomeDist_sd.append(0)
		F_ZtrayDist.append(0)
		M_ZtrayDist.append(0)
		F_ZtrayDist_sd.append(0)
		M_ZtrayDist_sd.append(0)
		F_StrayDist.append(0)
		M_StrayDist.append(0)
		F_StrayDist_sd.append(0)
		M_StrayDist_sd.append(0)
		
		# Track Subpopulation migration numbers here
		subpopmigration.append([]) # This adds a spot for next generation
		for isub in range(len(K)):
			subpopmigration[gen][isub]=0
			subpopmigration[gen+1].append([0]) # This adds spots for subpops in
		
		# Count ages for queue
		Ntself_pop = []
		for isub in range(0,len(SubpopIN)): # list of lists N by age for each patch
			SubpopIN_arr = np.array(SubpopIN[isub],dtype=dtype)
			thisone = SubpopIN_arr['age']
			thisone_ages = count_unique(thisone[np.where(thisone!=0)])
			# Create Nt from countages
			Nt = []
			# Exclude age 0s that were just added
			for iage in range(len(size_mean[0][0])): #Use len(size_mean) for # of ages
				if iage+1 in thisone_ages[0]:
					Nt.append(thisone_ages[1][thisone_ages[0].tolist().index(iage+1)])
				else:
					Nt.append(0)
			# Create list of Nt's, one for each patch
			Ntself_pop.append(Nt)	
		
		# Ignore queues if only one species
		if len(XQs) > 0:
			# Loop through queue spots, Putting Ntself_pop for grabbing for other species.
			for ispecies in range(len(XQs[spcNO])):
				XQs[spcNO][ispecies].put(Ntself_pop) 
		Ntother_pop = []
		popcount = 0
		# Ignore queues if only one species
		if len(XQs) > 0:
			for ispecies in range(len(XQs[spcNO])):			
				if ispecies != spcNO: 
					Ntother_pop.append(XQs[ispecies][spcNO].get(block = True))				
		# For continue processing of other species, put 0s in PutQ
		#N1_pop = [len(SubpopIN[x]) for x in range(0,len(SubpopIN))]
		#PutQ.put(N1_pop) # Do this for as many species as we have eventually loop (number of species -1)
		#N2_pop = GetQ.get(block = True) # Grabs the above to use, probably list it to house the species - 1
		
	# Return variables from this argument
	return SubpopIN
	
	# End::DoImmigration()