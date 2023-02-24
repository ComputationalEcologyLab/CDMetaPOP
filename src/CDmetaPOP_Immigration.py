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
def GetProbArray(Fxycdmatrix,Mxycdmatrix,offspring,answer,K,natal,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp):
	'''
	GetProbArray()
	This function gets indices for F and M specific cdmatrix values
	In direction cost, this is the row value; xycdmatrix[0] grabs column in original cdmatrix; xycdmatrix[:,0] grabs row vals in original cdmatrix or use xycdmatrix[currentpop][natalpop]
	'''
	
	# Index into offspring array
	currentoff = offspring
	
	# Subpopulation current offspring is in - index - 1
	currentsubpop = int(currentoff['EmiPop']) - 1
	natalsubpop = int(currentoff['NatalPop']) - 1
	indSex = currentoff['sex']
	
	# If straying process - use full cdmatrix
	if answer == 'strayer_emiPop':
		# Append the freegrid probabilities for the offspring choice
		if indSex == 'XX': # Female 
			probarray = copy.deepcopy(Fxycdmatrix[currentsubpop])
		elif indSex == 'XY': # Male 
			probarray = copy.deepcopy(Mxycdmatrix[currentsubpop])
		elif indSex == 'YY': # Male 
			probarray = copy.deepcopy(Mxycdmatrix[currentsubpop])
		else:
			print('Invalid offspring list.')
			sys.exit(-1)
		
		# Where natal grounds = 0, set prob to 0
		probarray[np.where(np.asarray(natal)==0)[0]] = 0.
		# Where K = 0, set prob to 0
		probarray[np.where(np.asarray(K)==0)[0]] = 0.
	
	# If straying process - use full cdmatrix
	elif answer == 'strayer_natalPop':
		
		# Append the freegrid probabilities for the offspring choice
		if indSex == 'XX': # Female 
			probarray = copy.deepcopy(Fxycdmatrix[natalsubpop])
		elif indSex == 'XY': # Male 
			probarray = copy.deepcopy(Mxycdmatrix[natalsubpop])
		elif indSex == 'YY': # Male 
			probarray = copy.deepcopy(Mxycdmatrix[natalsubpop])
		else:
			print('Invalid offspring list.')
			sys.exit(-1)
		
		# Where natal grounds = 0, set prob to 0
		probarray[np.where(np.asarray(natal)==0)[0]] = 0.
		# Where K = 0, set prob to 0
		probarray[np.where(np.asarray(K)==0)[0]] = 0.
		
		# Where currentpop to 0 probabilities set to zero
		if indSex == 'XX': # Female 
			nothere = copy.deepcopy(Fxycdmatrix[currentsubpop])
		elif indSex == 'XY': # Male 
			nothere = copy.deepcopy(Mxycdmatrix[currentsubpop])
		elif indSex == 'YY': # Male 
			nothere = copy.deepcopy(Mxycdmatrix[currentsubpop])
		probarray[np.where(nothere == 0.)[0]] = 0.
		
	# If local dispersal process
	elif answer == 'localDispersal':
		
		# Append the freegrid probabilities for the offspring choice
		if indSex == 'XX': # Female 
			probarray = copy.deepcopy(Fxycdmatrix[natalsubpop])
		elif indSex == 'XY': # Male 
			probarray = copy.deepcopy(Mxycdmatrix[natalsubpop])
		elif indSex == 'YY': # Male 
			probarray = copy.deepcopy(Mxycdmatrix[natalsubpop])
		else:
			print('Invalid offspring list.')
			sys.exit(-1)
		
		# Where natal grounds = 0, set prob to 0
		probarray[np.where(np.asarray(natal)==0)[0]] = 0.
		# Where K = 0, set prob to 0
		probarray[np.where(np.asarray(K)==0)[0]] = 0.
		
	# Return home attempt
	elif answer == 'immigrator':
		# If K is 0, then can not return home
		if K[natalsubpop] == 0:
			probarray = [0.0]
		else:
			# Append the freegrid probabilities for the offspring choice
			if indSex == 'XX': # Female 
				# Get probarray from current to natal - only one number!
				probarray_BA = [copy.deepcopy(Fxycdmatrix[currentsubpop][natalsubpop])]
				probarray_AB = [copy.deepcopy(Fxycdmatrix[natalsubpop][currentsubpop])]
			elif indSex == 'XY': # Male 
				# Get probarray from current to natal - only one number!
				probarray_BA = [copy.deepcopy(Mxycdmatrix[currentsubpop][natalsubpop])]
				probarray_AB = [copy.deepcopy(Mxycdmatrix[natalsubpop][currentsubpop])]
			elif indSex == 'YY': # Male 
				# Get probarray from current to natal - only one number!
				probarray_BA = [copy.deepcopy(Mxycdmatrix[currentsubpop][natalsubpop])]
				probarray_AB = [copy.deepcopy(Mxycdmatrix[natalsubpop][currentsubpop])]
			else:
				print('Invalid offspring list.')
				sys.exit(-1)
			# Check asymmetrical cost - get relative difference for probability
			if probarray_BA[0] == 0.0: # Absolutely can not make it back, e.g., complete barrier
				probarray = [0.0]
			if probarray_BA[0] < probarray_AB[0]: # Harder to make it back - higher cost lower prob, e.g., partial barrier
				probarray = [1. - ((probarray_AB[0] - probarray_BA[0]) / probarray_AB[0])]
			else: # Same probability/cost to make it back, then assume it migrates back
				probarray = [1.0]
	
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
def Immigration(SubpopIN,K,natal,Fxycdmatrix,Mxycdmatrix,gen,\
cdevolveans,fitvals,SelectionDeaths,DisperseDeaths,\
burningen_cdevolve,ProbPatch,ProbSuccess,\
cdmatrix_StrBack,ProbAge,Population,dtype,sizecall,size_mean,PackingDeaths,PopulationAge,packans,PackingDeathsAge,packpar1,homeattempt,timecdevolve,patchvals,PopTag,subpopmort_mat,Track_YYSelectionPackDeaths,Track_WildSelectionPackDeaths,cdmatrix_dispLocal,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,age_percmort,comp_coef,XQs,Kadj_track,Track_KadjEmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,f_leslie,f_leslie_std,ageDispProb):
	
	'''
	Immigration()
	This function settles individuals to a location through residency, straying, immigration.
	SubpopIN - (NatalPop,EmiPop,ImmiPop,age,sex,infection,name,genes)
	'''			
		
	# Store the keepers
	SubpopIN_keep = []
			
	# Get unique number of subpops
	nosubpops = len(K)
	
	# Add spot to track numbers of individuals
	DisperseDeaths.append([])
	SelectionDeaths.append([])
	ProbSuccess.append([]) # Prob switches
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
	
	# Decide where everybody moves - loop through subpop, doesn't need to be shuffled.
	for isub in range(len(SubpopIN)):
		
		# Loop through individuals in this subpopulation
		for iind in range(len(SubpopIN[isub])):
			
			# Get individual, easier for indexing
			outpool = SubpopIN[isub][iind]
			originalpop = outpool['NatalPop']
			emipop = outpool['EmiPop']
			outpool_genes = outpool['genes']
			Indsex = outpool['sex']

			# Get this individuals original ClassVars file and bins for indexing
			natalP = int(SubpopIN[isub][iind]['classfile'].split('_')[0].split('P')[1])
			theseclasspars = int(SubpopIN[isub][iind]['classfile'].split('_')[1].split('CV')[1])
				
			# --------------------------------
			# Check if Migrant
			# --------------------------------
			#if isub + 1 != int(emipop): # If current patch is NOT the same as the individuals 'emipop' == 'Resident'
			if outpool['name'][0] == 'E':
							
				# -----------------------------------
				# Check for strayer first
				# -----------------------------------
				# No genotype tied to stray rate
				if cdevolveans != 'stray':
				
					# Check for patch stray probability
					Str_Patch = ProbPatch[int(originalpop)-1] # index 1 minus
					
					# Check for age/size stray
					indexofProb = outpool[sizecall]
					# If size control then get size nearest to values in age file
					if sizecall == 'size':
						closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb))
						Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0]
					else:
						Find = indexofProb
					# Check for ages over last age
					if Find > len(ProbAge[natalP][theseclasspars]) - 1:
						Find = len(ProbAge[natalP][theseclasspars]) - 1 # Make last age
					Str_Class = ProbAge[natalP][theseclasspars][Find]
										
				# If genotype tied to stray rate
				else:
					if outpool_genes[0] == 2: #AA
						indProb = float(fitvals[isub][0][0])
					elif outpool_genes[0] == 1 and outpool_genes[1] == 1: #Aa
						indProb = float(fitvals[isub][1][0])
					elif outpool_genes[1] == 2: #aa
						indProb = float(fitvals[isub][2][0])
					else: #Something else
						# Check for patch stray probability
						Str_Patch = ProbPatch[int(originalpop)-1] # index 1 minus
						
						# Check for age/size stray
						indexofProb = outpool[sizecall]
						# If size control then get size nearest to values in age file
						if sizecall == 'size':
							closestval = min(size_mean[natalP][theseclasspars], key=lambda x:abs(x-indexofProb))
							Find = np.where(np.asarray(size_mean[natalP][theseclasspars])==closestval)[0][0]
						else:
							Find = indexofProb
						# Check for ages over last age
						if Find > len(ProbAge[natalP][theseclasspars]) - 1:
							Find = len(ProbAge[natalP][theseclasspars]) - 1 # Make last age
						Str_Class = ProbAge[natalP][theseclasspars][Find]

				# Get sex class options if given
				if Indsex == 'XX': # Female
					Str_Class = float(Str_Class.split('~')[0])
				elif Indsex == 'XY': #Male
					# Check if more than value here
					if len(Str_Class.split('~')) > 1:
						Str_Class = float(Str_Class.split('~')[1])
					else:
						Str_Class = float(Str_Class.split('~')[0])
				else: # YY male
					# Check if more than value here
					if len(Str_Class.split('~')) == 3:
						Str_Class = float(Str_Class.split('~')[2])
					elif len(Str_Class.split('~')) == 2:
						Str_Class = float(Str_Class.split('~')[1])
					else:
						Str_Class = float(Str_Class.split('~')[0])
				
				# Then multiply these together
				indProb = Str_Patch * Str_Class
				
				# Decide if strayer?
				randProb = np.random.uniform()	# Get a random number			
				# Flip the coin for patch stray
				if randProb < indProb:			
									
					# Then stray
					indProbans = 'Yes'
							
				# Patch stray not a success
				else:
					indProbans = 'No'
				
				# ---------------------------------
				# Straying occurred for the Migrant
				# ---------------------------------
				if indProbans == 'Yes':				
					
					# Then Stray success - go anywhere from current location and use full cdmatrix, minus where K = 0 and natal grounds are 0
					# ---------------------------------------------------------
					
					# Get the probability array
					probarray = GetProbArray(cdmatrix_StrBack,cdmatrix_StrBack,outpool,'strayer_emiPop',K,natal,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp)
											
					# If spot available to stray to
					if sum(probarray) != 0.0:
						
						# CDEVOLVE
						if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)
							
							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:							
									# for option 3 in which has to be mature
									if cdevolveans == '1_mat' and outpool['mature'] == 0:
										differentialmortality = 0.0
									else:						
										# Call 1-locus selection model
										differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
							else:							
								# for option 3 in which has to be mature
								if cdevolveans == '1_mat' and outpool['mature'] == 0:
									differentialmortality = 0.0
								else:						
									# Call 1-locus selection model
									differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
													
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue
														
						# CDEVOLVE - 2 loci
						elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)

							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# for option 3 in which has to be mature
									if cdevolveans == '2_mat' and outpool['mature'] == 0:
										differentialmortality = 0.0
									else:
										# Call 2-locus selection model
										differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
							else:
								# for option 3 in which has to be mature
								if cdevolveans == '2_mat' and outpool['mature'] == 0:
									differentialmortality = 0.0
								else:
									# Call 2-locus selection model
									differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
												
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue
								
						# CDEVOLVE - Hindex
						elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)

							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# Call Hindex selection model
									differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
							else:
								# Call Hindex selection model
								differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
													
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue					
						
						# CDEVOLVE - Inbreeding F
						elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)
							
							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# Call 2-locus selection model
									differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
							else:
								# Call 2-locus selection model
								differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
													
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue
								
						# CDEVOLVE - Inbreeding * Outbreeding
						elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
						
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)
							
							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# Call Hindex selection model
									differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
							else:
								# Call Hindex selection model
								differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
														
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue
								
						# CDEVOLVE - Multiple loci selection model
						elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)

							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# Call Hindex selection model
									differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
							else:
								# Call Hindex selection model
								differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
												
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
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
								randcheck = np.random.uniform()
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
						
						# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
						# -------------------------
						straypop = str(iteminlist+1)
						outpool_name = outpool['name']
						outpool_name = outpool_name.split('_')
						
						name = 'S'+str(straypop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
						# Note 'NatalPop' gets updated here to the new straypop now
						recd = (straypop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(straypop)-1],outpool['species'],outpool['genes'])
									
						# Record outpool disperse information	
						SubpopIN_keep[int(straypop)-1].append(recd)		
						
						# Record outpool disperse information	
						SelectionDeaths[gen][int(straypop)-1].append(0)
						DisperseDeaths[gen][int(straypop)-1].append(0)
						ProbSuccess[gen].append(1)
												
					# If spot was not available to stray to ('Move Deaths')
					elif sum(probarray) == 0.0:
						# Then Break from the loop and move to next outpool	
						SelectionDeaths[gen][int(emipop)-1].append(0)
						DisperseDeaths[gen][int(emipop)-1].append(1)
						ProbSuccess[gen].append(1)
						continue
													
				#End::Prob success of stray
				# -------------------------------
				
				# -----------------------
				# No straying for Migrant
				# -----------------------
				else:
					
					# -------------------------------------
					# Attempt to move back to natal grounds
					
					# - Use partial cdmatrix - only natal supopulation values - should be only one number - check K = 0 and natal ground indeed 1
					probarray = GetProbArray(Fxycdmatrix,Mxycdmatrix,outpool,'immigrator',K,natal,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp)
					
					# Here check if makes it back
					randback = np.random.uniform()
					if randback >= probarray[0]: # Does not make it back 
						probarray[0] = 0.0
					else: # Does make it back
						probarray[0] = 1.0 
						
					# If it was possible to get back to original location
					if sum(probarray) != 0.0:
						
						# Then it is assumed it made it back to it's natal patch area int(originalpop)-1	
						# ------------------------------------------------------------------------------
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
						indProb = float(ageDispProb[natalP][theseclasspars][Find])
						
						# Decide if disperse?
						randProb = np.random.uniform()	# Get a random number			
						# Flip the coin for patch stray
						if randProb < indProb:									
							# Then Disperse
							indProbans = True
							# Also check to make sure want to local disperse entire populaion
							if not isinstance(cdmatrix_dispLocal,str):
								indProbans = True
							else:
								indProbans = False
						# Patch stray not a success
						else:
							indProbans = False
						
						# ------------------------------------------------------------------------
						# Local disperse this individual from it's natal patch once it moved back
						# ------------------------------------------------------------------------
						if indProbans:
							# Then locally move this individual - go anywhere from natal patch and use full cdmatrix, minus where K = 0 and natal grounds are 0
							# ---------------------------------------------------------
							
							# Get the probability array - be careful not to overwrite previous if statement!
							probarray_LD = GetProbArray(cdmatrix_dispLocal,cdmatrix_dispLocal,outpool,'localDispersal',K,natal,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp)
							
							# If spots available to move to:
							if sum(probarray_LD) != 0.0:
								# CDEVOLVE
								if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									# Select the w_choice item
									iteminlist = w_choice_item(probarray_LD)
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# for option 3 in which has to be mature
											if cdevolveans == '1_mat' and outpool['mature'] == 0:
												differentialmortality = 0.0
											else:						
												# Call 1-locus selection model
												differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '1_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:						
											# Call 1-locus selection model
											differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
							
								# CDEVOLVE - 2 loci
								elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									# Select the w_choice item
									iteminlist = w_choice_item(probarray_LD)
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# for option 3 in which has to be mature
											if cdevolveans == '2_mat' and outpool['mature'] == 0:
												differentialmortality = 0.0
											else:
												# Call 2-locus selection model
												differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '2_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
										
								# CDEVOLVE - Hindex
								elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									# Select the w_choice item
									iteminlist = w_choice_item(probarray_LD)
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
									else:
										# Call Hindex selection model
										differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue	

								# CDEVOLVE - Inbreeding F
								elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									# Select the w_choice item
									iteminlist = w_choice_item(probarray_LD)
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
									else:
										# Call 2-locus selection model
										differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
															
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue

								# CDEVOLVE - Inbreeding * Outbreeding
								elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
								
									# Select the w_choice item
									iteminlist = w_choice_item(probarray_LD)
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
									else:
										# Call Hindex selection model
										differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
								
								# CDEVOLVE - Multiple loci selection model
								elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									# Select the w_choice item
									iteminlist = w_choice_item(probarray_LD)

									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
									else:
										# Call Hindex selection model
										differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
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
									iteminlist = w_choice_item(probarray_LD)
									
									# What subpatchID is individual coming from
									from_subpatch = PopTag[isub]
															
									# What subpatchID is individual proposing to go to
									to_subpatch = PopTag[iteminlist]
									
									# If it is dispersing to another subpatchID
									if from_subpatch != to_subpatch:
										
										# grab the differential mortality associated with moving into this new subpatchID - from subpatch TO subpatch - cols are TO, read row, then col for correct spot
										differentialmortality = subpopmort_mat[int(to_subpatch)-1][int(from_subpatch)-1]
										# Check if survives
										randcheck = np.random.uniform()
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
									iteminlist = w_choice_item(probarray_LD)	
								
								# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-get in DoCal,age,sex,size,infection,capture,name
								# ---------------------------------------
								disppop = str(iteminlist+1)
								#immipop = outpool['EmiPop']
								immipop = int(originalpop)
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								name = 'ID'+str(disppop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
								
								recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
											
								# Record outpool disperse information	
								SubpopIN_keep[int(disppop)-1].append(recd)				
							
								SelectionDeaths[gen][int(disppop)-1].append(0)
								DisperseDeaths[gen][int(disppop)-1].append(0)
								ProbSuccess[gen].append(0)
								continue
															
							# No spots available to move to, then stays in natal patch
							elif sum(probarray_LD) == 0.0:
								
								# CDEVOLVE
								if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
															
									# Then it makes it back to original pop
									iteminlist = int(originalpop)-1
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# for option 3 in which has to be mature
											if cdevolveans == '1_mat' and outpool['mature'] == 0:
												differentialmortality = 0.0
											else:
												# Call 1-locus selection model
												differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '1_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 1-locus selection model
											differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(0)
										continue
																
								# CDEVOLVE - 2 loci
								elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									# Then it makes it back to original pop
									iteminlist = int(originalpop)-1
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# for option 3 in which has to be mature
											if cdevolveans == '2_mat' and outpool['mature'] == 0:
												differentialmortality = 0.0
											else:
												# Call 2-locus selection model
												differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '2_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(0)
										continue
										
								# CDEVOLVE - Inbreeding F
								elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									# Select the w_choice item
									iteminlist = int(originalpop) - 1
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
									else:
										# Call 2-locus selection model
										differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
															
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(0)
										continue
											
								# CDEVOLVE - Hindex
								elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):

									# Then it makes it back to original pop
									iteminlist = int(originalpop)-1
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
									else:
										# Call Hindex selection model
										differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(0)
										continue	

								# CDEVOLVE - Inbreeding * Outbreeding
								elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
								
									# Select the w_choice item
									iteminlist = int(originalpop)-1
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
									else:
										# Call Hindex selection model
										differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(0)
										continue
								
								# CDEVOLVE - Multiple loci selection model
								elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									# Select the w_choice item
									iteminlist = int(originalpop)-1
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
									else:
										# Call Hindex selection model
										differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(0)
										continue
								
								# If subpopulation differential mortality is on
								elif not isinstance(subpopmort_mat,str):
									
									# Then it attemps to make it back to original pop
									iteminlist = int(originalpop)-1
									
									# What subpatchID is individual coming from
									from_subpatch = PopTag[isub]
															
									# What subpatchID is individual proposing to go to
									to_subpatch = PopTag[iteminlist]
									
									# If it is dispersing to another subpatchID
									if from_subpatch != to_subpatch:
										
										# grab the differential mortality associated with moving into this new subpatchID - from subpatch TO subpatch - cols are TO, read row, then col for correct spot
										differentialmortality = subpopmort_mat[int(to_subpatch)-1][int(from_subpatch)-1]
										# Check if survives
										randcheck = np.random.uniform()
										dispersingto = iteminlist
										# If it does not survive
										if randcheck < differentialmortality:
											SelectionDeaths[gen][dispersingto].append(1)
											DisperseDeaths[gen][dispersingto].append(0)
											ProbSuccess[gen].append(1)
											continue					
								
								#If not cdevolve or if cdevolve but it is in burn in gen
								else:
																
									# Then it makes it back to original pop
									iteminlist = int(originalpop)-1								
								
								# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-to get in DoCal(),age,sex,capture,name
								# ---------------------------------
								immipop = str(iteminlist+1)
								immipop = int(originalpop)
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = 'I'+str(immipop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
								
								recd = (originalpop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(immipop)-1],outpool['species'],outpool['genes'])
											
								# Record outpool disperse information	
								SubpopIN_keep[int(immipop)-1].append(recd)		
								
								# Record outpool disperse information	
								SelectionDeaths[gen][int(immipop)-1].append(0)
								DisperseDeaths[gen][int(immipop)-1].append(0)
								ProbSuccess[gen].append(0)
									
						# ---------------------------------------------------------------------------------
						# Local dispersal is turned off - stays in it's natal patch area int(originalpop)-1
						# ---------------------------------------------------------------------------------	
						else:
							# CDEVOLVE
							if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
														
								# Then it makes it back to original pop
								iteminlist = int(originalpop)-1
								
								if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
									# Check if individual's age matches user specified selection age
									if outpool['age'] != int(timecdevolve.split(':')[1]):
										differentialmortality = 0.0
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '1_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 1-locus selection model
											differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
								else:
									# for option 3 in which has to be mature
									if cdevolveans == '1_mat' and outpool['mature'] == 0:
										differentialmortality = 0.0
									else:
										# Call 1-locus selection model
										differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
																
								# Then flip the coin to see if outpool survives its location
								randcheck = np.random.uniform()
								dispersingto = iteminlist
								# If outpool did not survive: break from loop, move to next outpool
								if randcheck < differentialmortality:
									SelectionDeaths[gen][dispersingto].append(1)
									DisperseDeaths[gen][dispersingto].append(0)
									ProbSuccess[gen].append(0)
									continue
															
							# CDEVOLVE - 2 loci
							elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
								
								# Then it makes it back to original pop
								iteminlist = int(originalpop)-1
								
								if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
									# Check if individual's age matches user specified selection age
									if outpool['age'] != int(timecdevolve.split(':')[1]):
										differentialmortality = 0.0
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '2_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
								else:
									# for option 3 in which has to be mature
									if cdevolveans == '2_mat' and outpool['mature'] == 0:
										differentialmortality = 0.0
									else:
										# Call 2-locus selection model
										differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
											
								# Then flip the coin to see if outpool survives its location
								randcheck = np.random.uniform()
								dispersingto = iteminlist
								# If outpool did not survive: break from loop, move to next outpool
								if randcheck < differentialmortality:
									SelectionDeaths[gen][dispersingto].append(1)
									DisperseDeaths[gen][dispersingto].append(0)
									ProbSuccess[gen].append(0)
									continue
									
							# CDEVOLVE - Hindex
							elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):

								# Then it makes it back to original pop
								iteminlist = int(originalpop)-1
								
								if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
									# Check if individual's age matches user specified selection age
									if outpool['age'] != int(timecdevolve.split(':')[1]):
										differentialmortality = 0.0
									else:
										# Call Hindex selection model
										differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
								else:
									# Call Hindex selection model
									differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
															
								# Then flip the coin to see if outpool survives its location
								randcheck = np.random.uniform()
								dispersingto = iteminlist
								# If outpool did not survive: break from loop, move to next outpool
								if randcheck < differentialmortality:
									SelectionDeaths[gen][dispersingto].append(1)
									DisperseDeaths[gen][dispersingto].append(0)
									ProbSuccess[gen].append(0)
									continue

							# CDEVOLVE - Inbreeding F
							elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
								
								# Select the w_choice item
								iteminlist = int(originalpop)-1
								
								if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
									# Check if individual's age matches user specified selection age
									if outpool['age'] != int(timecdevolve.split(':')[1]):
										differentialmortality = 0.0
									else:
										# Call 2-locus selection model
										differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
								else:
									# Call 2-locus selection model
									differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
														
								# Then flip the coin to see if outpool survives its location
								randcheck = np.random.uniform()
								dispersingto = iteminlist
								# If outpool did not survive: break from loop, move to next outpool
								if randcheck < differentialmortality:
									SelectionDeaths[gen][dispersingto].append(1)
									DisperseDeaths[gen][dispersingto].append(0)
									ProbSuccess[gen].append(0)
									continue
									
							# CDEVOLVE - Inbreeding * Outbreeding
							elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
								# Select the w_choice item
								iteminlist = int(originalpop)-1
								
								if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
									# Check if individual's age matches user specified selection age
									if outpool['age'] != int(timecdevolve.split(':')[1]):
										differentialmortality = 0.0
									else:
										# Call Hindex selection model
										differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
								else:
									# Call Hindex selection model
									differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
															
								# Then flip the coin to see if outpool survives its location
								randcheck = np.random.uniform()
								dispersingto = iteminlist
								# If outpool did not survive: break from loop, move to next outpool
								if randcheck < differentialmortality:
									SelectionDeaths[gen][dispersingto].append(1)
									DisperseDeaths[gen][dispersingto].append(0)
									ProbSuccess[gen].append(0)
									continue
							
							# CDEVOLVE - Multiple loci selection model
							elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
								
								# Select the w_choice item
								iteminlist = int(originalpop)-1
								
								if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
									# Check if individual's age matches user specified selection age
									if outpool['age'] != int(timecdevolve.split(':')[1]):
										differentialmortality = 0.0
									else:
										# Call Hindex selection model
										differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
								else:
									# Call Hindex selection model
									differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
																
								# Then flip the coin to see if outpool survives its location
								randcheck = np.random.uniform()
								dispersingto = iteminlist
								# If outpool did not survive: break from loop, move to next outpool
								if randcheck < differentialmortality:
									SelectionDeaths[gen][dispersingto].append(1)
									DisperseDeaths[gen][dispersingto].append(0)
									ProbSuccess[gen].append(0)
									continue
							
							# If subpopulation differential mortality is on
							elif not isinstance(subpopmort_mat,str):
								
								# Then it attemps to make it back to original pop
								iteminlist = int(originalpop)-1
								
								# What subpatchID is individual coming from
								from_subpatch = PopTag[isub]
														
								# What subpatchID is individual proposing to go to
								to_subpatch = PopTag[iteminlist]
								
								# If it is dispersing to another subpatchID
								if from_subpatch != to_subpatch:
									
									# grab the differential mortality associated with moving into this new subpatchID - from subpatch TO subpatch - cols are TO, read row, then col for correct spot
									differentialmortality = subpopmort_mat[int(to_subpatch)-1][int(from_subpatch)-1]
									# Check if survives
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If it does not survive
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue					
							
							#If not cdevolve or if cdevolve but it is in burn in gen
							else:
															
								# Then it makes it back to original pop
								iteminlist = int(originalpop)-1	
											
							# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-to get in DoCal(),age,sex,capture,name
							# --------------------------
							immipop = str(iteminlist+1)
							outpool_name = outpool['name']
							outpool_name = outpool_name.split('_')
							
							name = 'I'+str(immipop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
							
							recd = (originalpop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(immipop)-1],outpool['species'],outpool['genes'])
										
							# Record outpool disperse information	
							SubpopIN_keep[int(immipop)-1].append(recd)		
							
							# Record outpool disperse information	
							SelectionDeaths[gen][int(immipop)-1].append(0)
							DisperseDeaths[gen][int(immipop)-1].append(0)
							ProbSuccess[gen].append(0)						
												
					# If statement to check if there were not spots to disperse to; did not make it back
					elif sum(probarray) == 0.0:
						
						# Could not make it home (either K = 0 or exceded threshold (complete or partial barrier)
						if homeattempt == 'mortality':
							# Store information
							SelectionDeaths[gen][int(emipop)-1].append(0)
							DisperseDeaths[gen][int(emipop)-1].append(1)
							ProbSuccess[gen].append(0)
							# Then Break from the loop and move to next outpool
							continue
						
						# Then attempt to stray the individual one more time
						elif homeattempt == 'stray_emiPop':

							# Get the probability array
							probarray = GetProbArray(cdmatrix_StrBack,cdmatrix_StrBack,outpool,'strayer_emiPop',K,natal,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp)
		
							# If statement to check if there are spots for offpsring to stray to
							if sum(probarray) != 0.0:
								
								# Select the w_choice item
								iteminlist = w_choice_item(probarray)						
								
								# CDEVOLVE
								if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# for option 3 in which has to be mature
											if cdevolveans == '1_mat' and outpool['mature'] == 0:
												differentialmortality = 0.0
											else:
												# Call 1-locus selection model
												differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '1_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 1-locus selection model
											differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
																
								# CDEVOLVE - 2 loci
								elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# for option 3 in which has to be mature
											if cdevolveans == '2_mat' and outpool['mature'] == 0:
												differentialmortality = 0.0
											else:
												# Call 2-locus selection model
												differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '2_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
																
								# CDEVOLVE - Hindex
								elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
									else:
										# Call Hindex selection model
										differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue

								# CDEVOLVE - Inbreeding F
								elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
									else:
										# Call 2-locus selection model
										differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
															
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
								
								# CDEVOLVE - Inbreeding * Outbreeding
								elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
									else:
										# Call Hindex selection model
										differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
								
								# CDEVOLVE - Multiple loci selection model
								elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
									else:
										# Call Hindex selection model
										differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
				
								# If subpopulation differential mortality is on
								elif not isinstance(subpopmort_mat,str):
									
									# What subpatchID is individual coming from
									from_subpatch = PopTag[isub]
															
									# What subpatchID is individual proposing to go to
									to_subpatch = PopTag[iteminlist]
									
									# If it is dispersing to another subpatchID
									if from_subpatch != to_subpatch:
										
										# grab the differential mortality associated with moving into this new subpatchID - from subpatch TO subpatch - cols are TO, read row, then col for correct spot
										differentialmortality = subpopmort_mat[int(to_subpatch)-1][int(from_subpatch)-1]
										# Check if survives
										randcheck = np.random.uniform()
										dispersingto = iteminlist
										# If it does not survive
										if randcheck < differentialmortality:
											SelectionDeaths[gen][dispersingto].append(1)
											DisperseDeaths[gen][dispersingto].append(0)
											ProbSuccess[gen].append(1)
											continue								
								
								# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
								straypop = str(iteminlist+1)
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = 'Z'+str(straypop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
								
								recd = (straypop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(straypop)-1],outpool['species'],outpool['genes'])
											
								# Record outpool disperse information	
								SubpopIN_keep[int(straypop)-1].append(recd)		
								
								# Record outpool disperse information	
								SelectionDeaths[gen][int(straypop)-1].append(0)
								DisperseDeaths[gen][int(straypop)-1].append(0)
								ProbSuccess[gen].append(1)
														
							# If statement to check if there were not spots to disperse to in straying
							elif sum(probarray) == 0.0:
								# Then Break from the loop and move to next outpool	
								SelectionDeaths[gen][int(emipop)-1].append(0)
								DisperseDeaths[gen][int(emipop)-1].append(1)
								ProbSuccess[gen].append(1)
								continue

						# Then attempt to stray the individual one more time
						elif homeattempt == 'stray_natalPop':
							
							# Get the probability array
							probarray = GetProbArray(cdmatrix_StrBack,cdmatrix_StrBack,outpool,'strayer_natalPop',K,natal,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp)
		
							# If statement to check if there are spots for offpsring to stray to
							if sum(probarray) != 0.0:
								
								# Select the w_choice item
								iteminlist = w_choice_item(probarray)						
								
								# CDEVOLVE
								if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# for option 3 in which has to be mature
											if cdevolveans == '1_mat' and outpool['mature'] == 0:
												differentialmortality = 0.0
											else:
												# Call 1-locus selection model
												differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '1_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 1-locus selection model
											differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
																
								# CDEVOLVE - 2 loci
								elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# for option 3 in which has to be mature
											if cdevolveans == '2_mat' and outpool['mature'] == 0:
												differentialmortality = 0.0
											else:
												# Call 2-locus selection model
												differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
									else:
										# for option 3 in which has to be mature
										if cdevolveans == '2_mat' and outpool['mature'] == 0:
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
																
								# CDEVOLVE - Hindex
								elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):

									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
									else:
										# Call Hindex selection model
										differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue

								# CDEVOLVE - Inbreeding F
								elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
									
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call 2-locus selection model
											differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
									else:
										# Call 2-locus selection model
										differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
															
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
										
								# CDEVOLVE - Inbreeding * Outbreeding
								elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
								
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
									else:
										# Call Hindex selection model
										differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
								
								# CDEVOLVE - Multiple loci selection model
								elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
								
									if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
										# Check if individual's age matches user specified selection age
										if outpool['age'] != int(timecdevolve.split(':')[1]):
											differentialmortality = 0.0
										else:
											# Call Hindex selection model
											differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
									else:
										# Call Hindex selection model
										differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
																
									# Then flip the coin to see if outpool survives its location
									randcheck = np.random.uniform()
									dispersingto = iteminlist
									# If outpool did not survive: break from loop, move to next outpool
									if randcheck < differentialmortality:
										SelectionDeaths[gen][dispersingto].append(1)
										DisperseDeaths[gen][dispersingto].append(0)
										ProbSuccess[gen].append(1)
										continue
				
								# If subpopulation differential mortality is on
								elif not isinstance(subpopmort_mat,str):
									
									# What subpatchID is individual coming from
									from_subpatch = PopTag[isub]
															
									# What subpatchID is individual proposing to go to
									to_subpatch = PopTag[iteminlist]
									
									# If it is dispersing to another subpatchID
									if from_subpatch != to_subpatch:
										
										# grab the differential mortality associated with moving into this new subpatchID - from subpatch TO subpatch - cols are TO, read row, then col for correct spot
										differentialmortality = subpopmort_mat[int(to_subpatch)-1][int(from_subpatch)-1]
										# Check if survives
										randcheck = np.random.uniform()
										dispersingto = iteminlist
										# If it does not survive
										if randcheck < differentialmortality:
											SelectionDeaths[gen][dispersingto].append(1)
											DisperseDeaths[gen][dispersingto].append(0)
											ProbSuccess[gen].append(1)
											continue
								
								# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
								straypop = str(iteminlist+1)
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = 'Z'+str(straypop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
								
								recd = (straypop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(straypop)-1],outpool['species'],outpool['genes'])
											
								# Record outpool disperse information	
								SubpopIN_keep[int(straypop)-1].append(recd)		
								
								# Record outpool disperse information	
								SelectionDeaths[gen][int(straypop)-1].append(0)
								DisperseDeaths[gen][int(straypop)-1].append(0)
								ProbSuccess[gen].append(1)
														
							# If statement to check if there were not spots to disperse to in straying
							elif sum(probarray) == 0.0:
								# Then Break from the loop and move to next outpool	
								SelectionDeaths[gen][int(emipop)-1].append(0)
								DisperseDeaths[gen][int(emipop)-1].append(1)
								ProbSuccess[gen].append(1)
								continue
						
						else:
							print('Home attempt must be either mortality or stray options. See user manual.')
							sys.exit(-1)
									
			#End::Migrant Checks
			# ------------------------
				
			# ------------------------------------------------------------------------------
			# Else not a migrant (resident) - 'emipop' == 'Resident' or isub + 1 == 'emipop'
			# ------------------------------------------------------------------------------
			else:				
				# First check if it decides to local disperse give population and age
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
				indProb = float(ageDispProb[natalP][theseclasspars][Find])
				
				# Decide if disperse?
				randProb = np.random.uniform()	# Get a random number			
				# Flip the coin for patch stray
				if randProb < indProb:									
					# Then Disperse
					indProbans = True
					# Also check to make sure want to local disperse entire populaion
					if not isinstance(cdmatrix_dispLocal,str):
						indProbans = True
					else:
						indProbans = False
				# Patch stray not a success
				else:
					indProbans = False
				
				# ------------------------------
				# Local disperse this individual
				# ------------------------------
				if indProbans:
					# Then locally move this individual - go anywhere from current location and use full cdmatrix, minus where K = 0 and natal grounds are 0
					# ---------------------------------------------------------
					
					# Get the probability array
					probarray = GetProbArray(cdmatrix_dispLocal,cdmatrix_dispLocal,outpool,'localDispersal',K,natal,patchvals,cdevolveans,gen,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp)					
				
					# If spots available to move to:
					if sum(probarray) != 0.0:
						
						# CDEVOLVE
						if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)
							
							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# for option 3 in which has to be mature
									if cdevolveans == '1_mat' and outpool['mature'] == 0:
										differentialmortality = 0.0
									else:						
										# Call 1-locus selection model
										differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
							else:
								# for option 3 in which has to be mature
								if cdevolveans == '1_mat' and outpool['mature'] == 0:
									differentialmortality = 0.0
								else:						
									# Call 1-locus selection model
									differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],iteminlist)
															
															
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue
					
						# CDEVOLVE - 2 loci
						elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)

							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# for option 3 in which has to be mature
									if cdevolveans == '2_mat' and outpool['mature'] == 0:
										differentialmortality = 0.0
									else:
										# Call 2-locus selection model
										differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
							else:
								# for option 3 in which has to be mature
								if cdevolveans == '2_mat' and outpool['mature'] == 0:
									differentialmortality = 0.0
								else:
									# Call 2-locus selection model
									differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],iteminlist)
														
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue
								
						# CDEVOLVE - Hindex
						elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)

							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# Call Hindex selection model
									differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
							else:
								# Call Hindex selection model
								differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
														
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue

						# CDEVOLVE - Inbreeding F
						elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)
							
							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# Call 2-locus selection model
									differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
							else:
								# Call 2-locus selection model
								differentialmortality = DoFSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans)
														
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue
						
						# CDEVOLVE - Inbreeding * Outbreeding
						elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
						
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)
							
							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# Call Hindex selection model
									differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
							else:
								# Call Hindex selection model
								differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],iteminlist,EHom,cdevolveans,outpool['hindex'],patchvals[iteminlist])
														
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
							dispersingto = iteminlist
							# If outpool did not survive: break from loop, move to next outpool
							if randcheck < differentialmortality:
								SelectionDeaths[gen][dispersingto].append(1)
								DisperseDeaths[gen][dispersingto].append(0)
								ProbSuccess[gen].append(1)
								continue
						
						# CDEVOLVE - Multiple loci selection model
						elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
							
							# Select the w_choice item
							iteminlist = w_choice_item(probarray)

							if len(timecdevolve.split(':')) > 1: # User indicated age check			
						
								# Check if individual's age matches user specified selection age
								if outpool['age'] != int(timecdevolve.split(':')[1]):
									differentialmortality = 0.0
								else:
									# Call Hindex selection model
									differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
							else:
								# Call Hindex selection model
								differentialmortality = DoMLocusSelection(outpool['genes'],iteminlist,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
														
							# Then flip the coin to see if outpool survives its location
							randcheck = np.random.uniform()
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
								randcheck = np.random.uniform()
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
						
						# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-get in DoCal,age,sex,size,infection,capture,name
						# ---------------------------------------
						disppop = str(iteminlist+1)
						immipop = outpool['EmiPop']
						outpool_name = outpool['name']
						outpool_name = outpool_name.split('_')
						name = 'RD'+str(disppop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
						
						recd = (disppop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(disppop)-1],outpool['species'],outpool['genes'])
									
						# Record outpool disperse information	
						SubpopIN_keep[int(disppop)-1].append(recd)				
					
						SelectionDeaths[gen][int(disppop)-1].append(0)
						DisperseDeaths[gen][int(disppop)-1].append(0)
						ProbSuccess[gen].append(0)
						continue
							
					# No spots available to move to: 
					elif sum(probarray) == 0.0:
						print('Nowhere to locally disperse to. All individuals will die.')
						# Then Break from the loop and move to next outpool	
						SelectionDeaths[gen][int(emipop)-1].append(0)
						DisperseDeaths[gen][int(emipop)-1].append(1)
						ProbSuccess[gen].append(1)
						continue
						
				# -------------------------------------------------------
				# Local dispersal is turned off - Stays in Resident patch
				# -------------------------------------------------------
				else:					
					# CDEVOLVE
					if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
						if len(timecdevolve.split(':')) > 1: # User indicated age check
						
							# Check if individual's age matches user specified selection age
							if outpool['age'] != int(timecdevolve.split(':')[1]):
								differentialmortality = 0.0
							else:
								# for option 3 in which has to be mature
								if cdevolveans == '1_mat' and outpool['mature'] == 0:
									differentialmortality = 0.0
								else:
									# Call 1-locus selection model
									differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],int(emipop)-1)
						else:
							# for option 3 in which has to be mature
							if cdevolveans == '1_mat' and outpool['mature'] == 0:
								differentialmortality = 0.0
							else:
								# Call 1-locus selection model
								differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],int(emipop)-1)
													
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						dispersingto = int(emipop)-1
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(0)
							continue
													
					# CDEVOLVE - 2 loci
					elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
						if len(timecdevolve.split(':')) > 1: # User indicated age check	
							# Check if individual's age matches user specified selection age
							if outpool['age'] != int(timecdevolve.split(':')[1]):
								differentialmortality = 0.0
							else:
								# for option 3 in which has to be mature
								if cdevolveans == '2_mat' and outpool['mature'] == 0:
									differentialmortality = 0.0
								else:
									# Call 2-locus selection model
									differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],int(emipop)-1)
						else:
							# for option 3 in which has to be mature
							if cdevolveans == '2_mat' and outpool['mature'] == 0:
								differentialmortality = 0.0
							else:
								# Call 2-locus selection model
								differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],int(emipop)-1)
													
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						dispersingto = int(emipop)-1
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(0)
							continue
					 
					# CDEVOLVE - Hindex
					elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
						if len(timecdevolve.split(':')) > 1: # User indicated age check						
							# Check if individual's age matches user specified selection age
							if outpool['age'] != int(timecdevolve.split(':')[1]):
								differentialmortality = 0.0
							else:
								# Call Hindex selection model
								differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[int(emipop)-1])						
						# Run differential mortality check for all ages
						else:
							differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[int(emipop)-1])
													
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						dispersingto = int(emipop)-1
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(0)
							continue

					# CDEVOLVE - Inbreeding F
					elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
						if len(timecdevolve.split(':')) > 1: # User indicated age check	
							# Check if individual's age matches user specified selection age
							if outpool['age'] != int(timecdevolve.split(':')[1]):
								differentialmortality = 0.0
							else:
								# Call 2-locus selection model
								differentialmortality = DoFSelection(fitvals,outpool['genes'],int(emipop)-1,EHom,cdevolveans)
						else:
							# Call 2-locus selection model
							differentialmortality = DoFSelection(fitvals,outpool['genes'],int(emipop)-1,EHom,cdevolveans)
													
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						dispersingto = int(emipop)-1
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(0)
							continue
							
					# CDEVOLVE - Inbreeding * Outbreeding
					elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
						if len(timecdevolve.split(':')) > 1: # User indicated age check							
							# Check if individual's age matches user specified selection age
							if outpool['age'] != int(timecdevolve.split(':')[1]):
								differentialmortality = 0.0
							else:
								# Call Hindex selection model
								differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],int(emipop)-1,EHom,cdevolveans,outpool['hindex'],patchvals[int(emipop)-1])
						else:
							# Call Hindex selection model
							differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],int(emipop)-1,EHom,cdevolveans,outpool['hindex'],patchvals[int(emipop)-1])
													
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						dispersingto = int(emipop)-1
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(0)
							continue
					
					# CDEVOLVE - Multiple loci selection model
					elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find('Back') != -1):
						if len(timecdevolve.split(':')) > 1: # User indicated age check			
							
							# Check if individual's age matches user specified selection age
							if outpool['age'] != int(timecdevolve.split(':')[1]):
								differentialmortality = 0.0
							else:
								# Call Hindex selection model
								differentialmortality = DoMLocusSelection(outpool['genes'],int(emipop)-1,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
						else:
								# Call Hindex selection model
								differentialmortality = DoMLocusSelection(outpool['genes'],int(emipop)-1,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)		
					
													
						# Then flip the coin to see if outpool survives its location
						randcheck = np.random.uniform()
						dispersingto = int(emipop)-1
						# If outpool did not survive: break from loop, move to next outpool
						if randcheck < differentialmortality:
							SelectionDeaths[gen][dispersingto].append(1)
							DisperseDeaths[gen][dispersingto].append(0)
							ProbSuccess[gen].append(0)
							continue
						
					# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-get in DoCal,age,sex,size,infection,capture,name
					# ---------------------------------------
					immipop = outpool['EmiPop']
					outpool_name = outpool['name']
					outpool_name = outpool_name.split('_')
					name = 'R'+str(immipop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
					
					recd = (originalpop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['MID'],outpool['FID'],outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(immipop)-1],outpool['species'],outpool['genes'])
								
					# Record outpool disperse information	
					SubpopIN_keep[int(immipop)-1].append(recd)				
				
					SelectionDeaths[gen][int(immipop)-1].append(0)
					DisperseDeaths[gen][int(immipop)-1].append(0)
					ProbSuccess[gen].append(0)
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
				if len(XQs) > 0:
					alphas = comp_coef[isub].split(';') # Extracting alphas
					if gen >= startcomp: # Check timing 
						tempspeciesSum = [] # Create sum list for  >= 2 species
						for ispecies in range(len(Nother_pop)):
							tempspeciesSum.append(Nother_pop[ispecies][isub]*float(alphas[ispecies]))
						Kpop = int(Kpop - sum(tempspeciesSum)) # Adjustment or K
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
					Nage_samp_ind = np.random.choice(np.arange(len(SubpopIN_arr)),Kpop,replace=False).tolist()
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
	# Logistic selected - 
	# -------------------------------------------- 
	elif packans.split('_')[0] == 'logistic':	
		
		# ------------------------------------------
		# Get other species Ns from all Patches here
		# ------------------------------------------		
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
			#if multiprocessing.current_process().name == "S0":
			#	ForkablePdb().set_trace()
			# ----------------------------
			# Parameters for Second species (comp_coef = alpha2_1)
			# ----------------------------
			if len(XQs) > 0:
				alphas = comp_coef[isub].split(';') # Relabling here, so clear for us right now
				if gen >= startcomp:
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
			# May not need adjusted K anymore but could replace with adjusted N
			#Kadj_track[gen].append(Kpop)
			if sum(Nt) != 0:
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
				
				fecund = np.asarray(f_leslie[isub][0], dtype=float) #Import from classvars
				#if multiprocessing.current_process().name == "S0":
				#	ForkablePdb().set_trace()	
								
				Mb = [] # Create N x N-1 matrix
				for idx in range(len(fecund)-1): 
					for idx2 in range(len(fecund)):
						if idx==idx2:
							Mb.append(1-float(age_percmort[isub][0][idx]))
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
				# CDay; fix for low population sizes compared to age classes
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
			Track_YYSelectionPackDeaths[gen][isub] = 0 # Zero for now
			Track_WildSelectionPackDeaths[gen][isub] = 0 # Zero fro now
			
			# Count up each uniages
			age_adjusted = SubpopIN_keepK[isub]['age']			
			# Tracking age N
			for iage in range(len(PopulationAge[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				PopulationAge[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			PopulationAge[gen][iage].append(len(sizeindex))	

			
	else:
		print('See user manual for population model options.')
		sys.exit(-1)
	#pdb.set_trace()	
	# Summary numbers
	SelectionDeaths[gen].insert(0,sum(SelectionDeaths[gen]))
	DisperseDeaths[gen].insert(0,sum(DisperseDeaths[gen]))
	PackingDeaths[gen].insert(0,sum(PackingDeaths[gen]))
	Kadj_track[gen].insert(0,sum(Kadj_track[gen]))
	Track_YYSelectionPackDeaths[gen].insert(0,sum(Track_YYSelectionPackDeaths[gen]))
	Track_WildSelectionPackDeaths[gen].insert(0,sum(Track_WildSelectionPackDeaths[gen]))
	ProbSuccess[gen] = sum(ProbSuccess[gen])
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
def CalculateDispersalMetrics(OffDisperseIN,Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,FDispDistCD,MDispDistCD,\
FDispDistCDstd,MDispDistCDstd,subpopmigration,name,gen,Fthreshold,Mthreshold,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC):
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
						
			# Store migration numbers
			if name != 'All':
				if indTo != indFrom:
					subpopmigration[gen][int(indTo)-1].append(1)
			
			# If female - CD distances
			if indSex == 'XX':
				Fcount = Fcount + 1
				probval = Fxycdmatrix[int(indFrom)-1][int(indTo)-1]
				
				# If panmictic
				if Fdispmoveno == '4' or Fdispmoveno == '6': 
					cdval = 0.0
				# If prob matrix or FIDIMO
				elif Fdispmoveno == '9' or Fdispmoveno == '11':
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
					
				elif Fdispmoveno == '7':
					cdval = float(FB) + np.sqrt(-2*float(FC)**2 * np.log((probval*(FScaleMax-FScaleMin)+FScaleMin)/float(FA)))
					
				elif Fdispmoveno == '8':
					cdval = (1.-probval)*(FScaleMax-FScaleMin)+FScaleMin
				# If pareto
				elif Fdispmoveno == '10':
					if probval == max(Fxycdmatrix[isub]):
						cdval = 0.0
					else:
						cdval = pow(((float(FA)*float(FB)**float(FA))/probval),(1/(float(FA)+1))) - float(FB)
				else:
					print('Movement function does not exist')
					sys.exit(-1)
					
				FtempAvgDispDistCD.append(cdval)				
					
			# Else if Male
			elif indSex == 'XY' or indSex == 'YY': # assume males same movement matrix 	
				Mcount = Mcount + 1
				probval = Mxycdmatrix[int(indFrom)-1][int(indTo)-1]
				
				# If panmictic
				if Mdispmoveno == '4' or Mdispmoveno == '6': 
					cdval = 0.0
				# If prob matrix or FIDIMO
				elif Mdispmoveno == '9' or Mdispmoveno == '11':
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
					if probval == max(Mxycdmatrix[isub]):
						cdval = 0.0
					else:
						cdval = pow(((float(MA)*float(MB)**float(MA))/probval),(1/(float(MA)+1))) - float(MB)
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
def DoImmigration(SubpopIN,K,natal,Fdispmoveno,Mdispmoveno,\
Fxycdmatrix,Mxycdmatrix,gen,cdevolveans,fitvals,subpopmigration,\
SelectionDeaths,DisperseDeaths,burningen_cdevolve,Prob,ProbSuccess,\
StrBackno,cdmatrix_StrBack,ProbAge,Fthreshold,Mthreshold,Strthreshold,Population,dtype,sizeans,size_mean,PackingDeaths,N_Immigration_age,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC,StrScaleMax,StrScaleMin,StrA,StrB,StrC,packans,PackingDeathsAge,packpar1,homeattempt,timecdevolve,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,patchvals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,cdmatrix_dispLocal,dispLocalparA,dispLocalparB,dispLocalparC,thresh_dispLocal,dispLocal_ScaleMin,dispLocal_ScaleMax,dispLocalno,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,age_percmort,comp_coef,XQs,Kadj_track,Track_KadjEmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,f_leslie,f_leslie_std,ageDispProb):

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
	
		# Get size or age control here - for Mg Prob value
		if sizeans == 'Y':
			sizecall = 'size'
		elif sizeans == 'N':
			sizecall = 'age'
		else:
			print('Specify Y or N for size control parameters.')
			sys.exit(-1)
		
		SubpopIN = Immigration(SubpopIN,K,natal,\
		Fxycdmatrix,Mxycdmatrix,gen,\
		cdevolveans,fitvals,\
		SelectionDeaths,DisperseDeaths,\
		burningen_cdevolve,Prob,ProbSuccess,cdmatrix_StrBack,ProbAge,Population,dtype,sizecall,size_mean,PackingDeaths,N_Immigration_age,packans,PackingDeathsAge,packpar1,homeattempt,timecdevolve,patchvals,PopTag,subpopmort_mat,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,cdmatrix_dispLocal,plasticans,burningen_plastic,timeplastic,plastic_behaviorresp,age_percmort,comp_coef,XQs,Kadj_track,Track_KadjEmi,startcomp,spcNO,implementcomp,betas_selection,xvars_betas,maxfit,minfit,f_leslie,f_leslie_std,ageDispProb)
		
		# Calculate Dispersal Metrics for strayers 'S' 
		tempStrayS = CalculateDispersalMetrics(SubpopIN,StrBackno,StrBackno,cdmatrix_StrBack,cdmatrix_StrBack,\
		F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,subpopmigration,'S',gen,Strthreshold,Strthreshold,StrScaleMax,StrScaleMin,StrScaleMax,StrScaleMin,StrA,StrB,StrC,StrA,StrB,StrC)
		
		# Calculate Dispersal Metrics for strayers 'Z' 
		tempStrayZ = CalculateDispersalMetrics(SubpopIN,StrBackno,StrBackno,cdmatrix_StrBack,cdmatrix_StrBack,\
		F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,subpopmigration,'Z',gen,Strthreshold,Strthreshold,StrScaleMax,StrScaleMin,StrScaleMax,StrScaleMin,StrA,StrB,StrC,StrA,StrB,StrC)
		'''
		# Calculate Dispersal Metrics for immigrators only 'I'
		tempImmiI = CalculateDispersalMetrics(SubpopIN,	Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,\
		F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,subpopmigration,'I',gen,Fthreshold,Mthreshold,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC)
		'''
		# Calculate Dispersal Metrics for local dispersers only 'ID' and 'RD'
		tempImmiD = CalculateDispersalMetrics(SubpopIN,dispLocalno,dispLocalno,cdmatrix_dispLocal,cdmatrix_dispLocal,\
		F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,subpopmigration,'D',gen,thresh_dispLocal,thresh_dispLocal,dispLocal_ScaleMax,dispLocal_ScaleMin,dispLocal_ScaleMax,dispLocal_ScaleMin,dispLocalparA,dispLocalparB,dispLocalparC,dispLocalparA,dispLocalparB,dispLocalparC)
				
		# Track Subpopulation migration numbers here
		subpopmigration.append([]) # This adds a spot for next generation
		for isub in range(len(K)):
			subpopmigration[gen][isub]=sum(subpopmigration[gen][isub])
			subpopmigration[gen+1].append([0]) # This adds spots for subpops in
			
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
		ProbSuccess.append( [0 for x in range(0,len(SubpopIN)+1)] )
		# Age tracking variables here
		PackingDeathsAge.append( [0 for x in range(0,len(size_mean[0][0]))] )
		N_Immigration_age.append( [0 for x in range(0,len(size_mean[0][0]))] )
		ProbAge.append( [0 for x in range(0,len(size_mean[0][0]))] )
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