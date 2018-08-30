# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Immigration.py
# Author: Erin L Landguth
# Created: January 2013
# Description: This is the function/module file for immigration processes.
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
from scipy.stats import truncnorm
from CDmetaPOP_Offspring import DoOffspringVars
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
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()
	
# ---------------------------------------------------------------------------------------------------	
def GetProbArray(Fxycdmatrix,Mxycdmatrix,offspring,answer,K,natal):
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
			probarray = Fxycdmatrix[currentsubpop]
		elif indSex == 'XY': # Male 
			probarray = Mxycdmatrix[currentsubpop]
		elif indSex == 'YY': # Male 
			probarray = Mxycdmatrix[currentsubpop]
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
			probarray = Fxycdmatrix[natalsubpop]
		elif indSex == 'XY': # Male 
			probarray = Mxycdmatrix[natalsubpop]
		elif indSex == 'YY': # Male 
			probarray = Mxycdmatrix[natalsubpop]
		else:
			print('Invalid offspring list.')
			sys.exit(-1)
		
		# Where natal grounds = 0, set prob to 0
		probarray[np.where(np.asarray(natal)==0)[0]] = 0.
		# Where K = 0, set prob to 0
		probarray[np.where(np.asarray(K)==0)[0]] = 0.
		
		# Where currentpop to 0 probabilities set to zero
		if indSex == 'XX': # Female 
			nothere = Fxycdmatrix[currentsubpop]
		elif indSex == 'XY': # Male 
			nothere = Mxycdmatrix[currentsubpop]
		elif indSex == 'YY': # Male 
			nothere = Mxycdmatrix[currentsubpop]
		probarray[np.where(nothere == 0.)[0]] = 0.
		
	# Return home attempt
	elif answer == 'immigrator':
		# If K is 0, then can not return home
		if K[natalsubpop] == 0:
			probarray = [0.0]
		else:
			# Append the freegrid probabilities for the offspring choice
			if indSex == 'XX': # Female 
				# Get probarray from current to natal - only one number!
				probarray_BA = [Fxycdmatrix[currentsubpop][natalsubpop]]
				probarray_AB = [Fxycdmatrix[natalsubpop][currentsubpop]]
			elif indSex == 'XY': # Male 
				# Get probarray from current to natal - only one number!
				probarray_BA = [Mxycdmatrix[currentsubpop][natalsubpop]]
				probarray_AB = [Mxycdmatrix[natalsubpop][currentsubpop]]
			elif indSex == 'YY': # Male 
				# Get probarray from current to natal - only one number!
				probarray_BA = [Mxycdmatrix[currentsubpop][natalsubpop]]
				probarray_AB = [Mxycdmatrix[natalsubpop][currentsubpop]]
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
			
	return probarray
	
	# End::GetProbArray()
	
# ---------------------------------------------------------------------------------------------------	
def Immigration(SubpopIN,K,N0,natal,Fxycdmatrix,Mxycdmatrix,gen,\
cdevolveans,fitvals,subpopmigration,SelectionDeaths,DisperseDeaths,\
burningen,ProbPatch,ProbSuccess,\
cdmatrix_StrBack,ProbAge,Population,dtype,sizecall,size_mean,PackingDeaths,PopulationAge,packans,PackingDeathsAge,ithmcrundir,packpar1,noOffspring,Bearpairs,size_std,Femalepercent,sourcePop,transmissionprob,age_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,loci,muterate,mtdna,mutationans,geneswap,allelst,homeattempt,timecdevolve,N_beforePack_pop,N_beforePack_age,SelectionDeathsImm_Age0s,patchvals,assortmateModel,inheritans_classfiles,PopTag,subpopmort_mat,eggFreq,sexans,YYmat_slope,YYmat_int,YYmat_set):
	
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
			originalpop = outpool['NatalPop']
			emipop = outpool['EmiPop']
			outpool_genes = literal_eval(outpool['genes'])
			Indsex = outpool['sex']

			# Get this individuals original ClassVars file and bins for indexing
			natalP = int(SubpopIN[isub][iind]['classfile'].split('_')[0].split('P')[1])
			theseclasspars = int(SubpopIN[isub][iind]['classfile'].split('_')[1].split('CV')[1])
					
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
				if outpool_genes[0][0] == 2: #AA
					indProb = float(fitvals[isub][0][0])
				elif outpool_genes[0][0] == 1 and outpool_genes[0][1] == 1: #As
					indProb = float(fitvals[isub][1][0])
				elif outpool_genes[0][1] == 2: #aa
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
			randProb = rand()	# Get a random number			
			# Flip the coin for patch stray
			if randProb < indProb:			
								
				# Then stray
				indProbans = 'Yes'
						
			# Patch stray not a success
			else:
				indProbans = 'No'
			
			# --------------------------
			# Straying occurred
			# --------------------------
			if indProbans == 'Yes':				
				
				# Then Stray success - go anywhere from current location and use full cdmatrix, minus where K = 0 and natal grounds are 0
				# ---------------------------------------------------------
				
				# Get the probability array
				probarray = GetProbArray(cdmatrix_StrBack,cdmatrix_StrBack,outpool,'strayer_emiPop',K,natal)
										
				# If statement to check if there are spots for offpsring to stray to
				if sum(probarray) != 0.0:
					
					# CDEVOLVE
					if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
						
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
					elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
						
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
					elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
						
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
					
					# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
					straypop = str(iteminlist+1)
					outpool_name = outpool['name']
					outpool_name = outpool_name.split('_')
				
					name = 'S'+str(straypop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
					
					recd = (originalpop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(straypop)-1],outpool['genes'])
								
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
												
			#End::Prob success of stray
			# -------------------------------
			
			# ------------
			# No straying
			# ------------
			else:
			
				# -------------------------------------------------
				# Check if already in natal pop: does not immigrate
				# -------------------------------------------------
				if originalpop == emipop:
					# CDEVOLVE
					if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
												
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
					elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
						
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
					elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen) and (timecdevolve.find('Back') != -1):

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
					
					# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-get in DoCal,age,sex,size,infection,capture,name
					immipop = outpool['EmiPop']
					outpool_name = outpool['name']
					outpool_name = outpool_name.split('_')
					name = 'R'+str(immipop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
					
					recd = (originalpop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(immipop)-1],outpool['genes'])
								
					# Record outpool disperse information	
					SubpopIN_keep[int(immipop)-1].append(recd)				
				
					SelectionDeaths[gen][int(immipop)-1].append(0)
					DisperseDeaths[gen][int(immipop)-1].append(0)
					ProbSuccess[gen].append(0)
					continue				
				
				# Immigrate back to natal population
				# ----------------------------------
				else:				
					# - Use partial cdmatrix - only natal supopulation values - should be only one number - check K = 0 and natal ground indeed 1
					probarray = GetProbArray(Fxycdmatrix,Mxycdmatrix,outpool,'immigrator',K,natal)
					
					# Here check if makes it back
					randback = rand()
					if randback >= probarray[0]: # Does not make it back 
						probarray[0] = 0.0
					else: # Does make it back
						probarray[0] = 1.0 
						
					# If statement to check if there are spots for offpsring to disperse to
					if sum(probarray) != 0.0:
						
						# CDEVOLVE
						if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
													
							# Then it makes it back to original pop
							iteminlist = int(originalpop)-1
							
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
								ProbSuccess[gen].append(0)
								continue
														
						# CDEVOLVE - 2 loci
						elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
							
							# Then it makes it back to original pop
							iteminlist = int(originalpop)-1
							
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
								ProbSuccess[gen].append(0)
								continue
								
						# CDEVOLVE - Hindex
						elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen) and (timecdevolve.find('Back') != -1):

							# Then it makes it back to original pop
							iteminlist = int(originalpop)-1
							
							# Call Hindex selection model
							differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[iteminlist])
													
							# Then flip the coin to see if outpool survives its location
							randcheck = rand()
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
								randcheck = rand()
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
						immipop = str(iteminlist+1)
						outpool_name = outpool['name']
						outpool_name = outpool_name.split('_')
						
						name = 'I'+str(immipop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
						
						recd = (originalpop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(immipop)-1],outpool['genes'])
									
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
							probarray = GetProbArray(cdmatrix_StrBack,cdmatrix_StrBack,outpool,'strayer_emiPop',K,natal)
		
							# If statement to check if there are spots for offpsring to stray to
							if sum(probarray) != 0.0:
								
								# Select the w_choice item
								iteminlist = w_choice_item(probarray)						
								
								# CDEVOLVE
								if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
									
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
								elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
									
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
								elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen) and (timecdevolve.find('Back') != -1):

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
								
								# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
								straypop = str(iteminlist+1)
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = 'Z'+str(straypop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
								
								recd = (originalpop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(straypop)-1],outpool['genes'])
											
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
							probarray = GetProbArray(cdmatrix_StrBack,cdmatrix_StrBack,outpool,'strayer_natalPop',K,natal)
		
							# If statement to check if there are spots for offpsring to stray to
							if sum(probarray) != 0.0:
								
								# Select the w_choice item
								iteminlist = w_choice_item(probarray)						
								
								# CDEVOLVE
								if (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
									
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
								elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen) and (timecdevolve.find('Back') != -1):
									
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
								elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen) and (timecdevolve.find('Back') != -1):

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
								
								# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
								straypop = str(iteminlist+1)
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = 'Z'+str(straypop)+'_F'+str(emipop)+'_'+outpool_name[2]+'_'+outpool_name[3]+'_'+outpool_name[4]+'_'+outpool_name[5]
								
								recd = (originalpop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],outpool['sex'],outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['hindex'],outpool['classfile'],PopTag[int(straypop)-1],outpool['genes'])
											
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
																
				#End::Prob migrate stay in else statement
				# ---------------------------------------
				
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
	
	# Organize type data in SubpopIN_keep
	
	# Loop through each subpop, sort, and grab Kage
	SubpopIN_keepAge1plus = []
	PackingDeathsAge.append([])
	PackingDeathsAge[gen] = [[] for x in xrange(0,len(size_mean[0][0]))]
	N_beforePack_pop.append([])
	N_beforePack_age.append([])
	N_beforePack_age[gen] = [[] for x in xrange(0,len(size_mean[0][0]))]
		
	# -------------------
	# Packing is selected
	# -------------------
	if packans == 'packing':
		# Loop through each subpopulation
		for isub in xrange(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# ----------------------------------------
			# Get the number of eggs/fry in this patch
			# ----------------------------------------
			# Only if there are offspring
			if len(noOffspring) != 0: 
				# Get mother's information
				mothers_patch = Bearpairs[:,0]
				mothers_patch_ind = np.where(mothers_patch['NatalPop']==str(isub+1))[0]
				mothers_patch_file = mothers_patch[mothers_patch_ind]['classfile']
				mothers_patch_hindex = mothers_patch[mothers_patch_ind]['hindex']
				# Father's information
				fathers_patch_hindex = Bearpairs[:,1][mothers_patch_ind]['hindex']
				fathers_patch_file = Bearpairs[:,1][mothers_patch_ind]['classfile']
			
				# Get the offspring in patch
				offspring_patch = noOffspring[mothers_patch_ind]
				Popoffspring = sum(offspring_patch)
				offspring_patch_hindex = mothers_patch_hindex/2. + fathers_patch_hindex/2.
				
				offspring_size = [] # sizes
				#offspring_file = [] # files assigned
				for ifile in xrange(len(offspring_patch)):
					# Total offspring for this female
					theseoffspring = offspring_patch[ifile]
					
					# Then get the classfile assignment for these offspring based on the new hindex - make sure this is stochastic for each offspring or binomial draw
					
					# Grab the mothers and fathers classfile index
					mothers_thisfile = mothers_patch_file[ifile]
					mothers_natalP = int(mothers_thisfile.split('_')[0].split('P')[1])
					mothers_theseclasspars = int(mothers_thisfile.split('_')[1].split('CV')[1])
					fathers_thisfile = fathers_patch_file[ifile]
					fathers_natalP = int(fathers_thisfile.split('_')[0].split('P')[1])
					fathers_theseclasspars = int(fathers_thisfile.split('_')[1].split('CV')[1])
					
					# Get random assignment for either mother or fathers classfile
					# -------------------------------------------------------------
					randnos = np.random.sample(theseoffspring) # Create a list of random numbers
					# Copy and make a int for writing over classfile assignment
					offspring_mu = copy.deepcopy(randnos)
					offspring_sigma = copy.deepcopy(randnos)
					if inheritans_classfiles == 'random': # random assignment
						# Where the numbers are < 0.5, assign fathers class file information
						offspring_mu[np.where(randnos < 0.5)[0]] = size_mean[fathers_natalP][fathers_theseclasspars][0]
						offspring_sigma[np.where(randnos < 0.5)[0]] = size_std[fathers_natalP][fathers_theseclasspars][0]
												
						# Else assign to mother
						offspring_mu[np.where(randnos >= 0.5)[0]] = size_mean[mothers_natalP][mothers_theseclasspars][0]
						offspring_sigma[np.where(randnos >= 0.5)[0]] = size_std[mothers_natalP][mothers_theseclasspars][0]
						#temp_offspring_files[np.where(randnos >= 0.5)[0]] = mothers_thisfile
					else: # Hindex assignment
						theseoffspring_hindex = offspring_patch_hindex[ifile]
						# If rand nos less than hinex of offspring, make 1.0 Hindex file	
						offspring_mu[np.where(randnos < theseoffspring_hindex)[0]] = size_mean[0][0][0]
						offspring_sigma[np.where(randnos < theseoffspring_hindex)[0]] = size_std[0][0][0]
						# Else, they are greater than hindex of offspring, make 0.0 Hindex file
						offspring_mu[np.where(randnos >= theseoffspring_hindex)[0]] = size_mean[0][1][0]
						offspring_sigma[np.where(randnos >= theseoffspring_hindex)[0]] = size_std[0][1][0]
								
					# Get sizes
					# Case for when sigma = 0, temp replace
					sigma0_index = np.where(offspring_sigma == 0)[0]
					if len(sigma0_index) != 0:
						# Temp replace those values
						offspring_sigma[sigma0_index] = 0.0000001
					sizesamp = np.random.normal(offspring_mu,offspring_sigma)
									
					# Check to see if negative values, set to 0
					if not isinstance(sizesamp,float):
						sizesamp[np.where(sizesamp < 0)[0]] = 0.
						# Append to list
						offspring_size.append(sizesamp.tolist())
					else:
						if sizesamp < 0:
							sizesamp = 0.
							# Append to list
							offspring_size.append(sizesamp)
																	
				# Add sizes to SubpopIN_arr['size']
				offspring_size = np.asarray(sum(offspring_size,[]))
				tempSizePatch = np.concatenate((SubpopIN_arr['size'],offspring_size))
				
				# Add age 0 numbers to SubpopIN_arr['age']
				tempNewAge = np.zeros(Popoffspring,dtype='int')
				tempAgePatch = np.concatenate((SubpopIN_arr['age'],tempNewAge))
				
				# Store the classfiles for later reference
				#offspring_file = np.asarray(sum(offspring_file,[]))
				#tempFilePatch = np.concatenate((SubpopIN_arr['classfile'],offspring_file))
							
			# If there are no offspring
			else:
				# Get tracking number
				Popoffspring = 0
				# Get size and age for rest of patch
				tempSizePatch = SubpopIN_arr['size']
				tempAgePatch = SubpopIN_arr['age']
				#tempFilePatch = SubpopIN_arr['classfile']
						
			# -----------------------------------
			# Get countages - for 'age' adjusted
			# -----------------------------------
			# Switch here for size or age control
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, tempSizePatch)
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)
			else:
				# Count up each uniages
				countages = count_unique(tempAgePatch)
				age_adjusted = tempAgePatch
				
			# ------------------------
			# Get packing parameters
			# ------------------------
			# K,N for this population
			Kpop = K[isub]
			Npop = sum(countages[1])
			
			# Tracking N before packing
			# -------------------------
			N_beforePack_pop[gen].append(Npop) 
			for iage in xrange(len(size_mean[0][0])):
				sizeindex = np.where(age_adjusted==iage)[0]
				N_beforePack_age[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			N_beforePack_age[gen][iage].append(len(sizeindex))
			
			# Continue packing
			# ----------------
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
						Nage_index = np.where(tempAgePatch==Ageclass)[0]
															
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
				
			# --------------------------------------------
			# Select out the packed current Age 1+ to keep
			#---------------------------------------------
						
			# Clean up index - grabs largest to less than 0 class
			Nage_samp_ind_all = sum(Nage_samp_ind,[])
			
			# Find adults in samp list that survived
			Nage_ind_adults = np.arange(len(SubpopIN_arr))
			index = np.in1d(Nage_samp_ind_all,Nage_ind_adults)
			#index = np.in1d(Nage_ind_adults,Nage_samp_ind_all)
			Nage_samp_ind_adults = np.asarray(Nage_samp_ind_all)[index]	
			
			# Temp for Age 1+ to keep
			if len(Nage_samp_ind_adults) == 0: # No one survived
				deleteall = Nage_ind_adults
				SubpopIN_arr = np.delete(SubpopIN_arr,deleteall)
				SubpopIN_keepAge1plus.append(SubpopIN_arr)			
			else:
				SubpopIN_keepAge1plus.append(SubpopIN_arr[Nage_samp_ind_adults])
						
			# -------------------------------------------------------------
			# Get the new adjusted offspring for each Bearpair - use later
			# -------------------------------------------------------------
			if Npop != 0 and Popoffspring != 0:
				Nage_samp_ind_off = len(Nage_samp_ind_all)-len(Nage_samp_ind_adults)		
				if Nage_samp_ind_off < Popoffspring:
					# Split based on proportion
					offspring_patch = np.around(offspring_patch * float(Nage_samp_ind_off)/sum(offspring_patch))
					offspring_patch = np.asarray(offspring_patch,dtype='int')
					noOffspring[mothers_patch_ind] = offspring_patch
						
			# ------------------
			# Tracking numbers
			# ------------------
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
			PackingDeaths[gen][isub] = Npop - len(Nage_samp_ind_all)	
			
	# ---------------------
	# Packing is turned off
	# ---------------------
	elif packans == 'N':
		# Loop through each patch
		for isub in xrange(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# ----------------------------------------
			# Get the number of eggs/fry in this patch
			# ----------------------------------------
			# Only if there are offspring
			if len(noOffspring) != 0:
				mothers_patch = Bearpairs[:,0]
				mothers_patch_ind = np.where(mothers_patch['NatalPop']==str(isub+1))[0]
				mothers_patch_file = mothers_patch[mothers_patch_ind]['classfile']
				fathers_patch_file = Bearpairs[:,1][mothers_patch_ind]['classfile']
			
				# Get the offspring in patch
				offspring_patch = noOffspring[mothers_patch_ind]
				Popoffspring = sum(offspring_patch)
				
				# Get fry sizes - from mothers classfile
				offspring_size = [] # sizes
				for ifile in xrange(len(offspring_patch)):
					# Total offspring for this female
					theseoffspring = offspring_patch[ifile]
					
					# Grab the mothers and fathers classfile index
					mothers_thisfile = mothers_patch_file[ifile]
					mothers_natalP = int(mothers_thisfile.split('_')[0].split('P')[1])
					mothers_theseclasspars = int(mothers_thisfile.split('_')[1].split('CV')[1])
					fathers_thisfile = fathers_patch_file[ifile]
					fathers_natalP = int(fathers_thisfile.split('_')[0].split('P')[1])
					fathers_theseclasspars = int(fathers_thisfile.split('_')[1].split('CV')[1])		
					
					# Get random assignment for either mother or fathers classfile
					# -------------------------------------------------------------
					randnos = np.random.sample(theseoffspring) # Create a list of random numbers
					# Copy and make a int for writing over classfile assignment
					offspring_mu = copy.deepcopy(randnos)
					offspring_sigma = copy.deepcopy(randnos)
					if inheritans_classfiles == 'random': # random assignment
						# Where the numbers are < 0.5, assign fathers class file information
						offspring_mu[np.where(randnos < 0.5)[0]] = size_mean[fathers_natalP][fathers_theseclasspars][0]
						offspring_sigma[np.where(randnos < 0.5)[0]] = size_std[fathers_natalP][fathers_theseclasspars][0]
												
						# Else assign to mother
						offspring_mu[np.where(randnos >= 0.5)[0]] = size_mean[mothers_natalP][mothers_theseclasspars][0]
						offspring_sigma[np.where(randnos >= 0.5)[0]] = size_std[mothers_natalP][mothers_theseclasspars][0]
						#temp_offspring_files[np.where(randnos >= 0.5)[0]] = mothers_thisfile
					else: # Hindex assignment
						theseoffspring_hindex = offspring_patch_hindex[ifile]
						# If rand nos less than hinex of offspring, make 1.0 Hindex file	
						offspring_mu[np.where(randnos < theseoffspring_hindex)[0]] = size_mean[0][0][0]
						offspring_sigma[np.where(randnos < theseoffspring_hindex)[0]] = size_std[0][0][0]
						# Else, they are greater than hindex of offspring, make 0.0 Hindex file
						offspring_mu[np.where(randnos >= theseoffspring_hindex)[0]] = size_mean[0][1][0]
						offspring_sigma[np.where(randnos >= theseoffspring_hindex)[0]] = size_std[0][1][0]
					
					# Get sizes
					# Case for when sigma = 0, temp replace
					sigma0_index = np.where(offspring_sigma == 0)[0]
					if len(sigma0_index) != 0:
						# Temp replace those values
						offspring_sigma[sigma0_index] = 0.0000001
					sizesamp = np.random.normal(offspring_mu,offspring_sigma)
					
					# Check to see if negative values, set to 0
					if not isinstance(sizesamp,float):
						sizesamp[np.where(sizesamp < 0)[0]] = 0.
						# Append to list
						offspring_size.append(sizesamp.tolist())
					else:
						if sizesamp < 0:
							sizesamp = 0.
							# Append to list
							offspring_size.append(sizesamp)
									
				# Add these to SubpopIN_arr['size']
				offspring_size = np.asarray(sum(offspring_size,[]))
				tempSizePatch = np.concatenate((SubpopIN_arr['size'],offspring_size))
				
				# Add age 0 numbers to SubpopIN_arr['age']
				tempNewAge = np.zeros(Popoffspring,dtype='int')
				tempAgePatch = np.concatenate((SubpopIN_arr['age'],tempNewAge))
			
			# If there are no offspring
			else:
				# Get tracking number
				Popoffspring = 0
				# Get size and age for rest of patch
				tempSizePatch = SubpopIN_arr['size']
				tempAgePatch = SubpopIN_arr['age']	
						
			# K,N for this population
			Kpop = K[isub]
			Npop = len(tempAgePatch)
						
			# -----------------------------------
			# Get countages - for 'age' adjusted
			# -----------------------------------
			# Switch here for size or age control
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, tempSizePatch)
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)
			else:
				# Count up each uniages
				countages = count_unique(tempAgePatch)
				age_adjusted = tempAgePatch
			
			# Tracking N before packing
			# -------------------------
			N_beforePack_pop[gen].append(Npop) 
			for iage in xrange(len(N_beforePack_age[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				N_beforePack_age[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			N_beforePack_age[gen][iage].append(len(sizeindex))
			
			# Get indexing numbers
			Nage_ind_all = np.arange(len(tempAgePatch))
			Nage_ind_adults = np.arange(len(SubpopIN_arr))
			
			# If nobody in this patch
			if Npop == 0 or Kpop == 0:
				# Use empty index
				Nage_samp_ind_all = Nage_ind_adults
				Nage_samp_ind_adults = Nage_ind_adults
				PackingDeaths[gen][isub] = 0			
			else:
				# else do nothing if Npop is below Kpop
				if Npop <= Kpop:	
					# Grab index for adult and all
					Nage_samp_ind_adults = Nage_ind_adults
					Nage_samp_ind_all = Nage_ind_all
					
				# If more than Kpop, grab random number
				else:
					# Get index for all
					Nage_samp_ind_all = np.asarray(random.sample(Nage_ind_all,Kpop))
					# Get just the adults
					index = np.in1d(Nage_samp_ind_all,Nage_ind_adults)
					Nage_samp_ind_adults = Nage_samp_ind_all[index]				
			
			# Select out the Age1+ adults
			# ---------------------------
			SubpopIN_keepAge1plus.append(SubpopIN_arr[Nage_samp_ind_adults])
			
			# Get the new adjusted offspring list for each Bearpair - use later
			# ------------------------------------
			if Npop != 0 and Popoffspring != 0:
				Nage_samp_ind_off = len(Nage_samp_ind_all)-len(Nage_samp_ind_adults)
				if Nage_samp_ind_off < Popoffspring:
					# Split based on proportion
					offspring_patch = np.around(offspring_patch * float(Nage_samp_ind_off)/sum(offspring_patch))
					offspring_patch = np.asarray(offspring_patch,dtype='int')
					noOffspring[mothers_patch_ind] = offspring_patch
			
			# Tracking numbers - no packing deaths
			# ------------------------------------
			PackingDeaths[gen][isub] = 0
			for iage in xrange(len(size_mean[0][0])):
				# Just store 0 for packing deaths age
				PackingDeathsAge[gen][iage].append(0)			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
			
	# ---------------------
	# Packing is turned off - this is a very special case, where eggs only enter into population and can exceed K
	# ---------------------
	elif packans == 'N_keepeggs':
		
		# Loop through each patch
		for isub in xrange(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			
			# ----------------------------------------
			# Get the number of eggs/fry in this patch
			# ----------------------------------------
			# Only if there are offspring and ONLY offspring
			if len(noOffspring) != 0 and len(SubpopIN_arr) == 0:
				mothers_patch = Bearpairs[:,0]
				mothers_patch_ind = np.where(mothers_patch['NatalPop']==str(isub+1))[0]
				mothers_patch_file = mothers_patch[mothers_patch_ind]['classfile']
				fathers_patch_file = Bearpairs[:,1][mothers_patch_ind]['classfile']
			
				# Get the offspring in patch
				offspring_patch = noOffspring[mothers_patch_ind]
				Popoffspring = sum(offspring_patch)
				
				# Get fry sizes - from mothers classfile
				offspring_size = [] # sizes
				for ifile in xrange(len(offspring_patch)):
					# Total offspring for this female
					theseoffspring = offspring_patch[ifile]
					
					# Grab the mothers and fathers classfile index
					mothers_thisfile = mothers_patch_file[ifile]
					mothers_natalP = int(mothers_thisfile.split('_')[0].split('P')[1])
					mothers_theseclasspars = int(mothers_thisfile.split('_')[1].split('CV')[1])
					fathers_thisfile = fathers_patch_file[ifile]
					fathers_natalP = int(fathers_thisfile.split('_')[0].split('P')[1])
					fathers_theseclasspars = int(fathers_thisfile.split('_')[1].split('CV')[1])		
					
					# Get random assignment for either mother or fathers classfile
					# -------------------------------------------------------------
					randnos = np.random.sample(theseoffspring) # Create a list of random numbers
					# Copy and make a int for writing over classfile assignment
					offspring_mu = copy.deepcopy(randnos)
					offspring_sigma = copy.deepcopy(randnos)
					if inheritans_classfiles == 'random': # random assignment
						# Where the numbers are < 0.5, assign fathers class file information
						offspring_mu[np.where(randnos < 0.5)[0]] = size_mean[fathers_natalP][fathers_theseclasspars][0]
						offspring_sigma[np.where(randnos < 0.5)[0]] = size_std[fathers_natalP][fathers_theseclasspars][0]
												
						# Else assign to mother
						offspring_mu[np.where(randnos >= 0.5)[0]] = size_mean[mothers_natalP][mothers_theseclasspars][0]
						offspring_sigma[np.where(randnos >= 0.5)[0]] = size_std[mothers_natalP][mothers_theseclasspars][0]
						#temp_offspring_files[np.where(randnos >= 0.5)[0]] = mothers_thisfile
					else: # Hindex assignment
						theseoffspring_hindex = offspring_patch_hindex[ifile]
						# If rand nos less than hinex of offspring, make 1.0 Hindex file	
						offspring_mu[np.where(randnos < theseoffspring_hindex)[0]] = size_mean[0][0][0]
						offspring_sigma[np.where(randnos < theseoffspring_hindex)[0]] = size_std[0][0][0]
						# Else, they are greater than hindex of offspring, make 0.0 Hindex file
						offspring_mu[np.where(randnos >= theseoffspring_hindex)[0]] = size_mean[0][1][0]
						offspring_sigma[np.where(randnos >= theseoffspring_hindex)[0]] = size_std[0][1][0]
					
					# Get sizes
					# Case for when sigma = 0, temp replace
					sigma0_index = np.where(offspring_sigma == 0)[0]
					if len(sigma0_index) != 0:
						# Temp replace those values
						offspring_sigma[sigma0_index] = 0.0000001
					sizesamp = np.random.normal(offspring_mu,offspring_sigma)
					
					# Check to see if negative values, set to 0
					if not isinstance(sizesamp,float):
						sizesamp[np.where(sizesamp < 0)[0]] = 0.
						# Append to list
						offspring_size.append(sizesamp.tolist())
					else:
						if sizesamp < 0:
							sizesamp = 0.
							# Append to list
							offspring_size.append(sizesamp)
									
				# Add these to SubpopIN_arr['size']
				offspring_size = np.asarray(sum(offspring_size,[]))
				tempSizePatch = np.concatenate((SubpopIN_arr['size'],offspring_size))
				
				# Add age 0 numbers to SubpopIN_arr['age']
				tempNewAge = np.zeros(Popoffspring,dtype='int')
				tempAgePatch = np.concatenate((SubpopIN_arr['age'],tempNewAge))
			
			# If there are no offspring, but other individuals here
			elif len(noOffspring) == 0:
				# Get tracking number
				Popoffspring = 0
				# Get size and age for rest of patch
				tempSizePatch = SubpopIN_arr['size']
				tempAgePatch = SubpopIN_arr['age']	
						
			# K,N for this population
			Kpop = K[isub]
			Npop = len(tempAgePatch)
						
			# -----------------------------------
			# Get countages - for 'age' adjusted
			# -----------------------------------
			# Switch here for size or age control
			if sizecall == 'size':
				age_adjusted = np.searchsorted(size_mean_middles, tempSizePatch)
				# Count up each unique 'sizes'
				countages = count_unique(age_adjusted)
			else:
				# Count up each uniages
				countages = count_unique(tempAgePatch)
				age_adjusted = tempAgePatch
			
			# Tracking N before packing
			# -------------------------
			N_beforePack_pop[gen].append(Npop) 
			for iage in xrange(len(N_beforePack_age[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				N_beforePack_age[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			N_beforePack_age[gen][iage].append(len(sizeindex))
			
			# Get indexing numbers
			Nage_ind_all = np.arange(len(tempAgePatch))
			Nage_ind_adults = np.arange(len(SubpopIN_arr))
			
			# If nobody in this patch
			if Npop == 0 or Kpop == 0:
				# Use empty index
				Nage_samp_ind_all = Nage_ind_adults
				Nage_samp_ind_adults = Nage_ind_adults
				PackingDeaths[gen][isub] = 0			
			else:
				# else do nothing if Npop is below Kpop
				if Npop <= Kpop:	
					# Grab index for adult and all
					Nage_samp_ind_adults = Nage_ind_adults
					Nage_samp_ind_all = Nage_ind_all
					
				# If more than Kpop, grab random number
				else:
					# But only if it was NOT Age0s here
					if len(Nage_ind_adults) != 0:					
						# Get index for all
						Nage_samp_ind_all = np.asarray(random.sample(Nage_ind_all,Kpop))
						# Get just the adults
						index = np.in1d(Nage_samp_ind_all,Nage_ind_adults)
						Nage_samp_ind_adults = Nage_samp_ind_all[index]	
					else:
						Nage_samp_ind_all = np.asarray(random.sample(Nage_ind_all,Npop))
						# Get just the adults
						index = np.in1d(Nage_samp_ind_all,Nage_ind_adults)
						Nage_samp_ind_adults = Nage_samp_ind_all[index]
			
			# Select out the Age1+ adults
			# ---------------------------
			SubpopIN_keepAge1plus.append(SubpopIN_arr[Nage_samp_ind_adults])
			
			# Get the new adjusted offspring list for each Bearpair - use later
			# ------------------------------------
			if Npop != 0 and Popoffspring != 0:
				Nage_samp_ind_off = len(Nage_samp_ind_all)-len(Nage_samp_ind_adults)
				if Nage_samp_ind_off < Popoffspring:
					# Split based on proportion
					offspring_patch = np.around(offspring_patch * float(Nage_samp_ind_off)/sum(offspring_patch))
					offspring_patch = np.asarray(offspring_patch,dtype='int')
					noOffspring[mothers_patch_ind] = offspring_patch
			
			# Tracking numbers - no packing deaths
			# ------------------------------------
			PackingDeaths[gen][isub] = 0
			for iage in xrange(len(size_mean[0][0])):
				# Just store 0 for packing deaths age
				PackingDeathsAge[gen][iage].append(0)			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
			
	# --------------------------------------------
	# Logistic selected - not applied at this step
	# -------------------------------------------- 
	elif packans == 'logistic':
		
		# Loop through each subpop
		for isub in xrange(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepAge1plus.append(SubpopIN_arr)
			
			# K,N for this population
			Kpop = K[isub]
			Npop = len(SubpopIN_arr)
			
			# Count up each uniages
			age_adjusted = SubpopIN_arr['age']
			
			# Tracking numbers - N before packing, packing deaths (0), 
			# -----------------------------------
			N_beforePack_pop[gen].append(Npop)
			PackingDeaths[gen][isub] = 0 # No packing deaths here			
			for iage in xrange(len(N_beforePack_age[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				N_beforePack_age[gen][iage].append(len(sizeindex))
				PackingDeathsAge[gen][iage].append(0) # Store 0 for packing deaths age
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			N_beforePack_age[gen][iage].append(len(sizeindex))
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
							
	# Error check
	else:
		print('See user manual for population model options.')
		sys.exit(-1)
	
	# ------------------------------
	# Get the survived SubpopIN_Age0
	# ------------------------------
	offspring = DoOffspringVars(Bearpairs,Femalepercent,sourcePop,size_mean,transmissionprob,gen,sizecall,age_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,noOffspring,size_std,inheritans_classfiles,eggFreq,sexans,YYmat_slope,YYmat_int,YYmat_set)
	
	# Get dtype for offspring - first check to see if any Bearpairs exist.
	#if Bearpairs[0][0] != -9999:
	#if len(Bearpairs[0][0]) != 1:
	if not isinstance(Bearpairs[0][0],int):	
		offdtype = [('Mother',(str,len(Bearpairs[0][0]['genes']))),('Father',(str,len(Bearpairs[0][0]['genes']))),('NatalPop',(str,len(SubpopIN)+1)),('EmiPop',(str,len(SubpopIN)+1)),('ImmiPop',(str,len(SubpopIN)+1)),('EmiCD',float),('ImmiCD',float),('age',int),('sex',(str,2)),('size',float),('mature',int),('newmature',int),('infection',int),('name',(str,100)),('capture',int),('recapture',int),('layeggs',float),('M_hindex',float),('F_hindex',float),('classfile',(str,100)),('popID',(str,100))]
	else:
		offdtype = [('Mother',(str,2)),('Father',(str,2)),('NatalPop',(str,len(SubpopIN)+1)),('age',int),('sex',(str,2)),('size',float),('mature',int),('newmature',int),('infection',int),('name',(str,100)),('capture',int),('recapture',int),('layeggs',float),('M_hindex',float),('F_hindex',float),('classfile',(str,100)),('popID',(str,100))]

	offspring = np.asarray(offspring,dtype=offdtype) # Convert to array with dytpe
	
	# Update the Wright Fisher case for sex here
	if Femalepercent == 'WrightFisher':
		# If the subpopulation number is not even then sys exit
		if np.mod(len(offspring),2) == 1:
			print("You have WrightFisher specified and the offspring births must be even.")
			sys.exit(-1)
		# Then create half males and females and shuffle
		offsex = np.append(np.zeros(len(offspring)/2,"int"),np.ones(len(offspring)/2,"int"))
		offsex = np.asarray(offsex,dtype=str)
		offsex[np.where(offsex == '0')[0]] = 'XX'
		offsex[np.where(offsex == '1')[0]] = 'XY'
		np.random.shuffle(offsex)
		# And reassign the sex to offspring list
		offspring['sex'] = offsex
	
	# delete Bearpairs
	del(Bearpairs)
	
	# ---------------------------------------------------
	# Call AddAge0() and InheritGenes() and track numbers
	# ---------------------------------------------------
	SubpopIN_keepK = AddAge0s(SubpopIN_keepAge1plus,K,offspring,gen,Population,loci,muterate,mtdna,mutationans,dtype,geneswap,allelst,PopulationAge,sizecall,size_mean,cdevolveans,burningen,timecdevolve,fitvals,SelectionDeathsImm_Age0s,assortmateModel,patchvals,packans)				
	del offspring
	
	# ---------------
	# Summary numbers
	# ---------------
	SelectionDeaths[gen].insert(0,sum(SelectionDeaths[gen]))
	DisperseDeaths[gen].insert(0,sum(DisperseDeaths[gen]))
	PackingDeaths[gen].insert(0,sum(PackingDeaths[gen]))
	N_beforePack_pop[gen].insert(0,sum(N_beforePack_pop[gen]))
	ProbSuccess[gen] = sum(ProbSuccess[gen])
	# Age tracking
	for iage in xrange(len(PackingDeathsAge[gen])):		
		PackingDeathsAge[gen][iage] = sum(PackingDeathsAge[gen][iage])
		N_beforePack_age[gen][iage] = sum(N_beforePack_age[gen][iage])
		
	return SubpopIN_keepK
	
	# End::DoImmigration()
	
# ---------------------------------------------------------------------------------------------------	
def CalculateDispersalMetrics(OffDisperseIN,xgridcopy,ygridcopy,\
Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,FDispDistCD,MDispDistCD,\
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
	for isub in xrange(len(OffDisperseIN)):
		
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
			
		else:
			# Get all movers, but not Age0s
			miIndex_temp = np.asarray([i for i, val in enumerate(OffDisperseIN[isub]['name']) if not 'Age' in val])
			if len(miIndex_temp) != 0:
				Ind_temp = OffDisperseIN[isub][miIndex_temp]
				miIndex = np.asarray([i for i, val in enumerate(Ind_temp['name']) if not 'R' in val])
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
		for ioff in xrange(len(Ind)):
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
					
				elif Fdispmoveno == '7':
					cdval = float(FB) + np.sqrt(-2*float(FC)**2 * np.log((probval*(FScaleMax-FScaleMin)+FScaleMin)/float(FA)))
					
				elif Fdispmoveno == '8':
					cdval = (1.-probval)*(FScaleMax-FScaleMin)+FScaleMin
				# If pareto
				elif Fdispmoveno == '10':
					if probval == 1.0:
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
def DoImmigration(SubpopIN,K,N0,natal,Fdispmoveno,Mdispmoveno,\
Fxycdmatrix,Mxycdmatrix,gen,xgridcopy,\
ygridcopy,cdevolveans,fitvals,subpopmigration,\
SelectionDeaths,DisperseDeaths,burningen,Prob,ProbSuccess,\
StrBackno,cdmatrix_StrBack,ProbAge,Fthreshold,Mthreshold,Strthreshold,Population,dtype,sizeans,size_mean,PackingDeaths,N_Immigration_age,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC,StrScaleMax,StrScaleMin,StrA,StrB,StrC,packans,PackingDeathsAge,ithmcrundir,packpar1,noOffspring,Bearpairs,size_std,Femalepercent,sourcePop,transmissionprob,age_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,loci,muterate,mtdna,mutationans,geneswap,allelst,homeattempt,timecdevolve,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,patchvals,assortmateModel,inheritans_classfiles,PopTag,subpopmort_mat,eggFreq,sexans,YYmat_slope,YYmat_int,YYmat_set):

	'''
	DoImmigration()
	Disperse the individual back to patch
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
	
	SubpopIN = Immigration(SubpopIN,K,N0,natal,\
	Fxycdmatrix,Mxycdmatrix,gen,\
	cdevolveans,fitvals,subpopmigration,\
	SelectionDeaths,DisperseDeaths,\
	burningen,Prob,ProbSuccess,cdmatrix_StrBack,ProbAge,Population,dtype,sizecall,size_mean,PackingDeaths,N_Immigration_age,packans,PackingDeathsAge,ithmcrundir,packpar1,noOffspring,Bearpairs,size_std,Femalepercent,sourcePop,transmissionprob,age_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,loci,muterate,mtdna,mutationans,geneswap,allelst,homeattempt,timecdevolve,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s,patchvals,assortmateModel,inheritans_classfiles,PopTag,subpopmort_mat,eggFreq,sexans,YYmat_slope,YYmat_int,YYmat_set)
	
	# Calculate Dispersal Metrics for strayers 'S' 
	tempStrayN = CalculateDispersalMetrics(SubpopIN,xgridcopy,ygridcopy,\
	StrBackno,StrBackno,cdmatrix_StrBack,cdmatrix_StrBack,\
	F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,subpopmigration,'S',gen,Strthreshold,Strthreshold,StrScaleMax,StrScaleMin,StrScaleMax,StrScaleMin,StrA,StrB,StrC,StrA,StrB,StrC)
	
	# Calculate Dispersal Metrics for strayers 'Z' 
	tempStrayN = CalculateDispersalMetrics(SubpopIN,xgridcopy,ygridcopy,\
	StrBackno,StrBackno,cdmatrix_StrBack,cdmatrix_StrBack,\
	F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,subpopmigration,'Z',gen,Strthreshold,Strthreshold,StrScaleMax,StrScaleMin,StrScaleMax,StrScaleMin,StrA,StrB,StrC,StrA,StrB,StrC)
	
	# Calculate Dispersal Metrics for how far from home - all
	tempImmiN = CalculateDispersalMetrics(SubpopIN,xgridcopy,ygridcopy,\
	Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,\
	F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,subpopmigration,'All',gen,Fthreshold,Mthreshold,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC)
		
	# Track Subpopulation migration numbers here
	subpopmigration.append([]) # This adds a spot for next generation
	for isub in xrange(len(K)):
		subpopmigration[gen][isub]=sum(subpopmigration[gen][isub])
		subpopmigration[gen+1].append([0]) # This adds spots for subpops in
	
	# Return variables from this argument
	return SubpopIN
	
	# End::DoImmigration()