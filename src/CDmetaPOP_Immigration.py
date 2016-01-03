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
from CDmetaPOP_Modules import AddAge0s
from CDmetaPOP_Modules import Do1LocusSelection
from CDmetaPOP_Modules import Do2LocusSelection

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
	indSex = int(currentoff['sex'])
		
	# If straying process - use full cdmatrix
	if answer == 'strayer_emiPop':
		# Append the freegrid probabilities for the offspring choice
		if indSex == 0: # Female offspring
			probarray = Fxycdmatrix[currentsubpop]
		elif indSex == 1: # Male offspring
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
		if indSex == 0: # Female offspring
			probarray = Fxycdmatrix[natalsubpop]
		elif indSex == 1: # Male offspring
			probarray = Mxycdmatrix[natalsubpop]
		else:
			print('Invalid offspring list.')
			sys.exit(-1)
		
		# Where natal grounds = 0, set prob to 0
		probarray[np.where(np.asarray(natal)==0)[0]] = 0.
		# Where K = 0, set prob to 0
		probarray[np.where(np.asarray(K)==0)[0]] = 0.
		
		# Where currentpop to 0 probabilities set to zero
		if indSex == 0: # Female offspring
			nothere = Fxycdmatrix[currentsubpop]
		elif indSex == 1: # Male offspring
			nothere = Mxycdmatrix[currentsubpop]
		probarray[np.where(nothere == 0.)[0]] = 0.
		
	# Return home attempt
	elif answer == 'immigrator':
		# If K is 0, then can not return home
		if K[natalsubpop] == 0:
			probarray = [0.0]
		else:
			# Append the freegrid probabilities for the offspring choice
			if indSex == 0: # Female offspring
				# Get probarray from current to natal - only one number!
				probarray_BA = [Fxycdmatrix[currentsubpop][natalsubpop]]
				probarray_AB = [Fxycdmatrix[natalsubpop][currentsubpop]]
			elif indSex == 1: # Male offspring
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
cdmatrix_StrBack,ProbAge,Population,dtype,sizecall,size_mean,PackingDeaths,PopulationAge,packans,PackingDeathsAge,ithmcrundir,packpar1,noOffspring,Bearpairs,size_std,Femalepercent,sourcePop,transmissionprob,M_mature,F_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,loci,muterate,mtdna,mutationans,geneswap,allelst,homeattempt,timecdevolve,N_beforePack_pop,N_beforePack_age,SelectionDeathsImm_Age0s):
	
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
						
			# -----------------------------------
			# Check for strayer first
			# -----------------------------------
			# Check for patch stray probability
			Str_Patch = ProbPatch[int(originalpop)-1] # index 1 minus
			
			# Check for age/size stray
			indexofProb = outpool[sizecall]
			# If size control then get size nearest to values in age file
			if sizecall == 'size':
				closestval = min(size_mean[isub], key=lambda x:abs(x-indexofProb))
				Find = np.where(np.asarray(size_mean[isub])==closestval)[0][0]
			else:
				Find = indexofProb
			# Check for ages over last age
			if Find > len(ProbAge[isub]) - 1:
				Find = len(ProbAge[isub]) - 1 # Make last age
			Str_Class = ProbAge[isub][Find]

			# Then multiply these together
			indProb = Str_Patch * Str_Class
			
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
					if (cdevolveans == '1' or cdevolveans == '1_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
						
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
					elif (cdevolveans == '2' or cdevolveans == '2_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
						
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
							
					# If not cdevolve or if cdevolve but it is in burn in gen
					else:
						
						# Select the w_choice item
						iteminlist = w_choice_item(probarray)	
						
					# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
					straypop = str(iteminlist+1)
					outpool_name = outpool['name']
					outpool_name = outpool_name.split('_')
					
					name = 'S'+str(straypop)+'_'+str(emipop)+'_'+outpool_name[1]+'_'+outpool_name[2]+'_'+outpool_name[3]
					
					recd = (originalpop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],int(outpool['sex']),outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['genes'])
								
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
					if (cdevolveans == '1' or cdevolveans == '1_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
												
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
					elif (cdevolveans == '2' or cdevolveans == '2_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
						
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
					 
					# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-get in DoCal,age,sex,size,infection,capture,name
					immipop = outpool['EmiPop']
					outpool_name = outpool['name']
					outpool_name = outpool_name.split('_')
					name = 'R'+str(immipop)+'_'+str(emipop)+'_'+outpool_name[1]+'_'+outpool_name[2]+'_'+outpool_name[3]
					
					recd = (originalpop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],int(outpool['sex']),outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['genes'])
								
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
						if (cdevolveans == '1' or cdevolveans == '1_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
													
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
						elif (cdevolveans == '2' or cdevolveans == '2_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
							
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
								
						# If not cdevolve or if cdevolve but it is in burn in gen
						else:
														
							# Then it makes it back to original pop
							iteminlist = int(originalpop)-1	
											
						# Record string name of outpool,OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-to get in DoCal(),age,sex,capture,name
						immipop = str(iteminlist+1)
						outpool_name = outpool['name']
						outpool_name = outpool_name.split('_')
						
						name = 'I'+str(immipop)+'_'+str(emipop)+'_'+outpool_name[1]+'_'+outpool_name[2]+'_'+outpool_name[3]
						
						recd = (originalpop,emipop,immipop,outpool['EmiCD'],-9999,outpool['age'],int(outpool['sex']),outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['genes'])
									
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
								if (cdevolveans == '1' or cdevolveans == '1_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both') :
									
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
								elif (cdevolveans == '2' or cdevolveans == '2_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
									
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
																
								# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
								straypop = str(iteminlist+1)
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = 'Z'+str(straypop)+'_'+str(emipop)+'_'+outpool_name[1]+'_'+outpool_name[2]+'_'+outpool_name[3]
								
								recd = (originalpop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],int(outpool['sex']),outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['genes'])
											
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
								if (cdevolveans == '1' or cdevolveans == '1_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both') :
									
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
								elif (cdevolveans == '2' or cdevolveans == '2_mat') and gen >= burningen and (timecdevolve == 'Back' or timecdevolve == 'Both'):
									
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
																
								# Record string name of OrigninalSubpop,EmiSubpop,ImmiSubpop,EmiCD,ImmiCD-tofill in DoCalc,age,sex,capture,name
								straypop = str(iteminlist+1)
								outpool_name = outpool['name']
								outpool_name = outpool_name.split('_')
								
								name = 'Z'+str(straypop)+'_'+str(emipop)+'_'+outpool_name[1]+'_'+outpool_name[2]+'_'+outpool_name[3]
								
								recd = (originalpop,emipop,straypop,outpool['EmiCD'],-9999,outpool['age'],int(outpool['sex']),outpool['size'],outpool['mature'],outpool['newmature'],int(outpool['infection']),name,outpool['capture'],outpool['recapture'],outpool['layeggs'],outpool['genes'])
											
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
	# Organize type data in SubpopIN_keep
	
	# Loop through each subpop, sort, and grab Kage
	SubpopIN_keepAge1plus = []
	PackingDeathsAge.append([])
	PackingDeathsAge[gen] = [[] for x in xrange(0,len(size_mean[0]))]
	N_beforePack_pop.append([])
	N_beforePack_age.append([])
	N_beforePack_age[gen] = [[] for x in xrange(0,len(size_mean[0]))]
		
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
				mothers_patch = Bearpairs[:,0]
				mothers_patch_ind = np.where(mothers_patch['NatalPop']==str(isub+1))[0]
			
				# Get the offspring in patch
				offspring_patch = noOffspring[mothers_patch_ind]
				Popoffspring = sum(offspring_patch)
				
				# Get fry sizes - initialize assume first size
				mu,sigma = size_mean[isub][0],size_std[isub][0]			
				# Case here for sigma == 0
				if sigma != 0:
					lower, upper = 0,np.inf
					sizesamp  = truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma,size=Popoffspring)			
				else:
					sizesamp = np.asarray([mu]*Popoffspring)
				
				# Add these to SubpopIN_arr['size']
				tempSizePatch = np.concatenate((SubpopIN_arr['size'],sizesamp))
				
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
				
			# -----------------------------------
			# Get countages - for 'age' adjusted
			# -----------------------------------
			# Switch here for size or age control
			if sizecall == 'size':
				size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
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
			for iage in xrange(len(N_beforePack_age[gen])):
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
					if Ageclass > len(size_mean[isub])-1:
						indexforAgeclass = len(size_mean[isub]) - 1
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
						Kage = int(round(Kage))
						Kused = Kage
						Nage_samp_ind.append(random.sample(Nage_index,Kage))
						PackingDeathsAge[gen][indexforAgeclass].append(Nage-Kage)
											
					# The adjust the habitat or reallocation for next class
					if Kage == 0:
						Kage_hab_adj = Kage_hab
					else:
						Kage_hab_adj = Kage_hab - (Kused * Kage_hab / Kage)
				
			# --------------------------------------------
			# Select out thes packed current Age 1+ to keep
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
			
				# Get the offspring in patch
				offspring_patch = noOffspring[mothers_patch_ind]
				Popoffspring = sum(offspring_patch)
				
				# Get fry sizes - initialize assume first size
				mu,sigma = size_mean[isub][0],size_std[isub][0]			
				# Case here for sigma == 0
				if sigma != 0:
					lower, upper = 0,np.inf
					sizesamp  = truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma,size=Popoffspring)			
				else:
					sizesamp = np.asarray([mu]*Popoffspring)
				
				# Add these to SubpopIN_arr['size']
				tempSizePatch = np.concatenate((SubpopIN_arr['size'],sizesamp))
				
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
				size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
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
			for iage in xrange(len(PackingDeathsAge[gen])):
				# Just store 0 for packing deaths age
				PackingDeathsAge[gen][iage].append(0)			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
			
	# --------------------------------------------
	# Logistic selected - not applied at this step
	# -------------------------------------------- 
	elif packans == 'logistic':
	
		for isub in xrange(len(K)):
			
			# Get each SubpopIN pop as array
			SubpopIN_arr = np.array(SubpopIN_keep[isub],dtype=dtype)
			# Append all information to temp SubpopKeep variable
			SubpopIN_keepAge1plus.append(SubpopIN_arr)
			
			# K,N for this population
			Kpop = K[isub]
			Npop = len(SubpopIN_arr)
			
			# Get size adjusted age for tracking		
			if sizecall == 'size':
				size_mean_middles = np.asarray(size_mean[isub])[1:] - np.diff(np.asarray(size_mean[isub]).astype('f'))/2
				age_adjusted = np.searchsorted(size_mean_middles, SubpopIN_arr['size'])
			else: # age call
				# Count up each uniages
				age_adjusted = SubpopIN_arr['age']
			
			# Tracking N before packing
			# -------------------------
			N_beforePack_pop[gen].append(Npop) 
			for iage in xrange(len(N_beforePack_age[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				N_beforePack_age[gen][iage].append(len(sizeindex))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			N_beforePack_age[gen][iage].append(len(sizeindex))
			
			# Tracking numbers - no packing deaths
			# ------------------------------------
			PackingDeaths[gen][isub] = 0
			for iage in xrange(len(PackingDeathsAge[gen])):
				# Just store 0 for packing deaths age
				PackingDeathsAge[gen][iage].append(0)			
			# Store some numbers in this loop too.
			SelectionDeaths[gen][isub] = sum(SelectionDeaths[gen][isub])
			DisperseDeaths[gen][isub] = sum(DisperseDeaths[gen][isub])
				
	# Error check
	else:
		print('Packing answer must be Y or N.')
		sys.exit(-1)
		
	# ------------------------------
	# Get the survived SubpopIN_Age0
	# ------------------------------
	offspring = DoOffspringVars(Bearpairs,Femalepercent,sourcePop,size_mean,transmissionprob,gen,sizecall,M_mature,F_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,noOffspring,size_std)
	
	# Get dtype for offspring - first check to see if any Bearpairs exist.
	if Bearpairs[0][0] != -9999:	
		offdtype = [('Mother',(str,len(Bearpairs[0][0]['genes']))),('Father',(str,len(Bearpairs[0][0]['genes']))),('NatalPop',(str,10)),('EmiPop',(str,2)),('ImmiPop',(str,2)),('EmiCD',float),('ImmiCD',float),('age',int),('sex',int),('size',float),('mature',int),('newmature',int),('infection',int),('name',(str,20)),('capture',int),('recapture',int),('layeggs',float)]
	else:
		offdtype = [('Mother',(str,2)),('Father',(str,2)),('NatalPop',(str,10)),('age',int),('sex',int),('size',float),('mature',int),('newmature',int),('infection',int),('name',(str,20)),('capture',int),('recapture',int),('layeggs',float)]
	
	offspring = np.asarray(offspring,dtype=offdtype) # Convert to array with dytpe
	
	# Update the Wright Fisher case for sex here
	if Femalepercent == 'WrightFisher':
		# If the subpopulation number is not even then sys exit
		if np.mod(len(offspring),2) == 1:
			print("You have WrightFisher specified and the offspring births must be even.")
			sys.exit(-1)
		# Then create half males and females and shuffle
		offsex = np.append(np.zeros(len(offspring)/2,"int"),np.ones(len(offspring)/2,"int"))
		np.random.shuffle(offsex)
		# And reassign the sex to offspring list
		offspring['sex'] = offsex
	
	# delete Bearpairs
	del(Bearpairs)
	
	# ---------------------------------------------------
	# Call AddAge0() and InheritGenes() and track numbers
	# ---------------------------------------------------
	SubpopIN_keepK = AddAge0s(SubpopIN_keepAge1plus,K,offspring,gen,Population,loci,muterate,mtdna,mutationans,dtype,geneswap,allelst,PopulationAge,sizecall,size_mean,cdevolveans,burningen,timecdevolve,fitvals,SelectionDeathsImm_Age0s)				
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
			if int(indSex) == 0:
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
					
				FtempAvgDispDistCD.append(cdval)				
					
			# Else if Male
			elif int(indSex) == 1: 
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
				MtempAvgDispDistCD.append(cdval)

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
StrBackno,cdmatrix_StrBack,ProbAge,Fthreshold,Mthreshold,Strthreshold,Population,dtype,sizeans,size_mean,PackingDeaths,N_Immigration_age,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC,StrScaleMax,StrScaleMin,StrA,StrB,StrC,packans,PackingDeathsAge,ithmcrundir,packpar1,noOffspring,Bearpairs,size_std,Femalepercent,sourcePop,transmissionprob,M_mature,F_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,loci,muterate,mtdna,mutationans,geneswap,allelst,homeattempt,timecdevolve,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd):

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
	burningen,Prob,ProbSuccess,cdmatrix_StrBack,ProbAge,Population,dtype,sizecall,size_mean,PackingDeaths,N_Immigration_age,packans,PackingDeathsAge,ithmcrundir,packpar1,noOffspring,Bearpairs,size_std,Femalepercent,sourcePop,transmissionprob,M_mature,F_mature,Mmat_slope,Mmat_int,Fmat_slope,Fmat_int,Mmat_set,Fmat_set,loci,muterate,mtdna,mutationans,geneswap,allelst,homeattempt,timecdevolve,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s)
	
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