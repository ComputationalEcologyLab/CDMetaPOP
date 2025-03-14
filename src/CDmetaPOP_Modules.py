# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Modules.py
# Author: Erin L Landguth
# Created: June 2010
# Description: This is the function/module file for CDmetaPOP vX
# --------------------------------------------------------------------------------------------------
	
# Python specific functions
import os, copy, pdb, sys
from ast import literal_eval
import numpy as np 
from CDmetaPOP_PostProcess import DoOutput
import scipy.stats
from inspect import currentframe, getframeinfo
import pandas as pd
from CDmetaPOP_Disease import *

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

# ---------------------------------------------------------------------------
class ForkablePdb(pdb.Pdb):

    _original_stdin_fd = sys.stdin.fileno()
    _original_stdin = None

    def __init__(self):
        pdb.Pdb.__init__(self, nosigint=True)

    def _cmdloop(self):
        current_stdin = sys.stdin
        try:
            if not self._original_stdin:
                self._original_stdin = os.fdopen(self._original_stdin_fd)
            sys.stdin = self._original_stdin
            self.cmdloop()
        finally:
            sys.stdin = current_stdin

# ---------------------------------------------------------------------------
def validate(condition, error_message, exit_code=-1):
    """
    Generic error-checking function.

    Args:
        condition (bool): The condition to evaluate.
        error_message (str): The error message to display if the condition fails.
        exit_code (int): The code to exit with if the condition fails. Default is -1.
    """
    if condition:
        print(error_message)
        sys.exit(exit_code)
	
	# End::validate()
	
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

# ---------------------------------------------------------------------------
def move_allele(allele, direction, possiblealleles, offgenes):
    """Move an allele forward or backward based on the direction."""
    offgenes[allele] -= 1
    if direction == 'forward':
        offgenes[allele + 1] += 1
    elif direction == 'backward':
        offgenes[allele - 1] += 1

# ---------------------------------------------------------------------------
def sample_allele(allele_pool, possiblealleles, offgenes, mothergenes, mtdna_ans, lastloci, plastic_or_not):
	# Helper function for sampling alleles in inheritgenes()
	if len(allele_pool) > 0:
		sampled_index = np.random.choice(allele_pool, 1)[0]
		allele_position = possiblealleles[sampled_index]
		offgenes[allele_position] += 1
        # Reset to 1 if allele becomes 2 and only if plastic True
		if plastic_or_not and offgenes[allele_position] == 2:
			offgenes[allele_position] = 1
		# mtDNA is turned on, then override with mothergenes
		if mtdna_ans == 'Y' and lastloci:
			# Force last locus to be mothergenes - possible alleles are from the last loop above
			offgenes[possiblealleles] = mothergenes[possiblealleles]
	

# ---------------------------------------------------------------------------
def split_or_get(value, icdtime,cdclimgentime):
	"""
	Helper function to handle splitting by '|' and getting the corresponding value based on icdtime.
	"""
	# PopVars will come in as tuple if |	
	if isinstance(value, (list, tuple)):
		# Put in a check here to make sure same number of climate time check in values
		if len(value) != len(cdclimgentime):
			print('If using cdclimate, then variables swapped must be the same length. PopVars')
			sys.exit(-1)
		else: # Returns the correct cdclimate value 
			return value[icdtime] 
	else: # Patchbased values will not be split as a tuple
		if len(value.split('|')) > 1: # There is a cdclimate var to swap here		
			# Put in check for same sizes
			if len(value.split('|')) != len(cdclimgentime):
				print('If using cdclimate, then variables swapped must be the same length. PatchVars')
				sys.exit(-1)
			else: # If there is a cdclimate call for value, then split and take the index,   
				return value.split('|')[icdtime]
		else:
			return value
				

# ---------------------------------------------------------------------------		
def process_disp(disp_param, icdtime, sexchromo, cdclimgentime):
	"""Helper function to handle dispersal parameter processing."""
	return sexsplit(split_or_get(disp_param, icdtime,cdclimgentime), sexchromo)


# --------------------------------------------------------------------------
def PrepTextFile(textpath):
	'''
	PrepTextFile() - Prepare the input files
	'''
	
	return textpath.strip('\n').strip('\r')
	
	# End::PrepTextFile()

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
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()	
	
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
	return item,count
	
	#End::w_choice_general()

# ---------------------------------------------------------------------------------------------------	 
def stochastic_update(mu,sigma,list_value,allow_negative=None):
	''' Generate random normal value with given mu,sigma'''
	if mu != 'N' and mu != 'E':		
		mu = float(mu)		
		sigma = float(sigma)
		# Case here for sigma == 0
		if sigma != 0:
			# Call a truncated normal here
			#lower, upper = 0,np.inf
			#X = truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
			X = np.random.normal(mu,sigma)
			if not allow_negative:
				if X < 0:
					X = 0.
			list_value.append(round(X,3))
		else:
			list_value.append(round(mu,3))
	else:
		list_value.append(mu) 
	
	#End::stochastic_update
	
# ---------------------------------------------------------------------------------------------------	 
def sexsplit(var, sexchromo,lineno=None):
	'''
	Checks fields for sex split by tilda, returns tuple length 4
	'''
	var_parts = var.split('~')
	if len(var_parts) > 1 and len(var_parts) < 5:
		if len(var_parts) != sexchromo:
			print(f'Default value must equal the number of sex chromosomes in variable {var}')
			sys.exit(-1)
		
		FXXvar, MXYvar = var_parts[0], var_parts[1]
		if sexchromo == 2:
			MYYvar = FYYvar = 'N'
		elif sexchromo == 3:
			MYYvar = var_parts[2]
			FYYvar = 'N'
		elif sexchromo == 4:
			MYYvar = var_parts[2]
			FYYvar = var_parts[3]

	elif len(var_parts) == 1:
		FXXvar = MXYvar = MYYvar = FYYvar = var

	else:
		print(f'Default value split either 1 or length of sexchromo field for variable {var} at line {lineno}')
		sys.exit(-1)
		
	return FXXvar, MXYvar, MYYvar, FYYvar
	#End::sexsplit()

# ---------------------------------------------------------------------------------------------------	 
def calc_EHom(SubpopIN):
	'''
	calc_Ehom()
	Calculate population level/global Expected Homozygo for Fhi individual calculation used later
	Plink citation (Purcell et al. 2007)
	'''	
	tempgenes = SubpopIN[0]['genes']
	for isub in range(len(SubpopIN)):
		if isub != 0:
			tempgenes = np.concatenate((tempgenes,SubpopIN[isub]['genes']),axis=0)
	
	# ----------------------Loci numbers
	no2s = np.array(tempgenes==2).sum(axis=0)
	no1s = np.array(tempgenes==1).sum(axis=0)
	no0s = np.array(tempgenes==0).sum(axis=0)
	ploc = ((2*no2s) + no1s ) / (2 * len(tempgenes))
	qloc = 1. - ploc
	oneminus2pq = 1. - 2*ploc*qloc
	EHom = sum(oneminus2pq[np.arange(0,len(qloc),2)])
	return EHom
	#End::calc_EHom()

# ---------------------------------------------------------------------------------------------------
def updatePlasticGenes(Ind,cdevolveans,gen,geneswap,burningen_plastic,patchTemp,plasticans,timeplastic, gridsample, patchHab,plastic_signalresp):
	'''
	This function will check and update the plastic gene region.
	'''
	
	# Get location in genes array for plastic region
	# ---------------------------					
	Indgenes = Ind['genes']
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
	plaloci = 1 # For now, 1, but make this general later
	# Get index for plastic region
	plaloci_index = range(selloci*2,selloci*2+plaloci*2)
	
	# If first generation, set alleles in plastic region to 0 or 'turned off', unless geneswap not started
	if gen == geneswap:
		# But ensure that this is only done once if Out;Back specified, just initialize at second update
		if (timeplastic == 'Out;Back' or timeplastic == 'Back;Out') and gridsample == 'Middle':
			if Indgenes[plaloci_index[0]] == 2:
				Indgenes[plaloci_index[0]] = 1
			if Indgenes[plaloci_index[1]] == 2:
				Indgenes[plaloci_index[1]] = 1
		else: # otherwise initialize with timeplastic choice here	
			if Indgenes[plaloci_index[0]] == 2:
				Indgenes[plaloci_index[0]] = 1
			if Indgenes[plaloci_index[1]] == 2:
				Indgenes[plaloci_index[1]] = 1
			
	# Skip if delayed start time
    
	if gen >= burningen_plastic and plasticans.split('_')[0] == 'Temp':

		# If patch temp values are greater than/equal to threshold and check to make sure the alleles are still 0 (not turned on)
		# If patch temp values are greater than/equal to threshold
		if (patchTemp >= plastic_signalresp):
			
			get_plaallele1_index = plaloci_index[0]
			if Indgenes[get_plaallele1_index] == 1:
				Indgenes[get_plaallele1_index] = Indgenes[get_plaallele1_index]+1 
			get_plaallele2_index = plaloci_index[1]
			if Indgenes[get_plaallele2_index] == 1:
				Indgenes[get_plaallele2_index] = Indgenes[get_plaallele2_index]+1
                
	if gen >= burningen_plastic and plasticans.split('_')[0] == 'Hab':
		
		if (float(patchHab) >= plastic_signalresp):
			get_plaallele1_index = plaloci_index[0]
			if Indgenes[get_plaallele1_index] == 1:
				Indgenes[get_plaallele1_index] = Indgenes[get_plaallele1_index]+1 
			get_plaallele2_index = plaloci_index[1]
			if Indgenes[get_plaallele2_index] == 1:
				Indgenes[get_plaallele2_index] = Indgenes[get_plaallele2_index]+1

	#End::updatePlasticGenes()

# ---------------------------------------------------------------------------------------------------	
def callDiffMortality(cdevolveans,gen,burningen_cdevolve,timecdevolve,OutorBack,outpool,fitvals,location,patchvals,betas_selection,xvars_betas,maxfit,minfit,subpopmort_mat,PopTag,isub,EHom=None):
	'''
	This function condenses the large blocks of calls to cdevolve, returning the differentialmortality value.
	It calculates both selection and spatial mortality togethers
	'''
	
	# CDEVOLVE - No
	# -------------
	if cdevolveans == 'N':
		differentialmortality = 0.0
	
	# CDEVOLVE - Do1LocusSelection
	# ----------------------------
	elif (cdevolveans == '1' or cdevolveans == '1_mat' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link') and (gen >= burningen_cdevolve) and (timecdevolve.find(OutorBack) != -1):
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
					differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],location)
		else:
			# for option 3 in which has to be mature
			if cdevolveans == '1_mat' and outpool['mature'] == 0:
				differentialmortality = 0.0
			else:
				# Call 1-locus selection model
				differentialmortality = Do1LocusSelection(fitvals,outpool['genes'][0:2],location)
	
	# CDEVOLVE - Do2LocusSelection
	# ----------------------------
	elif (cdevolveans == '2' or cdevolveans == '2_mat') and (gen >= burningen_cdevolve) and (timecdevolve.find(OutorBack) != -1):
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
					differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],location)
		else:
			# for option 3 in which has to be mature
			if cdevolveans == '2_mat' and outpool['mature'] == 0:
				differentialmortality = 0.0
			else:
				# Call 2-locus selection model
				differentialmortality = Do2LocusSelection(fitvals,outpool['genes'][0:4],location)
								
	# CDEVOLVE - DoHindexSelection
	# ----------------------------
	elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= burningen_cdevolve) and (timecdevolve.find(OutorBack) != -1):
		if len(timecdevolve.split(':')) > 1: # User indicated age check						
			# Check if individual's age matches user specified selection age
			if outpool['age'] != int(timecdevolve.split(':')[1]):
				differentialmortality = 0.0
			else:
				# Call Hindex selection model
				differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[location])						
		# Run differential mortality check for all ages
		else:
			differentialmortality = DoHindexSelection(cdevolveans,outpool['hindex'],patchvals[location])
											
	# CDEVOLVE - DoFSelection
	# -----------------------
	elif (cdevolveans.split('_')[0] == 'F') and (gen >= burningen_cdevolve) and (timecdevolve.find(OutorBack) != -1):
		if len(timecdevolve.split(':')) > 1: # User indicated age check	
			# Check if individual's age matches user specified selection age
			if outpool['age'] != int(timecdevolve.split(':')[1]):
				differentialmortality = 0.0
			else:
				# Call 2-locus selection model
				differentialmortality = DoFSelection(fitvals,outpool['genes'],location,EHom,cdevolveans)
		else:
			# Call 2-locus selection model
			differentialmortality = DoFSelection(fitvals,outpool['genes'],location,EHom,cdevolveans)
										
	# CDEVOLVE - DoFHindexSelection (inbreeding and outbreeding)
	# ----------------------------------------------------------
	elif (cdevolveans.split('_')[0] == 'FHindex') and (gen >= burningen_cdevolve) and (timecdevolve.find(OutorBack) != -1):
		if len(timecdevolve.split(':')) > 1: # User indicated age check							
			# Check if individual's age matches user specified selection age
			if outpool['age'] != int(timecdevolve.split(':')[1]):
				differentialmortality = 0.0
			else:
				# Call Hindex selection model
				differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],location,EHom,cdevolveans,outpool['hindex'],patchvals[location])
		else:
			# Call Hindex selection model
			differentialmortality = DoFHindexSelection(fitvals,outpool['genes'],location,EHom,cdevolveans,outpool['hindex'],patchvals[location])
	
	# CDEVOLVE - DoMLocusSelection
	# ----------------------------
	elif (cdevolveans.split('_')[0] == 'P') and (gen >= burningen_cdevolve) and (timecdevolve.find(OutorBack) != -1):
		if len(timecdevolve.split(':')) > 1: # User indicated age check			
			
			# Check if individual's age matches user specified selection age
			if outpool['age'] != int(timecdevolve.split(':')[1]):
				differentialmortality = 0.0
			else:
				# Call Hindex selection model
				differentialmortality = DoMLocusSelection(outpool['genes'],location,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)
		else:
				# Call Hindex selection model
				differentialmortality = DoMLocusSelection(outpool['genes'],location,cdevolveans,betas_selection,xvars_betas,maxfit,minfit)	
	
	#CDEVOLVE - is on (one of the above) but not equal to or aftern the burningen time
	else:
		differentialmortality = 0.0		
	
	# ------------------------------
	# Check spatial mortality Death
	# ------------------------------						
	# If subpopulation differential mortality is on
	if not isinstance(subpopmort_mat,str):
		# What subpatchID is individual coming from
		from_subpatch = PopTag[isub]													
		# What subpatchID is individual proposing to go to
		to_subpatch = PopTag[location]							
		# If it is dispersing to another subpatchID
		if from_subpatch != to_subpatch:								
			# grab the differential mortality associated with moving into this new subpatchID - from subpatch TO subpatch - cols are TO, read row, then col for correct spot
			differentialmortality_SpatialSubPopMort = subpopmort_mat[int(to_subpatch)-1][int(from_subpatch)-1]
	else: 
		differentialmortality_SpatialSubPopMort = 0.0			
		
	# Calculated and Check Differential Mortality for both selection and spatial mortality
	# ------------------------------------------------------------------------------------
	differentialmortality_Total = 1. - ((1. - differentialmortality) * (1. - differentialmortality_SpatialSubPopMort))
	
	return differentialmortality_Total
	# END::callCDEVOLVE()



# ---------------------------------------------------------------------------------------------------	
def Do1LocusSelection(fitvals,genes,location):
	'''
	Do1LocusSelection()
	This function calculates individual differential mortality, ie,
	viability selection, for the 1-locus selection model.
	'''
	
	# If L0A0|L0A0 -- loci under selection:
	if int(genes[0]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][0])
																
	# If L0A0|L0A1 -- loci under selection:
	elif int(genes[0]) == 1 and int(genes[1]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][1])
																															
	# If L0A1|L0A1 -- loci under selection
	elif int(genes[1]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals[location][2])
		
	# Another genotype
	else:
		differentialmortality = 0.0		
	
	return differentialmortality
	
	# End::Do1LocusSelection()
	
# ---------------------------------------------------------------------------------------------------	
def Do2LocusSelection(fitvals,genes,location):
	'''
	Do2LocusSelection()
	This function calculates individual differential mortality, ie,
	viability selection, for the 2-locus selection model.
	'''
	
	# If L0A0|L0A0|L1A0|L1A0 - AABB -- loci under selection:
	if int(genes[0]) == 2 and int(genes[2]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][0])
									
	# If L0A0|L0A1|L1A0|L1A0 - AaBB -- loci under selection:
	elif int(genes[0]) == 1 and int(genes[1]) == 1 and int(genes[2]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][1])																															
	# If L0A1|L0A1|L1A0|L1A0 - aaBB -- loci under selection
	elif int(genes[1]) == 2 and int(genes[2]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals[location][2])
									
	# If L0A0|L0A0|L1A0|L1A1 - AABb -- loci under selection:
	elif int(genes[0]) == 2 and int(genes[2]) == 1 and int(genes[3]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][3])
									
	# If L0A0|L0A1|L1A0|L1A1 - AaBb -- loci under selection:
	elif int(genes[0]) == 1 and int(genes[1]) == 1 and int(genes[2]) == 1 and int(genes[3]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][4])
																	
	# If L0A1|L0A1|L1A0|L1A1 - aaBb -- loci under selection
	elif int(genes[1]) == 2 and int(genes[2]) == 1 and int(genes[3]) == 1:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals[location][5])
	
	# If L0A0|L0A0|L1A1|L1A1 - AAbb -- loci under selection:
	elif int(genes[0]) == 2 and int(genes[3]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][6])
									
	# If L0A0|L0A1|L1A1|L1A1 - Aabb -- loci under selection:
	elif int(genes[0]) == 1 and int(genes[1]) == 1 and int(genes[3]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals[location][7])
																															
	# If L0A1|L0A1|L1A1|L1A1 - aabb -- loci under selection
	elif int(genes[1]) == 2 and int(genes[3]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals[location][8])
	
	# Another genotype
	else:
		differentialmortality = 0.0	
	'''
	# Define a mapping of gene combinations to fitness values
	gene_map = {
		(2, 2, 0, 0): 0,  # L0A0|L0A0|L1A0|L1A0 - AABB
		(1, 1, 2, 0): 1,  # L0A0|L0A1|L1A0|L1A0 - AaBB
		(0, 2, 2, 0): 2,  # L0A1|L0A1|L1A0|L1A0 - aaBB
		(2, 0, 1, 1): 3,  # L0A0|L0A0|L1A0|L1A1 - AABb
		(1, 1, 1, 1): 4,  # L0A0|L0A1|L1A0|L1A1 - AaBb
		(0, 2, 1, 1): 5,  # L0A1|L0A1|L1A0|L1A1 - aaBb
		(2, 0, 1, 2): 6,  # L0A0|L0A0|L1A1|L1A1 - AAbb
		(1, 1, 0, 2): 7,  # L0A0|L0A1|L1A1|L1A1 - Aabb
		(0, 2, 0, 2): 8   # L0A1|L0A1|L1A1|L1A1 - aabb
	}

	# Check the gene combination and grab the corresponding fitness value
	key = tuple(int(genes[i]) for i in range(4))  # Create a tuple of the first 4 gene values
	differentialmortality = float(fitvals[location][gene_map.get(key, 9)])  # Default to 0.0 if not found in the map
	'''
	return differentialmortality
	
	# End::Do2LocusSelection()

# ---------------------------------------------------------------------------------------------------	
def DoHindexSelection(cdevolveans,hindex,X):
	'''
	DoHindexSelection()
	This function calculates individual differential mortality, based on the individuals Hindex, temperature or environment at location based on a Gaussian.
	'''
	
	if hindex != -9999:
		# Gaussian
		# Example Hindex_linear_0.231:0.77
		# ---------
		if cdevolveans.split('_')[1] == 'Gauss' or cdevolveans.split('_')[1] == 'gauss':
			# Get parameters	
			pars = cdevolveans.split('_')[2].split(':')
			min_temp = float(pars[0])
			max_temp = float(pars[1])
			C = float(pars[2])
			min_ParentHindex = float(pars[4])
			max_ParentHindex = float(pars[5])
			
			# Check min and max parent Hindex and get p value
			if (hindex <= min_ParentHindex) or (hindex >= max_ParentHindex):
				p = 1.0
			else:
				p = float(pars[3])	
			
			# Get fitness value
			fitness = p * np.exp(-((X - (min_temp + (max_temp-min_temp)*hindex))**2/(2.*C**2)))
			
		# Parabolic
		# ---------
		elif cdevolveans.split('_')[1] == 'Para' or cdevolveans.split('_')[1] == 'para':
			# Get parameters	
			pars = cdevolveans.split('_')[2].split(':')
			p = float(pars[0])
			h = float(pars[1])
			k = float(pars[2])
			
			# Get fitness value
			fitness = k + ((hindex - h)**2 / (4 * p))
		
		# Step
		# ----
		elif cdevolveans.split('_')[1] == 'Step' or cdevolveans.split('_')[1] == 'step':
			# Get parameters
			pars = cdevolveans.split('_')[2].split(':')
			p = float(pars[0])
			h = float(pars[1])
			k = float(pars[2])

			# Get fitness value
			if hindex <= p:
				fitness = h
			else:
				fitness = k	

		# Linear
		# ---------
		elif cdevolveans.split('_')[1] == 'Linear' or cdevolveans.split('_')[1] == 'linear':
			# Get parameters	
			pars = cdevolveans.split('_')[2].split(':')
			p = float(pars[0])
			h = float(pars[1])
			#k = float(pars[2])
			
			# Get fitness value
			fitness = h + (p * hindex)
			
		# Error
		# -----
		else:
			print('CDEvolve answer Hindex specified, and either Gauss, Para, Linear, or Step must be specified.')
			sys.exit(-1)	
		
		# Get mortality value
		differentialmortality = 1. - fitness
		# Check for range 0-1
		if differentialmortality > 1:
			differentialmortality = 1.
		if differentialmortality < 0:
			differentialmortality = 0.
	# Hindex was -9999
	# ----------------
	else:
		differentialmortality = 0. # assume no differential mortality
	
	return differentialmortality
	
	# End::DoHindexSelection()
		
# ---------------------------------------------------------------------------------------------------	
def DoFSelection(fitvals,genes,location,EHom,cdevolveans): 
	'''
	DoFSelection()
	This function calculates individual FHi (Keller et al.), then differential mortality, given equation specified. 
	# From cdpop offspring,fitvals1,tempfreegrid,iteminlist,loci,cdevolveans
	'''
	
	no2s = np.array(genes==2).sum(axis=0)
	no1s = np.array(genes==1).sum(axis=0)
	no0s = np.array(genes==0).sum(axis=0)
	OHomi = (no2s + no0s) / 2.
	#Fhi = (OHomi - EHom) / ((len(genes)/2.) - EHom) # -1 and -1
	Fhi = OHomi/(len(genes)/2) # Calculating observed homozygosity
	
	# Slope/Intercept for Linear/Logistic functions
	m = float(cdevolveans.split('_')[2].split(':')[0])
	bint = float(cdevolveans.split('_')[2].split(':')[1])
	
	# If Linear function
	# Example F_logistic_6.71:-7.025
	
	if cdevolveans.split('_')[1] == 'linear' or cdevolveans.split('_')[1] == 'Linear':			
		differentialmortality = m * Fhi + bint
	# If Logistic function - exp(b + m*x)/ (1 + exp(b + m*x))
	if cdevolveans.split('_')[1] == 'logistic' or cdevolveans.split('_')[1] == 'Logistic':
		differentialmortality = np.exp(bint + m * Fhi) / (1. + np.exp(bint + m * Fhi))	
	
	#differentialmortality = 1. - survival
	# Check for survival/mortality values outside the 0-1 range
	if differentialmortality > 1:
		differentialmortality = 1.
	if differentialmortality < 0:
		differentialmortality = 0.
	return differentialmortality
	
	# End::DoFSelection()
	
# ---------------------------------------------------------------------------------------------------	
def DoFHindexSelection(fitvals,genes,location,EHom,cdevolveans,hindex,X): 
	'''
	DoFHindexSelection()
	This function calculates individual FHi (Keller et al.), then differential mortality, given equation specified. 
	This function calculates individual differential mortality, based on the individuals Hindex, temperature or environment at location based on a Gaussian.
	averages 2 differential mortalities
	'''
		
	# F - Selction Portion of Code
	no2s = np.array(genes==2).sum(axis=0)
	no1s = np.array(genes==1).sum(axis=0)
	no0s = np.array(genes==0).sum(axis=0)
	OHomi = (no2s + no0s) / 2.
	#Fhi = (OHomi - EHom) / ((len(genes)/2.) - EHom) # -1 and -1
	Fhi = OHomi/(len(genes)/2)
	
	# Slope/Intercept for Linear/Logistic functions
	m = float(cdevolveans.split('_')[2].split(':')[0])
	bint = float(cdevolveans.split('_')[2].split(':')[1])
	# If Linear function
	# Example FHindex_logistic_6.7:-7.025_linear_0.231:0.77
	if cdevolveans.split('_')[1] == 'linear' or cdevolveans.split('_')[1] == 'Linear':			
		differentialmortality_F = m * Fhi + bint
	# If Logistic function - exp(b + m*x)/ (1 + exp(b + m*x))
	if cdevolveans.split('_')[1] == 'logistic' or cdevolveans.split('_')[1] == 'Logistic':
		differentialmortality_F = np.exp(bint + m * Fhi) / (1. + np.exp(bint + m * Fhi))	
	# Check for survival/mortality values outside the 0-1 range
	if differentialmortality_F > 1:
		differentialmortality_F = 1.
	if differentialmortality_F < 0:
		differentialmortality_F = 0.
		
	# Hindex - Selection Portion of code
	if hindex != -9999:
		# Gaussian
		# ---------
		if cdevolveans.split('_')[3] == 'Gauss' or cdevolveans.split('_')[3] == 'gauss':
			# Get parameters	
			#pars = cdevolveans.split('_')[2].split(':')
			min_temp = float(cdevolveans.split('_')[4].split(':')[0])
			max_temp = float(cdevolveans.split('_')[4].split(':')[1])
			C = float(cdevolveans.split('_')[4].split(':')[2])
			min_ParentHindex = float(cdevolveans.split('_')[4].split(':')[4])
			max_ParentHindex = float(cdevolveans.split('_')[4].split(':')[5])
			
			# Check min and max parent Hindex and get p value
			if (hindex <= min_ParentHindex) or (hindex >= max_ParentHindex):
				p = 1.0
			else:
				p = float(cdevolveans.split('_')[4].split(':')[3])	
			
			# Get fitness value
			fitness = p * np.exp(-((X - (min_temp + (max_temp-min_temp)*hindex))**2/(2.*C**2)))
			
		# Parabolic
		# ---------
		elif cdevolveans.split('_')[3] == 'Para' or cdevolveans.split('_')[3] == 'para':
			# Get parameters	
			#pars = cdevolveans.split('_')[2].split(':')
			p = float(cdevolveans.split('_')[4].split(':')[0])
			h = float(cdevolveans.split('_')[4].split(':')[1])
			k = float(cdevolveans.split('_')[4].split(':')[2])
			
			# Get fitness value
			fitness = k + ((hindex - h)**2 / (4 * p))
		
		# Step
		# ----
		elif cdevolveans.split('_')[3] == 'Step' or cdevolveans.split('_')[3] == 'step':
			# Get parameters
			#pars = cdevolveans.split('_')[2].split(':')
			p = float(cdevolveans.split('_')[4].split(':')[0])
			h = float(cdevolveans.split('_')[4].split(':')[1])
			k = float(cdevolveans.split('_')[4].split(':')[2])

			# Get fitness value
			if hindex <= p:
				fitness = h
			else:
				fitness = k		
			
		# Linear
		# ---------
		elif cdevolveans.split('_')[3] == 'Linear' or cdevolveans.split('_')[3] == 'linear':
			# Get parameters	
			#pars = cdevolveans.split('_')[2].split(':')
			p = float(cdevolveans.split('_')[4].split(':')[0])#float(pars[0])
			h = float(cdevolveans.split('_')[4].split(':')[1])#float(pars[1])
			#k = float(pars[2])
			
			# Get fitness value
			fitness = h + (p * hindex)
					
		# Error
		# -----
		else:
			print('CDEvolve answer Hindex specified, and either Gauss, Para, or Step must be specified.')
			sys.exit(-1)	
		
		# Get mortality value
		differentialmortality_H = 1. - fitness
		# Check for survival/mortality values outside the 0-1 range
		if differentialmortality_H > 1:
			differentialmortality_H = 1.
		if differentialmortality_H < 0:
			differentialmortality_H = 0.
	# Hindex was -9999
	# ----------------
	else:
		differentialmortality_H = 0. # assume no differential mortality
	
	#return (differentialmortality_F + differentialmortality_H) / 2.
	return (1.-(1.-differentialmortality_F)*(1.-differentialmortality_H))
	
	# End::DoFHindexSelection()

# ---------------------------------------------------------------------------------------------------	
def DoMLocusSelection(genes,iteminlist,cdevolveans,betas,xvars,maxfit,minfit):
	'''
	DoMLocusSelection()
	This function calculates offsprings differential mortality, given the individuals genotype, betas, and xvariables supplied for the linear additive model. 
	'''
	
	# Get individuals genes under selection
	# -----------------------------
	# Assume the first N loci are under selection (and the first A alleles under selection)
	indexto = int(cdevolveans.split('_')[2].split('L')[1]) * int(cdevolveans.split('_')[3].split('A')[1])
	selgenes = genes[0:indexto]
	
	# Apply linear model
	atspot_xvars = xvars[iteminlist] # X vars at this grid spot
	
	linmodel = [] # [xvar][loci][allele]
	# Loop through each environment
	for ixvar in range(int(cdevolveans.split('_')[1].split('X')[1])):
		xvar = float(atspot_xvars[ixvar])
		thesebetas = sum(betas[ixvar],[])
		# Loop through the alleles
		for iall in range(len(selgenes)):
			thisbeta = thesebetas[iall]
			thisallele = selgenes[iall]
			if cdevolveans.split('_')[4] == 'ModelY': # Code 1,0
				if thisallele == 2:
					thisallele = 1 # Change second copy to 1
			# allele X environment X beta effect
			linmodel.append(xvar * thisbeta * thisallele)		
	
	# Add the beta not
	linmodel.append(betas[-1])
	
	# Get the max and min values for rescaling 
	thismaxfit = maxfit[0] 
	thisminfit = minfit[0]
	
	# For the case in which the population fixated and all the same genotypes
	if thismaxfit - thisminfit == 0:
		Fitness = 1.
	else:
		# Rescale GXE calculation by max/min 
		Fitness = (sum(linmodel) - thisminfit) / (thismaxfit - thisminfit)
	
	# Fitness could be less than 0 because of rescaling cases sum(linmodel) < thisminfit
	if Fitness < 0:
		Fitness = 0.
	# Fitness could be greater than 1 because of rescale cases sum(linmodel) > thismaxfit
	if Fitness > 1:
		Fitness = 1.
	
	# Add the linear model together and logit
	#Fitness = np.exp(sum(linmodel)) / (1. + np.exp(sum(linmodel)))
	
	# Convert fitness to differential mortality - 1 - Finess
	differentialmortality = 1. - Fitness
		
	return differentialmortality
	
	# End::DoMLocusSelection()


# ---------------------------------------------------------------------------------------------------	 
def GetMetrics(SubpopIN,K,Population,K_track,loci,alleles,gen,Ho,Alleles,He,p1,p2,q1,q2,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,N_Age,sizecall,size_mean,ClassSizes_Mean,ClassSizes_Std,N_Class,packans,RDispersers,IDispersers,xvars_betas,betas_selection,maxfit,minfit,cdevolveans,disease_vars,DiseaseStates_pop,DiseaseStates_EnvRes ):
	'''
	GetMetrics()
	This function summarizes the genotypes and
	produces genetic metrics.
	Ho - Observed heterozygosity per generation
	He - Expected heterozygoisty per generation
	Alleles - Total number of unique alleles in genotype*individuals
	'''
	
	# List for total, left, and right
	unique_alleles = Alleles
	
	# Get allele location as sequence from alleles array
	allele_numbers = np.asarray(list(range(alleles[0])) * loci) # Assumes same number of alleles per loci
	allele_numbers = []
	for iall in alleles:
		allele_numbers.append(list(range(iall)))
	allele_numbers = np.asarray(sum(allele_numbers,[]))
		
	# Get length of classes and size bins if more than one classfile, then bin from min to max
	# This was in previous < 1.17 versions, assume binning based on first classvars now.
	if sizecall == 'Y' and packans.split('_')[0] != 'logistic':
		# Get the middles for finding closest values
		#size_mean_middles = np.asarray(size_bin)[1:] - np.diff(np.asarray(size_bin).astype('f'))/2
		size_bin = size_mean[0][0]
		size_mean_middles = np.asarray(size_bin)[1:] - np.diff(np.asarray(size_bin).astype('f'))/2		
	classno = len(size_mean[0][0])
	
	# Add spots for Age tracking
	N_Age.append([])
	N_Class.append([])
	AgeSizes_Mean.append([])
	AgeSizes_Std.append([])
	ClassSizes_Mean.append([])
	ClassSizes_Std.append([])
	N_Age[gen] = [[] for x in range(0,classno)]
	N_Class[gen] = [[] for x in range(0,classno)]
	AgeSizes_Mean[gen] = [[] for x in range(0,classno)]
	AgeSizes_Std[gen] = [[] for x in range(0,classno)]
	ClassSizes_Mean[gen] = [[] for x in range(0,classno)]
	ClassSizes_Std[gen] = [[] for x in range(0,classno)]
			
	# Extract the genes information from SubpopIN
	tempgenes = SubpopIN[0]['genes'] # start 
	tempgenesPop = []
	Population.append([]) # Storage numbers - add spot for generation
	K_track.append([]) # Storage numbers - add spot for generation
	# Create a list to store the subpopulation grid number information
	subgrids = []
	all_freq_sub = []
	ho_count_sub = []
	ho_sub = []
	all_freq_sq_sub = []
	homozygosity_sub = []
	he_sub = []
	sumsubpopsHo = []
	sumsubpopsHe = []
	alleles_sub = []
	Residors.append([])
	Strayers1.append([])
	Strayers2.append([])
	Immigrators.append([])
	PopSizes_Mean.append([])
	PopSizes_Std.append([])	
	RDispersers.append([])
	IDispersers.append([])
	DiseaseStates_pop.append([])
	DiseaseStates_EnvRes.append([])
	
	# For each subopulation	
	for isub in range(len(K)):
		# -------------------------------------
		# Get total genes array for all subpops
		# -------------------------------------
		# Cast genes as an numpy array as byte type
		if isub != 0:
			tempgenes = np.concatenate((tempgenes,SubpopIN[isub]['genes']),axis=0)
		tempgenesPop.append(SubpopIN[isub]['genes'])
		
		# -------------------------------------
		# Add information to Population tracker
		# -------------------------------------
		Population[gen].append(len(SubpopIN[isub]))
		K_track[gen].append(K[isub])
		PopSizes_Mean[gen].append(np.mean(SubpopIN[isub]['size']))
		PopSizes_Std[gen].append(np.std(SubpopIN[isub]['size']))
		# Extract the disperser type
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'R' in val])
		Residors[gen].append(len(tempname))
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'I' in val])
		Immigrators[gen].append(len(tempname))
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'S' in val])
		Strayers1[gen].append(len(tempname))
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'Z' in val])
		Strayers2[gen].append(len(tempname))
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'ID' in val])
		IDispersers[gen].append(len(tempname))
		tempname = np.asarray([i for i, val in enumerate(SubpopIN[isub]['name']) if 'RD' in val])
		RDispersers[gen].append(len(tempname))
		# Extract count of disease states
		state_counts = np.bincount(SubpopIN[isub]['states'], minlength=disease_vars['noStates'][isub])
		DiseaseStates_pop[gen].append(state_counts)				
		#pdb.set_trace()
		if disease_vars['ImpDisease'] != 'N' and gen == 0:
			DiseaseStates_EnvRes[gen].append(float(disease_vars['PathLoad'][isub].split('|')[gen]))
		elif disease_vars['ImpDisease'] != 'N' and gen > 0:
			DiseaseStates_EnvRes[gen].append(round(disease_vars['PathLoad'][isub],2))
		else:
			DiseaseStates_EnvRes[gen].append(0)	
		#pdb.set_trace()								
		# Switch here for size or age control - # Note that first size classes used for binning
		if sizecall == 'Y' and packans.split('_')[0] != 'logistic': 
			age_adjusted = np.searchsorted(size_mean_middles, SubpopIN[isub]['size'])
		else:
			# Count up each uniages
			age_adjusted = SubpopIN[isub]['age']
		
		# Tracking age N
		for iage in range(len(AgeSizes_Mean[gen])):
			sizeindex = np.where(age_adjusted==iage)[0]
			ageindex = np.where(SubpopIN[isub]['age'] == iage)[0]
			if len(sizeindex) == 0:
				ClassSizes_Mean[gen][iage].append([0])
				N_Class[gen][iage].append(0)				
			else:
				ClassSizes_Mean[gen][iage].append(SubpopIN[isub][sizeindex]['size'].tolist())
				N_Class[gen][iage].append(len(SubpopIN[isub][sizeindex]['size'].tolist()))
			if len(ageindex) == 0:
				AgeSizes_Mean[gen][iage].append([0])
				N_Age[gen][iage].append(0)
			else:				
				AgeSizes_Mean[gen][iage].append(SubpopIN[isub][ageindex]['size'].tolist())	
				N_Age[gen][iage].append(len(SubpopIN[isub][ageindex]['size'].tolist()))
		
		# Special case where age class is greater than lastage
		sizeindex = np.where(age_adjusted > iage)[0]
		if len(sizeindex) == 0:
			ClassSizes_Mean[gen][iage].append([0])
			N_Class[gen][iage].append(0)
		else: # Add them to last class
			ClassSizes_Mean[gen][iage].append(SubpopIN[isub][sizeindex]['size'].tolist())
			N_Class[gen][iage].append(len(SubpopIN[isub][sizeindex]['size'].tolist()))
		ageindex = np.where(SubpopIN[isub]['age'] > iage)[0]
		if len(ageindex) == 0:
			AgeSizes_Mean[gen][iage].append([0])
			N_Age[gen][iage].append(0)
		else: # Add them to last class
			AgeSizes_Mean[gen][iage].append(SubpopIN[isub][ageindex]['size'].tolist())	
			N_Age[gen][iage].append(len(SubpopIN[isub][ageindex]['size'].tolist()))
		
		# Temp storage for GetMetrics()
		subgrids.append([])
		all_freq_sub.append([])
		ho_count_sub.append([])
		ho_sub.append([])
		all_freq_sq_sub.append([])
		homozygosity_sub.append([])
		he_sub.append([])
		alleles_sub.append([])
	# Add Population total
	Population[gen].insert(0,sum(Population[gen]))
	K_track[gen].insert(0,sum(K))
	DiseaseStates_pop[gen].insert(0,np.sum(np.asarray(DiseaseStates_pop[gen]),axis=0).tolist())
	DiseaseStates_EnvRes[gen].insert(0,np.sum(np.asarray(DiseaseStates_EnvRes[gen]),axis=0).tolist())
	
	# Age tracking	
	for iage in range(len(AgeSizes_Mean[gen])):
		tempagesize = np.asarray(sum(AgeSizes_Mean[gen][iage],[]))
		tempagesize = tempagesize[np.where(tempagesize != 0)[0]]
		tempsize = np.asarray(sum(ClassSizes_Mean[gen][iage],[]))
		tempsize = tempsize[np.where(tempsize != 0)[0]]
		N_Age[gen][iage] = sum(N_Age[gen][iage])
		N_Class[gen][iage] = sum(N_Class[gen][iage])		
		AgeSizes_Mean[gen][iage] = np.mean(tempagesize)
		AgeSizes_Std[gen][iage] = np.std(tempagesize)
		ClassSizes_Mean[gen][iage] = np.mean(tempsize)
		ClassSizes_Std[gen][iage] = np.std(tempsize)
	
	# ------------------------------------------------
	# Get Total information
	# ------------------------------------------------
	# Cast genes as an numpy array as byte type
	genes_array_woNA = np.asarray(tempgenes,dtype='float')
	
	# The total number of alleles
	total_alleles = len(allele_numbers)
	
	# Get unique number of subpops
	nosubpops = len(K)
	
	# And then get the number of filled grids
	filledgrids = Population[gen][0]		
		
	#Calculate the number of homogenous alleles for total
	ho_count_tot = np.array(genes_array_woNA==2).sum() # Note: this is not calculated for Ho
	
	# Get allele frequency for total # Calculate the observed het for total
	if filledgrids != 0:
		all_freq_tot = np.asarray(np.nansum(genes_array_woNA,axis=0),dtype = 'float')
		#all_freq_tot = np.asarray(np.nansum(genes_array_woNA,axis=0),dtype = 'float').reshape(total_alleles)
		all_freq_tot = all_freq_tot/(2*filledgrids)
		ho_tot = (float(filledgrids*loci - ho_count_tot)/(loci*filledgrids))		
	else:
		all_freq_tot = np.zeros(total_alleles,float)
		ho_tot = 0.0
	
	# Append Ho information (Observed Het)
	Ho.append([ho_tot])
	
	# Create an array to fill up with allele frequencies - only for total
	all_freq_list = np.zeros((total_alleles,2))		
	all_freq_list[:,0] = allele_numbers
	all_freq_list[:,1] = all_freq_tot
	# Get the sqare of the allele frequency for total
	all_freq_sq_tot = all_freq_tot**2
	# Get total number of alleles
	alleles_tot = np.array(all_freq_tot>0.).sum()
	# Append allele total information
	unique_alleles.append([alleles_tot])	
	# Calculate the homozygosity for total populations
	homozygosity_tot = sum(all_freq_sq_tot)/loci
	
	# Get allele frequency totals for selection section
	p1.append(all_freq_tot[0])
	p2.append(all_freq_tot[1])
	q1.append(all_freq_tot[2])
	q2.append(all_freq_tot[3])
		
	# Store He for [Total]
	if filledgrids != 0:
		he_tot = (1. - homozygosity_tot)
	else:
		he_tot = 0.0
	# Append He information (Expected Het)
	He.append([he_tot])		
	
	# -----------------Subpop numbers
	for isub in range(nosubpops):
		# Cast genes as an numpy array as byte type
		genes_array_subpop = np.asarray(tempgenesPop[isub],dtype='float')
		# Calculate the number of homogenous alleles in each subpop
		ho_count_sub[isub].append(np.array(genes_array_subpop==2).sum())
		
		if Population[gen][isub+1] != 0:			
			# Get allele frequency for subpopulations
			all_freq_sub[isub].append(np.asarray(np.sum(genes_array_subpop,axis=0),dtype = 'float'))
			all_freq_sub[isub] = all_freq_sub[isub][0]/(2*Population[gen][isub+1])
			# Calculate the observed het in each subpop
			ho_sub[isub].append((float(Population[gen][isub+1]*loci - ho_count_sub[isub][0])/(loci*Population[gen][isub+1])))
		else:
			# Get allele frequency for subpopulations
			all_freq_sub[isub].append(np.zeros(total_alleles,float))
			all_freq_sub[isub] = all_freq_sub[isub][0]
			# Calculate the observed het in each subpop
			ho_sub[isub].append(0.0)
		# Append Ho information (Observed Het)
		Ho[gen].append(ho_sub[isub][0])
		
		# Get the square of the allele frequency for subpops
		all_freq_sq_sub[isub].append(all_freq_sub[isub]**2)
		# Calculate the homozygosity for subpopulations
		homozygosity_sub[isub].append(sum(all_freq_sq_sub[isub][0])/loci)
		# Get the total number of alleles in each subpop
		alleles_sub[isub].append(np.array(all_freq_sub[isub]>0.).sum())
		# Append allele total information
		unique_alleles[gen].append(alleles_sub[isub][0])
		
		# Store He for subpopulations
		if Population[gen][isub+1] != 0:
			he_sub[isub].append(1. - homozygosity_sub[isub][0])
		else:
			he_sub[isub].append(0.0)
		# Append He information (Expected Het)
		He[gen].append(he_sub[isub][0])	
		
	# ----------------------------------------------
	# Max/min GXE local rescaling for fitness values
	# ----------------------------------------------
	if cdevolveans.split('_')[0] == 'P':
		
		# For Global option, at the first generation only, get all possible genotypes
		if gen == 0:
			# Calculate the total genotype space - combination replacement CR(possible alleles,2)^possible loci
			#total_genotypespace = ( math.factorial(int(cdevolveans.split('_')[3].split('A')[1]) + 2 - 1) / (math.factorial(2) * math.factorial(int(cdevolveans.split('_')[3].split('A')[1]) - 1)) )**int(cdevolveans.split('_')[2].split('L')[1])
			
			# Grab XVars as array - check for cdclimate bars
			if len(xvars_betas[0][0].split('|')) > 1:
				xvars_betas_temp = []
				# Loop through each isub
				for isub in range(len(xvars_betas)):
					xvars_betas_temp.append([xvars_betas[isub][0].split('|')[0]])
				xvars_betas_temp = np.asarray(xvars_betas_temp,dtype=float)
			else:
				xvars_betas_temp = np.asarray(xvars_betas,dtype=float)
			
			if cdevolveans.split('_')[4] == 'ModelY': # Code 1,0
				multiFactor = 1
			else:
				multiFactor = 2
			
			max_linmodel2 = []
			min_linmodel2 = []
			# For each patch spot, grab the Xvar values
			for igrid in range(nosubpops):
				grid_xvars = xvars_betas_temp[igrid] # This grids X variables
				max_linmodel2.append([])
				min_linmodel2.append([])
				# Loop through each X var and max/min each calculation
				for ivar in range(len(grid_xvars)):
					Xvar = grid_xvars[ivar]
					betas = np.asarray(betas_selection[ivar])
					
					# if the Xvar is Zero:
					if Xvar == 0:
						# zeros out this calculation
						max_linmodel2[igrid].append(0)
						min_linmodel2[igrid].append(0)
						
					# If Xvar is positive	
					elif Xvar > 0: 
						# For max, grab max beta value (positive) at each locus and multiply by 2 or 1
						max_linmodel2[igrid].append(Xvar*np.sum(np.max(betas,1)*multiFactor))
						# For min, grab min beta value at each locus and multiply by 2 or 1
						min_linmodel2[igrid].append(Xvar*np.sum(np.min(betas,1)*multiFactor))
						
					# if Xvar is negative	
					elif Xvar < 0: 
						# For max, grab min beta value at each locus and multiply by 2 or 1
						max_linmodel2[igrid].append(Xvar*np.sum(np.min(betas,1)*multiFactor))
						# For min, grab max beta value (positive) at each locus and multiply by 2 or 1 
						min_linmodel2[igrid].append(Xvar*np.sum(np.max(betas,1)*multiFactor))	
					
				max_linmodel2[igrid] = sum(max_linmodel2[igrid])
				min_linmodel2[igrid] = sum(min_linmodel2[igrid])
				
				maxfit.append(max(max_linmodel2))
				minfit.append(min(min_linmodel2))
					
	#End::GetMetrics()
	
# ---------------------------------------------------------------------------------------------------	 
def InheritGenes(gen,offspring,loci,muterate,mtdna,mutationans,K,dtype,geneswap,allelst,assortmateModel,noalleles,plasticans,cdevolveans,egg_add_call,egg_add):
	'''
	InheritGenes()
	Pass along gentic information to survived offspring from parents
	Input: offspring [femalegenes,malegenes,NatalPop,EmiPop,ImmiPop,age,sex,size,name]
	Output: SubpopIN_Age0	[NatalPop,EmiPop,ImmiPop,age,sex,size,name,genes]		
	'''		
	
	# Get Plastic loci region
	# -----------------------
	if plasticans != 'N':		
		# If cdevolve is on
		if cdevolveans != 'N':
			# Then the first l loci are for selection, next for plastic region
			if cdevolveans.split('_')[0] == 'P': # This is for multilocus selection, not currently implemented, to be moved over from cdpop
				selloci = int(cdevolveans.split('_')[2].split('L')[1])
			elif cdevolveans in {'1', 'M', 'G', '1_mat', '1_G_ind', '1_G_link', 'stray', 'Hindex', 'FHindex'}:
				selloci = 1
			elif cdevolveans in {'2', 'MG', 'MG_ind', 'MG_link', '2_mat'}:
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
		plaloci_index = list(range(selloci*2,selloci*2+plaloci*2))
	else:
		plaloci_index = [-9999]
	
	# Create list for appending
	Age0_keep = []
	
	# If there are offspring
	if len(offspring) != 0:
	
		# Begin loop through offspring
		for i in range(len(offspring)):
			
			# Check for geneswap time
			if gen >= geneswap and geneswap != 'N':
				
				# Temp storage for i's mother's and father's genes, then offspring genes
				mothergenes=offspring[i]['Mother']
				fathergenes=offspring[i]['Father']
				
				# Temp genes storage for offspring
				offgenes = np.zeros(len(fathergenes),dtype =int)
				# Allele indices
				alleles = np.arange(len(mothergenes)) 
				
				# Loop through each locus
				for iloci in range(loci):
					# Allele indices to sample from - index into offgenes
					possiblealleles = alleles[sum(noalleles[:iloci]):sum(noalleles[:iloci+1])]
					
					# Determine if this locus is in the plastic region
					is_plastic_region = (plasticans != 'N') and (plaloci_index == list(range(iloci*2, iloci*2 + plaloci*2)))
					
					# Get alleles for parents
					F2 = np.where(fathergenes[possiblealleles] == 2)[0] 
					F1 = np.where(fathergenes[possiblealleles] == 1)[0] 
					M2 = np.where(mothergenes[possiblealleles] == 2)[0]
					M1 = np.where(mothergenes[possiblealleles] == 1)[0]

					FALL = np.concatenate((F2, F2, F1))
					MALL = np.concatenate((M2, M2, M1))
					
					if len(FALL) == 0 or len(MALL) == 0: 
						print('FALL or MALL in inheritgenes() is equal to 0.')
						sys.exit(-1)
				
					# Process plastic or non-plastic regions
					lastloci = True if iloci == loci - 1 else False # Check for mtdna
					if not is_plastic_region:
						# Non-plastic region: diploid inheritance
						sample_allele(FALL, possiblealleles, offgenes, mothergenes, mtdna, lastloci, False)
						sample_allele(MALL, possiblealleles, offgenes, mothergenes, mtdna, lastloci, False)
					else:
						# Plastic region: offspring does not inherit 2 or plastic tracking from parents
						sample_allele(FALL, possiblealleles, offgenes, mothergenes, mtdna, lastloci, True)
						sample_allele(MALL, possiblealleles, offgenes, mothergenes, mtdna, lastloci, True)					
					
					# Additional genetic checks on this locus - mutation models
					#----------------------------------------------------------
					if muterate != 0.0:
						for iloci in range(loci): # Loop through loci
							mutationrandnos = np.random.sample(2) # Get two random numbers for checking
							
							# Get current alleles at this locus 2 or ones
							thisloci = possiblealleles[np.where(offgenes[possiblealleles] != 0)[0]]  
							
							# Handle homozygous case
							if len(thisloci) == 1:
								thisloci = np.concatenate((thisloci, thisloci), axis=0)
								
							# Process each allele in this locus - 
							for iall in range(2):
								if mutationrandnos[iall] < muterate:
									if mutationans == 'random':
										# Random mutation model
										offgenes[thisloci[iall]] -= 1 # Remove that allele
										new_allele = np.random.choice(possiblealleles[possiblealleles != thisloci[iall]])
										offgenes[new_allele] += 1						

									elif mutationans == 'forward':
										# Forward mutation model
										if iall == 0 and offgenes[possiblealleles][0] != 0:
											move_allele(thisloci[iall], 'forward', possiblealleles, offgenes)

									elif mutationans == 'backward':
										# Backward mutation model
										if iall == 1 and offgenes[possiblealleles][1] != 0:
											move_allele(thisloci[iall], 'backward', possiblealleles, offgenes)

									elif mutationans == 'forwardbackward':
										# Random forward/backward mutation
										direction = 'backward' if np.random.uniform() < 0.5 and iall == 1 else 'forward'
										move_allele(thisloci[iall], direction, possiblealleles, offgenes)

									elif mutationans == 'plastic':
										# Special mutation model for plastic loci
										if cdevolveans != 'N':
											print('Check mutation model, cdevolve, and plastic loci.')
											sys.exit(-1)
										if iloci == 0 and offgenes[possiblealleles][iall] != 0:
											if iall == 0:
												move_allele(thisloci[iall], 'forward', possiblealleles, offgenes)

									else:
										# Unknown mutation model
										print('The mutation model does not exist.')
										sys.exit(-1)
				
				# Calculate hybrid index
				hindex = offspring[i]['M_hindex']/2 + offspring[i]['F_hindex']/2
				
			# If not in geneswap time, then initialize with allelst
			# -----------------------------------------------------
			else:				
				
				# Get genes - For each loci:
				sourcepop = int(offspring[i][egg_add_call])-1
				thisgenefile = np.random.randint(len(allelst[sourcepop]))
				offgenes = [] # Storage
				
				# Loop through each locus
				for j in range(loci):					
					
					# Take a random draw from the w_choice function at jth locus
					rand1 = w_choice_general(allelst[sourcepop][thisgenefile][j])[0]
					rand2 = w_choice_general(allelst[sourcepop][thisgenefile][j])[0]
										
					# Determine alleles (0, 1, or 2) for each possible allele at the locus
					for k in range(len(allelst[sourcepop][thisgenefile][j])):
						if rand1 == rand2:
							tempindall = 2 if k == rand1 else 0  # Homozygous
						else:
							tempindall = 1 if k in {rand1, rand2} else 0  # Heterozygous
																		
						# Append to offgenes)
						offgenes.append(tempindall)
				
				# Special case for assortative mating option 2 (strict) and when hindex is being used - inherit the parents first locus, so that AA gives to AA and aa gives to aa, etc. 
				if assortmateModel == '2':
					offgenes[:noalleles[0]] = offspring[i]['Mother'][:noalleles[0]]
					
				# Special case MtDNA inheritance
				# ------------------
				if mtdna == 'Y':
					offgenes[-1] = offspring[i]['Mother'][-1]
				
				# Calculate hybrid index
				if offgenes[0] == 2: #AA
					hindex = 1.0
				elif offgenes[1] == 2: #aa
					hindex = 0.0
				elif offgenes[0] == 1 and offgenes[1] == 1: #Aa
					hindex = 0.5
				else:
					hindex = -9999			
			
			# Then record new offspring information to Subpop location ['NatalPOp - or mating pop,subpop of mother,NASubpop,EmiCD,ImmiCD,age,sex,size,mataure,name,capture,layeggs,hindex,classfile,speciesID,genes]
			id = offspring[i]['name']
			offpop = offspring[i][egg_add_call]
			if egg_add == 'mating': # Same as offpop
				matedpop = offpop
			else: # First pop, need to look up where mother mated
				matedpop = ''.join(str(num) for num in [int(idig) for idig in offspring[i]['MID'].split('_')[0] if idig.isdigit()])
			
			recd = (matedpop,offpop,'NA',0.0,-9999,offspring[i]['age'],offspring[i]['sex'],offspring[i]['size'],offspring[i]['mature'],offspring[i]['newmature'],offspring[i]['states'],id,offspring[i]['MID'],offspring[i]['FID'],0,0,offspring[i]['layeggs'],hindex,offspring[i]['classfile'],offspring[i]['popID'],offspring[i]['species'],offgenes)
			# Record offspring information to SubpopIN 
			Age0_keep.append(recd)
					
	# Delete offspring variable
	del offspring	
	
	# Return SubpopIN back to array with types
	# Organize type data in SubpopIN - read in from global
	tempN = []
	
	# array and add dtype
	Age0_keep = np.array(Age0_keep,dtype=dtype)
	
	# Return variables from this argument
	return Age0_keep
	
	# End::InheritGenes()

# ---------------------------------------------------------------------------------------------------	
def growInd(Indloc,SubpopIN,sizeLoo,sizeR0,size_1,size_2,size_3,size_4,sizevals,isub,iind,growans,size_mean,gridsample,cdevolveans,gen,burningen_cdevolve,alleles,fitvals,sexchromo):
	'''
	Growth options
	'''	
	# Get age
	Indage = SubpopIN[isub][iind]['age']
	# Get sex and split options if provided
	Indsex = SubpopIN[isub][iind]['sex']
	if Indsex == 'FXX':
		sxspot = 0
	elif Indsex == 'MXY':
		sxspot = 1
	elif Indsex == 'MYY':
		sxspot = 2
	else:
		sxspot = 3
	# Check for cdevolve growth option - get new growth parameters
	if gen >= burningen_cdevolve: # Skip if selection is not on
		# If MG independent
		if cdevolveans == 'MG_ind':
			Indgenes = SubpopIN[isub][iind]['genes']
			# BB
			if Indgenes[sum(alleles[0:1]) + 0 + 1] == 2:
				genespot = 3
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]							
			# Bb
			elif Indgenes[sum(alleles[0:1]) + 0 + 1] == 1 and Indgenes[sum(alleles[0:1]) + 1 + 1] == 1:
				genespot = 4
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]
			# bb
			elif Indgenes[sum(alleles[0:1]) + 1 + 1] == 2:
				genespot = 5
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]
			
		# If MG linked
		elif cdevolveans == 'MG_link':
			Indgenes = SubpopIN[isub][iind]['genes']
			# AA - use BB
			if Indgenes[0] == 2:
				genespot = 3							
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]							
			# Aa - use Bb
			elif Indgenes[0] == 1 and Indgenes[1] == 1:
				genespot = 4
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]
			# aa - use bb
			elif Indgenes[1] == 2:
				genespot = 5
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]
			
		# If just Locus B (Growth)
		elif cdevolveans == 'G':
			Indgenes = SubpopIN[isub][iind]['genes']
			# BB
			if Indgenes[sum(alleles[0:1]) + 0 + 1] == 2:
				genespot = 0
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]							
			# Bb
			elif Indgenes[sum(alleles[0:1]) + 0 + 1] == 1 and Indgenes[sum(alleles[0:1]) + 1 + 1] == 1:
				genespot = 1
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]
			# bb
			elif Indgenes[sum(alleles[0:1]) + 1 + 1] == 2:
				genespot = 2
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]
			
		# If 1_G independent
		elif cdevolveans == '1_G_ind':
			Indgenes = SubpopIN[isub][iind]['genes']
			# BB
			if Indgenes[sum(alleles[0:1]) + 0 + 1] == 2:
				genespot = 3
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]							
			# Bb
			elif Indgenes[sum(alleles[0:1]) + 0 + 1] == 1 and Indgenes[sum(alleles[0:1]) + 1 + 1] == 1:
				genespot = 4
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]
			# bb
			elif Indgenes[sum(alleles[0:1]) + 1 + 1] == 2:
				genespot = 5
				growans = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[2]
				size_1 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[3]
				size_2 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[4]
				size_3 = fitvals[int(Indloc)-1][genespot][sxspot].split(':')[5]
			
		# If 1_G linked
		elif cdevolveans == '1_G_link':
			Indgenes = SubpopIN[isub][iind]['genes']
			# AA - use BB
			if Indgenes[sum(alleles[0:0]) + 0 + 1] == 2:
				genespot = 3
				growans = fitvals[int(Indloc)-1][genespot][0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][2]
				size_1 = fitvals[int(Indloc)-1][genespot][3]
				size_2 = fitvals[int(Indloc)-1][genespot][4]
				size_3 = fitvals[int(Indloc)-1][genespot][5]							
			# Aa - use Bb
			elif Indgenes[sum(alleles[0:0]) + 0 + 1] == 1 and Indgenes[sum(alleles[0:1]) + 1 + 1] == 1:
				genespot = 4
				growans = fitvals[int(Indloc)-1][genespot][0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][2]
				size_1 = fitvals[int(Indloc)-1][genespot][3]
				size_2 = fitvals[int(Indloc)-1][genespot][4]
				size_3 = fitvals[int(Indloc)-1][genespot][5]
			# aa - use bb
			elif Indgenes[sum(alleles[0:0]) + 1 + 1] == 2:
				genespot = 5
				growans = fitvals[int(Indloc)-1][genespot][0]
				sizeLoo = fitvals[int(Indloc)-1][genespot][1]
				sizeR0 = fitvals[int(Indloc)-1][genespot][2]
				size_1 = fitvals[int(Indloc)-1][genespot][3]
				size_2 = fitvals[int(Indloc)-1][genespot][4]
				size_3 = fitvals[int(Indloc)-1][genespot][5]										
			
	# If CDEvolve Answer for maturation values is not operating, still check for splits from PopVars
	if cdevolveans != 'MG_ind' and cdevolveans != 'MG_link' and cdevolveans != 'G' and cdevolveans != '1_G_ind' and cdevolveans != '1_G_link':
		# ------------------sizeLoo		
		if len(sizeLoo.split('~')) == 1:
			sizeLoo = sizeLoo
		elif len(sizeLoo.split('~')) != sexchromo:
			print('Number of growth parameters Loo must match sex_chromo.')
			sys.exit(-1)
		else:
			sizeLoo = sizeLoo.split('~')[sxspot]
		
		# -------------------sizeR0
		if len(sizeR0.split('~')) == 1:
			sizeR0 = sizeR0
		elif len(sizeR0.split('~')) != sexchromo:
			print('Number of growth parameters R0 must match sex_chromo.')
			sys.exit(-1)
		else:
			sizeR0 = sizeR0.split('~')[sxspot]
		
		# ------------------size_1
		if len(size_1.split('~')) == 1:
			size_1 = size_1
		elif len(size_1.split('~')) != sexchromo:
			print('Number of growth parameters must match sex_chromo.')
			sys.exit(-1)
		else:
			size_1 = size_1.split('~')[sxspot]
		
		# ------------------size_2
		if len(size_2.split('~')) == 1:
			size_2 = size_2
		elif len(size_2.split('~')) != sexchromo:
			print('Number of growth parameters must match sex_chromo.')
			sys.exit(-1)
		else:
			size_2 = size_2.split('~')[sxspot]
		
		# ------------------size_3
		if len(size_3.split('~')) == 1:
			size_3 = size_3
		elif len(size_3.split('~')) != sexchromo:
			print('Number of growth parameters must match sex_chromo.')
			sys.exit(-1)
		else:
			size_3 = size_3.split('~')[sxspot]
		
	# -----------------------------
	# Grow based on von Bertalanffy
	# -----------------------------
	if growans == 'vonB':
		grow = float(size_4[int(Indloc) - 1]) / 365. # grow days or proportion of time
		t0 = float(size_3)
		K = float(sizeR0)
		if gridsample == 'Middle': 	# This is the second DoUpdate() - when they are 'Back' ind.csv before DoEmigration()
			newsize = float(sizeLoo) * (1. - np.exp(-K * ((Indage+1) - t0))) # 365 days
			#newsize = float(sizeLoo) * (1. - np.exp(-K * ((Indage+grow) - t0)))
		elif gridsample == 'Sample' or gridsample == 'N': # THis is the third DoUpdate() - when they are 'Out' indSample.csv before DoImmigration()
			newsize = float(sizeLoo) * (1. - np.exp(-K * ((Indage+grow) - t0)))
			#newsize = float(sizeLoo) * (1. - np.exp(-K * ((Indage+1) - t0)))
		if newsize <= 0.:
			#print('Warning: von Bertalanffy growth producing negative values.')
			newsize = 0.
		SubpopIN[isub][iind]['size'] = newsize	
			
	# -------------------------------
	# Grow based on temp fit len/size
	# -------------------------------
	elif growans == 'temperature':
		if sizevals[int(Indloc) - 1] != 'N':
			tempval = float(sizevals[int(Indloc) - 1])
			grow = float(size_4[int(Indloc) - 1]) / 365. # Note size_4 is growdays
			# size_2 = CV				
			int_R = -float(sizeR0) * ((scipy.stats.norm(float(size_1),float(size_2)*float(size_1)).pdf(tempval)) / (scipy.stats.norm(float(size_1),float(size_2)*float(size_1)).pdf(float(size_1))))
			
			L_inc = float(sizeLoo) * (1. - np.exp(int_R * (Indage+1-float(size_3)))) * ((scipy.stats.norm(float(size_1),float(size_2)*float(size_1)).pdf(tempval)) / (scipy.stats.norm(float(size_1),float(size_2)*float(size_1)).pdf(float(size_1))))
			# Get the incremental growth
			
			#L_inc_age = L_inc * np.exp((Indage+1) * int_R) #old ebt version
			L_inc_age = L_inc * np.exp((Indage+1) * (float(sizeR0)*-1))	# new version		
			
			# Update the new size for this individual		
			newsize = SubpopIN[isub][iind]['size'] + (L_inc_age * (grow))
			if newsize <= 0.:
				print('Warning: temperature growth producing negative values.')
				sys.exit(-1)
			SubpopIN[isub][iind]['size'] = newsize				
	
	# ----------------------------------------------
	# Grow based on temperature model but hindex too
	# ----------------------------------------------
	elif growans == 'temperature_hindex':
		if sizevals[int(Indloc) - 1] != 'N':			
			tempval = float(sizevals[int(Indloc) - 1])
			grow = float(size_4[int(Indloc) - 1]) / 365.
							
			int_R = -float(sizeR0) * ((scipy.stats.norm(float(size_1),float(size_2)*float(size_1)).pdf(tempval)) / (scipy.stats.norm(float(size_1),float(size_2)*float(size_1)).pdf(float(size_1))))
			
			# Get the Loo for this HIndex
			Indhindex = SubpopIN[isub][iind]['hindex']
			bothLoo = sizeLoo.split(';')
			if len(bothLoo) != 2: # Error check to make sure user entered correctly
				print('Growth option temperature_hindex specified; growth_Loo should have minimum and maximum Loo values given separated by ;. See user manual.')
				sys.exit(-1)
			sizeLoo_min = float(bothLoo[0])
			sizeLoo_max = float(bothLoo[1])
			sizeLoo_hindex = Indhindex * (sizeLoo_max - sizeLoo_min) + sizeLoo_min
			
			L_inc = float(sizeLoo_hindex) * (1. - np.exp(int_R * (Indage+1-float(size_3)))) * ((scipy.stats.norm(float(size_1),float(size_2)*float(size_1)).pdf(tempval)) / (scipy.stats.norm(float(size_1),float(size_2)*float(size_1)).pdf(float(size_1))))
			# Get the incremental growth
			#L_inc_age = L_inc * np.exp((Indage+1) * int_R)
			L_inc_age = L_inc * np.exp((Indage+1) * (float(sizeR0)*-1))
			# Update the new size for this individual		
			newsize = SubpopIN[isub][iind]['size'] + (L_inc_age * (grow))
			if newsize <= 0.:
				print('Warning: temperature growth producing negative values.')
				sys.exit(-1)
			SubpopIN[isub][iind]['size'] = newsize
	
	# -------------------
	# Grow based on known
	# -------------------
	elif growans == 'known':
		if gridsample == 'N': # Only apply at 3rd DoUpdate
			# Get individuals classfile to use
			natalP = int(SubpopIN[isub][iind]['classfile'].split('_')[0].split('P')[1])
			theseclasspars = int(SubpopIN[isub][iind]['classfile'].split('_')[1].split('CV')[1])
			# Get individuals current size
			currentsize = SubpopIN[isub][iind]['size']
			# Find size closest too
			size_mean_middles = np.asarray(size_mean[0][0])[1:] - np.diff(np.asarray(size_mean[0][0]).astype('f'))/2
			closest_size_index = np.searchsorted(size_mean_middles, currentsize)
			# Move individual to the next size class, first check last class case
			if closest_size_index == len(size_mean[0][0])-1:
				next_size_index = len(size_mean[0][0])-1
			else:
				next_size_index = closest_size_index + 1
			# Set new size
			newsize = size_mean[0][0][next_size_index]
			SubpopIN[isub][iind]['size'] = newsize
			
	# Error check
	else:
		print('Growth options include, vonB, temperature. Check that you have entered the correct formate in growth_option in Popvars.csv field.')
		sys.exit(-1)
	
	#End::growInd()
# ---------------------------------------------------------------------------------	
def matureInd(lastage,SubpopIN,isub,iind,sizeans,age_mature,eggFreq_mu,eggFreq_sd,cdevolveans,fitvals,burningen_cdevolve,gen,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo):
	'''
	Mature, and check egg frequency interval here
	'''
	
	# Get Sex - and index for splits
	Indsex = SubpopIN[isub][iind]['sex']
	if Indsex == 'FXX':
		sxspot = 0
	elif Indsex == 'MXY':
		sxspot = 1
	elif Indsex == 'MYY':
		sxspot = 2
	else:
		sxspot = 3
	
	# Check if the individual is lastage + , then use last index
	if SubpopIN[isub][iind]['age'] >= lastage:
		Indage = lastage-1
	else: # Use the individuals age
		Indage = SubpopIN[isub][iind]['age']
	sizesamp = SubpopIN[isub][iind]['size']
	# Check if becomes mature
	# -----------------------
	if SubpopIN[isub][iind]['mature'] == 0:		
		matval = 0.0 # initialize
		# Check default age/size for maturity
		if Indsex == 'FXX':
			if FXXmat_set != 'N':
				if len(FXXmat_set.split('age')) == 2: # Maturation value set for age
					AgeMature = int(FXXmat_set.split('age')[1])
					if Indage >= AgeMature: # If the age is > than default mature value, then becomes mature.
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
		elif Indsex == 'MXY':
			if MXYmat_set != 'N':
				if len(MXYmat_set.split('age')) == 2: # Maturation value set for age
					AgeMature = int(MXYmat_set.split('age')[1])
					if Indage >= AgeMature: # If the age is > than default mature value, then becomes mature.
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
		elif Indsex == 'MYY':
			if MYYmat_set != 'N':
				if len(MYYmat_set.split('age')) == 2: # Maturation value set for age
					AgeMature = int(MYYmat_set.split('age')[1])
					if Indage >= AgeMature: # If the age is > than default mature value, then becomes mature.
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
		elif Indsex == 'FYY':
			if FYYmat_set != 'N':
				if len(FYYmat_set.split('age')) == 2: # Maturation value set for age
					AgeMature = int(FYYmat_set.split('age')[1])
					if Indage >= AgeMature: # If the age is > than default mature value, then becomes mature.
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
			if sizeans == 'N': # Age control				
				if len(age_mature[Indage].split('~')) == 1:
					matval = float(age_mature[Indage].split('~')[0])
				elif len(age_mature[Indage].split('~')) == sexchromo:
					matval = float(age_mature[Indage].split('~')[sxspot])
				else:
					print('ClassVars age maturation probabilities must be length 1 or length of number of sex_chromo specified.')
					sys.exit(-1)				
					
			# If size control specified, grab slope/int values from PopVars	
			elif sizeans == 'Y': # Size control

				if (cdevolveans == 'M' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link') and (gen >= burningen_cdevolve):
					tempgenes = SubpopIN[isub][iind]['genes']
					if tempgenes[0] == 2: # AA
						tempvals = fitvals[isub][0] # First spot AA									
					elif tempgenes[0] == 1 and tempgenes[1] == 1: # Aa
						tempvals = fitvals[isub][1] # Second spot Aa
					elif tempgenes[1] == 2: # aa
						tempvals = fitvals[isub][2] # third spot aa
					else:
						print('2 alleles only with M options in cdevolveans.')
						sys.exit(-1)
					
					# Then Replace mat vals	0.1:1.0~0.2:0.5~0.1:1.0~0.2:0.5|0.1:1.0~0.2:0.5~0.1:1.0~0.2:0.5 based on sex split
					# tempvals ['0.1:1.0','0.2:0.5','0.1:1.0','0.2:0.5']								
					if len(tempvals) == 1:
						tempmat = tempvals[0].split(':')
					elif len(tempvals) == sexchromo:
						tempmat = tempvals[sxspot].split(':')
					else:
						print('PopVars fitness values for cdevolveans M or MG must be length 1 or length of number of sex_chromo specified.')
						sys.exit(-1)					
					mat_slope = float(tempmat[0])
					mat_int = float(tempmat[1])
					matval = np.exp(mat_int + mat_slope * SubpopIN[isub][iind]['size']) / (1 + np.exp(mat_int + mat_slope * SubpopIN[isub][iind]['size']))
				else:
					if Indsex == 'FXX': # Female		
						matval = np.exp(float(FXXmat_int) + float(FXXmat_slope) * SubpopIN[isub][iind]['size']) / (1 + np.exp(float(FXXmat_int) + float(FXXmat_slope) * SubpopIN[isub][iind]['size']))
					elif Indsex == 'MXY': # Male			
						matval = np.exp(float(MXYmat_int) + float(MXYmat_slope) * SubpopIN[isub][iind]['size']) / (1 + np.exp(float(MXYmat_int) + float(MXYmat_slope) * SubpopIN[isub][iind]['size']))
					elif Indsex == 'MYY': # Male			
						matval = np.exp(float(MYYmat_int) + float(MYYmat_slope) * SubpopIN[isub][iind]['size']) / (1 + np.exp(float(MYYmat_int) + float(MYYmat_slope) * SubpopIN[isub][iind]['size']))
					elif Indsex == 'FYY': # Male			
						matval = np.exp(float(FYYmat_int) + float(FYYmat_slope) * SubpopIN[isub][iind]['size']) / (1 + np.exp(float(FYYmat_int) + float(FYYmat_slope) * SubpopIN[isub][iind]['size']))				
			# Error check 	
			else:
				print('Size control option not correct, N or Y.')
				sys.exit(-1)
		
		# Check matval
		randmat = np.random.uniform()
		if randmat < matval:
			SubpopIN[isub][iind]['mature'] = 1 # Becomes mature	
			SubpopIN[isub][iind]['newmature'] = 1# Becomes new mature (Note after DoUpdate write out, this gets reset
		else:
			SubpopIN[isub][iind]['mature'] = 0 # Does not mature
			SubpopIN[isub][iind]['newmature'] = 0
			SubpopIN[isub][iind]['layeggs'] = 0			
		
	# Check if mature female, then chance it lays eggs
	if SubpopIN[isub][iind]['mature'] == 1 and (Indsex == 'FXX' or Indsex == 'FYY'):
		tempEggFreq=[] # Temp list value to store egg lay events
		stochastic_update(eggFreq_mu,eggFreq_sd,tempEggFreq)
		tempEggFreq = tempEggFreq[0] # Note indexing into first spot since list created above
		if tempEggFreq < 1: # If egg laying is less than 1 event per year
			randegglay = np.random.uniform()
			if randegglay < tempEggFreq:
				SubpopIN[isub][iind]['layeggs'] = 1 # Lays eggs next year
			else:
				SubpopIN[isub][iind]['layeggs'] = 0 # Does not Lays eggs next year
		else: # egg laying is greater than 1 event per year
			SubpopIN[isub][iind]['layeggs'] = tempEggFreq # Lays eggs next year	
		
	#End::matureInd()

# ---------------------------------------------------------------------------------	
def capInd(lastage,SubpopIN,isub,iind,sizecall,size_mean,ClasscapProb,PopcapProb,sexchromo):
	'''
	Capture individuals
	'''
	
	# Get sex
	Indsex = SubpopIN[isub][iind]['sex']
	if Indsex == 'FXX':
		sxspot = 0
	elif Indsex == 'MXY':
		sxspot = 1
	elif Indsex == 'MYY':
		sxspot = 2
	else:
		sxspot = 3
		
	# Get the age adjusted number for binning and indexing into Capture Age
	if sizecall == 'Y':
		size_mean_middles = np.asarray(size_mean)[1:] - np.diff(np.asarray(size_mean).astype('f'))/2
		age_adjusted = np.searchsorted(size_mean_middles, SubpopIN[isub][iind]['size'])
	else:
		age_adjusted = SubpopIN[isub][iind]['age']
	# If above last age class
	if age_adjusted >= lastage:
		age_adjusted = lastage - 1
	
	# Patch adjusted capture probability
	capval_pop = PopcapProb
	
	# Age adjusted capture probability			
	capval_age = ClasscapProb[age_adjusted]
	if len(capval_age.split('~')) == 1:
		capval_age = capval_age
	elif len(capval_age.split('~')) != sexchromo:
		print('Number of age-specific capture probability parameters must match sex_chromo.')
		sys.exit(-1)
	else:
		capval_age = capval_age.split('~')[sxspot]
	
	# Patch level check first
	if capval_pop != 'N':
		capval_pop = float(capval_pop)
		randcapno = np.random.uniform()
		if randcapno < capval_pop: # Successful patch capture
			# Class level check second
			if capval_age != 'N':					
				capval_age = float(capval_age) # Convert to float
				randcapno = np.random.uniform()
				if randcapno < capval_age: # Successful capture	
					SubpopIN[isub][iind]['capture'] = 1
					SubpopIN[isub][iind]['recapture'] = SubpopIN[isub][iind]['recapture']+1
	
	#End::capInd()
	
# ---------------------------------------------------------------------------------------------------	 
def DoUpdate(packans,SubpopIN,K,xgridpop,ygridpop,gen,nthfile,ithmcrundir,loci,alleles,logfHndl,gridsample,growans = None,cdevolveans = None,fitvals = None,burningen_cdevolve = None,ClasscapProb=None,PopcapProb=None,NCap=None,CapClass=None,sizecall=None,size_mean=None,Nclass=None,eggFreq_mu=None,eggFreq_sd=None,sizevals=None,sizeLoo=None,sizeR0=None,size_1=None,size_2=None,size_3=None,size_4=None,plasticans=None,burningen_plastic=None,timeplastic=None,plastic_signalresp=None,geneswap = None,habvals=None,sexchromo=None,Track_DiseaseStates=None,Track_DiseaseStates_AfterDeaths=None, disease_vars=None,age_mature=None,FXXmat_set=None,FXXmat_int=None,FXXmat_slope=None,MXYmat_set=None,MXYmat_int=None,MXYmat_slope=None,MYYmat_set=None,MYYmat_int=None,MYYmat_slope=None,FYYmat_set=None,FYYmat_int=None,FYYmat_slope=None):
	
	'''
	DoUpdate()
	Update Age, Size and some tracker variables.
	Write out information to file.
	'''	
	
	# --------------------------------------------------
	# Tracking numbers for capturing (Middle and Sample and N)
	if gridsample != 'Initial':
		
		classno = len(size_mean[0][0])
		
		# Patch tracking
		NCap.append([])	
		
		# Add spots for Age tracking
		CapClass.append([]) # If capture 0/1
		CapClass[gen] = [[] for x in range(0,classno)]
				
		# Add spots for Age tracking
		Nclass.append([])
		Nclass[gen] = [[] for x in range(0,classno)]
		
		# Add spots for disease tracking
		Track_DiseaseStates.append([])
		Track_DiseaseStates_AfterDeaths.append([])
		#pdb.set_trace()				
		# ---------------------------------------------------------------------
		# Begin loop through subpopulations updating tasks at appropriate times
		for isub in range(len(K)):
			
			# -------------------------------------------------------
			# Disease Process: Get patch level summary disease states
			# -------------------------------------------------------
			indstates_inthispatch = SubpopIN[isub]['states'] 
			countstates_inthispatch = np.bincount(indstates_inthispatch, minlength=disease_vars['noStates'][isub])					
			#pdb.set_trace()			
			# Begin looping through individuals in subpop
			# -------------------------------------------			
			for iind in range(len(SubpopIN[isub])):
								
				# -----------------------------------------------------
				# Get this individuals original ClassVars file and bins
				# -----------------------------------------------------				
				#start_time1 = datetime.datetime.now()
				natalP = int(SubpopIN[isub][iind]['classfile'].split('_')[0].split('P')[1])
				theseclasspars = int(SubpopIN[isub][iind]['classfile'].split('_')[1].split('CV')[1])
				#print(datetime.datetime.now() -start_time1, "Grab classpars")
				
				# -----------------------------
				# Grow here - middle and sample
				# -----------------------------
				if growans != 'N':
					#Indloc = SubpopIN[isub][iind][sourcePop] # Get location	
					Indloc = str(isub+1)						
					validate(Indloc == 'NA', 'Error in individual location DoUpdate()')
					growInd(Indloc,SubpopIN,sizeLoo,sizeR0,size_1,size_2,size_3,size_4,sizevals,isub,iind,growans,size_mean[natalP][theseclasspars],gridsample,cdevolveans,gen,burningen_cdevolve,alleles,fitvals,sexchromo)
										
				# --------------------------------------------
				# Age here - Middle or Second Update
				# --------------------------------------------
				lastage = classno
				if gridsample == 'Middle':					
					SubpopIN[isub][iind]['age'] = SubpopIN[isub][iind]['age'] + 1
										
				# -------------------------------------------------------
				# Mature, egg lay frequency here - Sample or Third Update
				# -------------------------------------------------------
				if gridsample != 'Middle':
					matureInd(lastage,SubpopIN,isub,iind,sizecall,age_mature[natalP][theseclasspars],eggFreq_mu,eggFreq_sd,cdevolveans,fitvals,burningen_cdevolve,gen,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo)	
				
				# -------------------------------------------------------
				# Check Plastic signal response
				# -------------------------------------------------------
				if (plasticans != 'N') and (((gridsample == 'Middle') and (timeplastic.find('Back') != -1)) or (((gridsample == 'Sample') or (gridsample == 'N')) and (timeplastic.find('Out') != -1))):					
					updatePlasticGenes(SubpopIN[isub][iind],cdevolveans,gen,geneswap,burningen_plastic,sizevals[isub],plasticans,timeplastic,gridsample,habvals[isub],plastic_signalresp)
									
				# ---------------------------------
				# Capture here - Middle and Sample
				# ---------------------------------				
				capInd(lastage,SubpopIN,isub,iind,sizecall,size_mean[natalP][theseclasspars],ClasscapProb[natalP][theseclasspars],PopcapProb[isub],sexchromo)
				
				# ----------------------------------------
				# Update Disease State for this individul
				# ----------------------------------------
				#pdb.set_trace()
				if ((gridsample == 'Middle' and disease_vars['ImpDisease'] in ['Both','Back']) or (gridsample == 'Sample' and disease_vars['ImpDisease'] in ['Both','Out'])):
					moveStates(SubpopIN,isub,iind,countstates_inthispatch,disease_vars,gen)
											
			## End:: iind Loop for this patch-tracking numbers on updates next			
			# -----------------------------------------------------------------
			# For tracking disease states updates before deaths
			# -----------------------------------------------------------------
			Track_DiseaseStates[gen].append([])
			indstates_inthispatch = SubpopIN[isub]['states'] 
			updated_countstates = np.bincount(indstates_inthispatch, minlength=disease_vars['noStates'][isub])
			Track_DiseaseStates[gen][isub].extend(updated_countstates)
			
			# ------------------------------------------------------------------------
			# If mortality comparment in Disease Module, then remove these individuals
			# ------------------------------------------------------------------------
			#pdb.set_trace()
			if ((gridsample == 'Middle' and disease_vars['ImpDisease'] in ['Both','Back']) or (gridsample == 'Sample' and disease_vars['ImpDisease'] in ['Both','Out'])) and (disease_vars['DComp'][isub] != 'N'):
				removeDStates(SubpopIN,isub,disease_vars,gen)
			#pdb.set_trace()
			
			# -----------------------------------------------------------------
			# For tracking disease states updates before after deaths
			# -----------------------------------------------------------------
			Track_DiseaseStates_AfterDeaths[gen].append([])
			indstates_inthispatch = SubpopIN[isub]['states'] 
			updated_countstates = np.bincount(indstates_inthispatch, minlength=disease_vars['noStates'][isub])
			Track_DiseaseStates_AfterDeaths[gen][isub].extend(updated_countstates)
			
			# -----------------------------------------------------------------------------------
			# Update environmental reservoir concentration of contaminant - Indirect Transmission
			# -----------------------------------------------------------------------------------
			if ((gridsample == 'Middle' and disease_vars['ImpDisease'] in ['Both','Back']) or (gridsample == 'Sample' and disease_vars['ImpDisease'] in ['Both','Out'])) and disease_vars['TransMode'][isub] == 'Indirect':
				updateEnvRes(disease_vars,isub,countstates_inthispatch)

			
			# -----------------------------------------------------------------
			# For tracking age/size numbers, use min and max for multiple files
			# -----------------------------------------------------------------
			if sizecall == 'Y' and packans.split('_')[0] != 'logistic':
				# Get the middles for finding closest values
				size_bin = size_mean[0][0]
				size_mean_middles_bin = np.asarray(size_bin)[1:] - np.diff(np.asarray(size_bin).astype('f'))/2
				age_adjusted = np.searchsorted(size_mean_middles_bin, SubpopIN[isub]['size'])
			else:
				# Count up each uniages
				age_adjusted = SubpopIN[isub]['age']
			
			# ---------------------------------------------------
			# Track numbers here for capture - Middle and Sample
			# ---------------------------------------------------			
			# Storage tracker for Capture N total
			NCap[gen].append(sum(SubpopIN[isub]['capture']))			
			# Tracking classes: CaptureProb
			for iage in range(len(CapClass[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				if len(sizeindex) == 0:
					CapClass[gen][iage].append(0)
				else:
					CapClass[gen][iage].append(sum(SubpopIN[isub][sizeindex]['capture']))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			if len(sizeindex) == 0:
				CapClass[gen][iage].append(0)
			else: # Add them to last class
				CapClass[gen][iage].append(sum(SubpopIN[isub][sizeindex]['capture']))
			
			# ---------------------------------------------------
			# Track numbers here for N - after growth, not Initial
			# ---------------------------------------------------			
			# Tracking classes: CaptureProb
			for iage in range(len(Nclass[gen])):
				sizeindex = np.where(age_adjusted==iage)[0]
				if len(sizeindex) == 0:
					Nclass[gen][iage].append(0)
				else:
					Nclass[gen][iage].append(len(SubpopIN[isub][sizeindex]['size']))
			# Special case where age class is greater than lastage
			sizeindex = np.where(age_adjusted > iage)[0]
			if len(sizeindex) == 0:
				Nclass[gen][iage].append(0)
			else: # Add them to last class
				Nclass[gen][iage].append(len(SubpopIN[isub][sizeindex]['size']))
				
		##End::Patch loop
		
		# --------------------------------------------------------
		# Track numbers for N and capture - Middle and Sample
		# --------------------------------------------------------
		# Age tracking
		for iage in range(len(CapClass[gen])):		
			CapClass[gen][iage] = sum(CapClass[gen][iage])
			Nclass[gen][iage] = sum(Nclass[gen][iage])			
			
		# ---------------------------------------------------------
		# Track numbers for disease state totals - Middle and Sample
		# ---------------------------------------------------------
		Track_DiseaseStates[gen].insert(0,np.sum(np.asarray(Track_DiseaseStates[gen]),axis=0).tolist())
		Track_DiseaseStates_AfterDeaths[gen].insert(0,np.sum(np.asarray(Track_DiseaseStates_AfterDeaths[gen]),axis=0).tolist())			

	# ------------------------------------------------------------
	# Write out text file for generations specified by nthfile
	# ------------------------------------------------------------	
	# Check if nthfile == generation or "Initial" or "Middle" or "Sample"
	nthfile = np.asarray(nthfile)
	if gridsample == 'Initial':
		DoOutput(SubpopIN,xgridpop,ygridpop,gen,ithmcrundir,loci,alleles,logfHndl,gridsample)
	elif gridsample == 'Middle' or gridsample == 'Sample':
		# If extinct don't do that
		getyear = np.where(gen == nthfile)[0] # Checks for nthfile output
		checkPopN = [len(SubpopIN[x]) for x in range(0,len(SubpopIN))] # Checks for population extinction
		if len(getyear) != 0 and sum(checkPopN) != 0: # Only write output it not extinct
			DoOutput(SubpopIN,xgridpop,ygridpop,gen,ithmcrundir,loci,alleles,logfHndl,gridsample)
	
	# -----------------------------------------------------
	# Release the captured individuals - Middle and Sample
	# Reset New maturers as well - ALL
	# -----------------------------------------------------	
	# Begin loop through subpopulations
	for isub in range(len(K)):
		if gridsample != 'Initial':
			SubpopIN[isub]['capture'] = 0
		SubpopIN[isub]['newmature'] = 0
	
	# Return variables only if updated age	
	return SubpopIN
	
	# End::DoUpdate()
