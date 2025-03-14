# -------------------------------------------------------------------------------------------------
# CDmetaPOP_Disease.py
# Author: Erin L Landguth
# Created: December 2024
# Description: This is the function/module file for Disease processes.
# --------------------------------------------------------------------------------------------------
	
# Python specific functions
import os, copy, pdb, sys,random
import numpy as np
import pandas as pd
#from CDmetaPOP_Modules import validate # circular import somewhere 

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

# ---------------------------------------------------------------------------
def read_disease_vars_files(filelist, datadir,implementdisease,pathogen_load,disease_fitvals):    
	''' Helper function to read in DiseaseVars file'''
	
	#pdb.set_trace()
	noStates,InitCond,TranRates_Files,SusCompartment,infectCompartment,deathCompartment,offspringTransmissionAns,TranRates,ImplementWhen,StartDiseaseWhen,TransMode,TranRates_copy,RDefense,TDefense = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
	if implementdisease != 'N':
		for isub in range(len(filelist)):
			# Read the CSV file into a pandas DataFrame
			df = pd.read_csv(datadir+filelist[isub], header=None)
			noStates.append(int(df.iloc[1,0]))
			InitCond.append(np.asarray(df.iloc[1,1].split(';'),float).tolist())
			TranRates_Files.append(df.iloc[1,2])
			SusCompartment.append(df.iloc[1,3])
			infectCompartment.append(df.iloc[1,4])
			deathCompartment.append(df.iloc[1,5])
			offspringTransmissionAns.append(df.iloc[1,6])
			StartDiseaseWhen.append(int(df.iloc[1,7]))
			TransMode.append(df.iloc[1,8])
			RDefense.append(df.iloc[1,9])
			TDefense.append(df.iloc[1,10])
			
			# Read in Transition Matrix File
			df1 = pd.read_csv(datadir+TranRates_Files[isub],header=None)
			TranRates.append(np.asarray(df1,dtype=str).tolist())
			TranRates_copy.append(np.asarray(df1,dtype=str).tolist())
			#pdb.set_trace()		
			# Error Checks
			if RDefense[isub] != 'N':
				validate(len(RDefense[isub].split('_')) != 2, 'Resistant disease defense option is in incorrect format. Specify N or state to state.')
			if TDefense[isub] != 'N':
				validate(len(TDefense[isub].split('_')) != 2, 'Tolerant disease defense option is in incorrect format. Specify N or state to state.')
			validate(TransMode[isub].lower() != 'direct' and TransMode[isub].lower() !='indirect', 'Transmission mode can be direct or indirect.')
			validate(noStates[isub] != len(InitCond[isub]), 'Number of Disease States must equal Initial Conditions specified.')
			validate(sum(InitCond[isub]) != 1.0, 'Initial conditions for disease states must sum to 1.')		
			validate((TransMode[isub].lower() == 'direct') and ((len(TranRates[isub]) or len(TranRates[isub][0])) != noStates[isub]), 'Transition matrix must be number of States x number of States if direct transmission is specified.')
			validate((TransMode[isub].lower() == 'indirect') and ((len(TranRates[isub]) or len(TranRates[isub][0])) != noStates[isub]+1), 'Transition matrix must be number of States + 1 x number of States + 1 if indirect transmission is specified.')
			validate(offspringTransmissionAns[isub].lower() != 'susceptible', 'Incorrect option used for offspring transmission rules.')
			if deathCompartment[isub] != 'N':
				validate(noStates[isub] <= int(deathCompartment[isub]), 'Specific compartments (death) must be less than total number of states.')
			validate(noStates[isub] <= int(infectCompartment[isub]), 'Specific compartments (infect) must be less than total number of states.')
			validate(noStates[isub] <= int(SusCompartment[isub]), 'Specific compartments (susceptible) must be less than total number of states.')

	else: # No disease, but initialize with one state
		for isub in range(len(filelist)):
			#noStates = [1 for x in range(len(filelist
			noStates.append(1)
			InitCond.append([])
			TranRates.append([])
			TranRates_copy.append([])
			SusCompartment.append([])
			infectCompartment.append([])
			deathCompartment.append([])
			offspringTransmissionAns.append([])
			StartDiseaseWhen.append([])
			TransMode.append([])
			RDefense.append([])
			TDefense.append([])
			#pathogen_load.append('0')
			#disease_fitvals.append([])
	
	tupDiseaseVars = {'noStates':noStates,'InitCond':InitCond,'TransRates':TranRates,'SusComp':SusCompartment,'InfComp':infectCompartment,'DComp':deathCompartment,'OffAns':offspringTransmissionAns, 'StartDisease': StartDiseaseWhen, 'TransRates_ForStochasticUpdates':TranRates_copy,'TransMode':TransMode, 'ImpDisease':implementdisease, 'PathLoad':pathogen_load, 'ResDefense':RDefense, 'TolDefense':TDefense, 'DefenseVals':disease_fitvals}
	
	return tupDiseaseVars
	#End::read_disease_rate_files()

# ---------------------------------------------------------------------------
def moveStates(SubpopIN,isub,iind,countstates_inthispatch,disease_vars,gen):
	'''
	Update the disease state for each individual
	'''
	#pdb.set_trace() 
	# After burnin
	if gen >= disease_vars['StartDisease'][isub]:				
				
		# This individuals state
		istate = SubpopIN[isub][iind]['states']	
		
		# States that are susceptible
		susceptible_state = int(disease_vars['SusComp'][isub])
		
		# Get the compartment(s) involved in infection - currently only 1 allowed
		infected_state = int(disease_vars['InfComp'][isub]) # index into infected state
							
		# Grab transition rates for this state (columns)
		iStates_transmission_rates = [row[istate] for row in disease_vars['TransRates'][isub]]
		
		# If indirect transmission, then grab the pathogen load and it's transition rates
		if disease_vars['TransMode'][isub] == 'Indirect':
			
			# Get the P state - assume this is the last state for the transition matrix
			p_state = disease_vars['noStates'][isub] 			
			del(iStates_transmission_rates[p_state]) # Remove this transition rate from matrix (it is environment not host)
			
			# Get the pathogen load in this patch
			pathogen_load_atthispatch = disease_vars['PathLoad'][isub]
			
			# Get the transition rates for p state (columns) - These rates are used in the force of infection
			FROMpState_TOiStates = [row[p_state] for row in disease_vars['TransRates'][isub]] 
			del(FROMpState_TOiStates[p_state]) # Remove this p state similar as above
						
		# If disease defense is on, grab its genes
		if disease_vars['ResDefense'][isub].lower() != 'n' or disease_vars['TolDefense'][isub].lower() != 'n':
			igenes = SubpopIN[isub][iind]['genes']
			validate(len(igenes) < 4*2, 'Specify 4 loci for defense defense strategies.')
			# Get the third and fourth locus used for defense
			Clocus = igenes[4:6]		
			Dlocus = igenes[6:8]
			defense_vals = disease_vars['DefenseVals'][isub]
		
		# Index values that this state could move to (non zero entries)
		moveTOIndex = np.where(np.asarray(iStates_transmission_rates) > 0)[0]
		
		# Check for multiple S transitions - need to add this in.
		validate(len(moveTOIndex) > 1 and istate == 0, 'S cannot transition to more than 1 compartment at this time.')
			
		# Then if there are potential states to move to
		if len(moveTOIndex) != 0:					
			
			# Shuffle the moveTOIndex and loop through these possible moves
			random.shuffle(moveTOIndex)
			
			# And loop through possibilities until this individual moves (or not)
			for imove in moveTOIndex:
				# Random number
				williStateMove_randNO = np.random.uniform()
				# This states tranistion rate TO imove
				willthisStateMove_ToThisState = disease_vars['TransRates'][isub][imove][istate]
				
				# Check for disease defense upgrades
				if disease_vars['ResDefense'][isub].lower() != 'n': 
				
					if int(disease_vars['ResDefense'][isub].split('_')[0]) == istate and int(disease_vars['ResDefense'][isub].split('_')[1]) == imove:
						# Get Resistant Clocus new values
						if Clocus[0] == 2:
							willthisStateMove_ToThisState = defense_vals[0]
							
						elif Clocus[0] == 1:
							willthisStateMove_ToThisState = defense_vals[1]
							
						elif Clocus[0] == 0:
							willthisStateMove_ToThisState = defense_vals[2]
				
				# If this is the death compartment and disease defense is tolerant, update death rate
				if disease_vars['TolDefense'][isub].lower() != 'n':
					
					if int(disease_vars['TolDefense'][isub].split('_')[0]) == istate and int(disease_vars['TolDefense'][isub].split('_')[1]) == imove:	
						# Get Tolerant Dlocus new values
						if Dlocus[0] == 2:
							willthisStateMove_ToThisState = defense_vals[3]
							
						elif Dlocus[0] == 1:
							willthisStateMove_ToThisState = defense_vals[4]
							
						elif Dlocus[0] == 0:
							willthisStateMove_ToThisState = defense_vals[5]
				
				# For the Susceptible state - need proportion of Is in the population
				if istate == susceptible_state:
					
					#Calculate force of infection as transmission rate * I / N
					proportion_ofInfected = (countstates_inthispatch[infected_state] / sum(countstates_inthispatch))					
					direct_tran = proportion_ofInfected * willthisStateMove_ToThisState
					indirect_tran = 0
					
					# If indirect transmission, then grab the environmental reservoirs addition force of infection
					if disease_vars['TransMode'][isub] == 'Indirect': 
						# Calculate the additional force of infection
						concentrationP_inEnv = pathogen_load_atthispatch / sum(countstates_inthispatch)
						indirect_tran = concentrationP_inEnv * FROMpState_TOiStates[imove]
					
				# For all other states
				else:					
					
					direct_tran = willthisStateMove_ToThisState 
					indirect_tran = 0 # Note this might be a funciton of P?
					
				# Then assess if it moves					
				if williStateMove_randNO <= (direct_tran + indirect_tran):
					
					# It transitions, update state
					SubpopIN[isub][iind]['states'] = imove
					break # From inner loop

	# End::moveStates()
	
# ---------------------------------------------------------------------------
def updateEnvRes(disease_vars,isub,countstates_inthispatch):
	'''If indirect transmission then update contaminant in environment'''
	
	# Get P variables
	p_state = disease_vars['noStates'][isub] # assume this is the last state
	pathogen_load_atthispatch = disease_vars['PathLoad'][isub]
	TOpState_FROMiStates = disease_vars['TransRates'][isub][p_state] # These rates are used to update the pathogen load
	FROMpState_TOiStates = [row[p_state] for row in disease_vars['TransRates'][isub]] # These rates are used in the force of infection
	
	# Get the compartment(s) involved in infection - currently only 1 allowed
	infected_state = int(disease_vars['InfComp'][isub]) # index into infected state
	#pdb.set_trace()
	transition_rates_pstate = np.where(np.asarray(TOpState_FROMiStates) != 0)[0]
	tplus1_pathogen_load_atthispatch = []
	for iP in transition_rates_pstate:
		
		# if iP < p_state
		if iP < p_state:
			addthesePs = TOpState_FROMiStates[iP] * countstates_inthispatch[iP]
		else:
			addthesePs = TOpState_FROMiStates[iP] * pathogen_load_atthispatch
			
		tplus1_pathogen_load_atthispatch.append(addthesePs)
	
	new_pathogen_load_atthispatch = pathogen_load_atthispatch + sum(tplus1_pathogen_load_atthispatch)
				
	disease_vars['PathLoad'][isub] = new_pathogen_load_atthispatch
			
	# End::updateEnvRes()



# ---------------------------------------------------------------------------
def removeDStates(SubpopIN,isub,disease_vars,gen):
	'''
	Removes the dead individuals from each patch.
	'''
	
	keep_these_inds = np.where(SubpopIN[isub]['states'] != int(disease_vars['DComp'][isub]))[0]
	
	# Append all information to temp SubpopKeep variable
	SubpopIN[isub] = SubpopIN[isub][keep_these_inds]
	
	# End::removeDStates()

# ---------------------------------------------------------------------------
def test_patchlevel_diseasestates_update(gen,startdisease,gridsample,implementdisease,disease_vars,temp_tracking_eqn,SubpopIN_thissub,isub,countstates_inthispatch,noStates,diseaserates_inthispatch):
	'''
	This function is only for testing and getting the patch level state transitions
	In DoUpdate, after first Disease Process call in patch loop add:
	# Testing patch-level disease state updates
	#temp_tracking_eqn = [] -> this goes before patch loop
	#temp_tracking_eqn.append([])			#test_patchlevel_diseasestates_update(gen,startdisease,gridsample,implementdisease,disease_vars,temp_tracking_eqn,SubpopIN[isub],isub,countstates_inthispatch,noStates,diseaserates_inthispatch)
	'''
	
	if gen >= startdisease and ((gridsample == 'Middle' and implementdisease in ['Both','Back']) or (gridsample == 'Sample' and implementdisease in ['Both','Out'])):	
				
		if sum(countstates_inthispatch) == 0: # if population is 0 in this patch
			temp_tracking_eqn[isub].append(countstates_inthispatch)
									
		else:				
			'''
			# S to I: S*I*a/N
			StoI = countstates_inthispatch[0]*countstates_inthispatch[1]*diseaserates_inthispatch[0]/sum(countstates_inthispatch)
			# I to R: I*b - same idea
			ItoR = countstates_inthispatch[1]*diseaserates_inthispatch[1]
			
			# Then update next states
			Stplus1 = countstates_inthispatch[0] - round(StoI)
			Itplus1 = round(StoI)+countstates_inthispatch[1] - round(ItoR)
			Rtplus1 = countstates_inthispatch[2] + round(ItoR)
			temp_tracking_eqn[isub].append([Stplus1,Itplus1,Rtplus1])
			'''
			temp_xij_thatmove_eqn = []					
			# Generalized code of the above lines (for any number of states)
			for irate in range(len(diseaserates_inthispatch)):						
				
				# Get this transition rate
				thisRate = diseaserates_inthispatch[irate]
				
				# Assume first one is the transmission process
				if irate == 0: 
					SitoSj = countstates_inthispatch[irate]*countstates_inthispatch[irate+1]*thisRate/sum(countstates_inthispatch)
				else:
					SitoSj = thisRate*countstates_inthispatch[irate]
				
				temp_xij_thatmove_eqn.append(round(SitoSj))
								
			# Then update new disease state counts
			for istate in range(len(countstates_inthispatch)):						
				
				if istate == 0: # First State is the S*I*a/N
					thisStateChange = temp_xij_thatmove_eqn[istate]
					Stplus1 = countstates_inthispatch[istate] - thisStateChange
				
				elif istate == len(countstates_inthispatch)-1: # Last one					
					Stplus1 = returnthisStateChange + countstates_inthispatch[istate]
				
				else:
					thisStateChange = temp_xij_thatmove_eqn[istate]
					Stplus1 = returnthisStateChange + countstates_inthispatch[istate] - thisStateChange
					
				temp_tracking_eqn[isub].append(Stplus1)
				returnthisStateChange = thisStateChange
			
	else:
		temp_tracking_eqn[isub].append(countstates_inthispatch)
		
	#End::test_patchlevel_diseasestates_update()

# ---------------------------------------------------------------------------
def DoOut_AllTimeDiseasePatch(K_track,N_Init,Track_DiseaseStates_pop,Track_DiseaseStates_SecondUpdate,N_Emigration,N_EmiMortality, Track_DiseaseStates_ThirdUpdate,N_Immigration,N_ImmiMortality,ithmcrundir,Track_DiseaseStates_AddAge0s,Track_DiseaseStates_AddedInds,Track_DiseaseStates_AfterDeaths_SecondUpdate,Track_DiseaseStates_AfterDeaths_ThirdUpdate,Track_DiseaseStates_EnvRes):
	'''Output tracking for Disease Counts'''
	
	# Create time array
	time = np.arange(0,len(K_track),1)
		
	# Get unique number of subpops -1 for total
	nosubpops = len(K_track[0])-1
	
	# Create file to write info to
	outputfile = open(ithmcrundir+'summary_popAllTime_DiseaseStates.csv','w')
	
	# Write out the titles
	# Add Titles from xypoints
	outputtitle = ['Year','N_Initial','States_GetMetrics','EnvResivoir_GetMetrics','States_AfterAnyAddedInds','States_SecondUpdate','States_SecondUpdate_AfterDeathsRemoved','N_AfterEmigration','States_AddedAge0s','N_AfterEmiMort','States_ThirdUpdate','States_ThirdUpdate_AfterDeathsRemoved','N_AfterImmigration','N_AfterImmiMort']
	
	# Write out the title
	for i in range(len(outputtitle)-1):
		outputfile.write(outputtitle[i]+',')
	# To get return character on the end
	outputfile.write(str(outputtitle[len(outputtitle)-1])+'\n')				
		
	# Write to file
	for i in range(len(time)-1):		
		outputfile.write(str(time[i])+',')
		for j in range(nosubpops+1):
			outputfile.write(str(N_Init[i][j])+'|')
		outputfile.write(',')
		
		for j in range(len(Track_DiseaseStates_pop[i])): # Subpop + 1
			for istate in range(len(Track_DiseaseStates_pop[i][j])): # Loop through states
				outputfile.write(str(Track_DiseaseStates_pop[i][j][istate]))
				if istate != len(Track_DiseaseStates_pop[i][j]) - 1:  # Add ';' if not the last element
					outputfile.write(';')
			outputfile.write('|')
		outputfile.write(',')
		
		for j in range(nosubpops+1):
			outputfile.write(str(Track_DiseaseStates_EnvRes[i][j])+'|')
		outputfile.write(',')
		
		for j in range(len(Track_DiseaseStates_AddedInds[i])): # Subpop + 1
			for istate in range(len(Track_DiseaseStates_AddedInds[i][j])): # Loop through states
				outputfile.write(str(Track_DiseaseStates_AddedInds[i][j][istate]))
				if istate != len(Track_DiseaseStates_AddedInds[i][j]) - 1:  # Add ';' if not the last element
					outputfile.write(';')
			outputfile.write('|')
		outputfile.write(',')
		
		for j in range(len(Track_DiseaseStates_SecondUpdate[i])): # Subpop + 1
			for istate in range(len(Track_DiseaseStates_SecondUpdate[i][j])): # Loop through states
				outputfile.write(str(Track_DiseaseStates_SecondUpdate[i][j][istate]))
				if istate != len(Track_DiseaseStates_SecondUpdate[i][j]) - 1:  # Add ';' if not the last element
					outputfile.write(';')
			outputfile.write('|')
		outputfile.write(',')
		
		for j in range(len(Track_DiseaseStates_AfterDeaths_SecondUpdate[i])): # Subpop + 1
			for istate in range(len(Track_DiseaseStates_AfterDeaths_SecondUpdate[i][j])): # Loop through states
				outputfile.write(str(Track_DiseaseStates_AfterDeaths_SecondUpdate[i][j][istate]))
				if istate != len(Track_DiseaseStates_AfterDeaths_SecondUpdate[i][j]) - 1:  # Add ';' if not the last element
					outputfile.write(';')
			outputfile.write('|')
		outputfile.write(',')
		
		for j in range(nosubpops+1):
			outputfile.write(str(N_Emigration[i][j])+'|')
		outputfile.write(',')
		
		for j in range(len(Track_DiseaseStates_AddAge0s[i])): # Subpop + 1
			for istate in range(len(Track_DiseaseStates_AddAge0s[i][j])): # Loop through states
				outputfile.write(str(Track_DiseaseStates_AddAge0s[i][j][istate]))
				if istate != len(Track_DiseaseStates_AddAge0s[i][j]) - 1:  # Add ';' if not the last element
					outputfile.write(';')
			outputfile.write('|')
		outputfile.write(',')
		
		for j in range(nosubpops+1):
			outputfile.write(str(N_EmiMortality[i][j])+'|')
		outputfile.write(',')		
		
		for j in range(len(Track_DiseaseStates_ThirdUpdate[i])): # Subpop + 1
			for istate in range(len(Track_DiseaseStates_ThirdUpdate[i][j])): # Loop through states
				outputfile.write(str(Track_DiseaseStates_ThirdUpdate[i][j][istate]))
				if istate != len(Track_DiseaseStates_ThirdUpdate[i][j]) - 1:  # Add ';' if not the last element
					outputfile.write(';')
			outputfile.write('|')
		outputfile.write(',')
		
		for j in range(len(Track_DiseaseStates_AfterDeaths_ThirdUpdate[i])): # Subpop + 1
			for istate in range(len(Track_DiseaseStates_AfterDeaths_ThirdUpdate[i][j])): # Loop through states
				outputfile.write(str(Track_DiseaseStates_AfterDeaths_ThirdUpdate[i][j][istate]))
				if istate != len(Track_DiseaseStates_AfterDeaths_ThirdUpdate[i][j]) - 1:  # Add ';' if not the last element
					outputfile.write(';')
			outputfile.write('|')
		outputfile.write(',')
		
		for j in range(nosubpops+1):
			outputfile.write(str(N_Immigration[i][j])+'|')
		outputfile.write(',')				
					
		for j in range(nosubpops+1):
			outputfile.write(str(N_ImmiMortality[i][j])+'|')
		outputfile.write('\n')		
		
	# Logging message
	#stringout = 'The file summary_popAllTime_DiseaseStates.csv has been created'
	#logMsg(logfHndl,stringout)	
	
	# Close file
	outputfile.close()
	
	#End::DoOut_AllTimeDiseasePatch()



