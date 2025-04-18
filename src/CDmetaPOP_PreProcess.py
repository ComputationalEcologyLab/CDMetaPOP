# -------------------------------------------------------------------------------------------------
# CDmetaPOP_PreProcess.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file pre processing.
# --------------------------------------------------------------------------------------------------

from scipy.stats import truncnorm
from scipy.stats import norm
import numpy as np 
		
# General python imports
import os,sys,pdb,copy
from ast import literal_eval 
from CDmetaPOP_Modules import *
from inspect import currentframe, getframeinfo

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
	# The case where all of the values in lst are the same
	if len(lst) == count:
		count = count-1
	return item,count
	#End::w_choice_general()		
	
# ----------------------------------------------------------------------------------
def loadFile(filename, header_lines=0, delimiter=None, cdpop_inputvars=False): ###
	'''
	Used to load file hearders according to current UNICOR standards
	as of 5/19/2011
	'''
	try:
		inputfile = open(filename)
	except (IOError,OSError) as e:
		print(("Load file: %s the file (%s) is not available!"%(e,filename)))
		sys.exit(-1)
	header_dict = {}
	data_list = []
	index_list = [] ###
	
	if delimiter != None:
		lines = [ln.rstrip().split(delimiter) for ln in inputfile]
	else:
		lines = [ln.rstrip().split() for ln in inputfile]
	# Close file
	inputfile.close()
		
	for i,line in enumerate(lines):
		if i < header_lines:
			if len(line) <= 1:
				print("Only one header value in line, skipping line...")
				continue
			#else:	
			elif len(line) == 2: ###
				header_dict[line[0]] = line[1]
			### working with a file where the first line is all header keys
			### and the following lines are their data
			elif cdpop_inputvars: ###
				for j in range(len(line)): ###
					header_dict[line[j]] = [] ###
					index_list.append(line[j]) ###
		else:
			#This is a fix to remove empty entries from from a line if they
			#exist, this is to fix a problem with reading in cdmatrices
			
			for i in range(line.count('')):
				line.remove('')
			data_list.append(line)
			if cdpop_inputvars: ###
				#tempTuple = ()
				for j in range(len(line)): ###
					# remove the following lines should the tuple representing bar-delimited values break anything -TJJ
					if line[j].find('|') != -1:
						tempList = line[j].split('|')
						line[j] = tuple(tempList)
					#---
					
					header_dict[index_list[j]].append(line[j]) ###
	
	if not cdpop_inputvars:
		return header_dict, data_list
	else:
		n_jobs = len(lines) - header_lines
		return header_dict, index_list, n_jobs
	#End::loadFile
	
# ---------------------------------------------------------------------------------------------------	 
def GetMaxCDValue(threshold,cdmatrix):	
	'''
	GetMaxCDValue()
	This function calculates the maximum movement thresholds.
	'''	
	
	# movement threshold if max specified
	if str(threshold).endswith('max'):
		# If max
		if len(threshold.strip('max')) == 0:
			threshold = np.amax(cdmatrix)
		else:
			threshold = (float(threshold.strip('max'))/100.)*np.amax(cdmatrix)
	else:
		threshold = float(threshold)
	
	return threshold
	#End::GetMaxCDValue()
	
# ---------------------------------------------------------------------------------------------------	 
def ReadCDMatrix(cdmatrixfilename,function,threshold,A,B,C):
	'''
	ReadMateCDMatrix()
	This function reads in the mating cost distance matrix.
	'''	
	
	# Check statements
	if os.path.exists(cdmatrixfilename):
		# Open file for reading
		inputfile = open(cdmatrixfilename,'r')
	else:
		print(("CDmetaPOP ReadCDMatrix() error: open failed, could not open %s"%(cdmatrixfilename)))
		sys.exit(-1)
	
	# Read lines from the file
	lines = inputfile.readlines()
	
	# Close the file
	inputfile.close()
	
	# Create an empty matrix to append to 
	bigCD = []
	
	# Split up each line in file and append to empty matrix, x
	for spot in lines:
		thisline = spot.strip('\n').split(',')
		bigCD.append(thisline[0:len(lines)])
	bigCD = np.asarray(bigCD,dtype='float')
	
	# Delete lines from earlier
	del(lines)
		
	# Store number of files
	nofiles = len(bigCD)
	
	# Calculate max and min of bigCD matrix
	minbigCD = np.amin(bigCD)
	maxbigCD = np.amax(bigCD)
	
	# Error checks on these values
	if minbigCD < 0:
		print('Cost matrix values should not be negative.')
		sys.exit(-1)
	if maxbigCD != 0:
		if maxbigCD < 1 and (function == '2' or function == '5'):
			print('Careful use of cost distance values less than 1 when using options 2 or 5.')
			sys.exit(-1)
	
	# Get maximum cdvalue to use for movethreshold if specified
	threshold = GetMaxCDValue(threshold,bigCD)
	
	# Create a matrix of to be filled 
	cdmatrix = []
	if function != '9':
		# Fill up matrix with float value of array x
		for j in range(nofiles):
			cdmatrix.append([])
			for k in range(nofiles):
				
				# For the linear function
				if function == '1':
					scale_min = 0.
					scale_max = threshold
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) > threshold:
						cdmatrix[j].append(0.0)
					# If threshold is 0 (philopatry) set to 1 - can't dived by 0
					elif float(bigCD[j][k]) <= threshold and threshold == 0.0:
						cdmatrix[j].append(1.0)
					# Else calculated function value and if not philopatry
					elif float(bigCD[j][k]) <= threshold and threshold != 0.0:
						cdmatrix[j].append(-(1./threshold)*float(bigCD[j][k]) + 1)
					else:
						print('Something off in linear function values.')
						sys.exit(-1)
							
				# For the inverse square function
				elif function == '2':			
					# This function gets rescale: calculate here
					if threshold == 0:
						scale_min = 0.
					else:
						scale_min = 1./(pow(threshold,2))
					scale_max = 1.
					
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) > threshold:
						cdmatrix[j].append(0.0)
					# If threshold is 0 (philopatry) set to 1 - can't dived by 0
					elif float(bigCD[j][k]) <= threshold and threshold == 0.0:
						cdmatrix[j].append(1.0)
					# If cd mat is 0. 
					elif float(bigCD[j][k]) <= threshold and threshold != 0.0 and float(bigCD[j][k]) == 0.0 or (minbigCD == maxbigCD or int(maxbigCD) == 0):
						cdmatrix[j].append(1.0)
					# Else calculated function value
					elif float(bigCD[j][k]) <= threshold and threshold != 0.0 and float(bigCD[j][k]) != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
						invsq_val = 1./(pow(float(bigCD[j][k]),2))
						invsq_val = (invsq_val - scale_min) / (scale_max - scale_min)
						cdmatrix[j].append(invsq_val)# Else something else.
					else:
						print('Something off in inv squ function values.')
						sys.exit(-1)
						
				# Nearest neighbor function here
				elif function == '3':
					print('Nearest neighbor function is not currently implemented.')
					print('You can use Linear function with neighbor threshold for approximation.')
					sys.exit(-1)
				
				# Random function here
				elif function == '4':
					scale_min = 0.
					scale_max = threshold
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) > threshold:
						cdmatrix[j].append(0.0)
					# Else calculated function value
					else:
						cdmatrix[j].append(1.0)
				
				# For the negative binomial function
				elif function == '5':
					
					# This function gets rescale: calculate here
					scale_min = A*pow(10,-B*float(threshold))
					scale_max = A*pow(10,-B*float(minbigCD))
				
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) > threshold:
						cdmatrix[j].append(0.0)
					# Rescaled value divide by zero check cases
					elif float(bigCD[j][k]) <= threshold and threshold == 0.0 and (minbigCD == maxbigCD or int(maxbigCD) == 0):
						cdmatrix[j].append(1.0)
					# Else calculated function value
					elif float(bigCD[j][k]) <= threshold and threshold != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
						negexp = A*pow(10,-B*float(bigCD[j][k]))
						negexp = (negexp - scale_min) / (scale_max - scale_min)
						cdmatrix[j].append(negexp)
					# Else something else.
					else:
						print('Something off in neg exp function values.')
						sys.exit(-1)
						
				# For in a subpopulation only
				elif function == '6':
					scale_min = 0.
					scale_max = 1.
					# Check if within the same subpopulation/patch
					if j == k:					
						cdmatrix[j].append(1.0)
					else:
						cdmatrix[j].append(0.0)
				
				# For Gaussian function 
				elif function == '7':
				
					# This function gets rescale: calculate here
					scale_min = A*np.exp(-((float(threshold)-B)**2)/(2*C**2))
					scale_max = A*np.exp(-((float(minbigCD)-B)**2)/(2*C**2))
				
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) > threshold:
						cdmatrix[j].append(0.0)
					# Rescaled value divide by zero check cases
					elif float(bigCD[j][k]) <= threshold and threshold == 0.0 and (minbigCD == maxbigCD or int(maxbigCD) == 0):
						cdmatrix[j].append(1.0)
					# Else calculated function value
					elif float(bigCD[j][k]) <= threshold and threshold != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
						gauss_val = A*np.exp(-((float(bigCD[j][k])-B)**2)/(2*C**2))
						gauss_val = (gauss_val - scale_min) / (scale_max - scale_min)
						cdmatrix[j].append(gauss_val)
					# Else something else.
					else:
						print('Something off in gauss function values.')
						sys.exit(-1)
						
				# For cost distance matrix only function 
				elif function == '8':
				
					scale_min = minbigCD
					scale_max = threshold
					
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) > threshold:
						cdmatrix[j].append(0.0) 
					# Rescaled value divide by zero check cases - philopatry
					elif (float(bigCD[j][k]) <= threshold and threshold == 0.0) and (minbigCD == maxbigCD or int(maxbigCD) == 0 or threshold == minbigCD):
						cdmatrix[j].append(1.0)
					# If cd mat is 0. 
					elif (float(bigCD[j][k]) <= threshold and threshold != 0.0 and float(bigCD[j][k]) == 0.0) or (minbigCD == maxbigCD or int(maxbigCD) == 0 or threshold == minbigCD):
						cdmatrix[j].append(1.0)
					# Else calculated function value
					elif (float(bigCD[j][k]) <= threshold and threshold != 0.0 and float(bigCD[j][k]) != 0.0) and (minbigCD != maxbigCD and int(maxbigCD) != 0 and threshold != minbigCD):
						cd_val = (float(bigCD[j][k]) - scale_min) / (scale_max - scale_min)
						cdmatrix[j].append(1. - cd_val)
					# Else something else.
					else:
						print('Something off in 8 function values.')
						sys.exit(-1)
					
				# For the Pareto function: (a*(b^a))/(cost distance^(a+1))
				elif function == '10':
					scale_min = None
					scale_max = None
					if B == 0 and minbigCD == 0:
						print('Warning: B = 0 in Pareto function and zeros in cost distance matrix. Divide by zero issue.')
						sys.exit(-1)
					if threshold == 0:
						print('Threshold can not be 0 in Pareto function.')
						sys.exit(-1)
					cdmatval = float(bigCD[j][k]) + B
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) > threshold:
						cdmatrix[j].append(0.0)
					# Else calculated function value
					else:
						pareto = (A*(B**A))/(cdmatval**(A+1))
						cdmatrix[j].append(pareto)
						
				# For FIDIMO equation f(x) = p X 1/sqrt(2*pi*sigstat^2) x exp(-(cost distance -mu)^2/2*sigstat^2) + (1-p) x 1/sqrt(2*pi*sigmob^2) x exp(-(cost distance-mu)^2/2*sigmob^2)
				# p = A; sigstat = B; sigmob = C
				elif function == '11':
					scale_min = None
					scale_max = None
					if B == 0 or C == 0:
						print('Warning: B = 0 or C = 0 in FIDIMO function. Divide by zero issue.')
						sys.exit(-1)
					if threshold == 0:
						print('Threshold can not be 0 in FIDIMO function.')
						sys.exit(-1)
					cdmatval = float(bigCD[j][k])
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) > threshold:
						cdmatrix[j].append(0.0)
					# Else calculated function value
					else:
						fidimo = A * (1./np.sqrt(2*np.pi*B**2)) * np.exp(-((cdmatval - 0.)**2/(2*B**2))) + (1-A) * (1./np.sqrt(2*np.pi*C**2)) * np.exp(-((cdmatval - 0.)**2/(2*C**2)))
						cdmatrix[j].append(fidimo)
					
				# error
				else:
					print('This movement function option does not exist.')
					sys.exit(-1)
	else: # For option 9
		cdmatrix = bigCD
		scale_min = minbigCD
		scale_max = maxbigCD
		
	# Transpose all matrices here (this was because gdistance asymetrical matrices read by column
	cdmatrix = np.transpose(np.asarray(cdmatrix))
	
	# Delete variables
	del(bigCD)
		
	# Return variables
	tupReadMat = cdmatrix,threshold,scale_min,scale_max
	return tupReadMat
	#End::ReadCDMatrix
	
# ---------------------------------------------------------------------------------------------------	
def CreateAlleleList(loci,alleles,xgenes):
	'''
	CreateAlleleList()
	This function creates a list for the allele storage.
	'''		
	
	# Store all information in a list [loci][allele#,probability]
	allelst = []
	for i in range(loci):
		allelst.append([])
		for k in range(alleles[i]):
			allspot = sum(alleles[0:i]) + k + 1
			#allelst[i].append([int(k),float(xgenes[alleles[i]*i+1+k][1])])
			allelst[i].append([int(k),float(xgenes[allspot][1])])
	
	# Return variables
	return allelst
	#End::CreateAlleleList()
	
# ---------------------------------------------------------------------------------------------------	 
def InitializeGenes(datadir,allefreqfilename,loci,alleles):	
	
	'''
	InitializeGenes()	
	'''
	
	allelst = []
	# Loop through allelefrequency files
	for ifile in range(len(allefreqfilename)):
		allelst.append([]) # add a spot for this patch 
		
		fileans = allefreqfilename[ifile]
		fileans = fileans.split(';') # Here for a spatial separation use this deliminator check.
				
		# Then loop through each file
		for i_splitpatch in range(len(fileans)):
			# If genetic structure intialized by a file...
			if not (fileans[i_splitpatch] == 'random' or fileans[i_splitpatch] == 'random_var'):
				
				# Check statements
				if os.path.exists(datadir+fileans[i_splitpatch]):
					# Open file for reading
					inputfile = open(datadir+fileans[i_splitpatch],'r')
				else:
					print(("CDmetaPOP InitializeGenes() error: open failed, could not open %s"%(fileans[i_splitpatch])))
					sys.exit(-1)
					
				# Read lines from the file
				lines = inputfile.readlines()
				
				#Close the file
				inputfile.close()
				
				# Create an empty matrix to append to
				xgenes = []
				
				# Split up each line in file and append to empty matrix, x
				for i in lines:
					thisline = i.replace('"','')
					thisline = thisline.strip('\n').strip('\r').strip(' ').split(',')
					xgenes.append(thisline)
				
				# Error check here
				if (len(xgenes)-1) != sum(alleles):
					print('Allele frequency file is not the specified number of loci and alleles as in in input file.')
					sys.exit(-1)
				
				# Delete lines from earlier
				del(lines)
				
				# Call CreateAlleleList()
				allelst[ifile].append(CreateAlleleList(loci,alleles,xgenes))
			
			# If genetic structure is to be initialize by random
			elif fileans[i_splitpatch] == 'random' or fileans[i_splitpatch] == 'random_var':
				
				# Create even distribution
				xgenes = []
				xgenes.append(['Allele List','Frequency'])
				for iloci in range(loci):
					for iall in range(alleles[iloci]):
						xgenes.append(['L'+str(iloci)+'A'+str(iall),str(1.0/alleles[iloci])])
				
				# Call CreateAlleleList()
				allelst[ifile].append(CreateAlleleList(loci,alleles,xgenes))
	
	# Delete x variable
	del(xgenes)
	
	# Return variables
	return allelst
	#End::InitializeGenes()
	
# ---------------------------------------------------------------------------------------------------	 
def InitializeAge(K,agefilename,datadir):
	'''
	InitializeAge()
	This function initializes the age of each population
	with an age distribution list.
	'''
	
	# Store all information in a list [age,probability]
	agelst = [] # [age,probability] for age distribution
	ageclass = []
	ageno = []
	sexratio = []
	age_percmort_out = []
	age_percmort_back = []
	size_percmort_out = []
	size_percmort_back = []
	age_percmort_out_sd = []
	age_percmort_back_sd = []
	size_percmort_out_sd = []
	size_percmort_back_sd = []
	age_MgOUT = []
	age_MgBACK = []
	age_S = []
	age_mu = []
	age_size_mean = []
	age_size_std = []
	age_mature = []
	age_sigma = []
	age_cap_out = []
	age_cap_back = []	
	lencheck = [] # temp error check
	f_leslie = []
	f_leslie_std = []
	age_DispProb = []
	
	# Then loop through each file
	for isub in range(len(agefilename)):
		
		# Add spots for this patch
		agelst.append([])
		ageclass.append([])
		ageno.append([])
		sexratio.append([])
		age_percmort_out.append([])
		age_percmort_back.append([])
		size_percmort_out.append([])
		size_percmort_back.append([])
		age_percmort_out_sd.append([])
		age_percmort_back_sd.append([])
		size_percmort_out_sd.append([])
		size_percmort_back_sd.append([])
		age_MgOUT.append([])
		age_MgBACK.append([])
		age_S.append([])
		age_mu.append([])
		age_size_mean.append([])
		age_size_std.append([])
		age_mature.append([])
		age_sigma.append([])
		age_cap_out.append([])
		age_cap_back.append([])
		lencheck.append([])
		f_leslie.append([])
		f_leslie_std.append([])
		age_DispProb.append([])
		
		fileans = agefilename[isub]
		fileans = fileans.split(';') # Patch split by 'species' - or hindex applications
		
		# The loop through each file
		for i_splitpatch in range(len(fileans)):
		
			# Add spot for this file
			agelst[isub].append([]) # [age,probability] for age distribution
			ageclass[isub].append([])
			ageno[isub].append([])
			sexratio[isub].append([])
			age_percmort_out[isub].append([])
			age_percmort_back[isub].append([])
			size_percmort_out[isub].append([])
			size_percmort_back[isub].append([])
			age_percmort_out_sd[isub].append([])
			age_percmort_back_sd[isub].append([])
			size_percmort_out_sd[isub].append([])
			size_percmort_back_sd[isub].append([])
			age_MgOUT[isub].append([])
			age_MgBACK[isub].append([])
			age_S[isub].append([])
			age_mu[isub].append([])
			age_size_mean[isub].append([])
			age_size_std[isub].append([])
			age_mature[isub].append([])
			age_sigma[isub].append([])
			age_cap_out[isub].append([])
			age_cap_back[isub].append([])
			lencheck[isub].append([])
			f_leslie[isub].append([])
			f_leslie_std[isub].append([])
			age_DispProb[isub].append([])
			
			# Check statements
			if os.path.exists(datadir+fileans[i_splitpatch]):
				# Open file for reading
				inputfile = open(datadir+fileans[i_splitpatch],'r')
			else:
				print(("CDmetaPOP InitializeAge() error: open failed, could not open %s"%(datadir+fileans[i_splitpatch])))
				sys.exit(-1)
			
			# Read lines from the file
			lines = inputfile.readlines()
			
			#Close the file
			inputfile.close()
			
			# Create an empty matrix to append to
			xage = []

			# Split up each line in file and append to empty matrix, x
			for i in lines:
				thisline = i.strip('\n').strip('\r').strip(' ').split('\r')
				for j in range(len(thisline)):
					xage.append(thisline[j].split(','))
			lencheck[isub][i_splitpatch].append(len(xage)-1)
			# Error on length of file
			if len(xage[0]) != 24:
				print('ClassVars file is not in the correct format (24 fields needed)')
				sys.exit(-1)
			# Delete lines from earlier
			del(lines)
			
			# Loop through file appending values
			for i in range(len(xage)-1):	
				ageclass[isub][i_splitpatch].append(int(xage[i+1][0]))
				age_size_mean[isub][i_splitpatch].append(float(xage[i+1][1]))
				age_size_std[isub][i_splitpatch].append(float(xage[i+1][2]))
				ageno[isub][i_splitpatch].append(float(xage[i+1][3]))
				sexratio[isub][i_splitpatch].append(xage[i+1][4])
				age_percmort_out[isub][i_splitpatch].append(xage[i+1][5])
				age_percmort_out_sd[isub][i_splitpatch].append(xage[i+1][6])
				age_percmort_back[isub][i_splitpatch].append(xage[i+1][7])
				age_percmort_back_sd[isub][i_splitpatch].append(xage[i+1][8])
				size_percmort_out[isub][i_splitpatch].append(xage[i+1][9])
				size_percmort_out_sd[isub][i_splitpatch].append(xage[i+1][10])
				size_percmort_back[isub][i_splitpatch].append(xage[i+1][11])				
				size_percmort_back_sd[isub][i_splitpatch].append(xage[i+1][12])
				age_MgOUT[isub][i_splitpatch].append(xage[i+1][13])
				age_MgBACK[isub][i_splitpatch].append(xage[i+1][14])
				age_S[isub][i_splitpatch].append(xage[i+1][15])
				age_DispProb[isub][i_splitpatch].append(xage[i+1][16])
				age_mature[isub][i_splitpatch].append(xage[i+1][17])
				age_mu[isub][i_splitpatch].append(xage[i+1][18])
				age_sigma[isub][i_splitpatch].append(xage[i+1][19])
				f_leslie[isub][i_splitpatch].append(xage[i+1][20])
				f_leslie_std[isub][i_splitpatch].append(xage[i+1][21])
				age_cap_out[isub][i_splitpatch].append(xage[i+1][22])
				age_cap_back[isub][i_splitpatch].append(xage[i+1][23])	
			
			# Get age distribution list
			for i in range(len(ageno[isub][i_splitpatch])):
				agelst[isub][i_splitpatch].append([ageclass[isub][i_splitpatch][i],ageno[isub][i_splitpatch][i]/sum(ageno[isub][i_splitpatch])])
			
			# Error check on sexratio
			if 'WrightFisher' in sexratio[isub][i_splitpatch]:
				if sexratio[isub][i_splitpatch][1:] != sexratio[isub][i_splitpatch][-1:]:
					print('Wright Fisher specified in ClassVars.csv file. All age classes must be WrightFisher.')
					sys.exit(-1)
			
			# Deletes
			del(xage)
	
	# Error check, all patches must have the same length of classes
	check = sum(sum(sum(lencheck,[]),[]))
	check = np.mod(check,len(agelst[0][0]))
	if check != 0:
		print('ClassVars all must have the same number of classes.')
		sys.exit(-1)
	del(lencheck)	
	
	# Return variables
	tupAgeFile = agelst,age_size_mean,age_size_std,sexratio,age_percmort_out,age_percmort_out_sd,age_percmort_back,age_percmort_back_sd,size_percmort_out,size_percmort_out_sd,size_percmort_back,size_percmort_back_sd,age_MgOUT,age_MgBACK,age_S,age_DispProb,age_mature,age_mu,age_sigma,f_leslie,f_leslie_std,age_cap_out,age_cap_back	
		
	return tupAgeFile
	#End::InitializeAge()

# ---------------------------------------------------------------------------------------------------	 
def InitializeID(K,N):
	'''
	InitializeID()
	This function initializes the location of each individuals for the id varialbe
	{Initial,Residor,Immigrant,Emigrant,Stayor}_{Year born}_{Natal Pop}_{Numeric ID}_{species ID} and produces the subpop temp list with corresponding PopTag ID and Species ID
	'''
	id = []
	subpop = []
	speciesID = []
	for isub in range(len(K)):
		if sum(np.asarray(np.asarray(K,str)[isub].split(';'),int)) > 0: # split checks for AddIndividual call; tempN0 is read in
		#if int(K[isub]) > 0:
			Nspecies = N[isub].split(';')
			for ispec in range(len(Nspecies)): # Loop through each species group
				for iind in range(int(Nspecies[ispec])): # Loop through each ind
					'''
					# See if spot fills based on Nvals
					probfill = float(N[isub]/float(K[isub]))
					randno = np.random.uniform()
					if randno <= probfill:
						# Get name
						name = 'R'+str(isub+1)+'_F'+str(isub+1)+'_m-1f-1'+'_P'+str(isub+1)+'_Y-1_U'+str(iind)
						id.append(name)
						subpop.append(isub+1)
					'''
					
					# Get name
					name = 'R'+str(isub+1)+'_F'+str(isub+1)+'_m-1f-1'+'_P'+str(isub+1)+'_Y-1_U'+str(iind)
					id.append(name)
					subpop.append(isub+1)
					speciesID.append(ispec+1)
	id = np.asarray(id)
	subpop = np.asarray(subpop,dtype = '|U6')
	speciesID = np.asarray(speciesID)
	
	return id,subpop,speciesID
	#End::InitializeID()

# ---------------------------------------------------------------------------------------------------	 
def InitializeVars(sexratio,agelst,cdinfect,loci,alleles,allelst,age_size_mean,age_size_std,subpop,age_mature,sizeans,cdevolveans,fitvals,burningen_cdevolve,addans,sexans,speciesID,N0,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,eggFreq_mu,eggFreq_sd):

	'''
	InitializeVars()
	This function initializes the age,sex,infection,genes of each individual based for the id variable
	'''
	age = []
	sex = []
	size = []
	infection = []
	genes = []
	mature = []
	capture = []
	recapture = []
	layEggs = []
	hindex = []
	whichClassFile = []
	#pdb.set_trace()		
	# Just loop through actual individuals, else this can take a long while - carful of indexing
	for iind in range(len(subpop)):
		
		# ---------------
		# Get patch number (minus 1 for indexing)
		# ----------------
		isub = int(subpop[iind]) - 1
		
		# Error checks on ; splits for N0, classvars file and allele frequency file
		if len(N0[isub].split(';')) != len(allelst[isub]):
			print('gene files split by species with ; length must equal gene init file split by ;')
			sys.exit(-1)
		if len(N0[isub].split(';')) != len(agelst[isub]):
			print('N0 split by species with ; length must equal classvars file split by ;')
			sys.exit(-1)
		
		# --------------------------
		# Get genes - For each loci:
		# --------------------------		
		genes.append([]) # And store genes information
		# First check to see if there is more than one file that can be used for this patch and then randomly choose which one to initialize this individual; make sure this file is stored for ClassVars selection later
		#thisgenefile = np.random.randint(len(allelst[isub]))
		thisgenefile = speciesID[iind]-1 # Match the speciesID (subtract 1 to start from 0) to the gene file order		
		for j in range(loci):							
			# Take a random draw from the w_choice function at jth locus
			rand1 = w_choice_general(allelst[isub][thisgenefile][j])[0]
			rand2 = w_choice_general(allelst[isub][thisgenefile][j])[0]			
			# Append assignment onto indall array - run through each condition for assignment of 1s or 2s or 0s
			# 	1s = heterozygous at that locus
			#	2s = homozygous at that locus
			#	0s = absence of allele
			for k in range(alleles[j]):									
				# Assignment of 2, the rest 0
				if rand1 == rand2: 
					if k < rand1 or k > rand1:
						tempindall = 0
					elif k == rand1:
						tempindall = 2					
				# Assignment of 1s, the rest 0
				if rand1 != rand2:
					if k < min(rand1,rand2) or k > max(rand1,rand2):
						tempindall = 0
					elif k == rand1 or k == rand2:
						tempindall = 1
					else:
						tempindall = 0							
				# Add to genes list
				genes[iind].append(tempindall)
		#pdb.set_trace()
		# ---------------------------------------------
		# Get AA / aa p value for genetag
		# ---------------------------------------------
		if genes[iind][0] == 2: # 1/0 in gene init file
			hindex.append(1.0)
		elif genes[iind][1] == 2: # 0/1
			hindex.append(0.0)
		elif genes[iind][0] == 1 and genes[iind][1] == 1:# 0.5/0.5 gene init file
			hindex.append(0.5)		
		else:
			hindex.append(-9999)		
		
		# -------------------------
		# Store classvars file used - passed in speciesID splits from N0, so match N0 to thisfile here
		# -------------------------
		# First make sure to grab the correct ClassVars file - Then check cases
		if len(allelst[isub]) == len(agelst[isub]):
			thisfile = thisgenefile # match files
		elif len(allelst[isub]) > len(agelst[isub]):
			if len(agelst[isub]) == 1:
				thisfile = 0 # only 1 class file use it
			else:
				thisfile = np.random.randint(len(agelst[isub])) # not ideal, but eg 3 genes and 2 classvars, randomly pick a class vars
		else:
			thisfile = np.random.randint(len(agelst[isub])) # not ideal, but eg 2 genes and 3 classvars, randomly pick a classvars
		# Append the original patch and class vars file location for indexing to later.
		whichClassFile.append('P'+str(isub)+'_CV'+str(thisfile))
		# ---------------
		# Get age
		# ---------------		
		agetemp = w_choice_general(agelst[isub][thisfile])[0]
		age.append(agetemp)
		
		# ---------------------------------
		# Set Size here
		# ---------------------------------
		# Set the size - add in Todd & Ng method
		mu,sigma = age_size_mean[isub][thisfile][agetemp],age_size_std[isub][thisfile][agetemp]			
		# Case here for sigma == 0
		if sigma != 0:
			#lower, upper = 0,np.inf
			#sizesamp  = truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)			
			sizesamp = np.random.normal(mu,sigma)
			if sizesamp < 0:
				sizesamp = 0.
			
		else:
			sizesamp = mu
		sizesamp = round(sizesamp,3)
		size.append(sizesamp)
		
		# ------------------------------
		# Sex set here
		# ------------------------------		
		# Case for Wright Fisher or not
		if sexratio[isub][thisfile][agetemp] != 'WrightFisher':
			# Get each sex prob, add index, convert to prob
			tempratios = sexratio[isub][thisfile][agetemp].split('~')
			# Error check, must add to 100 and must be 3 values
			if len(tempratios) != sexchromo:
				print('Must specify percentages for sex ratio in ClassVars equal to field sex_chromo in PopVars.')
				sys.exit(-1)
			if sum(np.asarray(tempratios,dtype=float)) != 1.:
				print('Sex ratio classes must sum to 1 in ClassVars file.')
				sys.exit(-1)
			ratiolst = []
			tempstrlst = ['FXX','MXY','MYY','FYY']
			for iratiolst in range(len(tempratios)):
				ratiolst.append([tempstrlst[iratiolst],float(tempratios[iratiolst])])
			#pdb.set_trace()
			offsex = w_choice_general(ratiolst)[0] 										
			sex.append(offsex)
		# Special case for WrightFisher
		else: 
			# Error check, WrightFisher can not include YYs
			if float(sexratio[isub][thisfile][agetemp].split(';')[2]) != 0:
				print('Wright Fisher option specified for sex ratios. YY individuals can not be considered.')
				sys.exit(-1)			
			offsex = int(2*np.random.uniform())
			if offsex == 0:
				offsex == 'FXX'
			else:
				offsex == 'MXY'
			sex.append(offsex) # temporary fill
		# Indexing used later for sex splits
		if offsex == 'FXX':
			sxspot = 0
		elif offsex == 'MXY':
			sxspot = 1
		elif offsex == 'MYY':
			sxspot = 2
		else:
			sxspot = 3
		# ------------------------
		# Set the infection status
		# ------------------------
		if cdinfect == 'Y':
			# Append a random number 0 or 1
			infection.append(int(2*np.random.uniform()))
		else:
			infection.append(0)	
			
		# -------------------
		# Capture probability
		# -------------------
		if addans == 'N':
			capture.append(0)
			recapture.append(0)
		else: # For adding individuals 
			capture.append(1)
			recapture.append(0)
		
		# ---------------------------------------------
		# Set maturity Y or N and get egg lay last year
		# ---------------------------------------------		
		#pdb.set_trace()
		matval = 0.0 # Initialize
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
			if sizeans == 'N': # Age control
				
				if len(age_mature[isub][thisfile][agetemp].split('~')) == 1:
					matval = float(age_mature[isub][thisfile][agetemp].split('~')[0])
				elif len(age_mature[isub][thisfile][agetemp].split('~')) != sexchromo:
					print('ClassVars age maturation probabilities must be length 1 or length of number of sex_chromo specified.')
					sys.exit(-1)
				else:
					matval = float(age_mature[isub][thisfile][agetemp].split('~')[sxspot])
					
			# If size control specified, grab slope/int values from PopVars	
			elif sizeans == 'Y': # Size control

				if (cdevolveans == 'M' or cdevolveans == 'MG_ind' or cdevolveans == 'MG_link') and burningen_cdevolve <= 0: # cdevolve answer mature
					tempgenes = genes[iind]
					if tempgenes[0] == 2: # AA
						tempvals = fitvals[isub][0] # First spot AA
						# In case cdclimate is on, grab first split
						tempvals = tempvals.split('|')[0]
						# Then split by sex ~
						tempvals = tempvals.split('~')					
					elif tempgenes[0] == 1 and tempgenes[1] == 1: # Aa
						tempvals = fitvals[isub][1] # Second spot Aa
						# In case cdclimate is on, grab first split
						tempvals = tempvals.split('|')[0]
						# Then split by sex ~
						tempvals = tempvals.split('~')
					elif tempgenes[1] == 2: # aa
						tempvals = fitvals[isub][2] # third spot aa
						# In case cdclimate is on, grab first split
						tempvals = tempvals.split('|')[0]
						# Then split by sex ~
						tempvals = tempvals.split('~')
					else:
						print('2 alleles only with M options in cdevolveans.')
						sys.exit(-1)
					
					# Then Replace mat vals	0.1:1.0~0.2:0.5~0.1:1.0~0.2:0.5|0.1:1.0~0.2:0.5~0.1:1.0~0.2:0.5
					# tempvals ['0.1:1.0','0.2:0.5','0.1:1.0','0.2:0.5']					
					if len(tempvals) == 1:
						tempmat = tempvals[0].split(':')
					elif len(tempvals) != sexchromo:
						print('PopVars fitness values for cdevolveans M or MG must be length 1 or length of number of sex_chromo specified.')
						sys.exit(-1)							
					else:
						tempmat = tempvals[sxspot].split(':')
					if offsex == 'FXX':
						FXXmat_slope = float(tempmat[0])
						FXXmat_int = float(tempmat[1])						
					elif offsex == 'MXY':
						MXYmat_slope = float(tempmat[0])
						MXYmat_int = float(tempmat[1])					
					elif offsex == 'MYY':
						MYYmat_slope = float(tempmat[0])
						MYYmat_int = float(tempmat[1])
					elif offsex == 'FYY':
						FYYmat_slope = float(tempmat[0])
						FYYmat_int = float(tempmat[1])				
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
		
		# Check probability mature and egg laying
		randmat = np.random.uniform()		
		if randmat < matval:
			mature.append(1)
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
							layEggs.append(1)
						else:
							layEggs.append(0)
					else: # egg laying is greater than 1 event per year
						layEggs.append(np.round(tempEggFreq,0))
				else:
					layEggs.append(0)
			else:
				if tempEggFreq < 1: # If egg laying is less than 1 event per year
					randegglay = np.random.uniform()
					if randegglay < eggFreq:
						layEggs.append(1)
					else:
						layEggs.append(0)
				else: # egg laying is greater than 1 event per year
					layEggs.append(np.round(tempEggFreq,0))
		else:
			mature.append(0)
			layEggs.append(0)
			
	# Return Vars
	return age,sex,size,infection,genes,mature,capture,layEggs,recapture,hindex,whichClassFile
	#End::InitializeVars()
	
# ---------------------------------------------------------------------------------------------------	 
def ReadXY(xyfilename):
	'''
	ReadMateXYCDMatrix()
	This function reads in the xy values for the cost distance matrix.
	'''		
	
	# Check statements
	if os.path.exists(xyfilename):
		# Open file for reading
		inputfile = open(xyfilename,'r')
	else:
		print(("CDmetaPOP ReadXY() error: open failed, could not open %s"%(xyfilename)))
		sys.exit(-1)
	
	# Read lines from the file
	lines = inputfile.readlines()
	
	#Close the file
	inputfile.close()
	
	# Create an empty matrix to append to
	xy = []
	
	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.rstrip('\n').rstrip('\r').split(',')
		xy.append(thisline)
		
	# Delete lines from earlier
	del(lines)
	
	# Return variables
	return xy
	#End::ReadXY()

# ---------------------------------------------------------------------------------------------------	 
def DoCDClimate(datadir,icdtime,cdclimgentime,matecdmatfile,dispOutcdmatfile,dispBackcdmatfile,straycdmatfile,matemoveno,dispmoveOutno,dispmoveBackno,StrBackno,matemovethresh,dispmoveOutthresh,dispmoveBackthresh,StrBackthresh,matemoveparA,matemoveparB,matemoveparC,dispmoveOutparA,dispmoveOutparB,dispmoveOutparC,dispmoveBackparA,dispmoveBackparB,dispmoveBackparC,StrBackparA,StrBackparB,StrBackparC,MgOut_patch_prob,Str_patch_prob,K,outsizevals,backsizevals,outgrowdays,backgrowdays,fitvals,popmort_back,popmort_out,eggmort,Kstd,popmort_back_sd,popmort_out_sd,eggmort_sd,outsizevals_sd,backsizevals_sd,outgrowdays_sd,backgrowdays_sd,pop_capture_back,pop_capture_out,cdevolveans,N0_pass,allefreqfiles_pass,classvarsfiles_pass,assortmateModel_pass,assortmateC_pass,subpopmort_pass,PopTag,dispLocalcdmatfile,dispLocalno,dispLocalparA,dispLocalparB,dispLocalparC,dispLocalthresh,comp_coef_pass,betaFile_selection,xvars_betas_pass,outhabvals_pass,backhabvals_pass,plastic_signalresp_pass,plastic_behaviorresp_pass,plasticans,muterate_pass,sexchromo,MgBack_patch_prob,Disperse_patch_prob):
	'''
	DoCDCliamte()
	Reads in cost distance matrices and converts to probabilities.
	'''
	
	# -------------------------------
	# Extract cdclimate values here
	# -------------------------------
	# Store cdmat file information - header file (loadFile()) passes tuple or string if only 1
	if not isinstance(cdclimgentime, (list,tuple)):
		# Error checks
		if cdclimgentime != ['0']:
			print('If not using CDClimate option, set begin time loop with cdclimgentime at 0.')
			sys.exit(-1)
	
	#---------------------------------------
	# Population-based parameters - PopVars
	# --------------------------------------
	
	# ------------------------------------------------
	# Mating ----------------------------------
	if isinstance(matecdmatfile, (list,tuple)):
		matecdmatfile = datadir+matecdmatfile[icdtime]
	else:	
		matecdmatfile = datadir+matecdmatfile
	# Function numbers next ----------------
	if isinstance(matemoveno, (list,tuple)):
		matemoveno = matemoveno[icdtime]
	# A, B, C pars next --------------------
	if isinstance(matemoveparA, (list,tuple)):
		matemoveparA = float(matemoveparA[icdtime])
	else:
		matemoveparA = float(matemoveparA)
	if isinstance(matemoveparB, (list,tuple)):
		matemoveparB = float(matemoveparB[icdtime])
	else:
		matemoveparB = float(matemoveparB)
	if isinstance(matemoveparC, (list,tuple)):
		matemoveparC = float(matemoveparC[icdtime])
	else:
		matemoveparC = float(matemoveparC)
	# Threshold ----------------------------------------
	if isinstance(matemovethresh, (list,tuple)):
		matemovethreshpass = matemovethresh[icdtime]
	else:
		matemovethreshpass = matemovethresh
	# ---------------------------------------
	# Read in cdmatrix.csv - For Mating
	# ---------------------------------------	
	tupReadMat = ReadCDMatrix(matecdmatfile,matemoveno,\
	matemovethreshpass,matemoveparA,matemoveparB,matemoveparC)
	matecdmatrix = np.asarray(tupReadMat[0])
	matemovethresh = tupReadMat[1]
	mate_ScaleMin = tupReadMat[2]
	mate_ScaleMax = tupReadMat[3]
	
	# ------------------------------------------------
	# Dispersal out ----------------------------------
	if isinstance(dispOutcdmatfile, (list,tuple)): # tuple and cdclimate call
		tupVal = sexsplit(dispOutcdmatfile[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispOutcdmatfile,sexchromo) # single string
	FXXdispOutcdmatfile = datadir+tupVal[0]
	MXYdispOutcdmatfile = datadir+tupVal[1]
	MYYdispOutcdmatfile = datadir+tupVal[2]
	FYYdispOutcdmatfile = datadir+tupVal[3]
	# Function numbers next ----------------
	if isinstance(dispmoveOutno, (list,tuple)):
		tupVal = sexsplit(dispmoveOutno[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveOutno,sexchromo) # single string
	FXXdispmoveOutno = tupVal[0]
	MXYdispmoveOutno = tupVal[1]
	MYYdispmoveOutno = tupVal[2]
	FYYdispmoveOutno = tupVal[3]	
	# A, B, C next ---------------------------------------
	if isinstance(dispmoveOutparA, (list,tuple)):
		tupVal = sexsplit(dispmoveOutparA[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveOutparA,sexchromo) # single string
	FXXdispmoveOutparA = tupVal[0]
	MXYdispmoveOutparA = tupVal[1]
	MYYdispmoveOutparA = tupVal[2]
	FYYdispmoveOutparA = tupVal[3]
	if isinstance(dispmoveOutparB, (list,tuple)):
		tupVal = sexsplit(dispmoveOutparB[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveOutparB,sexchromo) # single string
	FXXdispmoveOutparB = tupVal[0]
	MXYdispmoveOutparB = tupVal[1]
	MYYdispmoveOutparB = tupVal[2]
	FYYdispmoveOutparB = tupVal[3]
	if isinstance(dispmoveOutparC, (list,tuple)):
		tupVal = sexsplit(dispmoveOutparC[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveOutparC,sexchromo) # single string
	FXXdispmoveOutparC = tupVal[0]
	MXYdispmoveOutparC = tupVal[1]
	MYYdispmoveOutparC = tupVal[2]
	FYYdispmoveOutparC = tupVal[3]
	# Threshold -------------------------------------------
	if isinstance(dispmoveOutthresh, (list,tuple)):
		tupVal = sexsplit(dispmoveOutthresh[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveOutthresh,sexchromo) # single string
	FXXdispmoveOutthreshpass = tupVal[0]
	MXYdispmoveOutthreshpass = tupVal[1]
	MYYdispmoveOutthreshpass = tupVal[2]
	FYYdispmoveOutthreshpass = tupVal[3]	
	# --------------------------------------
	# Read in cdmatrix.csv - Dispersal Out
	# --------------------------------------
	if sexchromo == 2 or sexchromo == 3 or sexchromo == 4:
		# For FXX
		tupReadMat = ReadCDMatrix(FXXdispOutcdmatfile,FXXdispmoveOutno,FXXdispmoveOutthreshpass,float(FXXdispmoveOutparA),float(FXXdispmoveOutparB),float(FXXdispmoveOutparC))
		FXXdispOutcdmatrix = np.asarray(tupReadMat[0])
		FXXdispOutthresh = tupReadMat[1]
		FXXdispOut_ScaleMin = tupReadMat[2]
		FXXdispOut_ScaleMax = tupReadMat[3]
		# For MXY
		tupReadMat = ReadCDMatrix(MXYdispOutcdmatfile,MXYdispmoveOutno,MXYdispmoveOutthreshpass,float(MXYdispmoveOutparA),float(MXYdispmoveOutparB),float(MXYdispmoveOutparC))
		MXYdispOutcdmatrix = np.asarray(tupReadMat[0])
		MXYdispOutthresh = tupReadMat[1]
		MXYdispOut_ScaleMin = tupReadMat[2]
		MXYdispOut_ScaleMax = tupReadMat[3]	
		# For MYY
		MYYdispOutcdmatrix = 'N'
		MYYdispOutthresh = 'N'
		MYYdispOut_ScaleMin = 'N'
		MYYdispOut_ScaleMax = 'N'
		# For FYY		
		FYYdispOutcdmatrix = 'N'
		FYYdispOutthresh = 'N'
		FYYdispOut_ScaleMin = 'N'
		FYYdispOut_ScaleMax = 'N'
	if sexchromo == 3 or sexchromo == 4:
		# For MYY
		tupReadMat = ReadCDMatrix(MYYdispOutcdmatfile,MYYdispmoveOutno,MYYdispmoveOutthreshpass,float(MYYdispmoveOutparA),float(MYYdispmoveOutparB),float(MYYdispmoveOutparC))
		MYYdispOutcdmatrix = np.asarray(tupReadMat[0])
		MYYdispOutthresh = tupReadMat[1]
		MYYdispOut_ScaleMin = tupReadMat[2]
		MYYdispOut_ScaleMax = tupReadMat[3]	
	if sexchromo == 4:
		# For FYY
		tupReadMat = ReadCDMatrix(FYYdispOutcdmatfile,FYYdispmoveOutno,FYYdispmoveOutthreshpass,float(FYYdispmoveOutparA),float(FYYdispmoveOutparB),float(FYYdispmoveOutparC))
		FYYdispOutcdmatrix = np.asarray(tupReadMat[0])
		FYYdispOutthresh = tupReadMat[1]
		FYYdispOut_ScaleMin = tupReadMat[2]
		FYYdispOut_ScaleMax = tupReadMat[3]	
	
	# -------------------------------------------------
	# Dispersal back ----------------------------------
	if isinstance(dispBackcdmatfile, (list,tuple)): # tuple and cdclimate call
		tupVal = sexsplit(dispBackcdmatfile[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispBackcdmatfile,sexchromo) # single string
	FXXdispBackcdmatfile = datadir+tupVal[0]
	MXYdispBackcdmatfile = datadir+tupVal[1]
	MYYdispBackcdmatfile = datadir+tupVal[2]
	FYYdispBackcdmatfile = datadir+tupVal[3]	
	# Function numbers next ----------------
	if isinstance(dispmoveBackno, (list,tuple)):
		tupVal = sexsplit(dispmoveBackno[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveBackno,sexchromo) # single string
	FXXdispmoveBackno = tupVal[0]
	MXYdispmoveBackno = tupVal[1]
	MYYdispmoveBackno = tupVal[2]
	FYYdispmoveBackno = tupVal[3]
	# A, B, C next ---------------------------------------
	if isinstance(dispmoveBackparA, (list,tuple)):
		tupVal = sexsplit(dispmoveBackparA[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveBackparA,sexchromo) # single string
	FXXdispmoveBackparA = tupVal[0]
	MXYdispmoveBackparA = tupVal[1]
	MYYdispmoveBackparA = tupVal[2]
	FYYdispmoveBackparA = tupVal[3]
	if isinstance(dispmoveBackparB, (list,tuple)):
		tupVal = sexsplit(dispmoveBackparB[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveBackparB,sexchromo) # single string
	FXXdispmoveBackparB = tupVal[0]
	MXYdispmoveBackparB = tupVal[1]
	MYYdispmoveBackparB = tupVal[2]
	FYYdispmoveBackparB = tupVal[3]
	if isinstance(dispmoveBackparC, (list,tuple)):
		tupVal = sexsplit(dispmoveBackparC[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveBackparC,sexchromo) # single string
	FXXdispmoveBackparC = tupVal[0]
	MXYdispmoveBackparC = tupVal[1]
	MYYdispmoveBackparC = tupVal[2]
	FYYdispmoveBackparC = tupVal[3]
	# Threshold -------------------------------------------
	if isinstance(dispmoveBackthresh, (list,tuple)):
		tupVal = sexsplit(dispmoveBackthresh[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispmoveBackthresh,sexchromo) # single string
	FXXdispmoveBackthreshpass = tupVal[0]
	MXYdispmoveBackthreshpass = tupVal[1]
	MYYdispmoveBackthreshpass = tupVal[2]
	FYYdispmoveBackthreshpass = tupVal[3]
	# --------------------------------------
	# Read in cdmatrix.csv - Dispersal Back
	# --------------------------------------
	if sexchromo == 2 or sexchromo == 3 or sexchromo == 4:
		# For FXX
		tupReadMat = ReadCDMatrix(FXXdispBackcdmatfile,FXXdispmoveBackno,FXXdispmoveBackthreshpass,float(FXXdispmoveBackparA),float(FXXdispmoveBackparB),float(FXXdispmoveBackparC))
		FXXdispBackcdmatrix = np.asarray(tupReadMat[0])
		FXXdispBackthresh = tupReadMat[1]
		FXXdispBack_ScaleMin = tupReadMat[2]
		FXXdispBack_ScaleMax = tupReadMat[3]
		# For MXY
		tupReadMat = ReadCDMatrix(MXYdispBackcdmatfile,MXYdispmoveBackno,MXYdispmoveBackthreshpass,float(MXYdispmoveBackparA),float(MXYdispmoveBackparB),float(MXYdispmoveBackparC))
		MXYdispBackcdmatrix = np.asarray(tupReadMat[0])
		MXYdispBackthresh = tupReadMat[1]
		MXYdispBack_ScaleMin = tupReadMat[2]
		MXYdispBack_ScaleMax = tupReadMat[3]	
		# For MYY
		MYYdispBackcdmatrix = 'N'
		MYYdispBackthresh = 'N'
		MYYdispBack_ScaleMin = 'N'
		MYYdispBack_ScaleMax = 'N'
		# For FYY
		FYYdispBackcdmatrix = 'N'
		FYYdispBackthresh = 'N'
		FYYdispBack_ScaleMin = 'N'
		FYYdispBack_ScaleMax = 'N'
	if sexchromo == 3 or sexchromo == 4:
		# For MYY
		tupReadMat = ReadCDMatrix(MYYdispBackcdmatfile,MYYdispmoveBackno,MYYdispmoveBackthreshpass,float(MYYdispmoveBackparA),float(MYYdispmoveBackparB),float(MYYdispmoveBackparC))
		MYYdispBackcdmatrix = np.asarray(tupReadMat[0])
		MYYdispBackthresh = tupReadMat[1]
		MYYdispBack_ScaleMin = tupReadMat[2]
		MYYdispBack_ScaleMax = tupReadMat[3]	
	if sexchromo == 4:
		# For FYY
		tupReadMat = ReadCDMatrix(FYYdispBackcdmatfile,FYYdispmoveBackno,FYYdispmoveBackthreshpass,float(FYYdispmoveBackparA),float(FYYdispmoveBackparB),float(FYYdispmoveBackparC))
		FYYdispBackcdmatrix = np.asarray(tupReadMat[0])
		FYYdispBackthresh = tupReadMat[1]
		FYYdispBack_ScaleMin = tupReadMat[2]
		FYYdispBack_ScaleMax = tupReadMat[3]	
	
	# -------------------------------------------
	# Straying ----------------------------------
	if isinstance(straycdmatfile, (list,tuple)): # tuple and cdclimate call
		tupVal = sexsplit(straycdmatfile[icdtime],sexchromo)
	else:
		tupVal = sexsplit(straycdmatfile,sexchromo) # single string
	FXXStrcdmatfile = datadir+tupVal[0]
	MXYStrcdmatfile = datadir+tupVal[1]
	MYYStrcdmatfile = datadir+tupVal[2]
	FYYStrcdmatfile = datadir+tupVal[3]		
	# Function numbers next ----------------
	if isinstance(StrBackno, (list,tuple)):
		tupVal = sexsplit(StrBackno[icdtime],sexchromo)
	else:
		tupVal = sexsplit(StrBackno,sexchromo) # single string
	FXXStrno = tupVal[0]
	MXYStrno = tupVal[1]
	MYYStrno = tupVal[2]
	FYYStrno = tupVal[3]
	# A, B, C next ---------------------------------------
	if isinstance(StrBackparA, (list,tuple)):
		tupVal = sexsplit(StrBackparA[icdtime],sexchromo)
	else:
		tupVal = sexsplit(StrBackparA,sexchromo) # single string
	FXXStrparA = tupVal[0]
	MXYStrparA = tupVal[1]
	MYYStrparA = tupVal[2]
	FYYStrparA = tupVal[3]
	if isinstance(StrBackparB, (list,tuple)):
		tupVal = sexsplit(StrBackparB[icdtime],sexchromo)
	else:
		tupVal = sexsplit(StrBackparB,sexchromo) # single string
	FXXStrparB = tupVal[0]
	MXYStrparB = tupVal[1]
	MYYStrparB = tupVal[2]
	FYYStrparB = tupVal[3]
	if isinstance(StrBackparC, (list,tuple)):
		tupVal = sexsplit(StrBackparC[icdtime],sexchromo)
	else:
		tupVal = sexsplit(StrBackparC,sexchromo) # single string
	FXXStrparC = tupVal[0]
	MXYStrparC = tupVal[1]
	MYYStrparC = tupVal[2]
	FYYStrparC = tupVal[3]
	# Threshold -------------------------------------------
	if isinstance(StrBackthresh, (list,tuple)):
		tupVal = sexsplit(StrBackthresh[icdtime],sexchromo)
	else:
		tupVal = sexsplit(StrBackthresh,sexchromo) # single string
	FXXStrthreshpass = tupVal[0]
	MXYStrthreshpass = tupVal[1]
	MYYStrthreshpass = tupVal[2]
	FYYStrthreshpass = tupVal[3]	
	# --------------------------------------
	# Read in cdmatrix.csv - Straying
	# --------------------------------------
	if sexchromo == 2 or sexchromo == 3 or sexchromo == 4:
		# For FXX
		tupReadMat = ReadCDMatrix(FXXStrcdmatfile,FXXStrno,FXXStrthreshpass,float(FXXStrparA),float(FXXStrparB),float(FXXStrparC))
		FXXStrcdmatrix = np.asarray(tupReadMat[0])
		FXXStrthresh = tupReadMat[1]
		FXXStr_ScaleMin = tupReadMat[2]
		FXXStr_ScaleMax = tupReadMat[3]
		# For MXY
		tupReadMat = ReadCDMatrix(MXYStrcdmatfile,MXYStrno,MXYStrthreshpass,float(MXYStrparA),float(MXYStrparB),float(MXYStrparC))
		MXYStrcdmatrix = np.asarray(tupReadMat[0])
		MXYStrthresh = tupReadMat[1]
		MXYStr_ScaleMin = tupReadMat[2]
		MXYStr_ScaleMax = tupReadMat[3]	
		# For MYY
		MYYStrcdmatrix = 'N'
		MYYStrthresh = 'N'
		MYYStr_ScaleMin = 'N'
		MYYStr_ScaleMax = 'N'
		# For FYY
		FYYStrcdmatrix = 'N'
		FYYStrthresh = 'N'
		FYYStr_ScaleMin = 'N'
		FYYStr_ScaleMax = 'N'
	if sexchromo == 3 or sexchromo == 4:
		# For MYY
		tupReadMat = ReadCDMatrix(MYYStrcdmatfile,MYYStrno,MYYStrthreshpass,float(MYYStrparA),float(MYYStrparB),float(MYYStrparC))
		MYYStrcdmatrix = np.asarray(tupReadMat[0])
		MYYStrthresh = tupReadMat[1]
		MYYStr_ScaleMin = tupReadMat[2]
		MYYStr_ScaleMax = tupReadMat[3]	
	if sexchromo == 4:
		# For FYY
		tupReadMat = ReadCDMatrix(FYYStrcdmatfile,FYYStrno,FYYStrthreshpass,float(FYYStrparA),float(FYYStrparB),float(FYYStrparC))
		FYYStrcdmatrix = np.asarray(tupReadMat[0])
		FYYStrthresh = tupReadMat[1]
		FYYStr_ScaleMin = tupReadMat[2]
		FYYStr_ScaleMax = tupReadMat[3]	
	
	#---------------------------------------------------	
	# Local Dispersal ----------------------------------
	if isinstance(dispLocalcdmatfile, (list,tuple)): # tuple and cdclimate call
		tupVal = sexsplit(dispLocalcdmatfile[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispLocalcdmatfile,sexchromo) # single string
	FXXdispLocalcdmatfile = datadir+tupVal[0]
	MXYdispLocalcdmatfile = datadir+tupVal[1]
	MYYdispLocalcdmatfile = datadir+tupVal[2]
	FYYdispLocalcdmatfile = datadir+tupVal[3]	
	# Function numbers next ----------------
	if isinstance(dispLocalno, (list,tuple)):
		tupVal = sexsplit(dispLocalno[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispLocalno,sexchromo) # single string
	FXXdispLocalno = tupVal[0]
	MXYdispLocalno = tupVal[1]
	MYYdispLocalno = tupVal[2]
	FYYdispLocalno = tupVal[3]
	# A, B, C next ---------------------------------------
	if isinstance(dispLocalparA, (list,tuple)):
		tupVal = sexsplit(dispLocalparA[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispLocalparA,sexchromo) # single string
	FXXdispLocalparA = tupVal[0]
	MXYdispLocalparA = tupVal[1]
	MYYdispLocalparA = tupVal[2]
	FYYdispLocalparA = tupVal[3]
	if isinstance(dispLocalparB, (list,tuple)):
		tupVal = sexsplit(dispLocalparB[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispLocalparB,sexchromo) # single string
	FXXdispLocalparB = tupVal[0]
	MXYdispLocalparB = tupVal[1]
	MYYdispLocalparB = tupVal[2]
	FYYdispLocalparB = tupVal[3]
	if isinstance(dispLocalparC, (list,tuple)):
		tupVal = sexsplit(dispLocalparC[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispLocalparC,sexchromo) # single string
	FXXdispLocalparC = tupVal[0]
	MXYdispLocalparC = tupVal[1]
	MYYdispLocalparC = tupVal[2]
	FYYdispLocalparC = tupVal[3]
	# Threshold -------------------------------------------
	if isinstance(dispLocalthresh, (list,tuple)):
		tupVal = sexsplit(dispLocalthresh[icdtime],sexchromo)
	else:
		tupVal = sexsplit(dispLocalthresh,sexchromo) # single string
	FXXdispLocalthreshpass = tupVal[0]
	MXYdispLocalthreshpass = tupVal[1]
	MYYdispLocalthreshpass = tupVal[2]
	FYYdispLocalthreshpass = tupVal[3]
	# --------------------------------------
	# Read in cdmatrix.csv - Local Dispersal
	# --------------------------------------
	if sexchromo == 2 or sexchromo == 3 or sexchromo == 4:
		# For FXX
		tupReadMat = ReadCDMatrix(FXXdispLocalcdmatfile,FXXdispLocalno,FXXdispLocalthreshpass,float(FXXdispLocalparA),float(FXXdispLocalparB),float(FXXdispLocalparC))
		FXXdispLocalcdmatrix = np.asarray(tupReadMat[0])
		FXXdispLocalthresh = tupReadMat[1]
		FXXdispLocal_ScaleMin = tupReadMat[2]
		FXXdispLocal_ScaleMax = tupReadMat[3]
		# For MXY
		tupReadMat = ReadCDMatrix(MXYdispLocalcdmatfile,MXYdispLocalno,MXYdispLocalthreshpass,float(MXYdispLocalparA),float(MXYdispLocalparB),float(MXYdispLocalparC))
		MXYdispLocalcdmatrix = np.asarray(tupReadMat[0])
		MXYdispLocalthresh = tupReadMat[1]
		MXYdispLocal_ScaleMin = tupReadMat[2]
		MXYdispLocal_ScaleMax = tupReadMat[3]	
		# For MYY
		MYYdispLocalcdmatrix = 'N'
		MYYdispLocalthresh = 'N'
		MYYdispLocal_ScaleMin = 'N'
		MYYdispLocal_ScaleMax = 'N'
		# For FYY
		FYYdispLocalcdmatrix = 'N'
		FYYdispLocalthresh = 'N'
		FYYdispLocal_ScaleMin = 'N'
		FYYdispLocal_ScaleMax = 'N'
	if sexchromo == 3 or sexchromo == 4:
		# For MYY
		tupReadMat = ReadCDMatrix(MYYdispLocalcdmatfile,MYYdispLocalno,MYYdispLocalthreshpass,float(MYYdispLocalparA),float(MYYdispLocalparB),float(MYYdispLocalparC))
		MYYdispLocalcdmatrix = np.asarray(tupReadMat[0])
		MYYdispLocalthresh = tupReadMat[1]
		MYYdispLocal_ScaleMin = tupReadMat[2]
		MYYdispLocal_ScaleMax = tupReadMat[3]	
	if sexchromo == 4:
		# For FYY
		tupReadMat = ReadCDMatrix(FYYdispLocalcdmatfile,FYYdispLocalno,FYYdispLocalthreshpass,float(FYYdispLocalparA),float(FYYdispLocalparB),float(FYYdispLocalparC))
		FYYdispLocalcdmatrix = np.asarray(tupReadMat[0])
		FYYdispLocalthresh = tupReadMat[1]
		FYYdispLocal_ScaleMin = tupReadMat[2]
		FYYdispLocal_ScaleMax = tupReadMat[3]	
	
	#----------------------------------------------------------
	# Continue PopVars variable
	# ---------------------------------------------------------
	# Assortative mating model pars ---------------------------
	if isinstance(assortmateModel_pass, (list,tuple)):
		assortmateModel = str(assortmateModel_pass[icdtime])
	else:
		assortmateModel = str(assortmateModel_pass)
	if isinstance(assortmateC_pass, (list,tuple)):
		assortmateC = float(assortmateC_pass[icdtime])
	else:
		assortmateC = float(assortmateC_pass)		
	
	# Subpopmort percent matrix file, check and read in 
	if isinstance(subpopmort_pass,(list,tuple)):
		if subpopmort_pass[icdtime] != 'N':
			subpopmort_mat = ReadXY(datadir + subpopmort_pass[icdtime])
			subpopmort_mat = np.asarray(np.asarray(subpopmort_mat)[1:,1:],dtype='float')
			if (len(np.unique(PopTag)) != len(subpopmort_mat[0])):
				print('Subpatch id total does not equal the number of subpopulation mortality percentages given. Also make sure Subpatch id is consecutive.')
				sys.exit(-1)
		else:
			subpopmort_mat = subpopmort_pass[icdtime]
	else:
		if subpopmort_pass != 'N':
			subpopmort_mat = ReadXY(datadir + subpopmort_pass)
			subpopmort_mat = np.asarray(np.asarray(subpopmort_mat)[1:,1:],dtype='float')
			if (len(np.unique(PopTag)) != len(subpopmort_mat[0])):
				print('Subpatch id total does not equal the number of subpopulation mortality percentages given. Also make sure Subpatch id is consecutive.')
				sys.exit(-1)
		else:
			subpopmort_mat = subpopmort_pass
			
	# Plastic Temp/Hab values - plastic_signalres_pass,plastic_behaviorres_pass
	if plasticans != 'N':
		if isinstance(plastic_signalresp_pass, (list,tuple)):
			plastic_signalresp = float(plastic_signalresp_pass[icdtime])
		else:
			plastic_signalresp = float(plastic_signalresp_pass)
		if isinstance(plastic_behaviorresp_pass, (list,tuple)):
			plastic_behaviorresp = float(plastic_behaviorresp_pass[icdtime])
		else:
			plastic_behaviorresp = float(plastic_behaviorresp_pass)
	else:
		plastic_signalresp = plastic_signalresp_pass
		plastic_behaviorresp  = plastic_behaviorresp_pass
			
	# Mutation Rate
	if isinstance(muterate_pass, (list,tuple)):
		muterate = float(muterate_pass[icdtime])
	else:
		muterate = float(muterate_pass)	
			
	# ---------------------------------------------
	# Read in Beta File for Multiple loci selection
	# ---------------------------------------------
	if isinstance(betaFile_selection,(list,tuple)):
		tempbetaFile_selection = datadir+betaFile_selection[icdtime]
	else:
		tempbetaFile_selection = datadir+betaFile_selection
	tempbetas_selection = []	
	if cdevolveans.split('_')[0] == 'P':
		# Read in Beta File
		betavals = ReadXY(tempbetaFile_selection)
	
		# Error check on beta file - should be x number of betas alleles * number of xvars * number of loci under selection
		if (len(betavals)-1)*(len(betavals[0])-1) != int(cdevolveans.split('_')[3].split('A')[1])*int(cdevolveans.split('_')[1].split('X')[1])*int(cdevolveans.split('_')[2].split('L')[1]):
			print('Beta file for selection is incorrect. Specify the beta for each variableXlociXallele combination.')
			sys.exit(-1)
		# Then extract betavals - in order of variable, loci, allele [var][loci][allele]
		for ixvar in range(int(cdevolveans.split('_')[1].split('X')[1])):
			tempbetas_selection.append([])
			for iloci in range(int(cdevolveans.split('_')[2].split('L')[1])):
				tempbetas_selection[ixvar].append([])
				for iall in range(int(cdevolveans.split('_')[3].split('A')[1])):
					tempbetas_selection[ixvar][iloci].append(float(betavals[iall+1][ixvar*(int(cdevolveans.split('_')[2].split('L')[1]))+iloci+1]))
		# Add beta0 - will be the last spot in betas vars 
		if len(betavals[0][0]) == 0:
			tempbetas_selection.append(0.0)
		else:
			tempbetas_selection.append(float(betavals[0][0])) 
	
	# ----------------------
	# Patch based parameters
	# ----------------------
	tempMgBack_patch_prob = []
	tempDisperse_patch_prob = []
	tempStr_patch_prob = []
	tempMgOut_patch_prob = []
	tempoutsize = []
	tempbacksize = []
	tempoutgrow = []
	tempbackgrow = []
	tempoutsize_sd = []
	tempbacksize_sd = []
	tempoutgrow_sd = []
	tempbackgrow_sd = []
	tempfitvals = []
	tempK = []
	tempKstd = []
	temppopmort_back = []
	temppopmort_out = []
	tempeggmort = []
	temppopmort_back_sd = []
	temppopmort_out_sd = []
	tempeggmort_sd = []
	temppopCapOut = []
	temppopCapBack = []
	tempN0 = []
	tempAllelefile = []
	tempClassVarsfile = []
	tempcompcoef = []
	tempxvars_betas = []
	tempouthabvals = []
	tempbackhabvals = []
	for isub in range(len(K)):

		# For addindividual applications N0 split must have matching classvars and genefile splits
		if len(allefreqfiles_pass[isub].split('|')) != len(N0_pass[isub].split('|')) != len(classvarsfiles_pass[isub].split('|')):
			print('N0 split by | for Add Individual applications, must have matching Gene Initialization files and ClassVars files split by | - cdclimate.')
			sys.exit(-1)
		if len(MgBack_patch_prob[isub].split('|')) > 1:
			tempMgBack_patch_prob.append(float(MgBack_patch_prob[isub].split('|')[icdtime]))
		else:
			tempMgBack_patch_prob.append(float(MgBack_patch_prob[isub]))
		if len(Disperse_patch_prob[isub].split('|')) > 1:
			tempDisperse_patch_prob.append(float(Disperse_patch_prob[isub].split('|')[icdtime]))
		else:
			tempDisperse_patch_prob.append(float(Disperse_patch_prob[isub]))	
		if len(Str_patch_prob[isub].split('|')) > 1:
			tempStr_patch_prob.append(float(Str_patch_prob[isub].split('|')[icdtime]))
		else:
			tempStr_patch_prob.append(float(Str_patch_prob[isub]))
		if len(MgOut_patch_prob[isub].split('|')) > 1:
			tempMgOut_patch_prob.append(float(MgOut_patch_prob[isub].split('|')[icdtime]))
		else:
			tempMgOut_patch_prob.append(float(MgOut_patch_prob[isub]))
		if len(outsizevals[isub].split('|')) > 1:
			tempoutsize.append(outsizevals[isub].split('|')[icdtime])
		else:
			tempoutsize.append(outsizevals[isub])
		if len(backsizevals[isub].split('|')) > 1:
			tempbacksize.append(backsizevals[isub].split('|')[icdtime])
		else:
			tempbacksize.append(backsizevals[isub])		
		if len(outgrowdays[isub].split('|')) > 1:
			tempoutgrow.append(outgrowdays[isub].split('|')[icdtime])
		else:
			tempoutgrow.append(outgrowdays[isub])
		if len(backgrowdays[isub].split('|')) > 1:
			tempbackgrow.append(backgrowdays[isub].split('|')[icdtime])
		else:
			tempbackgrow.append(backgrowdays[isub])
		
		if len(outsizevals_sd[isub].split('|')) > 1:
			tempoutsize_sd.append(outsizevals_sd[isub].split('|')[icdtime])
		else:
			tempoutsize_sd.append(outsizevals_sd[isub])
		
		if len(backsizevals_sd[isub].split('|')) > 1:
			tempbacksize_sd.append(backsizevals_sd[isub].split('|')[icdtime])
		else:
			tempbacksize_sd.append(backsizevals_sd[isub])
		
		if len(outgrowdays_sd[isub].split('|')) > 1:
			tempoutgrow_sd.append(outgrowdays_sd[isub].split('|')[icdtime])
		else:
			tempoutgrow_sd.append(outgrowdays_sd[isub])
		
		if len(backgrowdays_sd[isub].split('|')) > 1:
			tempbackgrow_sd.append(backgrowdays_sd[isub].split('|')[icdtime])
		else:
			tempbackgrow_sd.append(backgrowdays_sd[isub])
		
		if len(K[isub].split('|')) > 1:
			tempK.append(int(K[isub].split('|')[icdtime]))
		else:
			tempK.append(int(K[isub]))
			
		if len(Kstd[isub].split('|')) > 1:
			tempKstd.append(int(Kstd[isub].split('|')[icdtime]))
		else:
			tempKstd.append(int(Kstd[isub]))
		
		if len(popmort_back[isub].split('|')) > 1:
			temppopmort_back.append(popmort_back[isub].split('|')[icdtime])
		else:
			temppopmort_back.append(popmort_back[isub])
		
		if len(popmort_out[isub].split('|')) > 1:
			temppopmort_out.append(popmort_out[isub].split('|')[icdtime])
		else:
			temppopmort_out.append(popmort_out[isub])
		
		if len(eggmort[isub].split('|')) > 1:
			tempeggmort.append(eggmort[isub].split('|')[icdtime])
		else:
			tempeggmort.append(eggmort[isub])
		if len(popmort_back_sd[isub].split('|')) > 1:
			temppopmort_back_sd.append(float(popmort_back_sd[isub].split('|')[icdtime]))
		else:
			temppopmort_back_sd.append(float(popmort_back_sd[isub]))
					
		if len(popmort_out_sd[isub].split('|')) > 1:		
			temppopmort_out_sd.append(float(popmort_out_sd[isub].split('|')[icdtime]))
		else:
			temppopmort_out_sd.append(float(popmort_out_sd[isub]))
		
		if len(eggmort_sd[isub].split('|')) > 1:
			tempeggmort_sd.append(eggmort_sd[isub].split('|')[icdtime])
		else:
			tempeggmort_sd.append(eggmort_sd[isub])
		
		if len(pop_capture_back[isub].split('|')) > 1:
			temppopCapBack.append(pop_capture_back[isub].split('|')[icdtime])
		else:
			temppopCapBack.append(pop_capture_back[isub])
			
		if len(pop_capture_out[isub].split('|')) > 1:		
			temppopCapOut.append(pop_capture_out[isub].split('|')[icdtime])
		else:
			temppopCapOut.append(pop_capture_out[isub])			
		
		if len(N0_pass[isub].split('|')) > 1:
			
			tempN0.append(N0_pass[isub].split('|')[icdtime])
		else:
			tempN0.append(N0_pass[isub])
			
		if len(allefreqfiles_pass[isub].split('|')) > 1:
			tempAllelefile.append(allefreqfiles_pass[isub].split('|')[icdtime])
		else:
			tempAllelefile.append(allefreqfiles_pass[isub])
			
		if len(classvarsfiles_pass[isub].split('|')) > 1: # More than one time value
			tempClassVarsfile.append(classvarsfiles_pass[isub].split('|')[icdtime])
		else: # No bar split, just one value
			tempClassVarsfile.append(classvarsfiles_pass[isub])
			
		if len(comp_coef_pass[isub].split('|')) > 1: # More than one time value
			tempcompcoef.append(comp_coef_pass[isub].split('|')[icdtime])
		else: # No bar split, just one value
			tempcompcoef.append(comp_coef_pass[isub])
		
		# 1 and 2 Locus Fitness Variables
		if len(fitvals) > 0:			
			tempfitvals.append([])
			for i in range(len(fitvals[isub])): # loop through 9 genotype combinations
				if len(fitvals[isub][i].split('|')) > 1: # More than one time
					# Error Check to make sure the times match
					if len(fitvals[isub][i].split('|')) != len(cdclimgentime):
						print('CDCLIMATE specified times - length of - does not match - length of - fitness values given.')
						sys.exit(-1)
					if len(fitvals[isub][i].split('|')[icdtime].split('~')) > 1: # G or M parameters []
						tempfitvals[isub].append(fitvals[isub][i].split('|')[icdtime].split('~'))
					else: # just fitness values 1 or 2
						tempfitvals[isub].append(fitvals[isub][i].split('|')[icdtime])
				else: # Just one value
					if len(fitvals[isub][i].split('~')) > 1: # G or M parameters []
						tempfitvals[isub].append(fitvals[isub][i].split('~'))
					else: # just fitness values
						tempfitvals[isub].append(fitvals[isub][i])					
			# Error checks
			if cdevolveans == 'G':
				if len(tempfitvals[isub][0][0].split(':')) != 6: 
					print('CDEVOLVE answer is G, 6 parameter values must be entered for growth equation, see user manual.')
					sys.exit(-1)
			if cdevolveans == 'MG_ind' or cdevolveans == 'MG_link':
				if len(tempfitvals[isub][0][0].split(':')) != 2 and len(len(tempfitvals[isub][3][0].split(':'))) != 6:
					print('CDEVOLVE answer is MG, 6 parameter values must be entered for growth equation and 2 for maturation equation, see user manual.')
					sys.exit(-1)					
		# MLocus Selection Variables
		if len(xvars_betas_pass) > 0:		
			tempxvars_betas.append([])
			for i in range(len(xvars_betas_pass[isub])):
				if len(xvars_betas_pass[isub][i].split('|')) > 1:
					tempxvars_betas[isub].append(xvars_betas_pass[isub][i].split('|')[icdtime])
				else:
					tempxvars_betas[isub].append(xvars_betas_pass[isub][i])		
		
		# Error check on grow days, must be equal to 365 if both entered
		if tempoutsize[isub] != 'N' and tempbacksize[isub] != 'N':
			if float(tempoutgrow[isub]) + float(tempbackgrow[isub]) > 365.:
				print('Grow days back and out must be <= 365.')
				sys.exit(-1)
				
		if len(outhabvals_pass[isub].split('|')) > 1:
			tempouthabvals.append(outhabvals_pass[isub].split('|')[icdtime])
		else:
			tempouthabvals.append(outhabvals_pass[isub])
		if len(backhabvals_pass[isub].split('|')) > 1:
			tempbackhabvals.append(backhabvals_pass[isub].split('|')[icdtime])
		else:
			tempbackhabvals.append(backhabvals_pass[isub])
			
	# -----------------------------------------------------
	# Class Vars
	# -----------------------------------------------------
	tupAgeFile = InitializeAge(K,tempClassVarsfile,datadir)	
	#agelst = tupAgeFile[0]
	#age_size_mean = tupAgeFile[1]
	#age_size_std = tupAgeFile[2]
	#sexratio = tupAgeFile[3]	
	age_percmort_out = tupAgeFile[4]
	age_percmort_out_sd = tupAgeFile[5]
	age_percmort_back = tupAgeFile[6]
	age_percmort_back_sd = tupAgeFile[7]
	size_percmort_out = tupAgeFile[8]
	size_percmort_out_sd = tupAgeFile[9]
	size_percmort_back = tupAgeFile[10]
	size_percmort_back_sd = tupAgeFile[11]
	age_MgOUT = tupAgeFile[12]
	age_MgBACK = tupAgeFile[13]
	age_S = tupAgeFile[14]
	age_DispProb = tupAgeFile[15]
	age_mature = tupAgeFile[16]
	age_mu = tupAgeFile[17]
	age_sigma = tupAgeFile[18]
	f_leslie = tupAgeFile[19]
	f_leslie_std = tupAgeFile[20]
	age_cap_out = tupAgeFile[21]
	age_cap_back = tupAgeFile[22]	
			
	# Return this functions variables - mate, dispOut, dispBack, Str_patch_prob, dispLocal (order 1) FXX, MXY, MYY, FYY (order 2) matrix,thresh,scalemin,scalemax,a,b,c,no
	tupClimate = matecdmatrix,FXXdispOutcdmatrix,MXYdispOutcdmatrix,MYYdispOutcdmatrix,FYYdispOutcdmatrix,\
	FXXdispBackcdmatrix,MXYdispBackcdmatrix,MYYdispBackcdmatrix,FYYdispBackcdmatrix,\
	FXXStrcdmatrix, MXYStrcdmatrix,MYYStrcdmatrix, FYYStrcdmatrix,\
	FXXdispLocalcdmatrix, MXYdispLocalcdmatrix,MYYdispLocalcdmatrix,FYYdispLocalcdmatrix,\
	matemovethresh,FXXdispOutthresh,MXYdispOutthresh,MYYdispOutthresh,FYYdispOutthresh,\
	FXXdispOutthresh,MXYdispOutthresh,MYYdispOutthresh,FYYdispOutthresh,\
	FXXStrthresh,MXYStrthresh,MYYStrthresh,FYYStrthresh,\
	FXXdispLocalthresh,MXYdispLocalthresh,MYYdispLocalthresh,FYYdispLocalthresh,\
	mate_ScaleMin,FXXdispOut_ScaleMin,MXYdispOut_ScaleMin,MYYdispOut_ScaleMin,FYYdispOut_ScaleMin,\
	FXXdispBack_ScaleMin,MXYdispBack_ScaleMin,MYYdispBack_ScaleMin,FYYdispBack_ScaleMin,\
	FXXStr_ScaleMin,MXYStr_ScaleMin,MYYStr_ScaleMin,FYYStr_ScaleMin,\
	FXXdispLocal_ScaleMin,MXYdispBack_ScaleMin,MYYdispBack_ScaleMin,FYYdispBack_ScaleMin,\
	mate_ScaleMax,FXXdispOut_ScaleMax,MXYdispOut_ScaleMax,MYYdispOut_ScaleMax,FYYdispOut_ScaleMax,\
	FXXdispBack_ScaleMax,MXYdispBack_ScaleMax,MYYdispBack_ScaleMax,FYYdispBack_ScaleMax,\
	FXXStr_ScaleMax,MXYStr_ScaleMax,MYYStr_ScaleMax,FYYStr_ScaleMax,\
	FXXdispLocal_ScaleMax,MXYdispBack_ScaleMax,MYYdispBack_ScaleMax,FYYdispBack_ScaleMax,\
	matemoveparA,FXXdispmoveOutparA,MXYdispmoveOutparA,MYYdispmoveOutparA,FYYdispmoveOutparA,\
	FXXdispmoveBackparA,MXYdispmoveBackparA,MYYdispmoveBackparA,FYYdispmoveBackparA,\
	FXXStrparA,MXYStrparA,MYYStrparA,FYYStrparA,\
	FXXdispLocalparA,MXYdispLocalparA,MYYdispLocalparA,FYYdispLocalparA,\
	matemoveparB,FXXdispmoveOutparB,MXYdispmoveOutparB,MYYdispmoveOutparB,FYYdispmoveOutparB,\
	FXXdispmoveBackparB,MXYdispmoveBackparB,MYYdispmoveBackparB,FYYdispmoveBackparB,\
	FXXStrparB,MXYStrparB,MYYStrparB,FYYStrparB,\
	FXXdispLocalparB,MXYdispLocalparB,MYYdispLocalparB,FYYdispLocalparB,\
	matemoveparC,FXXdispmoveOutparC,MXYdispmoveOutparC,MYYdispmoveOutparC,FYYdispmoveOutparC,\
	FXXdispmoveBackparC,MXYdispmoveBackparC,MYYdispmoveBackparC,FYYdispmoveBackparC,\
	FXXStrparC,MXYStrparC,MYYStrparC,FYYStrparC,\
	FXXdispLocalparC,MXYdispLocalparC,MYYdispLocalparC,FYYdispLocalparC,\
	matemoveno,FXXdispmoveOutno,MXYdispmoveOutno,MYYdispmoveOutno,FYYdispmoveOutno,\
	FXXdispmoveBackno,MXYdispmoveBackno,MYYdispmoveBackno,FYYdispmoveBackno,\
	FXXStrno,MXYStrno,MYYStrno,FYYStrno,\
	FXXdispLocalno,MXYdispLocalno,MYYdispLocalno,FYYdispLocalno,\
	tempMgOut_patch_prob,tempStr_patch_prob,tempoutsize,tempbacksize,tempoutgrow,tempbackgrow,tempfitvals,tempK,temppopmort_back,temppopmort_out,tempeggmort,tempKstd,temppopmort_back_sd,temppopmort_out_sd,tempeggmort_sd,tempoutsize_sd,tempbacksize_sd,tempoutgrow_sd,tempbackgrow_sd,temppopCapBack,temppopCapOut,tempN0,tempAllelefile,tempClassVarsfile,assortmateModel, assortmateC,subpopmort_mat,tempcompcoef,tempbetas_selection,tempxvars_betas,tempouthabvals,tempbackhabvals,plastic_signalresp,plastic_behaviorresp,muterate,age_percmort_out,age_percmort_out_sd,age_percmort_back,age_percmort_back_sd,size_percmort_out,size_percmort_out_sd,size_percmort_back,size_percmort_back_sd,age_MgOUT,age_MgBACK,age_S,age_DispProb,age_mature,age_mu,age_sigma,f_leslie,f_leslie_std,age_cap_out,age_cap_back,tempMgBack_patch_prob,tempDisperse_patch_prob	
	return tupClimate
	#End::DoCDClimate()

# ---------------------------------------------------------------------------------------------------	
def DoStochasticUpdate(K_mu,K_std,popmort_back_mu,popmort_back_sd,popmort_out_mu,popmort_out_sd,eggmort_mu,eggmort_sd,outsizevals_mu,outsizevals_sd,backsizevals_mu,backsizevals_sd,outgrowdays_mu,outgrowdays_sd,backgrowdays_mu,backgrowdays_sd,age_percmort_out_mu,age_percmort_out_sd,age_percmort_back_mu,age_percmort_back_sd,size_percmort_out_mu,size_percmort_out_sd,size_percmort_back_mu,size_percmort_back_sd,age_percmort_back_mu_egg,age_percmort_back_sd_egg,cor_mat,age_mu, age_sigma,f_leslie_mu,f_leslie_std,sexchromo):	
	'''
	Here update any stochastic variables. Add in Todd and Ng method for unbias draw.
	Generate correlated deviates
	'''
	
	# --------------------------------
	# For the patch specific parameters
	# Get correlated means
	# -------------------------------
	K = []
	popmort_back = []
	popmort_out = []
	eggmort_patch = []
	outsizevals = []
	backsizevals = []
	outgrowdays = []
	backgrowdays = []
	
	# For no cor_mat answer 'N'
	if isinstance(cor_mat,str):			
		for isub in range(len(K_mu)):
			# K ------------------
			stochastic_update(K_mu[isub],K_std[isub],K)
			# mort out --------------
			stochastic_update(popmort_out_mu[isub],popmort_out_sd[isub],popmort_out)
			# mort back ---------------
			stochastic_update(popmort_back_mu[isub],popmort_back_sd[isub],popmort_back)
			# egg mort ------------------
			stochastic_update(eggmort_mu[isub],eggmort_sd[isub],eggmort_patch)
			# temp vals out ----------------
			stochastic_update(outsizevals_mu[isub],outsizevals_sd[isub],outsizevals, True)
			# temp vals back ----------------
			stochastic_update(backsizevals_mu[isub],backsizevals_sd[isub],backsizevals, True)
			# grow days out ----------------
			mu = outgrowdays_mu[isub]
			sigma = outgrowdays_sd[isub]
			if mu != 'N':
				mu = float(mu)
				sigma = float(sigma)
				# Case here for sigma == 0
				if sigma != 0:
					# Call a truncated normal here
					#lower, upper = 0,365
					#X = truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
					X = np.random.normal(mu,sigma)
					if X < 0:
						X = 0.
					elif X > 365:
						X = 365.
					outgrowdays.append(round(X,3))
				else:
					outgrowdays.append(mu)				
			else:
				outgrowdays.append(mu)
			# grow days back ----------------
			mu = backgrowdays_mu[isub]
			sigma = backgrowdays_sd[isub]
			if mu != 'N':
				mu = float(mu)
				sigma = float(sigma)
				# Case here for sigma == 0
				if sigma != 0:
					# Call a truncated normal here
					#lower, upper = 0,365
					#X = truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
					X = np.random.normal(mu,sigma)
					if X < 0:
						X = 0.
					elif X > 365:
						X = 365.
					backgrowdays.append(round(X,3))
				else:
					backgrowdays.append(mu)
			else:
				backgrowdays.append(mu)
	
	# If cor_mat used
	else:
		# Generate deviates, use same deviates for each patch
		# Generate 8 independent normally distributed random
		# variables (with mean 0 and std. dev. 1). Use same for each patch
		#x = norm.rvs(size=(len(cor_mat[0]), patches))?
		x = norm.rvs(size=(len(cor_mat[0]), 1))
		
		# Combine patch mu, patch sd
		patch_mu = []
		patch_mu.append(K_mu)
		patch_mu.append(popmort_out_mu)
		patch_mu.append(outsizevals_mu)
		patch_mu.append(outgrowdays_mu)
		patch_mu.append(popmort_back_mu)
		patch_mu.append(backsizevals_mu)
		patch_mu.append(backgrowdays_mu)
		patch_mu.append(eggmort_mu)
		patch_sd = []
		patch_sd.append(K_std)
		patch_sd.append(popmort_out_sd)
		patch_sd.append(outsizevals_sd)
		patch_sd.append(outgrowdays_sd)
		patch_sd.append(popmort_back_sd)
		patch_sd.append(backsizevals_sd)
		patch_sd.append(backgrowdays_sd)
		patch_sd.append(eggmort_sd)
		patch_mu = np.asarray(patch_mu)
		patch_sd = np.asarray(patch_sd,dtype = 'float')
		
		# Create covariance matrix for each patch: D*R*D
		for isub in range(len(K_mu)):
			D = np.diag(patch_sd[:,isub])
			
			cov_mat = np.dot(D,cor_mat)
			cov_mat = np.dot(cov_mat,D)
					
			# Get this patch values to be correlated.
			thispatch = patch_mu[:,isub]
			# Get the parameters that are 'N' or 'E'
			Nindex = np.where(thispatch == 'N')[0]
			Eindex = np.where(thispatch == 'E')[0]
			# Turn these to 0 for now and convert to array
			thispatch[Nindex] = 0.0
			thispatch[Eindex] = 0.0
			thispatch = np.asarray(thispatch,dtype = 'float')
			
			# Generate the random samples.
			y = np.random.multivariate_normal(thispatch, cov_mat, size=1)	
			
			# Clean up data - turn values back to N or E
			y = np.asarray(y,dtype='str')[0]
			y[Nindex] = 'N'
			y[Eindex] = 'E'
			
			# Then append to each variable, int, checking < 0 cases
			# K ---------------------
			K.append(int(float(y[0])))
			if K[isub] < 0:
				K[isub] = 0
			# mort out----------------
			if y[1] == 'N' or y[1] == 'E':
				popmort_out.append(y[1])
			else:
				popmort_out.append(round(float(y[1]),3))
				if popmort_out[isub] < 0:
					popmort_out[isub] = 0		
			# temperautre out---------
			if y[2] == 'N':
				outsizevals.append(y[2])
			else:
				outsizevals.append(round(float(y[2]),3))
			if outsizevals[isub] < 0:
				outsizevals[isub] = 0
			# grow days out-----------
			y[3]
			if y[3] == 'N':
				outgrowdays.append(y[3])
			else:
				outgrowdays.append(round(float(y[3]),3))
			if outgrowdays[isub] < 0:
				outgrowdays[isub] = 0	
			# mort back---------------
			if y[4] == 'N' or y[4] == 'E':
				popmort_back.append(y[4])
			else:
				popmort_back.append(round(float(y[4]),3))
				if popmort_back[isub] < 0:
					popmort_back[isub] = 0	
			# temperature back--------
			if y[5] == 'N':
				backsizevals.append(y[5])
			else:
				backsizevals.append(round(float(y[5]),3))
			if backsizevals[isub] < 0:
				backsizevals[isub] = 0
			# grow days back ---------
			y[6]
			if y[6] == 'N':
				backgrowdays.append(y[6])
			else:
				backgrowdays.append(round(float(y[6]),3))
			if backgrowdays[isub] < 0:
				backgrowdays[isub] = 0
			# eggmort-----------------
			y[7]
			eggmort_patch.append(round(float(y[7]),3))
			if eggmort_patch[isub] < 0:
				eggmort_patch[isub] = 0
	
	# --------------------------------
	# For the age specific parameters
	# -------------------------------
	age_percmort_out = []
	age_percmort_back = []
	size_percmort_out = []
	size_percmort_back = []
	f_ind = []
	f_leslie = []
	# Split up into the subpops
	for isub in range(len(age_percmort_back_mu)):
		age_percmort_out.append([])
		age_percmort_back.append([])
		size_percmort_out.append([])
		size_percmort_back.append([])
		f_ind.append([])
		f_leslie.append([])
		# Split up for patch values
		for ifile in range(len(age_percmort_back_mu[isub])):
			age_percmort_out[isub].append([])
			age_percmort_back[isub].append([])
			size_percmort_out[isub].append([])
			size_percmort_back[isub].append([])
			f_ind[isub].append([])
			f_leslie[isub].append([])
			# Then loop through each class value
			for iage in range(len(age_percmort_back_mu[isub][ifile])):				
				# age mu or fecundity for each individual (mean eggs) -----------------no sex split
				# ---------------------------------------------------------------------------------
				stochastic_update(age_mu[isub][ifile][iage],age_sigma[isub][ifile][iage],f_ind[isub][ifile])
								
				# f leslie --- -----------------no sex split
				# ------------------------------------------
				stochastic_update(f_leslie_mu[isub][ifile][iage],f_leslie_std[isub][ifile][iage],f_leslie[isub][ifile])			
				
				# age mort back ----------------
				# Split if sex ratios given
				temp_age_mu = sexsplit(age_percmort_back_mu[isub][ifile][iage],sexchromo)
				temp_age_sd = sexsplit(age_percmort_back_sd[isub][ifile][iage],sexchromo)
				temp_store_mu = []
				for isex in range(len(temp_age_mu)):
					stochastic_update(temp_age_mu[isex],temp_age_sd[isex],temp_store_mu)
				# Merge back with ~ 	
				age_percmort_back[isub][ifile].append('~'.join(np.asarray(temp_store_mu,str)))
				
				# age mort out ----------------
				# Split if sex ratios given
				temp_age_mu = sexsplit(age_percmort_out_mu[isub][ifile][iage],sexchromo)
				temp_age_sd = sexsplit(age_percmort_out_sd[isub][ifile][iage],sexchromo)
				temp_store_mu = []
				for isex in range(len(temp_age_mu)):
					stochastic_update(temp_age_mu[isex],temp_age_sd[isex],temp_store_mu)
				# Merge back with ~ 	
				age_percmort_out[isub][ifile].append('~'.join(np.asarray(temp_store_mu,str)))
				
				# size mort back  ----------------
				# Split if sex ratios given
				temp_age_mu = sexsplit(size_percmort_back_mu[isub][ifile][iage],sexchromo)
				temp_age_sd = sexsplit(size_percmort_back_sd[isub][ifile][iage],sexchromo)
				temp_store_mu = []				
				for isex in range(len(temp_age_mu)):
					stochastic_update(temp_age_mu[isex],temp_age_sd[isex],temp_store_mu)
				# Merge back with ~ 	
				size_percmort_back[isub][ifile].append('~'.join(np.asarray(temp_store_mu,str)))
								
				# size mort out ----------------
				# Split if sex ratios given
				temp_age_mu = sexsplit(size_percmort_out_mu[isub][ifile][iage],sexchromo)
				temp_age_sd = sexsplit(size_percmort_out_sd[isub][ifile][iage],sexchromo)
				temp_store_mu = []
				for isex in range(len(temp_age_mu)):
					stochastic_update(temp_age_mu[isex],temp_age_sd[isex],temp_store_mu)
				# Merge back with ~ 	
				size_percmort_out[isub][ifile].append('~'.join(np.asarray(temp_store_mu,str)))
			
	# ----------------------------------
	# For one numbers - PopVars
	# ----------------------------------
	temp_store_mu = []		
	stochastic_update(age_percmort_back_mu_egg,age_percmort_back_sd_egg,temp_store_mu)
	eggmort_age = temp_store_mu[0]
	if eggmort_age != 'N':
		if eggmort_age > 1.:
			eggmort_age = 1.
			
	return K,popmort_back,popmort_out,eggmort_patch,outsizevals,backsizevals,outgrowdays,backgrowdays,age_percmort_out,age_percmort_back,	size_percmort_out,size_percmort_back,eggmort_age,f_ind,f_leslie
	#End::DoStochasticUpdate()
	
# ---------------------------------------------------------------------------------------------------	 
def DoPreProcess(outdir,datadir,irun,ithmcrun,xyfilename,loci,alleles,gen,logfHndl,cdevolveans,cdinfect,subpopemigration,subpopimmigration,sizeans,burningen_cdevolve,cor_mat_ans,inheritans_classfiles,sexans,spcNO,ibatch,betaFile_selection,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,eggFreq_mu,eggFreq_sd):

	'''
	DoPreProcess()
	This function does all the pre-processing work before
	CDPOP begins its time loops.
	'''
	# ----------------------------
	# Create directory
	# ----------------------------	
	ithmcrundir = outdir+'run'+str(irun)+'batch'+str(ibatch)+'mc'+str(ithmcrun)+'species'+str(spcNO)+'/'
	os.mkdir(ithmcrundir)
	
	# ------------------------------------------------------------------
	# Read in xy points file and store info in list
	# ------------------------------------------------------------------
	xy = ReadXY(xyfilename)
	
	# Error statement for column data
	if cdevolveans.split('_')[0] != 'P':
		if len(xy[1]) != 50:
			print('PatchVars input file is not correct version, see example input files.')
			sys.exit(-1)
	else: # MLoci selection is on, then extra columns possible in PatchVars
		if (len(xy[1]) - int(cdevolveans.split('_')[1].split('X')[1])) != 50:
			print('PatchVars input file must be 48 columns plus the specified number of variables operating in the multiple loci selection model; see example input files.')
			sys.exit(-1)
		
	# Store all information in lists by variable name and store
	Pop = []
	xgridpop = []
	ygridpop = []
	PopTag = []
	K_temp = []
	Kstd_temp = []
	N0_temp = []	
	natal_patches = []
	migrate_patches = []
	#disperse_patches = []
	#stray_patches = []
	allefreqfiles_temp = []
	classvarsfiles_temp = []
	popmort_out = []
	popmort_out_sd = []
	popmort_back = []
	popmort_back_sd = []
	MgOut_patch_prob = []
	MgBack_patch_prob = []
	Str_patch_prob = []	
	Disperse_patch_prob = []
	newmortperc = []
	newmortperc_sd = []
	setmigrate = []
	outsizevals = [] # growth values out
	outgrowdays = [] # grow days out
	backsizevals = [] # growth values back
	backgrowdays = [] # grow days back
	outsizevals_sd = [] # growth values out
	outgrowdays_sd = [] # grow days out
	backsizevals_sd = [] # growth values back
	backgrowdays_sd = [] # grow days back
	pop_capture_back_pass = [] # Grab first one to go into first DoUpdate
	pop_capture_out = []
	fitvals = [] # selection values
	outhabvals = []
	backhabvals = []
	comp_coef_temp = [] # Competition coef by patch
	xvars = [] # selection for multilocus/hindex options
	for i in range(len(xy)-1):
		# Add line number reference here for future versions
		#refLine = getframeinfo(inspect.currentframe()).lineno
		Pop.append(xy[i+1][0])
		xgridpop.append(float(xy[i+1][1]))
		ygridpop.append(float(xy[i+1][2]))
		PopTag.append(xy[i+1][3])
		K_temp.append(xy[i+1][4])
		Kstd_temp.append(xy[i+1][5])
		N0_temp.append(xy[i+1][6])
		natal_patches.append(int(xy[i+1][7]))
		migrate_patches.append(int(xy[i+1][8]))
		#disperse_patches.append(int(xy[i+1][9]))
		#stray_patches.append(int(xy[i+1][10]))		
		allefreqfiles_temp.append(xy[i+1][9])
		classvarsfiles_temp.append(xy[i+1][10])
		popmort_out.append(xy[i+1][11])
		popmort_out_sd.append(xy[i+1][12])
		popmort_back.append(xy[i+1][13])
		popmort_back_sd.append(xy[i+1][14])
		newmortperc.append(xy[i+1][15])
		newmortperc_sd.append(xy[i+1][16])
		MgOut_patch_prob.append(xy[i+1][17])
		setmigrate.append(xy[i+1][18])
		MgBack_patch_prob.append(xy[i+1][19])
		Str_patch_prob.append(xy[i+1][20])			
		Disperse_patch_prob.append(xy[i+1][21])		
		outsizevals.append(xy[i+1][22])
		outsizevals_sd.append(xy[i+1][23])
		outgrowdays.append(xy[i+1][24])
		outgrowdays_sd.append(xy[i+1][25])
		backsizevals.append(xy[i+1][26])
		backsizevals_sd.append(xy[i+1][27])
		backgrowdays.append(xy[i+1][28])		
		backgrowdays_sd.append(xy[i+1][29])
		pop_capture_out.append(xy[i+1][30])
		pop_capture_back_pass.append(xy[i+1][31])
		outhabvals.append(xy[i+1][32])
		backhabvals.append(xy[i+1][33])
		comp_coef_temp.append(xy[i+1][49]) # Bar allowed here
		xvars_indexspot = 50 # for getting the xvars spot for indexing
		if cdevolveans == '1' or cdevolveans == 'M' or cdevolveans == '1_mat' or cdevolveans == 'stray':
			fitvals.append([xy[i+1][34],xy[i+1][35],xy[i+1][36]])
		elif cdevolveans == 'G':
			fitvals.append([xy[i+1][37],xy[i+1][38],xy[i+1][39]])
		elif cdevolveans == '2' or cdevolveans == '2_mat':
			fitvals.append([xy[i+1][40],xy[i+1][41],xy[i+1][42],xy[i+1][43],xy[i+1][44],xy[i+1][45],xy[i+1][46],xy[i+1][47],xy[i+1][48]])
		elif cdevolveans == 'MG_ind' or cdevolveans == 'MG_link' or cdevolveans == '1_G_ind' or cdevolveans == '1_G_link':
			fitvals.append([xy[i+1][34],xy[i+1][35],xy[i+1][36],xy[i+1][37],xy[i+1][38],xy[i+1][39]])
		elif cdevolveans.split('_')[0] == 'P':
			xvars.append([])
			for ixvars in range(int(cdevolveans.split('_')[1].split('X')[1])):
				xvars[i].append(xy[i+1][xvars_indexspot+ixvars])		
			
	# Delete x variable
	del(xy)
	
	# --------------------------------------------
	# Read in correlation matrix
	# --------------------------------------------
	if cor_mat_ans == 'N':
		cor_mat = 'N'
	else:
		cor_mat = ReadXY(datadir+cor_mat_ans)
		cor_mat = np.asarray(np.asarray(cor_mat)[1:,1:],dtype='float')
		
	# -----------------------------------------------------------------------
	# Extract variables needed for initialization that vary with cdclimategen
	# -----------------------------------------------------------------------	
	# Get K for the first generation, but return K_temp to be read into CDClimate module, also get first capture probability back.
	# one check on N0 > 0 and natal_patches
	K = []
	Kstd = []
	pop_capture_back = []
	N0 = []
	allefreqfiles = []
	classvarsfiles = []
	for isub in range(len(K_temp)):
		mu = int(K_temp[isub].split('|')[0])
		sigma = int(Kstd_temp[isub].split('|')[0])
		K.append(mu)
		Kstd.append(sigma)
		pop_capture_back.append(pop_capture_back_pass[isub].split('|')[0])
		#N0.append(int(N0_temp[isub].split('|')[0]))
		N0.append(N0_temp[isub].split('|')[0])
		allefreqfiles.append(allefreqfiles_temp[isub].split('|')[0])
		classvarsfiles.append(classvarsfiles_temp[isub].split('|')[0])
		#temp0check = len(np.where(np.asarray(N0[0].split(';'))=='0')[0]) # Check for N values in split string
		if ';' in N0[isub]:
			if natal_patches[isub] == 0:
				N0[isub] = '0'
		elif int(N0[isub]) > 0 and natal_patches[isub] == 0:
			stringout = 'N0 specified in nonnatal grounds. Initializing N0 at patch '+str(isub+1)+' to 0.'
			logMsg(logfHndl,stringout)
			N0[isub] = '0'		
	
	# --------------------------------
	# Initialize subpop and ID field
	# --------------------------------
	id,subpop,speciesID = InitializeID(K,N0)
	
	# --------------------------------------------
	# Initialize genetic structure - distribution 
	# --------------------------------------------
	allelst = InitializeGenes(datadir,allefreqfiles,loci,alleles)
	
	# ------------------------------------------------
	# Initialize age structure - file and distribution
	# ------------------------------------------------ 	
	tupAgeFile = InitializeAge(K,classvarsfiles,datadir)	
	agelst = tupAgeFile[0]
	age_size_mean = tupAgeFile[1]
	age_size_std = tupAgeFile[2]
	sexratio = tupAgeFile[3]	
	age_mature = tupAgeFile[16]
	
	# ------------------------------------------------------------------
	# Initialize rest of variables: age,sex,infection,genes,size,mature...
	# ------------------------------------------------------------------
	age,sex,size,infection,genes,mature,capture,layEggs,recapture,hindex,whichClassFile = InitializeVars(sexratio,agelst,cdinfect,loci,alleles,allelst,age_size_mean,age_size_std,subpop,age_mature,sizeans,cdevolveans,fitvals,burningen_cdevolve,'N',sexans,speciesID,N0,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,eggFreq_mu,eggFreq_sd)
	
	# ------------------------------------
	# For multiple files, error check here
	# ------------------------------------	
	if inheritans_classfiles == 'Hindex':
		# Assume first patch has all of the information
		if allelst[0][0][0][0][1] != 1.0:
			print('First allele frequency file must be Hindex = 1.0 for offspring inherit answer Hindex.')
			sys.exit(-1)
		if allelst[0][1][0][1][1] != 1.0:
			print('Second allele frequency file must be Hindex = 0.0 for offspring inherit answer Hindex.')
			sys.exit(-1)
	# ----------------------------------------------
	# Store class variable SubpopIN_Init
	# ----------------------------------------------
	SubpopIN = []
	
	# Get unique patches
	unisubpops = len(Pop)
	
	# Organize type data in SubpopIN - here return this and also update dynamically.
	dtype = [('NatalPop',(str,len(str(unisubpops))+1)),('EmiPop',(str,len(str(unisubpops))+1)),('ImmiPop',(str,len(str(unisubpops))+1)),('EmiCD',float),('ImmiCD',float),('age',int),('sex',(str,3)),('size',float),('mature',int),('newmature',int),('infection',int),('name',(str,100)),('MID',(str,100)),('FID',(str,100)),('capture',int),('recapture',int),('layeggs',float),('hindex',float),('classfile',(str,100)),('popID',(str,100)),('species',int),('genes',('i',sum(alleles)))]
	
	# Get N here - N maybe slighlty different then specified due to random draws
	N = []
	subpopemigration.append([]) # These are tracking variables init here
	subpopimmigration.append([]) # There are tracking variables init here
	
	# Set class variable by list of populations
	for isub in range(unisubpops):
		# Stick error statement here for multiple files of genes and classvars
		if len(allelst[isub]) != len(agelst[isub]):
			stringout = 'Warning: Number of gene files and ClassVars files specified for patch, '+str(isub)+'are not the same. Will continue with initialization using files given using first files.'
			logMsg(logfHndl,stringout)	
	
		# Storage variables
		subpopemigration[0].append([0])
		subpopimmigration[0].append([0])
		SubpopIN.append([])
		N.append([])
		
		# Get number in this patch
		noinsub = len(np.where(subpop == Pop[isub])[0])
		
		# If K is 0
		if noinsub == 0:
			N[isub].append(0) # Track N
			
		# If K does not equal 0
		else:
			N[isub].append(noinsub) # Track N
			# Update the Wright Fisher case for sex here
			# ------------------------------------------
			if sexratio[isub][0] == 'WrightFisher':
				# If the subpopulation number is not even then sys exit
				if np.mod(noinsub,2) == 1:
					print("You have WrightFisher turned and this population must be even.")
					sys.exit(-1)
				# Then create half males and females and shuffle
				sex = np.append(np.zeros(noinsub/2,"int"),np.ones(noinsub/2,"int"))
				sex = np.asarray(sex,dtype=str)
				sex[np.where(sex == '0')[0]] = 'FXX'
				sex[np.where(sex == '1')[0]] = 'MXY'
				np.random.shuffle(sex)
			
			# Loop through individuals in subpop
			# ----------------------------------
			
			for iind in range(noinsub):
				# Check if it is an NA spot
				indspot = np.where(subpop == Pop[isub])[0][iind]
							
				# Record individual to subpopulation
				# ---------------------------------				
				# Update the Wright Fisher case for sex here
				if sexratio[isub][0] == 'WrightFisher':				
					# Subpop,EmiPop(NA),ImmiPop(NA),EmiCD,ImmiCD,age,sex,infection,name/id,motherid,fatherid,capture,recapture,layeggs,genes,mature,newmature, speciesID					
					recd = (subpop[indspot],'NA','NA',-9999,-9999,age[indspot],sex[iind],size[indspot],mature[indspot],mature[indspot],infection[indspot],id[indspot],-9999,-9999,capture[indspot],recapture[indspot],layEggs[indspot],hindex[indspot],whichClassFile[indspot],PopTag[isub],speciesID[indspot],np.asarray(genes[indspot]))
					SubpopIN[isub].append(recd)
				
				# Not special Wright Fisher case
				else:			
					# Subpop,EmiPop(NA),ImmiPop(NA),EmiCD,ImmiCD,age,sex,infection,name/id,motherid,fatherid,capture,recapture,layeggs,genes,mature, newmature speciesID
					recd = (subpop[indspot],'NA','NA',-9999,-9999,age[indspot],sex[indspot],size[indspot],mature[indspot],mature[indspot],infection[indspot],id[indspot],-9999,-9999,capture[indspot],recapture[indspot],layEggs[indspot],hindex[indspot],whichClassFile[indspot],PopTag[isub],speciesID[indspot],np.asarray(genes[indspot]))
					SubpopIN[isub].append(recd)
		# Convert to array with dytpe	
		#pdb.set_trace()
		SubpopIN[isub] = np.asarray(SubpopIN[isub],dtype=dtype)
	
	# Clean up N
	N = sum(N,[])
	
	# --------------------------
	# Error Checks
	# --------------------------	
	# For now, subpops need to be ordered 1 to N and not skipping, no 0s
	if len(np.where(np.asarray(Pop)=='0')[0]) != 0:
		print('Subpopulation identification field can not have 0 values.')
		sys.exit(-1)
	tempcheck = []
	for i in range(len(np.unique(Pop))):
		tempcheck.append(int(np.unique(Pop)[i]))
	tempcheck = np.sort(tempcheck)
	if len(tempcheck) > 1:
		for i in range(len(tempcheck)-1):
			if tempcheck[i+1]-tempcheck[i] > 1:
				print('Subpopulation identification field must be labeled sequentially or a single value.')
				sys.exit(-1)

	# ---------------------------------------------
	# Read in Beta File for Multiple loci selection
	# ---------------------------------------------
	if isinstance(betaFile_selection,(list,tuple)):
		tempbetaFile_selection = datadir+betaFile_selection[0] # Here read in first one in case multiple cdclimate
	else:
		tempbetaFile_selection = datadir+betaFile_selection
	tempbetas_selection = []	
	if cdevolveans.split('_')[0] == 'P':
		# Read in Beta File
		betavals = ReadXY(tempbetaFile_selection)
	
		# Error check on beta file - should be x number of betas alleles * number of xvars * number of loci under selection
		if (len(betavals)-1)*(len(betavals[0])-1) != int(cdevolveans.split('_')[3].split('A')[1])*int(cdevolveans.split('_')[1].split('X')[1])*int(cdevolveans.split('_')[2].split('L')[1]):
			print('Beta file for selection is incorrect. Specify the beta for each variableXlociXallele combination.')
			sys.exit(-1)
		# Then extract betavals - in order of variable, loci, allele [var][loci][allele]
		for ixvar in range(int(cdevolveans.split('_')[1].split('X')[1])):
			tempbetas_selection.append([])
			for iloci in range(int(cdevolveans.split('_')[2].split('L')[1])):
				tempbetas_selection[ixvar].append([])
				for iall in range(int(cdevolveans.split('_')[3].split('A')[1])):
					tempbetas_selection[ixvar][iloci].append(float(betavals[iall+1][ixvar*(int(cdevolveans.split('_')[2].split('L')[1]))+iloci+1]))
		# Add beta0 - will be the last spot in betas vars 
		if len(betavals[0][0]) == 0:
			tempbetas_selection.append(0.0)
		else:
			tempbetas_selection.append(float(betavals[0][0]))
	
	# Delete other storage variables
	del(size)
	del(age)
	del(id)
	del(genes)
	del(sex)
	del(infection)
	del(subpop)
	del(mature)
	del(capture)
	del(recapture)
	del(layEggs)
	del(hindex)
	del(whichClassFile)
	del(speciesID)
	
	# Return this functions variables
	tupPreProcess = ithmcrundir,\
	fitvals,allelst,subpopemigration,subpopimmigration,age_size_mean,age_size_std,xgridpop,ygridpop,\
	SubpopIN,N,K,dtype,outsizevals,backsizevals,\
	popmort_out,popmort_back,MgOut_patch_prob,Str_patch_prob,newmortperc,setmigrate,outgrowdays,backgrowdays,K_temp,Kstd_temp,Kstd,popmort_out_sd,popmort_back_sd,newmortperc_sd,outsizevals_sd,backsizevals_sd,outgrowdays_sd,backgrowdays_sd,pop_capture_back_pass,pop_capture_out,pop_capture_back,natal_patches,cor_mat,migrate_patches,N0_temp,allefreqfiles_temp,classvarsfiles_temp,PopTag,comp_coef_temp,xvars,tempbetas_selection,outhabvals,backhabvals,MgBack_patch_prob,Disperse_patch_prob
	
	return tupPreProcess	
	#End::DoPreProcess()

# ---------------------------------------------------------------------------------------------------	 		
def DoUserInput(fileans):
	
	'''
	DoUserInput()
	This function reads in the user input and 
	stores the variables.
	'''
	
	# Open file for reading
	inputfile = open(fileans,'r')

	# Read lines from the file
	lines = inputfile.readlines()

	#Close the file
	inputfile.close()

	# Create an empty matrix to append to
	inputvariables = []

	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.split(',')
		inputvariables.append(thisline)
		
	# Delete lines
	del(lines)

	return inputvariables
	
	#End::DoUserInput()

# -------------------------------------------------------------------------	
def AddIndividuals(SubpopIN,tempN0,tempAllelefile,tempClassVarsfile,datadir,loci,alleles,sizeans,cdinfect,cdevolveans,burningen_cdevolve,fitvals,dtype,N,natal_patches,gen,PopTag,sexans,logfHndl,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,eggFreq_mu,eggFreq_sd):
	'''
	AddIndividuals()
	This function adds more individuals with given classvars 
	allele frequency file information.
	'''
	
	# ---------------------------------------------
	# First error check on natal_patches grounds and N0 > 0
	# ---------------------------------------------
	for isub in range(len(tempN0)):
		if sum(np.asarray(tempN0[isub].split(';'),int)) > 0 and natal_patches[isub] == 0:
			stringout = 'N0 specified in nonnatal grounds. Initializing N0 at patch '+str(isub+1)+' to 0.'
			logMsg(logfHndl,stringout)
			tempN0[isub] = 0
	
	# --------------------------------
	# Initialize subpop field
	# --------------------------------
	id,subpop,speciesID = InitializeID(tempN0,tempN0)
	
	# --------------------------------------------
	# Initialize genetic structure - distribution 
	# --------------------------------------------
	allelst = InitializeGenes(datadir,tempAllelefile,loci,alleles)
	
	# ------------------------------------------------
	# Initialize age structure - file and distribution
	# ------------------------------------------------
	tupAgeFile = InitializeAge(tempN0,tempClassVarsfile,datadir)
		
	agelst = tupAgeFile[0]
	age_size_mean = tupAgeFile[1]
	age_size_std = tupAgeFile[2]
	sexratio = tupAgeFile[3]	
	age_mature = tupAgeFile[16]
	
	# ------------------------------------------------------------------
	# Initialize rest of variables: age,sex,infection,genes,size,mature
	# ------------------------------------------------------------------
	age,sex,size,infection,genes,mature,capture,layEggs,recapture,hindex,whichClassFile = InitializeVars(sexratio,agelst,cdinfect,loci,alleles,allelst,\
	age_size_mean,age_size_std,subpop,age_mature,sizeans,cdevolveans,fitvals,burningen_cdevolve,'Y',sexans,speciesID,tempN0,FXXmat_set,FXXmat_int,FXXmat_slope,MXYmat_set,MXYmat_int,MXYmat_slope,MYYmat_set,MYYmat_int,MYYmat_slope,FYYmat_set,FYYmat_int,FYYmat_slope,sexchromo,eggFreq_mu,eggFreq_sd)
	
	# ---------------------------------------------
	# Store class variable SubpopIN_add
	# ---------------------------------------------
	SubpopIN_keep = [] # Pass this one on
	
	# Get unique patches
	unisubpops = len(tempN0)
	
	# Set class variable by list of populations
	for isub in range(unisubpops):
		SubpopIN_add = [] # Temp array to concatenate
		
		# Get each SubpopIN pop as array
		SubpopIN_arr = np.array(SubpopIN[isub],dtype=dtype)		
		
		# Get number in this patch
		noinsub = len(np.where(subpop == str(isub+1))[0])
		
		# For tracking, N
		N[isub] = N[isub]+noinsub
		
		# If Nj does not equal 0
		if noinsub != 0:
			
			# Update the Wright Fisher case for sex here
			# ------------------------------------------
			if sexratio[isub][0] == 'WrightFisher':
				# If the subpopulation number is not even then sys exit
				if np.mod(noinsub,2) == 1:
					print("You have WrightFisher turned and this population must be even.")
					sys.exit(-1)
				# Then create half males and females and shuffle
				sex = np.append(np.zeros(noinsub/2,"int"),np.ones(noinsub/2,"int"))
				np.random.shuffle(sex)
		
			# Loop through individuals in subpop
			# ----------------------------------
			for iind in range(noinsub):				
				# Grab this index location
				indspot = np.where(subpop == str(isub+1))[0][iind]
				
				# Create ID for this individual
				# Get name
				name = 'R'+str(isub+1)+'_F'+str(isub+1)+'_m-1f-1'+'_P'+str(isub+1)+'_Y'+str(gen)+'_UN'+str(iind)
				
				# Record individual to subpopulation
				# ---------------------------------				
				# Update the Wright Fisher case for sex here
				if sexratio[isub][0] == 'WrightFisher':				
					# Subpop,EmiPop(NA),ImmiPop(NA),EmiCD,ImmiCD,age,sex,infection,name/id,motherid,fatherid,capture,recapture,layeggs,genes,mature,newmature
					recd = (subpop[indspot],subpop[indspot],subpop[indspot],-9999,-9999,age[indspot],sex[iind],size[indspot],mature[indspot],mature[indspot],infection[indspot],name,-9999,-9999,capture[indspot],recapture[indspot],layEggs[indspot],hindex[indspot],whichClassFile[indspot],PopTag[isub],speciesID[indspot],np.asarray(genes[indspot]))
				
				# Not special Wright Fisher case
				else:			
					# Subpop,EmiPop(NA),ImmiPop(NA),EmiCD,ImmiCD,age,sex,infection,name/id,motherid,fatherid,capture,recapture,layeggs,genes,mature, newmature
					recd = (subpop[indspot],subpop[indspot],subpop[indspot],-9999,-9999,age[indspot],sex[indspot],size[indspot],mature[indspot],mature[indspot],infection[indspot],name,-9999,-9999,capture[indspot],recapture[indspot],layEggs[indspot],hindex[indspot],whichClassFile[indspot],PopTag[isub],speciesID[indspot],np.asarray(genes[indspot]))
				SubpopIN_add.append(recd)
		
		# Convert to array with dytpe		
		SubpopIN_add = np.asarray(SubpopIN_add,dtype=dtype)
		
		# Append all information to temp SubpopKeep variable
		SubpopIN_keep.append(np.concatenate([SubpopIN_arr,SubpopIN_add]))
	
	# Delete storage variables
	del(size)
	del(age)
	del(id)
	del(genes)
	del(sex)
	del(infection)
	del(subpop)
	del(mature)
	del(capture)
	del(recapture)
	del(layEggs)
	del(hindex)
	del(whichClassFile)
	del(speciesID)
	
	return SubpopIN_keep
	#End::AddIndividuals()