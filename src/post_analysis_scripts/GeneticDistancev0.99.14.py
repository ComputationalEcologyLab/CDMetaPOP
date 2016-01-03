# ----------------------------------------------------------------------------
# v0.99.14 - 2015July - Create version for CDmetaPOP
# v7.0 - 2014Sept - NAs ignored and calculate a ED matrix at every time step. Speed up with numby intersection
# v6.0 - 2014April -- NAs ignored.
# v5.0 - 2012April5 -- Will calculate on specified generations.
# v4.0 - 2012March29 -- Will list folders in directory. Will not
#	work for NA.
# v3.0 - 2011June27 -- Assumes files in grid*.csv, but for new version of
#	CDPOP v1.0. Just Dps.
# v2.0 - 2008Dec22 -- Assumes files are in the gridNNN format, but reads in 
#	*.csv, not grid*.csv
# v1.0 - March 2008 -- Assumes files are named grid*.csv
# GeneticDistance.py
# Author: Erin L Landguth
# Created: March 2008 
# Description: This program calculates a genetic distance matrix using:
#	Bray-Curtis Method:
#	formula 1 - 2W/(A+B), where W is the minimum value between the two comp-
#	arison's A and B.  The specific application is to calculate this distance
#	matrix for genotypes of a population with n individuals: Genetic Distance.
#	Proportion of Shared Alleles:
#	Nei's:
#	1 - sqrt(ithfreq*jthfreq)/loci
#	Proportion of Shared alleles: 
#	1 - proportion of shared alleles between individuals.	
# Program Input: directory of *.csv files
# Program Output: oldname+Gdmatrix.csv 
# Program Steps:
#	1. User input information.
#	2. fileList all of the *.csv in directory
#	3. Run genetic distance method
# ----------------------------------------------------------------------------

import glob							# The power of glob
# Import statements
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
import time, datetime,os,pdb		# Other libraries 

# Timing events: start
start_time = datetime.datetime.now()

# ------------------------------------------
# Step 1: Get user information
# ------------------------------------------

# Store directory path name
directory = 'D:/projects/CDmetaPOP/BullTrout/Results/Gen-1_indfiles_riverine/test/'

# Sample the loci to run analysis on (for loci under selection)
loci = 14
#selloci = xrange(1,20)
#selloci = xrange(0,1)
selloci = xrange(loci)

# Number of alleles per locus
noalleperlocus = 34
alleles = int(noalleperlocus)*np.ones(loci,int)
# If variable alleles per locus
#alleles = array([6,10])

# The number of individuals
nogrids = 1587

# The generations to run
gen = range(21,101,1)
gen = [-1]

# -------------------------------	
# Step 2: List files in directory
# -------------------------------
# List folders in this dir
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(directory)

# Loop through folderList
for ifold in xrange(len(folderList)):
	
	# -----------------------------------	
	# Step 3: Run genetic distance method
	# -----------------------------------
		
	# ------------ Genetic Distance Matrix: Proportion of shared alleles -----------------------
	# List all files with .csv extensions (but only the grid ones)

	print '\n'
	print 'Creating the proportion of shared alleles genetic distance matrices...'

	# Get the first globbed file read in
	for i in xrange(len(gen)):
		
		# Open file for reading
		inputfile = open(folderList[ifold]+'/ind'+str(gen[i])+'.csv','r')
		
		# Read lines from the file
		lines = inputfile.readlines()
				
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		x = []
		
		# Split up each line in file and append to empty matrix, x
		for l in lines:
			thisline = l.split(',')
			x.append(thisline)
				
		# Store genetic information: genes[bear], but need them as float values
		Xvals = [] # Read in x and y
		Yvals = []		
		genes = []
		tempcount = 0
		for k in range(len(x)-1):			
			# Check for NA values
			if x[1+k][5] != 'NA':				
				# Read in X,Y
				Xvals.append(float(x[k+1][1]))
				Yvals.append(float(x[k+1][2]))
				# Create spot in genes
				#genes.append([])
				# Get list from read in file
				genes.append(x[k+1][13:int(13+np.sum(alleles))])
				tempcount = tempcount + 1
			
		# Create a matrix of zeros to be filled
		gendmatrix = np.zeros((tempcount,tempcount),float)
		geodmatrix = np.zeros((tempcount,tempcount),float)
	
		# Loop through each individual k
		for k in range(tempcount):
			tempK = np.asarray(genes[k],dtype=int)
			# Compare individual k to every other inidividual j
			for j in range(tempcount):
				tempJ = np.asarray(genes[j],dtype=int)
				# Create a tempvariable to be written over for each comparison
				tempmin=[]
				
				tempK_2 = np.where(tempK==2)[0]
				tempK_1 = np.where(tempK==1)[0]
				tempJ_2 = np.where(tempJ==2)[0]
				tempJ_1 = np.where(tempJ==1)[0]
				
				# Find the shared alleles between k and j checking the 4 conditions
				
				# Condition 2, 2
				tempmin.append(sum(np.in1d(tempK_2,tempJ_2)*2))
				# Condition 2,1
				tempmin.append(sum(np.in1d(tempK_2,tempJ_1)*1))
				# Condition 1,2
				tempmin.append(sum(np.in1d(tempK_1,tempJ_2)*1))
				# Condition 1,1
				tempmin.append(sum(np.in1d(tempK_1,tempJ_1)*1))
				
				# Write the Dps value to gendmatrix
				gendmatrix[k][j] = 1-float(np.nansum(tempmin))/(2*loci)
				
				geodmatrix[k][j] = float(np.sqrt((Xvals[k]-Xvals[j])**2+(Yvals[k]-Yvals[j])**2))
		
		# Strip directory/filename of grid and add 'Gdmatrix.csv'
		gdpathname = folderList[ifold]+'/Gdmatrix'+str(gen[i])+'.csv'
		geopathname = folderList[ifold]+'/Edmatrix'+str(gen[i])+'.csv'
		
		# Create file to write matrix to
		outputfile = open(gdpathname,'w')		
		# Sequence each row in the matrix
		for seqrow in gendmatrix:		
			# Grab each element in each row and write element to outputfile
			for ele in range(len(seqrow)):
				outputfile.write(str(seqrow[ele]))
				# Add comma
				outputfile.write(',')			
			# Return line
			outputfile.write('\n')		
		# Close file
		outputfile.close()
		# Create file to write matrix to
		outputfile = open(geopathname,'w')		
		# Sequence each row in the matrix
		for seqrow in geodmatrix:		
			# Grab each element in each row and write element to outputfile
			for ele in range(len(seqrow)):
				outputfile.write(str(seqrow[ele]))
				# Add comma
				outputfile.write(',')			
			# Return line
			outputfile.write('\n')		
		# Close file
		outputfile.close()
		
		print '\n'
		print 'The genetic distance matrix '+gdpathname+' and geodistance matrix '+geopathname+' have been created in '+str(datetime.datetime.now() -start_time)
		
