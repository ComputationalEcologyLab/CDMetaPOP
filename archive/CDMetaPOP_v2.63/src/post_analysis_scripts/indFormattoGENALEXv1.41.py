# -------------------------------------------------------------------------------------------------
# indFormattoGENALEX.py
# Author: Erin L Landguth
# Created: Feb 6, 2012
# Updated: Oct 1, 2015
# Description: This program will take the grid{nthfile}.csv genetic information 
#	files from CDPOP and convert them to the genalex.csv format.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
import pdb,random,os,sys,glob,datetime,time
from ast import literal_eval

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
#End::count_unique()

# Timing events: start
start_time = datetime.datetime.now()

# ------------------------------------	
# Step 1: Get user information
# ------------------------------------

# Directory location monte carlo folders - assume all folders to analyze
directory = 'D:/UM/test/'

# Number of loci
loci = 208

# Number of maximum possible alleles per loci
maxalleles = 36

# Number of populations
subpopno = 1

#Location in ind file for indexing
locus1spot = 17

# -------------------------------	
# Step 2: Prelimary work
# -------------------------------

# List folders in this dir
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(directory)

# Create a storage vector for alleles
alleles = maxalleles*np.ones(loci,int)

#Updated for new format with flexible alleles
alleles = np.asarray([17,29,30,12,36,20,25,28,18,11,16,13,13,9,14,13,14,32,13,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2])
	
# Create a genes vector, appending loci information with alleles to it
genes_genform = []
for iloci in xrange(loci):
	locitemp = np.arange(1,alleles[iloci]+1,1)
	genes_genform.append(list(locitemp))
	
# Begin loop through folders if specificied
for ifolder in xrange(len(folderList)):
	
	# --------------------------------------	
	# Step 3: List files in directory
	# --------------------------------------
	datfileList = glob.glob(folderList[ifolder]+'/'+'ind*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# ---------------------------	
	# Step 4: Read in file 
	# ---------------------------
	# Loop through each grid
	for igrid in xrange(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# Get list of subpop numbers for unique counts
		subpop_cdpop = []
		for ispot in xrange(nogrids): 
			subpop_cdpop.append(x[ispot+1][0])
								
		# Get the number of individuals in each pop
		unipatches = count_unique(subpop_cdpop)
		patchN = []
		# loop through all the patches
		for ipatch in xrange(subpopno):
				
			# Where is this patch in the ones that got read in, unipathes, careful of strings and integers.
			findno = np.where(unipatches[0]==str(ipatch+1))[0]
			
			# Check if it was there
			if len(findno) != 0:	
				# Then index into how many patches there
				patchN.append(unipatches[1][findno[0]])
			elif len(findno) == 0:			
				# Then just add in 0
				patchN.append(0)			
			else:
				print('Error in genepop file and patch nos.')
				sys.exit(-1)
		
		# Grab the rest of cdpop information appending to patch list
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		sex_cdpop = []
		subpop_cdpop = []
		GenFormgenes = []
		# Add spot for each patch
		for ipop in xrange(subpopno):
			id_cdpop.append([])
			x_cdpop.append([])
			y_cdpop.append([])
			age_cdpop.append([])
			sex_cdpop.append([])
			subpop_cdpop.append([])
			GenFormgenes.append([])
			for ispot in xrange(patchN[ipop]):
				# Counter into list x
				counter = ispot+sum(patchN[0:ipop])
				subpop_cdpop[ipop].append(x[1+counter][0])			
				x_cdpop[ipop].append(float(x[1+counter][1]))
				y_cdpop[ipop].append(float(x[1+counter][2]))
				id_cdpop[ipop].append(x[1+counter][3])
				sex_cdpop[ipop].append(x[1+counter][4])
				age_cdpop[ipop].append(x[1+counter][5])
				# Genes				
				GenFormgenes[ipop].append([])
				# Get each genotype
				for jspot in xrange(loci):					
					# Cdpop genes
					genes_cdpop = x[1+counter][int(locus1spot+sum(alleles[0:jspot])):int(locus1spot+sum(alleles[0:jspot+1]))]
					
					# Add gene individual spot 
					GenFormgenes[ipop][ispot].append([])
					# Loop through each allele spot at that locus
					for ithallele in xrange(alleles[jspot]):
						
						# Check if allele spot is 1
						if int(genes_cdpop[ithallele]) == 1:						
							# Then store that unique allele number
							GenFormgenes[ipop][ispot][jspot].append(genes_genform[jspot][ithallele])
						
						# Check if allele spot is 2
						elif int(genes_cdpop[ithallele]) == 2:						
							# Then store that unique allele number
							GenFormgenes[ipop][ispot][jspot].append(genes_genform[jspot][ithallele])
							GenFormgenes[ipop][ispot][jspot].append(genes_genform[jspot][ithallele])
							
						# Error check
						elif int(genes_cdpop[ithallele]) != 0:							
							print('Something wrong in gene genepop format. Email Erin.')
							sys.exit(-1)					
				
		# Delete x variable
		del(x)
			
		# Create file to write matrix to
		outputfilename = filename.split('ind')
		outputfile = open(outputfilename[0]+'/genalex_ind'+outputfilename[1].strip('.csv')+'.csv','w')
			
		# Write out the first and second line of GENALEX format
		outputfile.write(str(loci)+',')
		outputfile.write(str(nogrids)+',')
		outputfile.write(str(subpopno)+',')
		outputfile.write(str(nogrids)+',')
		for i in xrange(loci*alleles[0]):
			outputfile.write(',')
		outputfile.write('\n')
		outputfile.write(filename+',,,,')
		for i in xrange(loci*alleles[0]):
			outputfile.write(',')
		outputfile.write('\n')
		
		# Write out the third line of GENALEX format
		outputfile.write('Individual ID,Population,')
		for i in range(loci):
			outputfile.write('locus'+str(i+1)+'a,')
			outputfile.write('locus'+str(i+1)+'b,')
		# The trailing white space and XY information
		outputfile.write(',X,Y\n')
		
		# Loop through each ind spot and output		
		for ipop in xrange(subpopno):
			for ispot in xrange(patchN[ipop]):		
		
				outputfile.write('indiv'+str(ispot)+',')
				outputfile.write(str(subpop_cdpop[ipop][ispot])+',')
				
				# Loop through each locus
				for ithloci in xrange(loci):
				
					# Loop through each allele spot at that locus
					for ithallele in xrange(2):
					
						outputfile.write(str(GenFormgenes[ipop][ispot][ithloci][ithallele])+',')
					
				# Write out trailing information - white space and x,y information
				outputfile.write(',')
				outputfile.write(str(x_cdpop[ipop][ispot]).strip('[').strip(']')+',')
				outputfile.write(str(y_cdpop[ipop][ispot]).strip('[').strip(']')+'\n')	
						
		# Close file
		outputfile.close()
	print('GENALEX grid format file conversion complete.')
					
	print '\n'
	print 'Total conversion time for folder',folderList[ifolder],': ',datetime.datetime.now() - start_time