# -------------------------------------------------------------------------------------------------
# CDmetaPOP_PostProcess.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file post processing.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
import pdb,random,os,sys,glob
from ast import literal_eval 

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
	return item,count
	
	#End::w_choice_general()

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_general(loci,alleles,ithmcrundir,logfHndl,subgridtotal):
	'''
	DoGridOut_general()
	Output ind.csv in general genotype format	
	'''	
	
	subpopno = len(subgridtotal[0])	- 1
	
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in xrange(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'ind*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each ind
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
					genes_cdpop = x[1+counter][int(17+sum(alleles[0:jspot])):int(17+sum(alleles[0:jspot+1]))]
					
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
		outputfile = open(outputfilename[0]+'/general_ind'+outputfilename[1],'w')
				
		# Write out the titles that match general ind format
		title = ['Subpopulation','X','Y','ID','sex','age']
			
		# Write out the title
		for ititle in range(len(title)):
			outputfile.write(title[ititle])
			outputfile.write(',')
			
		# Write out the loci title 
		for i in range(loci-1):
			outputfile.write('Locus'+str(i+1)+'a')
			outputfile.write(',')
			outputfile.write('Locus'+str(i+1)+'b')
			outputfile.write(',')
		# To get a return character on the end of the title
		outputfile.write('Locus'+str(loci-1+1)+'a')
		outputfile.write(',')
		outputfile.write('Locus'+str(loci-1+1)+'b')
		outputfile.write('\n')	
					
		# Loop through each ind spot and output
		for ipop in xrange(subpopno):
			for ispot in xrange(patchN[ipop]):
					
				outputfile.write(subpop_cdpop[ipop][ispot]+',')
				outputfile.write(str(float(x_cdpop[ipop][ispot]))+',')
				outputfile.write(str(float(y_cdpop[ipop][ispot]))+',')
				outputfile.write(str(id_cdpop[ipop][ispot])+',')
				outputfile.write(str(sex_cdpop[ipop][ispot])+',')
				outputfile.write(str(age_cdpop[ipop][ispot])+',')				
				
				# Loop through each locus
				for ithloci in xrange(loci-1):
				
					# Loop through each allele spot at that locus
					for ithallele in xrange(2):
					
						outputfile.write(str(GenFormgenes[ipop][ispot][ithloci][ithallele])+',')
					
				# Return charater on end
				outputfile.write(str(GenFormgenes[ipop][ispot][loci-1][0])+',')
				outputfile.write(str(GenFormgenes[ipop][ispot][loci-1][1])+'\n')
											
		# Logging message
		stringout = 'The file ind'+outputfilename[0]+'/general'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	
	print('General ind format file conversion complete.')
	# End::DoGridOut_general()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_genalex(loci,alleles,ithmcrundir,logfHndl,subgridtotal):
	'''
	DoGridOut_genalex()
	Output ind.csv in genalex genotype format	
	'''	
	subpopno = len(subgridtotal[0])	- 1
		
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in xrange(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'ind*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each ind
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
			#genes[patch][indivdual][locus][allele]
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
					genes_cdpop = x[1+counter][int(17+sum(alleles[0:jspot])):int(17+sum(alleles[0:jspot+1]))]
					
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
		outputfile = open(outputfilename[0]+'/genalex_ind'+outputfilename[1].strip('.csv')+'.gen','w')
			
		# Write out the first and second line of GENALEX format
		outputfile.write(str(loci)+',')
		outputfile.write(str(nogrids)+',')
		outputfile.write(str(subpopno)+',')
		outputfile.write(str(nogrids)+'\n')
		outputfile.write(filename+'\n')		
			
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
					for ithallele in xrange(len(GenFormgenes[ipop][ispot][ithloci])):
					
						outputfile.write(str(GenFormgenes[ipop][ispot][ithloci][ithallele])+',')
					
				# Write out trailing information - white space and x,y information
				outputfile.write(',')
				outputfile.write(str(x_cdpop[ipop][ispot]).strip('[').strip(']')+',')
				outputfile.write(str(y_cdpop[ipop][ispot]).strip('[').strip(']')+'\n')	
											
		# Logging message
		stringout = 'The file ind'+outputfilename[0]+'/genalex'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	
	print('GENALEX ind format file conversion complete.')
	# End::DoGridOut_genalex()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_structure(loci,alleles,ithmcrundir,logfHndl,subgridtotal):
	'''
	DoGridOut_structure()
	Output ind.csv in structure genotype format	
	'''	
	
	subpopno = len(subgridtotal[0])	- 1	
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in xrange(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'ind*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each ind
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
		genes_cdpop = []
		GenFormgenes = []
		# Add spot for each patch
		for ipop in xrange(subpopno):
			id_cdpop.append([])
			x_cdpop.append([])
			y_cdpop.append([])
			age_cdpop.append([])
			sex_cdpop.append([])
			subpop_cdpop.append([])
			#genes[patch][indivdual][locus][allele]
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
					genes_cdpop = x[1+counter][int(17+sum(alleles[0:jspot])):int(17+sum(alleles[0:jspot+1]))]
					
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
				
		# Create file to write matrix to
		outputfilename = filename.split('ind')
		outputfile = open(outputfilename[0]+'/structure_ind'+outputfilename[1].strip('.csv')+'.stru','wb')	
		
		# Write out the first line of structure format
		for i in range(loci - 1):
			outputfile.write('locus'+str(i+1)+'_1\t')
		outputfile.write('locus'+str(loci - 1)+'_1\n')
					
		# Write out the genes of each individual 
		for ipop in xrange(subpopno):
						
			# Loop through each ind spot and output
			for ispot in xrange(patchN[ipop]):
			
				# ID write out
				outputfile.write(id_cdpop[ipop][ispot]+'\t')
				# Pop number and then a 1?
				outputfile.write(subpop_cdpop[ipop][ispot]+'\t1\t')
				
				# Loop through each locus
				for ithloci in xrange(loci - 1):
					outputfile.write(str(GenFormgenes[ipop][ispot][ithloci][0])+'\t')
					outputfile.write(str(GenFormgenes[ipop][ispot][ithloci][1])+'\t')
				outputfile.write(str(GenFormgenes[ipop][ispot][loci-1][0])+'\t')
				outputfile.write(str(GenFormgenes[ipop][ispot][loci-1][1])+'\n')
						
		# Logging message
		stringout = 'The file ind'+outputfilename[0]+'/structure'+outputfilename[1]+'.stru has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
		
	print('STRUCTURE ind format file conversion complete.')
	# End::DoGridOut_structure()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_genepop(loci,alleles,ithmcrundir,logfHndl,subgridtotal):
	'''
	DoGridOut_genalex()
	Output ind.csv in genalex genotype format	
	'''	
	subpopno = len(subgridtotal[0])	- 1
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in xrange(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'ind*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
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
					genes_cdpop = x[1+counter][int(17+sum(alleles[0:jspot])):int(17+sum(alleles[0:jspot+1]))]
					
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
		outputfile = open(outputfilename[0]+'/genepop_ind'+outputfilename[1].strip('.csv')+'.gen','w')
			
		# Write out the first and second line of GENEPOP format
		outputfile.write(outputfilename[0]+'ind'+outputfilename[1]+'\n')
		for i in xrange(loci):
			outputfile.write('LOCUS-'+str(i+1)+'\n')
		
		# Write out the genes of each individual by population
		for ipop in xrange(subpopno):
			#outputfile.write('POP'+str(ipop+1)+'\n')
			outputfile.write('POP\n')
			
			# Loop through each ind spot and output
			for ispot in xrange(patchN[ipop]):
				
				# Individual ID
				outputfile.write(id_cdpop[ipop][ispot]+', ')

				# Loop through each locus
				for ithloci in xrange(loci):
					templociname = ''
					# Loop through each allele spot at that locus
					for ithallele in xrange(len(GenFormgenes[ipop][ispot][ithloci])):
						# Add 100
						templociname = templociname + str(GenFormgenes[ipop][ispot][ithloci][ithallele]+100)
						
					outputfile.write(templociname+' ')
				# Return character
				outputfile.write('\n')				
													
		# Logging message
		stringout = 'The file ind'+outputfilename[0]+'/genepop'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	print('GENEPOP ind format file conversion complete.')
	# End::DoGridOut_genepop()	

# ---------------------------------------------------------------------------------------------------	 
def DoOutput(SubpopIN,xgridpop,ygridpop,gen,ithmcrundir,loci,alleles,logfHndl,gridsample):
	'''
	DoOutput()
	Generate file of individuals
	'''	
	
	if gridsample == 'Initial' or gridsample == 'Middle':
		# Create file to write info to
		outputfile = open(ithmcrundir+'ind'+str(gen)+'.csv','w')
	elif gridsample == 'Sample':
		# Create file to write info to
		outputfile = open(ithmcrundir+'indSample'+str(gen)+'.csv','w')
	else:
		print('gridsample not correct format label. N or Sample.')
		sys.exit(-1)
	
	# Write out the titles - Add Titles from xypoints
	title = ['PatchID','XCOORD','YCOORD','ID','sex','age','size','mature','newmature','layeggs','capture','recapture','infection','CDist','Hindex','Species','ClassFile','SubPatchID']
	
	# Write out the title from xy points
	for i in xrange(len(title)):
		# Write out FID number
		outputfile.write(title[i]+',')
			
	# Write out the loci title info
	# Loop for loci length
	for i in xrange(loci-1):
		# Loop for allele length
		for j in xrange(alleles[i]):
			outputfile.write('L'+str(i)+'A'+str(j)+',')
	# To get a return character on the end of the title
	for i in xrange(alleles[loci-1]-1):
		outputfile.write('L'+str(loci-1)+'A'+str(i)+',')
	outputfile.write('L'+str(loci-1)+'A'+str(alleles[loci-1]-1)+'\n')
	
	# Then loop through each subpopulation
	for isub in xrange(len(SubpopIN)):
		Npop = len(SubpopIN[isub])
				
		# Loop through each N spot in this subpopultion
		for iK in xrange(Npop):			
			# Output the Subpopulation,X,Y - always assume in order
			outputfile.write(str(isub+1)+',')
			outputfile.write(str(xgridpop[isub])+',')
			outputfile.write(str(ygridpop[isub])+',')
			# Grab this individual and genes
			Ind = SubpopIN[isub][iK]
			#Ind_genes = literal_eval(Ind['genes']) - when 'str' dtype
			#Ind_genes = sum(Ind_genes,[])
			Ind_genes = np.asarray(Ind['genes'],'str')			
			# Output this individuals information
			outputfile.write(Ind['name']+',')
			outputfile.write(str(Ind['sex'])+',')
			outputfile.write(str(Ind['age'])+',')
			outputfile.write(str(Ind['size'])+',')
			outputfile.write(str(Ind['mature'])+',')
			outputfile.write(str(Ind['newmature'])+',')
			outputfile.write(str(Ind['layeggs'])+',')
			outputfile.write(str(Ind['capture'])+',')
			outputfile.write(str(Ind['recapture'])+',')
			outputfile.write(str(Ind['infection'])+',')
			if gridsample == 'Sample':
				outputfile.write(str(Ind['EmiCD'])+',')
			elif gridsample == 'Initial' or gridsample == 'Middle':
				outputfile.write(str(Ind['ImmiCD'])+',')
			outputfile.write(str(Ind['hindex'])+',')
			outputfile.write(str(Ind['species'])+',')
			outputfile.write(str(Ind['classfile'])+',')
			outputfile.write(str(Ind['popID'])+',')
			#outputfile.write(repr(Ind_genes).strip('[]')+'\n')
			for iall in xrange(sum(alleles)-1):
				outputfile.write(str(Ind_genes[iall])+',')
			outputfile.write(str(Ind_genes[sum(alleles)-1])+'\n')
																
	if gridsample == 'Initial' or 'Middle':
		# Logging message
		stringout = 'The file ind'+str(gen)+'.csv has been created'
		logMsg(logfHndl,stringout)
	elif gridsample == 'Sample':
		# Logging message
		stringout = 'The file indSample'+str(gen)+'.csv has been created'
		logMsg(logfHndl,stringout)
	
	# Close file
	outputfile.close()
	
	# End::DoOutput()

# ---------------------------------------------------------------------------------------------------	
def DoOut_AllTimeClass(K_track,ithmcrundir,logfHndl,N_Init_Age,N_back_age,PackingDeathsEmiAge,N_Emigration_age,AgeDeathsOUT,N_out_age,PackingDeathsImmAge,N_Immigration_age,AgeDeathsIN,AgeSizes_Mean,AgeSizes_Std,Capture_Back,Capture_Out,ClassSizes_Mean,ClassSizes_Std,N_Init_Class,size_mean,SizeDeathsOUT,SizeDeathsIN,N_beforePack_Immi_age):
	'''
	Create summary_classAllTime.csv file.
	'''
	
	# Create time array
	time = np.arange(0,len(K_track),1)
		
	# Get unique number of subpops
	nosubpops = len(K_track[0])-1
	
	# Create file to write info to
	outputfile = open(ithmcrundir+'summary_classAllTime.csv','w')
	
	# Write out the titles
	# Add Titles from xypoints
	outputtitle = ['Year','Ages','N_Initial_Age','AgeSize_Mean','AgeSize_Std','Size_Classes','N_Initial_Class','ClassSize_Mean','ClassSize_Std','N_GrowthBack','Capture_Back','PackingDeaths_Emigration','N_AfterEmigration','Deaths_AfterEmiMort_Age','Deaths_AfterEmiMort_Size','N_GrowthOut','Capture_Out','N_BeforePacking_AddAge0s','PackingDeaths_Immigration','N_Immigration','Deaths_AfterImmiMort_Age','Deaths_AfterImmiMort_Size']
	
	# Write out the title
	for i in xrange(len(outputtitle)-1):
		outputfile.write(outputtitle[i])
		outputfile.write(',')
	# To get return character on the end
	outputfile.write(str(outputtitle[len(outputtitle)-1]))				
	outputfile.write('\n')
	
	# Write to file
	for i in xrange(len(time)-1):		
		outputfile.write(str(time[i])+',')	
		for j in xrange(len(size_mean[0][0])): # Ages
			outputfile.write(str(j)+'|')
		outputfile.write(',')
		for j in xrange(len(N_Init_Age[i])):
			outputfile.write(str(N_Init_Age[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(AgeSizes_Mean[i])):
			outputfile.write(str(AgeSizes_Mean[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(AgeSizes_Std[i])):
			outputfile.write(str(AgeSizes_Std[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(size_mean[0][0])): #Sizes
			outputfile.write(str(size_mean[0][0][j])+'|')
		outputfile.write(',')
		for j in xrange(len(N_Init_Class[i])):
			outputfile.write(str(N_Init_Class[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(ClassSizes_Mean[i])):
			outputfile.write(str(ClassSizes_Mean[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(ClassSizes_Std[i])):
			outputfile.write(str(ClassSizes_Std[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(N_back_age[i])):
			outputfile.write(str(N_back_age[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(Capture_Back[i])):
			outputfile.write(str(Capture_Back[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(PackingDeathsEmiAge[i])):
			outputfile.write(str(PackingDeathsEmiAge[i][j])+'|')
		outputfile.write(',')				
		for j in xrange(len(N_Emigration_age[i])):
			outputfile.write(str(N_Emigration_age[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(len(AgeDeathsOUT[i])):
			outputfile.write(str(AgeDeathsOUT[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(SizeDeathsOUT[i])):
			outputfile.write(str(SizeDeathsOUT[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(len(N_out_age[i])):
			outputfile.write(str(N_out_age[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(Capture_Out[i])):
			outputfile.write(str(Capture_Out[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(N_beforePack_Immi_age[i])):
			outputfile.write(str(N_beforePack_Immi_age[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(len(PackingDeathsImmAge[i])):
			outputfile.write(str(PackingDeathsImmAge[i][j])+'|')
		outputfile.write(',')				
		for j in xrange(len(N_Immigration_age[i])):
			outputfile.write(str(N_Immigration_age[i][j])+'|')
		outputfile.write(',')	
		for j in xrange(len(AgeDeathsIN[i])):
			outputfile.write(str(AgeDeathsIN[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(SizeDeathsIN[i])):
			outputfile.write(str(SizeDeathsIN[i][j])+'|')
		outputfile.write('\n')						
	
	# For the tracking numbers that only have time counting
	i = len(time)-1
	outputfile.write(str(time[i])+',')
	for j in xrange(len(size_mean[0][0])):
		outputfile.write(str(j)+'|')
	outputfile.write(',')
	for j in xrange(len(N_Init_Age[i])):
		outputfile.write(str(N_Init_Age[i][j])+'|')
	outputfile.write(',')
	for j in xrange(len(AgeSizes_Mean[0])):
		outputfile.write(str(AgeSizes_Mean[i][j])+'|')
	outputfile.write(',')
	for j in xrange(len(AgeSizes_Std[0])):
		outputfile.write(str(AgeSizes_Std[i][j])+'|')
	outputfile.write(',')
	for j in xrange(len(size_mean[0][0])):
		outputfile.write(str(size_mean[0][0][j])+'|')
	outputfile.write(',')
	for j in xrange(len(N_Init_Class[i])):
		outputfile.write(str(N_Init_Class[i][j])+'|')
	outputfile.write(',')	
	for j in xrange(len(ClassSizes_Mean[0])):
		outputfile.write(str(ClassSizes_Mean[i][j])+'|')
	outputfile.write(',')
	for j in xrange(len(ClassSizes_Std[0])):
		outputfile.write(str(ClassSizes_Std[i][j])+'|')
	outputfile.write(',')
	for j in xrange(len(N_back_age[0])):
		outputfile.write('NA|')
	outputfile.write(',')
	for j in xrange(len(Capture_Back[0])):
		outputfile.write('NA|')
	outputfile.write(',')
	for j in xrange(len(PackingDeathsEmiAge[0])):
		outputfile.write('NA|')
	outputfile.write(',')		
	for j in xrange(len(N_Emigration_age[0])):
		outputfile.write('NA|')
	outputfile.write(',')		
	for j in xrange(len(AgeDeathsOUT[0])):
		outputfile.write('NA|')
	outputfile.write(',')
	for j in xrange(len(SizeDeathsOUT[0])):
		outputfile.write('NA|')
	outputfile.write(',')	
	for j in xrange(len(N_out_age[0])):
		outputfile.write('NA|')
	outputfile.write(',')
	for j in xrange(len(Capture_Out[0])):
		outputfile.write('NA|')
	outputfile.write(',')
	for j in xrange(len(N_beforePack_Immi_age[0])):
		outputfile.write('NA|')
	outputfile.write(',')		
	for j in xrange(len(PackingDeathsImmAge[0])):
		outputfile.write('NA|')
	outputfile.write(',')	
	for j in xrange(len(N_Immigration_age[0])):
		outputfile.write('NA|')
	outputfile.write(',')		
	for j in xrange(len(AgeDeathsIN[0])):
		outputfile.write('NA|')
	outputfile.write(',')	
	for j in xrange(len(SizeDeathsIN[0])):
		outputfile.write('NA|')
	outputfile.write('\n')
	
	
	# Logging message
	stringout = 'The file summary_classAllTime.csv has been created'
	logMsg(logfHndl,stringout)	
	
	# Close file
	outputfile.close()
	# End::DoOut_Class()
# ---------------------------------------------------------------------------------------------------	
def DoOut_AllTimePatch(K_track,ithmcrundir,logfHndl,N_Init,ToTFemales,ToTMales,BreedFemales,BreedMales,Female_BreedEvents,Births,EggDeaths,SelectionDeathsEmi,DisperseDeathsEmi,PackingDeathsEmi,N_Emigration,PopDeathsOUT,N_EmiMortality,SelectionDeathsImm,DisperseDeathsImm,PackingDeathsImm,N_Immigration,PopDeathsIN,N_ImmiMortality,Alleles,He,Ho,p1,p2,q1,q2,MateDistCD,MateDistCDstd,F_EmiDist,F_EmiDist_sd,M_EmiDist,M_EmiDist_sd,F_HomeDist,F_HomeDist_sd,M_HomeDist,M_HomeDist_sd,F_StrayDist,F_StrayDist_sd,M_StrayDist,M_StrayDist_sd,Infected,subpopemigration,subpopimmigration,MgSuccess,AdultNoMg,StrSuccess,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,MatureCount,ImmatureCount,Capture_Back,Capture_Out,N_beforePack_Immi_pop,SelectionDeathsImm_Age0s,F_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist,M_ZtrayDist_sd,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,ToTYYMales,BreedYYMales,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,RDispersers,IDispersers,BirthsYY):

	'''
	OutputPatch.csv file created
	'''
	
	# Get growth rate here
	tempPop = np.asarray(N_Init,dtype='float')[:,0]
	growthPop = tempPop[1:]/tempPop[0:(len(tempPop)-1)]
	
	# Create time array
	time = np.arange(0,len(K_track),1)
		
	# Get unique number of subpops -1 for total
	nosubpops = len(K_track[0])-1
	
	# Create file to write info to
	outputfile = open(ithmcrundir+'summary_popAllTime.csv','w')
	
	# Write out the titles
	# Add Titles from xypoints
	outputtitle = ['Year','K','GrowthRate','N_Initial','PopSizes_Mean','PopSizes_Std','N_Females','N_Males','N_YYMales','N_BreedFemales','N_BreedMales','N_BreedYYMales','MatureCount','ImmatureCount','EggLayEvents','Births','EggDeaths','MyyProgeny','Capture_Back','SelectionDeaths_Emigration','MoveDeaths_Emigration','PackingDeaths_Emigration','YYSelectionPackingDeaths_Emi','WildSelectionPackingDeaths_Emi','N_Emigration','Deaths_EmiMort','N_EmiMortality','Capture_Out','SelectionDeaths_Immigration','MoveDeaths_Immigration','N_beforePacking_AddAge0s','PackingDeaths_Immigration','YYSelectionPackingDeaths_Immi','WildSelectionPackingDeaths_Immi','SelectionDeaths_Age0s_Immigration','N_Immigration','Deaths_ImmiMort','N_ImmiMortality','Alleles','He','Ho','p1','p2','q1','q2','MateDist','MateDist_SD','Female_EmigrationDist','Female_EmigrationDist_SD','Male_EmigrationDist','Male_EmigrationDist_SD','Female_FromHomeDist','Female_FromHomeDist_SD','Male_FromHomeDist','Male_FromHomeDist_SD','Female_StrayerDist','Female_StrayerDist_SD','Male_StrayerDist','Male_StrayerDist_SD','Female_HomeAttemptStrayDist','Female_HomeAttemptStrayDist_SD','Male_HomeAttemptStrayDist','Male_HomeAttemptStrayDist_SD','Infected','Residors','Strayers_1','Strayers_2','Immigrators','ResidentDispersers','ImmigrantDispersers','AA_aa_Mates','AA_AA_Mates','aa_aa_Mates','AA_Aa_Mates','aa_Aa_Mates','Aa_Aa_Mates']
	
	# Write out the title
	for i in xrange(len(outputtitle)-1):
		outputfile.write(outputtitle[i])
		outputfile.write(',')
	# To get return character on the end
	outputfile.write(str(outputtitle[len(outputtitle)-1]))				
	outputfile.write('\n')
	
	# Write to file
	for i in xrange(len(time)-1):		
		outputfile.write(str(time[i])+',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(K_track[i+1][j])+'|') # Don't use first K
		outputfile.write(',')
		outputfile.write(str(growthPop[i])+',')		
		for j in xrange(nosubpops+1):
			outputfile.write(str(N_Init[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops):
			outputfile.write(str(PopSizes_Mean[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops):
			outputfile.write(str(PopSizes_Std[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(nosubpops+1):
			outputfile.write(str(ToTFemales[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(ToTMales[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(ToTYYMales[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(nosubpops+1):
			outputfile.write(str(BreedFemales[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(BreedMales[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(BreedYYMales[i][j])+'|')
		outputfile.write(',')		
		outputfile.write(str(MatureCount[i])+',')
		outputfile.write(str(ImmatureCount[i])+',')
		outputfile.write(str(Female_BreedEvents[i])+',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Births[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(len(EggDeaths[i])):
			outputfile.write(str(EggDeaths[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(BirthsYY[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(Capture_Back[i])):
			outputfile.write(str(Capture_Back[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(len(SelectionDeathsEmi[i])):
			outputfile.write(str(SelectionDeathsEmi[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(DisperseDeathsEmi[i])):
			outputfile.write(str(DisperseDeathsEmi[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(PackingDeathsEmi[i])):
			outputfile.write(str(PackingDeathsEmi[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(Track_YYSelectionPackDeathsEmi[i])):
			outputfile.write(str(Track_YYSelectionPackDeathsEmi[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(Track_WildSelectionPackDeathsEmi[i])):
			outputfile.write(str(Track_WildSelectionPackDeathsEmi[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(nosubpops+1):
			outputfile.write(str(N_Emigration[i][j])+'|')
		outputfile.write(',')			
		for j in xrange(len(PopDeathsOUT[i])):
			outputfile.write(str(PopDeathsOUT[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(N_EmiMortality[i][j])+'|')
		outputfile.write(',')
		
		for j in xrange(len(Capture_Out[i])):
			outputfile.write(str(Capture_Out[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(len(SelectionDeathsImm[i])):
			outputfile.write(str(SelectionDeathsImm[i][j])+'|')
		outputfile.write(',')				
		for j in xrange(len(DisperseDeathsImm[i])):
			outputfile.write(str(DisperseDeathsImm[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(N_beforePack_Immi_pop[i])):
			outputfile.write(str(N_beforePack_Immi_pop[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(PackingDeathsImm[i])):
			outputfile.write(str(PackingDeathsImm[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(Track_YYSelectionPackDeathsImmi[i])):
			outputfile.write(str(Track_YYSelectionPackDeathsImmi[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(Track_WildSelectionPackDeathsImmi[i])):
			outputfile.write(str(Track_WildSelectionPackDeathsImmi[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(SelectionDeathsImm_Age0s[i])):
			outputfile.write(str(SelectionDeathsImm_Age0s[i][j])+'|')
		outputfile.write(',')		
		for j in xrange(nosubpops+1):
			outputfile.write(str(N_Immigration[i][j])+'|')
		outputfile.write(',')				
		for j in xrange(len(PopDeathsIN[i])):
			outputfile.write(str(PopDeathsIN[i][j])+'|')
		outputfile.write(',')			
		for j in xrange(nosubpops+1):
			outputfile.write(str(N_ImmiMortality[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Alleles[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(He[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Ho[i][j])+'|')
		outputfile.write(',')
		outputfile.write(str(p1[i])+',')
		outputfile.write(str(p2[i])+',')
		outputfile.write(str(q1[i])+',')
		outputfile.write(str(q2[i])+',')
		outputfile.write(str(MateDistCD[i])+',')
		outputfile.write(str(MateDistCDstd[i])+',')
		outputfile.write(str(F_EmiDist[i])+',')
		outputfile.write(str(F_EmiDist_sd[i])+',')
		outputfile.write(str(M_EmiDist[i])+',')
		outputfile.write(str(M_EmiDist_sd[i])+',')
		outputfile.write(str(F_HomeDist[i])+',')
		outputfile.write(str(F_HomeDist_sd[i])+',')
		outputfile.write(str(M_HomeDist[i])+',')
		outputfile.write(str(M_HomeDist_sd[i])+',')
		outputfile.write(str(F_StrayDist[i])+',')
		outputfile.write(str(F_StrayDist_sd[i])+',')
		outputfile.write(str(M_StrayDist[i])+',')
		outputfile.write(str(M_StrayDist_sd[i])+',')
		outputfile.write(str(F_ZtrayDist[i])+',')
		outputfile.write(str(F_ZtrayDist_sd[i])+',')
		outputfile.write(str(M_ZtrayDist[i])+',')
		outputfile.write(str(M_ZtrayDist_sd[i])+',')		
		outputfile.write(str(Infected[i])+',')
		for j in xrange(nosubpops):
			outputfile.write(str(Residors[i+1][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops):
			outputfile.write(str(Strayers1[i+1][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops):
			outputfile.write(str(Strayers2[i+1][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops):
			outputfile.write(str(Immigrators[i+1][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops):
			outputfile.write(str(RDispersers[i+1][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops):
			outputfile.write(str(IDispersers[i+1][j])+'|')
		outputfile.write(',')
		outputfile.write(str(Track_AAaaMates[i])+',')
		outputfile.write(str(Track_AAAAMates[i])+',')
		outputfile.write(str(Track_aaaaMates[i])+',')
		outputfile.write(str(Track_AAAaMates[i])+',')
		outputfile.write(str(Track_aaAaMates[i])+',')
		outputfile.write(str(Track_AaAaMates[i])+',')		
		outputfile.write('\n')		
		
	# Logging message
	stringout = 'The file summary_popAllTime.csv has been created'
	logMsg(logfHndl,stringout)	
	
	# Close file
	outputfile.close()
	# End::DoOut_Patch

# ---------------------------------------------------------------------------------------------------	
def DoOut_Patch(K_track,ithmcrundir,logfHndl,N_Init,ToTFemales,ToTMales,BreedFemales,BreedMales,Births,EggDeaths,SelectionDeathsEmi,DisperseDeathsEmi,PackingDeathsEmi,N_Emigration,PopDeathsOUT,N_EmiMortality,SelectionDeathsImm,DisperseDeathsImm,PackingDeathsImm,N_Immigration,PopDeathsIN,N_ImmiMortality,Alleles,He,Ho,subpopemigration,subpopimmigration,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,nthfile,Capture_Back,Capture_Out,ToTYYMales,BreedYYMales,RDispersers,IDispersers,BirthsYY):
	'''
	summary_pop{year}_foldername.csv file created
	'''
	
	# Create patch array
	nopops = np.arange(1,len(K_track[0]),1)
	
	# Field headers 
	outputtitle = ['Subpopulation','K','N_Initial','PopSizes_Mean','PopSizes_Std','ToTFemales','ToTMales','ToTYYMales','BreedFemales','BreedMales','BreedYYMales','Births','EggDeaths','MyyProgeny','Capture_Back','SelectionDeaths_Emigration','MoveDeaths_Emigration','PackingDeaths_Emigration','N_AfterEmigration','Deaths_AfterEmiMort','N_AfterEmiMortality','Capture_Out','SelectionDeaths_Immigration','MoveDeaths_Immigration','PackingDeaths_Immigration','N_AfterImmigration','Deaths_AfterImmiMort','N_AfterImmiMortality','SubpopEmigration','SubpopImmigration','Residors','Strayers_1','Strayers_2','Immigrators','ResidentDispersers','ImmigrantDispersers','Alleles','He','Ho']
	
	# Begin loop through nthfile
	for itime in nthfile:
	
		# Create file to write info to
		outputfile = open(ithmcrundir+'summary_pop'+str(itime)+'.csv','w')
		
		# Write out the title
		for i in xrange(len(outputtitle)-1):
			outputfile.write(outputtitle[i])
			outputfile.write(',')
		# To get return character on the end
		outputfile.write(str(outputtitle[len(outputtitle)-1]))				
		outputfile.write('\n')
		
		# Write to file - loop through patches
		for ipop in xrange(len(nopops)):		
			outputfile.write(str(nopops[ipop])+',')
			outputfile.write(str(K_track[itime+1][ipop+1])+',') # Don't use first K
			outputfile.write(str(N_Init[itime][ipop+1])+',')
			outputfile.write(str(PopSizes_Mean[itime][ipop])+',')
			outputfile.write(str(PopSizes_Std[itime][ipop])+',')
			outputfile.write(str(ToTFemales[itime][ipop+1])+',')
			outputfile.write(str(ToTMales[itime][ipop+1])+',')
			outputfile.write(str(ToTYYMales[itime][ipop+1])+',')
			outputfile.write(str(BreedFemales[itime][ipop+1])+',')
			outputfile.write(str(BreedMales[itime][ipop+1])+',')
			outputfile.write(str(BreedYYMales[itime][ipop+1])+',')
			outputfile.write(str(Births[itime][ipop+1])+',')
			outputfile.write(str(EggDeaths[itime][ipop+1])+',')
			outputfile.write(str(BirthsYY[itime][ipop+1])+',')
			outputfile.write(str(Capture_Back[itime][ipop])+',')
			outputfile.write(str(SelectionDeathsEmi[itime][ipop+1])+',')
			outputfile.write(str(DisperseDeathsEmi[itime][ipop+1])+',')
			outputfile.write(str(PackingDeathsEmi[itime][ipop+1])+',')
			outputfile.write(str(N_Emigration[itime][ipop+1])+',')
			outputfile.write(str(PopDeathsOUT[itime][ipop+1])+',')
			outputfile.write(str(N_EmiMortality[itime][ipop+1])+',')
			outputfile.write(str(Capture_Out[itime][ipop])+',')
			outputfile.write(str(SelectionDeathsImm[itime][ipop+1])+',')
			outputfile.write(str(DisperseDeathsImm[itime][ipop+1])+',')
			outputfile.write(str(PackingDeathsImm[itime][ipop+1])+',')
			outputfile.write(str(N_Immigration[itime][ipop+1])+',')
			outputfile.write(str(PopDeathsIN[itime][ipop+1])+',')
			outputfile.write(str(N_ImmiMortality[itime][ipop+1])+',')				
			outputfile.write(str(subpopemigration[itime][ipop])+',')
			outputfile.write(str(subpopimmigration[itime][ipop])+',')
			outputfile.write(str(Residors[itime][ipop])+',')
			outputfile.write(str(Strayers1[itime][ipop])+',')
			outputfile.write(str(Strayers2[itime][ipop])+',')
			outputfile.write(str(Immigrators[itime][ipop])+',')
			outputfile.write(str(RDispersers[itime][ipop])+',')
			outputfile.write(str(IDispersers[itime][ipop])+',')
			outputfile.write(str(Alleles[itime][ipop+1])+',')
			outputfile.write(str(He[itime][ipop+1])+',')
			outputfile.write(str(Ho[itime][ipop+1])+'\n')
			
		# For Total numbers
		outputfile.write('Totals,')
		outputfile.write(str(K_track[itime+1][0])+',')
		outputfile.write(str(N_Init[itime][0])+',')
		outputfile.write('NA,')
		outputfile.write('NA,')
		outputfile.write(str(ToTFemales[itime][0])+',')
		outputfile.write(str(ToTMales[itime][0])+',')
		outputfile.write(str(BreedFemales[itime][0])+',')
		outputfile.write(str(BreedMales[itime][0])+',')
		outputfile.write(str(Births[itime][0])+',')
		outputfile.write(str(EggDeaths[itime][0])+',')
		outputfile.write(str(BirthsYY[itime][0])+',')
		outputfile.write(str(sum(Capture_Back[itime]))+',')		
		outputfile.write(str(SelectionDeathsEmi[itime][0])+',')		
		outputfile.write(str(DisperseDeathsEmi[itime][0])+',')	
		outputfile.write(str(PackingDeathsEmi[itime][0])+',')
		outputfile.write(str(N_Emigration[itime][0])+',')
		outputfile.write(str(PopDeathsOUT[itime][0])+',')
		outputfile.write(str(N_EmiMortality[itime][0])+',')
		outputfile.write(str(sum(Capture_Out[itime]))+',')
		outputfile.write(str(SelectionDeathsImm[itime][0])+',')
		outputfile.write(str(DisperseDeathsImm[itime][0])+',')
		outputfile.write(str(PackingDeathsImm[itime][0])+',')
		outputfile.write(str(N_Immigration[itime][0])+',')
		outputfile.write(str(PopDeathsIN[itime][0])+',')
		outputfile.write(str(N_ImmiMortality[itime][0])+',')
		outputfile.write(str(sum(subpopemigration[itime]))+',')
		outputfile.write(str(sum(subpopimmigration[itime]))+',')
		outputfile.write(str(sum(Residors[itime]))+',')
		outputfile.write(str(sum(Strayers1[itime]))+',')
		outputfile.write(str(sum(Strayers2[itime]))+',')
		outputfile.write(str(sum(Immigrators[itime]))+',')
		outputfile.write(str(sum(RDispersers[itime]))+',')
		outputfile.write(str(sum(IDispersers[itime]))+',')
		outputfile.write(str(Alleles[itime][0])+',')
		outputfile.write(str(He[itime][0])+',')
		outputfile.write(str(Ho[itime][0])+'\n')
		
		# Logging message
		stringout = 'The files summary_pop.csv have been created'
		logMsg(logfHndl,stringout)	
		
		# Close file
		outputfile.close()

	# End::DoOut_Patch()
	
# ---------------------------------------------------------------------------------------------------
def DoOut_Class(ithmcrundir,logfHndl,N_Init_Age,N_back_age,PackingDeathsEmiAge,N_Emigration_age,AgeDeathsOUT,N_out_age,PackingDeathsImmAge,N_Immigration_age,AgeDeathsIN,AgeSizes_Mean,AgeSizes_Std,nthfile,Capture_Back,Capture_Out,size_mean,sizeans,ClassSizes_Mean,ClassSizes_Std,N_Init_Class,SizeDeathsOUT,SizeDeathsIN):

	'''
	Create summary_class{year}.csv file.
	'''
	
	# Create class array
	noclass = np.arange(0,len(size_mean[0][0]),1)
	
	# Get class values
	if sizeans == 'Y':
		classvals = size_mean[0][0]
	else:
		classvals = range(0,len(size_mean[0][0]),1)
	
	# Write out the titles
	# Add Titles from xypoints
	outputtitle = ['Age','Class','N_Initial_Age','AgeSize_Mean','AgeSize_Std','N_Initial_Class','ClassSize_Mean','ClassSize_Std','N_GrowthBack_Class','Capture_Back','PackingDeaths_Emigration_Class','N_AfterEmigration_Class','Deaths_AfterEmiMort_Age','Deaths_AfterEmiMort_Size','N_GrowthOut_Class','Capture_Out','PackingDeaths_Immigration_Class','N_Immigration_Class','Deaths_AfterImmiMort_Age','Deaths_AfterImmiMort_Size']
			
	# Begin loop through nthfile
	for itime in nthfile:
	
		# Create file to write info to
		outputfile = open(ithmcrundir+'summary_class'+str(itime)+'.csv','w')			
		
		# Write out the title
		for i in xrange(len(outputtitle)-1):
			outputfile.write(outputtitle[i])
			outputfile.write(',')
		# To get return character on the end
		outputfile.write(str(outputtitle[len(outputtitle)-1]))				
		outputfile.write('\n')
		
		# Write to file - loop through patches
		for iage in xrange(len(noclass)):
			outputfile.write(str(iage)+',')
			outputfile.write(str(classvals[iage])+',')
			outputfile.write(str(N_Init_Age[itime][iage])+',')
			outputfile.write(str(AgeSizes_Mean[itime][iage])+',')
			outputfile.write(str(AgeSizes_Std[itime][iage])+',')
			outputfile.write(str(N_Init_Class[itime][iage])+',')
			outputfile.write(str(ClassSizes_Mean[itime][iage])+',')
			outputfile.write(str(ClassSizes_Std[itime][iage])+',')
			outputfile.write(str(N_back_age[itime][iage])+',')
			outputfile.write(str(Capture_Back[itime][iage])+',')
			outputfile.write(str(PackingDeathsEmiAge[itime][iage])+',')
			outputfile.write(str(N_Emigration_age[itime][iage])+',')
			outputfile.write(str(AgeDeathsOUT[itime][iage])+',')
			outputfile.write(str(SizeDeathsOUT[itime][iage])+',')
			outputfile.write(str(N_out_age[itime][iage])+',')
			outputfile.write(str(Capture_Out[itime][iage])+',')
			outputfile.write(str(PackingDeathsImmAge[itime][iage])+',')
			outputfile.write(str(N_Immigration_age[itime][iage])+',')
			outputfile.write(str(AgeDeathsIN[itime][iage])+',')
			outputfile.write(str(SizeDeathsIN[itime][iage])+'\n')
			
		# For totals
		outputfile.write('Totals,')
		outputfile.write(str(len(noclass))+',')
		outputfile.write(str(sum(N_Init_Age[itime]))+',')
		outputfile.write('NA,')
		outputfile.write('NA,')
		outputfile.write(str(sum(N_Init_Class[itime]))+',')
		outputfile.write('NA,')
		outputfile.write('NA,')
		if N_back_age[itime][0] != 'NA':
			outputfile.write(str(sum(N_back_age[itime]))+',')
		else:
			outputfile.write('NA,')
		outputfile.write(str(sum(Capture_Back[itime]))+',')
		outputfile.write(str(sum(PackingDeathsEmiAge[itime]))+',')
		outputfile.write(str(sum(N_Emigration_age[itime]))+',')
		outputfile.write(str(sum(AgeDeathsOUT[itime]))+',')
		outputfile.write(str(sum(SizeDeathsOUT[itime]))+',')
		if N_out_age[itime][0] != 'NA':
			outputfile.write(str(sum(N_out_age[itime]))+',')
		else:
			outputfile.write('NA,')
		outputfile.write(str(sum(Capture_Out[itime]))+',')
		outputfile.write(str(sum(PackingDeathsImmAge[itime]))+',')
		outputfile.write(str(sum(N_Immigration_age[itime]))+',')
		outputfile.write(str(sum(AgeDeathsIN[itime]))+',')
		outputfile.write(str(sum(SizeDeathsIN[itime]))+'\n')		
		
		# Logging message
		stringout = 'The files summary_class{year}.csv has been created'
		logMsg(logfHndl,stringout)	
		
		# Close file
		outputfile.close()
	
	# End::DoOut_Class()
	
# ---------------------------------------------------------------------------------------------------	 
def DoPostProcess(ithmcrundir,loci,alleles,looptime,ToTFemales,ToTMales,\
BreedFemales,BreedMales,Births,PopDeathsIN,PopDeathsOUT,Alleles,He,Ho,\
MateDistCD,MateDistCDstd,nthfile,logfHndl,p1,p2,q1,q2,\
subpopemigration,subpopimmigration,FAvgMate,MAvgMate,FSDMate,\
MSDMate,SelectionDeathsEmi,SelectionDeathsImm,\
DisperseDeathsEmi,DisperseDeathsImm,Female_BreedEvents,\
gridformat,MgSuccess,AdultNoMg,\
StrSuccess,EggDeaths,K_track,N_Init,N_Emigration,N_EmiMortality,N_Immigration,N_ImmiMortality,Infected,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,AgeSizes_Mean,AgeSizes_Std,PackingDeathsEmi,PackingDeathsImm,N_Init_Age,N_Emigration_age,N_Immigration_age,AgeDeathsOUT,AgeDeathsIN,PackingDeathsEmiAge,PackingDeathsImmAge,MatureCount,ImmatureCount,N_back_age,N_out_age,outputans,gen,CaptureCount_Back,CaptureCount_ClassBack,CaptureCount_Out,CaptureCount_ClassOut,size_mean,sizeans,ClassSizes_Mean,ClassSizes_Std,N_Init_Class,SizeDeathsOUT,SizeDeathsIN,N_beforePack_Immi_pop,N_beforePack_Immi_age,SelectionDeathsImm_Age0s,F_StrayDist,M_StrayDist,F_StrayDist_sd,M_StrayDist_sd,F_ZtrayDist,M_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist_sd,F_HomeDist,M_HomeDist,F_HomeDist_sd,M_HomeDist_sd,F_EmiDist,M_EmiDist,F_EmiDist_sd,M_EmiDist_sd,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,ToTYYMales,BreedYYMales,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,RDispersers,IDispersers,BirthsYY):
	'''
	DoPostProcess()
	Create Distance Matrices - Geographic, Genetic, and Cost
	and output.csv file.
	'''	
	
	# If extinct before nthfile end check here
	if gen < max(nthfile):
		nthfile = np.where(np.asarray(nthfile) <= gen)[0]
	
	# ------------------------
	# Grid format options
	# ------------------------
	# General format
	if gridformat == 'general':		
		DoGridOut_general(loci,alleles,ithmcrundir,logfHndl,K_track)
		
	# GENALEX format
	elif gridformat == 'genalex':
		DoGridOut_genalex(loci,alleles,ithmcrundir,logfHndl,K_track)
		
	# STRUCTURE format
	elif gridformat == 'structure':
		DoGridOut_structure(loci,alleles,ithmcrundir,logfHndl,K_track)
	
	# GENEPOP format
	elif gridformat == 'genepop':
		DoGridOut_genepop(loci,alleles,ithmcrundir,logfHndl,K_track)
			
	# -------------------------------------------------
	# output for each gen and patch vars
	# -------------------------------------------------
	DoOut_AllTimePatch(K_track,ithmcrundir,logfHndl,N_Init,ToTFemales,ToTMales,BreedFemales,BreedMales,Female_BreedEvents,Births,EggDeaths,SelectionDeathsEmi,DisperseDeathsEmi,PackingDeathsEmi,N_Emigration,PopDeathsOUT,N_EmiMortality,SelectionDeathsImm,DisperseDeathsImm,PackingDeathsImm,N_Immigration,PopDeathsIN,N_ImmiMortality,Alleles,He,Ho,p1,p2,q1,q2,MateDistCD,MateDistCDstd,F_EmiDist,F_EmiDist_sd,M_EmiDist,M_EmiDist_sd,F_HomeDist,F_HomeDist_sd,M_HomeDist,M_HomeDist_sd,F_StrayDist,F_StrayDist_sd,M_StrayDist,M_StrayDist_sd,Infected,subpopemigration,subpopimmigration,MgSuccess,AdultNoMg,StrSuccess,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,MatureCount,ImmatureCount,CaptureCount_Back,CaptureCount_Out,N_beforePack_Immi_pop,SelectionDeathsImm_Age0s,F_ZtrayDist,F_ZtrayDist_sd,M_ZtrayDist,M_ZtrayDist_sd,Track_AAaaMates,Track_AAAAMates,Track_aaaaMates,Track_AAAaMates,Track_aaAaMates,Track_AaAaMates,ToTYYMales,BreedYYMales,Track_YYSelectionPackDeathsEmi,Track_WildSelectionPackDeathsEmi,Track_YYSelectionPackDeathsImmi,Track_WildSelectionPackDeathsImmi,RDispersers,IDispersers,BirthsYY)	
	
	# -------------------------------------------------
	# output for each generation and class vars
	# -------------------------------------------------
	DoOut_AllTimeClass(K_track,ithmcrundir,logfHndl,N_Init_Age,N_back_age,PackingDeathsEmiAge,N_Emigration_age,AgeDeathsOUT,N_out_age,PackingDeathsImmAge,N_Immigration_age,AgeDeathsIN,AgeSizes_Mean,AgeSizes_Std,CaptureCount_ClassBack,CaptureCount_ClassOut,ClassSizes_Mean,ClassSizes_Std,N_Init_Class,size_mean,SizeDeathsOUT,SizeDeathsIN,N_beforePack_Immi_age)

	if outputans == 'Y':
			
		# ---------------------------------------------
		# Output Patch values at given nthfile
		# ---------------------------------------------
		
		DoOut_Patch(K_track,ithmcrundir,logfHndl,N_Init,ToTFemales,ToTMales,BreedFemales,BreedMales,Births,EggDeaths,SelectionDeathsEmi,DisperseDeathsEmi,PackingDeathsEmi,N_Emigration,PopDeathsOUT,N_EmiMortality,SelectionDeathsImm,DisperseDeathsImm,PackingDeathsImm,N_Immigration,PopDeathsIN,N_ImmiMortality,Alleles,He,Ho,subpopemigration,subpopimmigration,Residors,Strayers1,Strayers2,Immigrators,PopSizes_Mean,PopSizes_Std,nthfile,CaptureCount_Back,CaptureCount_Out,ToTYYMales,BreedYYMales,RDispersers,IDispersers,Births)
		
		# ---------------------------------------------
		# Output Class values at given nthfile
		# ---------------------------------------------
		
		DoOut_Class(ithmcrundir,logfHndl,N_Init_Age,N_back_age,PackingDeathsEmiAge,N_Emigration_age,AgeDeathsOUT,N_out_age,PackingDeathsImmAge,N_Immigration_age,AgeDeathsIN,AgeSizes_Mean,AgeSizes_Std,nthfile,CaptureCount_ClassBack,CaptureCount_ClassOut,size_mean,sizeans,ClassSizes_Mean,ClassSizes_Std,N_Init_Class,SizeDeathsOUT,SizeDeathsIN)
		
	#End:DoPostProcess()