# -----------------------------------------------------------------------------
# genesvalidate.py
# 2014 10 09: Erin Landguth
# This script grabs output.csv files, splits bars, and grabs diagnostics for checking CDPOP/CDFISH/CDmetaPOP
# v0 - Initial script for He/Ho check
# v1 - Allele plots
# ----------------------------------------------------------------------------- 

# Load modules
import os,pdb,pandas
from pylab import *	
import scipy as sp			

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."

# ---------
# User info
# ---------
# Directory locations of folders
dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT373_testPackingInit/testPacking50gens/age0/"
plottitle = 'Age 0'
savename = 'Age0_'
batchno = 0 # batch number to run
dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT373_testPackingInit/testPacking50gens/age1/"
plottitle = 'Age 1'
savename = 'Age1_'
batchno = 1 # batch number to run
dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT373_testPackingInit/testPacking50gens/age2/"
plottitle = 'Age 2'
savename = 'Age2_'
batchno = 2 # batch number to run
dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT373_testPackingInit/testPacking50gens/allage/"
plottitle = 'Age 0-5'
savename = 'Age0to5_'
batchno = 3 # batch number to run

outdir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT373_testPackingInit/testPacking50gens/"

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 50 # Number of years 
startgenes = 0
nthfile = range(startgenes,gen,1)
maxA = 100 # loci * alleles
mcno = 2 # Number of MCs

# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# ---------------------------------
# Storage variables
# ---------------------------------

alleles = []
He = []
Ho = []

# Loop through batches
for ibatch in xrange(1):

	alleles.append([])
	He.append([])
	Ho.append([])
	
	# Loop through MCs
	for imc in xrange(mcno):
	
		# Open output.csv file in folder
		inputfile = open(dir+'batchrun'+str(batchno)+'mcrun'+str(imc)+'/output.csv')
		
		# Read lines from the file
		lines = inputfile.readlines()

		#Close the file
		inputfile.close()

		# Create an empty matrix to append to
		values = []

		# Split up each line in file and append to empty matrix for generation specified
		for i in xrange(len(lines)):
			thisline = lines[i].split(',')
			values.append(thisline)

		# Delete lines
		del(lines)
		
		alleles[ibatch].append([])
		He[ibatch].append([])
		Ho[ibatch].append([])
		
		# Then Loop through generations/time
		for iout in xrange(gen):
		
			alleles[ibatch][imc].append([])
			He[ibatch][imc].append([])
			Ho[ibatch][imc].append([])
			
			# Grab all patch values - patch values with total
			for j in xrange(len(values[1+iout][30].split('|'))-1):
				alleles[ibatch][imc][iout].append(int(values[1+iout][29].split('|')[j]))
				He[ibatch][imc][iout].append(float(values[1+iout][30].split('|')[j]))
				Ho[ibatch][imc][iout].append(float(values[1+iout][31].split('|')[j]))
			
			
alleles = np.asarray(alleles)
He = np.asarray(He)
Ho = np.asarray(Ho)
allD = (alleles - 1) / (maxA - 1.)

# --------------------------------------------
# Get mean over Monte Carlosfor each batch run
# --------------------------------------------
all_m = np.nansum(alleles[:][:],axis=1)/mcno
all_sd = np.std(alleles[:][:],axis=1)	
error = qnorm*all_sd/(mcno)
all_l = all_m-error
all_r = 	all_m+error
all_min = np.min(alleles[:][:],axis=1)
all_max = np.max(alleles[:][:],axis=1)

He_m = np.nansum(He[:][:],axis=1)/mcno
He_sd = np.std(He[:][:],axis=1)	
error = qnorm*He_sd/(mcno)
He_l = He_m-error
He_r = 	He_m+error
He_min = np.min(He[:][:],axis=1)
He_max = np.max(He[:][:],axis=1)

Ho_m = np.nansum(Ho[:][:],axis=1)/mcno
Ho_sd = np.std(Ho[:][:],axis=1)	
error = qnorm*Ho_sd/(mcno)
Ho_l = Ho_m-error
Ho_r = 	Ho_m+error
Ho_min = np.min(Ho[:][:],axis=1)
Ho_max = np.max(Ho[:][:],axis=1)

allD_m = np.nansum(allD[:][:],axis=1)/mcno
allD_sd = np.std(allD[:][:],axis=1)	
error = qnorm*allD_sd/(mcno)
allD_l = allD_m-error
allD_r = 	allD_m+error

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------

# Calculate Non-fisherian rate of loss + noselfing+1
rateofloss_nonfisherian=[]
equ1_Tot=[]
equ1_1 = []
N = 1000
for i in range(len(nthfile)):
	rateofloss_nonfisherian.append(He_m[0][0][0]*exp((-1./(2*N+1))*nthfile[i]))
	equ1_Tot.append(He_m[0][0][0]*((1-(1./(2*N+1)))**nthfile[i]))

# Plot Total He,Ho vs. time 
figure()
plot(nthfile,He_m[0][:,0][nthfile],'k',label='Expected',ms=10)
#plot(time,HeTot_Left,'-.r',ms=10)
#plot(time,HeTot_Right,'-.k',ms=10)
plot(nthfile,Ho_m[0][:,0][nthfile],':k',label='Observed',ms=10)
#plot(time,HoTot_Left,'-.k',ms=10)
#plot(time,HoTot_Right,'-.k',ms=10)
#plot(nthfile,equ1_Tot,'-k',ms=10,label='Equation 1')
xlabel('generations',fontsize=18)
ylabel('Heterozygosity',fontsize=18)
axis([startgenes,gen,0.0,1.0])
title(plottitle)
legend(loc=3)
# Updating fontsize on axes
fontsize=16
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
savefig(dir+savename+'HeHo.png',dpi=savedpi)
show()