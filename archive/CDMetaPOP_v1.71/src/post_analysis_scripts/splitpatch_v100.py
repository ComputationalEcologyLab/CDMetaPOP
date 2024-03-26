# -----------------------------------------------------------------------------
# splitpatch.py
# 2015 05 25: Erin Landguth
# This script grabs summary.csv files, splits bars for given parameter, writes 
# seperate file.
# v0.99 - Initial script for version 0.99 files.
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
dir = 'D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v1.00_20151130/All_MovementDistances/'
savename = "WCT_100Max_NInit_"
# batch number - batchrun(batchno)mcrun{mcrun}
batchno = 2

dir = 'D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v1.00_20151130/All_StrayRates/'
savename = "NInit_WCT_50max_str0.001_"
savename = "NInit_WCT_50max_str0.01_"
savename = "NInit_WCT_50max_str0.05_"
batchno = 2

dir = 'D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v1.00_20151130/All_Landscapes/'
savename = "NInit_WCT_50max_str0.05_Riverine_"
#savename = "NInit_WCT_50max_str0.05_ExBarr"
#savename = "NInit_WCT_50max_str0.05_FutBarr"
#savename = "NInit_WCT_50max_str0.05_ChangeBarr"
# batch number - batchrun(batchno)mcrun{mcrun}
batchno = 0

# Direction location of patchvars.csv to get X,Y values
xydir = 'D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v1.00_20151130/PatchVars_WCT1384.csv'

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 125 # Number of years 
nthfile = range(0,gen,1)
mcno = 2 # Number of MCs

# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# ---------------------------------
# Read in XY file and store XY
# ---------------------------------
Patch = []
X = []
Y = []
# Open output.csv file in folder
inputfile = open(xydir)

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

# Grab values from file
for i in xrange(len(values)-1):
	Patch.append(int(values[1+i][0]))
	X.append(float(values[1+i][1]))
	Y.append(float(values[1+i][2]))	

# ---------------------------------
# Read in summary file information
# ---------------------------------

N_Init = [] 
R = []
S = []
I = []
He = []

# Loop through MCs
for imc in xrange(mcno):

	N_Init.append([])
	R.append([])
	S.append([])
	I.append([])

	# Open output.csv file in folder
	inputfile = open(dir+'batchrun'+str(batchno)+'mcrun'+str(imc)+'/summary_popAllTime.csv')
	
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
	
	# Then Loop through generations/time
	for iout in xrange(gen):
	
		N_Init[imc].append([])
		R[imc].append([])
		S[imc].append([])
		I[imc].append([])
		
				
		# Grab all patch values - patch values with total 
		for j in xrange(1,len(values[1+iout][3].split('|'))-1):
			N_Init[imc][iout].append(int(values[1+iout][3].split('|')[j]))
			#R[imc][iout].append(int(values[1+iout][57].split('|')[j]))
			#S[imc][iout].append(int(values[1+iout][58].split('|')[j]))
			#I[imc][iout].append(int(values[1+iout][60].split('|')[j]))
		
N_Init = np.asarray(N_Init)
#R = np.asarray(R)
#S = np.asarray(S)
#I = np.asarray(I)

# --------------------------------------------
# Get mean over Monte Carlos
# --------------------------------------------
N_Init_m = np.nansum(N_Init,axis=0)/mcno
N_Init_sd = np.std(N_Init,axis=0)	
error = qnorm*N_Init_sd/(mcno)
N_Init_l = N_Init_m-error
N_Init_r = 	N_Init_m+error
N_Init_min = np.min(N_Init,axis=0)
N_Init_max = np.max(N_Init,axis=0)
'''
R_m = np.nansum(R,axis=0)/mcno
R_sd = np.std(R,axis=0)	
error = qnorm*R_sd/(mcno)
R_l = R_m-error
R_r = 	R_m+error
R_min = np.min(R,axis=0)
R_max = np.max(R,axis=0)

S_m = np.nansum(S,axis=0)/mcno
S_sd = np.std(S,axis=0)	
error = qnorm*S_sd/(mcno)
S_l = S_m-error
S_r = 	S_m+error
S_min = np.min(S,axis=0)
S_max = np.max(S,axis=0)

I_m = np.nansum(I,axis=0)/mcno
I_sd = np.std(I,axis=0)	
error = qnorm*I_sd/(mcno)
I_l = I_m-error
I_r = 	I_m+error
I_min = np.min(I,axis=0)
I_max = np.max(I,axis=0)
'''
# -----------------------------
# Write to file 1
# -----------------------------
# And write out title at end of file
outputfile = open(dir+savename+'1.csv','w')

# Create year string
timeloop = range(2000,2000+gen,1)

# Write out title
outputfile.write('Patch,X,Y,')
for iout in xrange(gen-1):
	outputfile.write('T'+str(timeloop[iout])+',')
outputfile.write('T'+str(timeloop[iout+1])+'\n')

# Write out values
for ipatch in xrange(len(N_Init_m[0])):
	outputfile.write(str(Patch[ipatch])+','+str(X[ipatch])+','+str(Y[ipatch])+',')
	for itime in xrange(len(N_Init_m)-1):
		outputfile.write(str(N_Init_m[itime][ipatch])+',')
	outputfile.write(str(N_Init_m[itime+1][ipatch])+'\n')

#Close the file
outputfile.close()	

# -----------------------------
# Write to file 2
# -----------------------------

# And write out title at end of file
outputfile = open(dir+savename+'2.csv','w')

# Create year string
timeloop = range(1995,1995+gen,1)

# Write out title
outputfile.write('Patch,X,Y,Value,Time\n')
for itime in xrange(len(timeloop)):
	timeout = '1/1/'+str(timeloop[itime])
	for ipatch in xrange(len(N_Init_m[0])):
		outputfile.write(str(Patch[ipatch])+','+str(X[ipatch])+','+str(Y[ipatch])+','+str(N_Init_m[itime][ipatch])+','+timeout+'\n')

#Close the file
outputfile.close()


