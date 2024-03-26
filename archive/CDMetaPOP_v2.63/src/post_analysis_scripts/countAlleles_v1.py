# -------------------------------------------------
# CountL0.py
# Erin Landguth
# June 8, 2015
# Description: This script reads in ind files (or indSample files)
# and counts the number of 2, 1, 0 at first locus (this is for 
# the selection module and identification of different strains)
# v0 - Initial script
# v1 - Add in N by L0 counts and plot.
# ---------------------------------------------------

# Import statements
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."				
import pdb,os,copy		
from random import *
from pylab import *	

# Functions needed
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
#-----------------
# User Input
#-----------------
# Directory location of output.csv
#dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1378_v09912_Hybrid_20150625/Selection_0/"
# Outputname 
#outname = "Patch1378_0Selection_"

#dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1378_v09912_Hybrid_20150625/Selection_25/"
# Outputname 
#outname = "Patch1378_25Selection_"

#dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1378_v09912_Hybrid_20150625/Selection_50/"
# Outputname 
#outname = "Patch1378_50Selection_"

#dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1378_v09912_Hybrid_20150625/Selection_100/"
# Outputname 
#outname = "Patch1378_100Selection_"

dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/ExtBarrier_Selection_0/"
outname = "Patch1384_NoSelection_ExtPartBarr_"

dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/ExtBarrier_Selection_0_ReducedRT/"
outname = "Patch1384_NoSelection_ExtPartBarr_ReducedRT_"

dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/ExtBarrier_Selection_100_ReducedRT/"
outname = "Patch1384_100Selection_ExtPartBarr_ReducedRT_"

dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/ExtBarrier_Selection_100/"
outname = "Patch1384_100Selection_ExtPartBarr_"

dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/Riverine_Selection_0/"
outname = "Patch1384_NoSelection_Riverine_"

dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/Riverine_Selection_100/"
outname = "Patch1384_100Selection_Riverine_"

dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/Riverine_Selection_0_ReducedRT/"
outname = "Patch1384_NoSelection_Riverine_ReducedRT_"

dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/Riverine_Selection_100_ReducedRT/"
outname = "Patch1384_100Selection_Riverine_ReducedRT_"


# batch number - batchrun(batchno)mcrun{mcrun}
batchno = 0

# Location and name for xy, costdistance to subpop location
xyfilename = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v100_HybridBarriers_20150816/PatchWCT1384_Selection100.csv"

# Number of monte carlo runs - used to average
mcruns = 1

# Generation to extract summary from
gen = range(9,110,1)
gen = [0,13,50,99]
gen = [0,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,50,60,80,99]
gen = range(0,130,5)
#gen = [0,13,99]
# Number of patches
n = 1384

# Which grid files to do ind{}.csv or indSample{}.csv (ind or indSample)
gridformat = 'ind'

# Calculate each statistic - Error (use qnorm(0.975)=1.959964 - couldn't find equivalent qnorm) 
qnorm = 1.959964

# ------------------
# End User Input
# ------------------

# List folders in this dir
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

#-----------------------------------------------
# Initial File IO - read in cost distance values
#-----------------------------------------------	
# Open file to extract XY values
xyinputfile = open(xyfilename,'r')

# Read lines from the file
lines = xyinputfile.readlines()

#Close the file
xyinputfile.close()

# Create an empty matrix to append to
xyvalues = []

# Split up each line in file and append to empty matrix for generation specified
for i in xrange(len(lines)):
	thisline = lines[i].split(',')
	xyvalues.append(thisline)

# Get x and y values
X = []
Y = []
K = []
Patch = []
for i in xrange(len(xyvalues)-1):
	X.append(float(xyvalues[i+1][1]))
	Y.append(float(xyvalues[i+1][2]))
	K.append(int(xyvalues[i+1][3]))
	Patch.append(int(xyvalues[i+1][0]))

# Delete lines
del(lines)

# All information
L0_2 = []
L0_1 = []
L0_0 = []
N = []
N_L0_2 = []
N_L0_1 = []
N_L0_0 = []

# --------------------
# Begin gen loop
# --------------------
for igen in xrange(len(gen)):
	
	# ----------------------------	
	#  Preliminary vector storage
	# ----------------------------
	# All information
	L0_2.append([])
	L0_1.append([])
	L0_0.append([])
	N.append([])
	N_L0_2.append([])
	N_L0_1.append([])
	N_L0_0.append([])
	
	# -----------------------------------	
	# Read in and store metrics
	# -----------------------------------
	# Loop over folders
	for imc in xrange(mcruns): 
		
		# Add a batch spot to the vectors
		L0_2[igen].append([])
		L0_1[igen].append([])
		L0_0[igen].append([])
		N[igen].append([])
		
		# --------------------------
		#  Read in information 
		# --------------------------		
		# Open file to extract number of migrants
		xyinputfile = open(dir+"batchrun"+str(batchno)+"mcrun"+str(imc)+"/"+gridformat+str(gen[igen])+".csv",'r')

		# Read lines from the file
		lines = xyinputfile.readlines()

		#Close the file
		xyinputfile.close()

		# Create an empty matrix to append to
		x = []

		# Split up each line in file and append to empty matrix for generation specified
		for i in xrange(len(lines)):
			thisline = lines[i].split(',')
			x.append(thisline)
		
		N_gen = []
		L0 = []
		# Extract values from this grid	
		for i in xrange(len(x)-1):
			N_gen.append(int(x[i+1][0]))
			L0.append(int(x[i+1][14]))			
		# Get unique L0
		L0 = np.asarray(L0)
		uniL0 = count_unique(L0)
		
		# Get the number in each patch
		for i in range(1,n+1):
			indexHere = np.where(np.asarray(N_gen) == i)[0]
			N[igen][imc].append(len(indexHere))
			# Index into L0
			countL0 = count_unique(L0[indexHere])
			L0_0[igen][imc].append(len(np.where(L0[indexHere] == 0)[0]))
			L0_1[igen][imc].append(len(np.where(L0[indexHere] == 1)[0]))
			L0_2[igen][imc].append(len(np.where(L0[indexHere] == 2)[0]))
		
		N_L0_2[igen].append(sum(L0_2[igen][imc]))
		N_L0_1[igen].append(sum(L0_1[igen][imc]))
		N_L0_0[igen].append(sum(L0_0[igen][imc]))
		
# --------------------------------
# Summary MC
# --------------------------------
N = np.asarray(N)
L0_0 = np.asarray(L0_0)
L0_1 = np.asarray(L0_1)
L0_2 = np.asarray(L0_2)
N_L0_0 = np.asarray(N_L0_0)
N_L0_1 = np.asarray(N_L0_1)
N_L0_2 = np.asarray(N_L0_2)

N_m = np.nansum(N,axis = 1)/mcruns
L0_0 = np.nansum(L0_0,axis = 1)/mcruns
L0_1 = np.nansum(L0_1,axis = 1)/mcruns
L0_2 = np.nansum(L0_2,axis = 1)/mcruns

N_L0_0_m = np.nansum(N_L0_0,axis = 1)/mcruns
N_L0_1_m = np.nansum(N_L0_1,axis = 1)/mcruns
N_L0_2_m = np.nansum(N_L0_2,axis = 1)/mcruns
N_L0_0_sd = np.std(N_L0_0,axis = 1)/mcruns
N_L0_1_sd = np.std(N_L0_1,axis = 1)/mcruns
N_L0_2_sd = np.std(N_L0_2,axis = 1)/mcruns	
error = N_L0_0_sd*qnorm/mcruns
N_L0_0_l = N_L0_0_m-error
N_L0_0_r = 	N_L0_0_m+error
error = N_L0_1_sd*qnorm/mcruns
N_L0_1_l = N_L0_1_m-error
N_L0_1_r = 	N_L0_1_m+error
error = N_L0_2_sd*qnorm/mcruns
N_L0_2_l = N_L0_2_m-error
N_L0_2_r = 	N_L0_2_m+error
	
# -------------------------
# Write to file
# -------------------------		
for igen in xrange(len(gen)):
	outputfile_1 = open(dir+'L0summary_Time'+str(gen[igen])+'_'+outname+'.csv','w')
	outputfile_2 = open(dir+'L0summary_Time'+str(gen[igen])+'_'+outname+'_forARC.csv','w')

	# Write title
	outputfile_1.write('Patch,X,Y,K,N,L0_2,L0_1,L0_0\n')
	outputfile_2.write('Patch,X,Y,K,N,L0_2,L0_1,L0_0\n')

	# Write each info
	for i in range(len(K)):
		outputfile_1.write(str(Patch[i])+',')
		outputfile_1.write(str(X[i])+',')
		outputfile_1.write(str(Y[i])+',')
		outputfile_1.write(str(K[i])+',')
		outputfile_1.write(str(N_m[igen][i])+',')
		outputfile_1.write(str(L0_2[igen][i])+',')
		outputfile_1.write(str(L0_1[igen][i])+',')
		outputfile_1.write(str(L0_0[igen][i])+'\n')
		
		outputfile_2.write(str(Patch[i])+',')
		outputfile_2.write(str(X[i])+',')
		outputfile_2.write(str(Y[i])+',')
		outputfile_2.write(str(K[i])+',')
		if N_m[igen][i] <= 1:
			outputfile_2.write('-9999,-9999,-9999,-9999\n')
		else:
			outputfile_2.write(str(N_m[igen][i])+',')
			outputfile_2.write(str(L0_2[igen][i])+',')
			outputfile_2.write(str(L0_1[igen][i])+',')
			outputfile_2.write(str(L0_0[igen][i])+'\n')
		
	outputfile_1.close()
	outputfile_2.close()	

# --------------------------------------------------
# Plotting N L0 2,1,0 totals through time
# --------------------------------------------------
figure()
plot(gen,N_L0_2_m,'-b',label='WCT',linewidth=2)
if mcruns > 1:
	plot(gen,N_L0_2_l,'-.b',label='',linewidth=2)
	plot(gen,N_L0_2_r,'-.b',label='',linewidth=2)
	plot(gen,N_L0_1_l,'-.g',label='',linewidth=2)
	plot(gen,N_L0_1_r,'-.g',label='',linewidth=2)
	plot(gen,N_L0_0_l,'-.r',label='',linewidth=2)
	plot(gen,N_L0_0_r,'-.r',label='',linewidth=2)
plot(gen,N_L0_1_m,'-g',label='Cutbows',linewidth=2)
plot(gen,N_L0_0_m,'-r',label='RT',linewidth=2)

xlabel('Time',fontsize=18)
ylabel('N',fontsize=18)
#title(plottitle,fontsize=21)
#axis([-0.1,gen,np.min(N_init_pop_m),np.max(N_init_pop_m)])
#axis([-0.1,gen,50000,150000])
legend(loc=0)
savefig(dir+outname+'_N_L0.png',dpi=400)

show()