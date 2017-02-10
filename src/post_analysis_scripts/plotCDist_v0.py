# ----------------------------------------------------------
# plotCDist_v0.py
# Erin Landguth/Casey Day
# June 15, 2015
# Description: This script reads in CDDist movement column and plots 
# histogram to compare movement scenarios
# v0 - Initial script
# ----------------------------------------------------------

# Import statements
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."				
import pdb,os,copy		
from random import *
import pylab as P

# Functions needed
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
#-----------------
# User Input
#-----------------
# Directory location of output.csv
dir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2016Sampling/WCT2016Nov22_noGenes2x2_6kmMove_modKv3/"

# Outputname 
outname = "_Diagnostics_WCT2016Nov22_noGenes2x2_6kmMove_modKv3_OUT_"
outdir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2016Sampling/WCT2016Nov22_noGenes2x2_6kmMove_modKv3/summary/"
xlabel = 'Probability'

# batch number - batchrun(batchno)mcrun{mcrun}
batchno = 0

# MCRun number
mcrun = 0

# Generations to look at
gen = [100]

# Which grid files to do ind{}.csv or indSample{}.csv (ind or indSample)
gridformat = 'indSample'

# ------------------
# End User Input
# ------------------

# List folders in this dir
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

cdist = []
N = []
Age0s = []

# --------------------
# Begin gen loop
# --------------------
for igen in xrange(len(gen)):

	cdist.append([])
	N.append([])
	Age0s.append([])
	
	# --------------------------
	#  Read in information 
	# --------------------------		
	# Open file to extract number of migrants
	xyinputfile = open(dir+"batchrun"+str(batchno)+"mcrun"+str(mcrun)+"/"+gridformat+str(gen[igen])+".csv",'r')

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
	
	# Extract values from this grid	
	for i in xrange(len(x)-1):
		cdist[igen].append(float(x[i+1][13]))
		N[igen].append(int(x[i+1][0]))
	

# Get N numbers
for igen in xrange(len(gen)):	
	N[igen] = len(N[igen])
	cdist_temp = np.asarray(cdist[igen],dtype=float)
	Age0s[igen].append(len(np.where(cdist_temp==-9999.0)[0]))
	
# ---------------------
# Plotting
# ---------------------
for igen in xrange(len(gen)):
	cdist_temp = np.asarray(cdist[igen],dtype=float)
	cdist_noAge0s = cdist_temp[np.where(cdist_temp!=-9999.0)[0]]
	if len(cdist_noAge0s) != 0:
		P.figure(igen)
		P.hist(cdist_noAge0s, histtype='bar')
		P.xlabel(xlabel)
		P.title("Year "+str(gen[igen])+"; N_move = "+str(N[igen])+ "; N_nomove = "+str(Age0s[igen][0]))
		P.savefig(outdir+outname+'CostDistance_'+str(gen[igen])+'.png',dpi=200)
		
P.show()	
	












