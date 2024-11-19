# -------------------------------------------------
# plotSizeAge.py
# Erin Landguth
# July 27, 2015
# Description: This script reads in size/age information
# for given year and plots.
# v0 - Initial script
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
dir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2016Sampling/WCT2016Nov22_noGenes2x2/"
# Outputname 
outname = "_Diagnostics_WCT2016Nov22_noGenes2x2_OUT_"
outdir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2016Sampling/WCT2016Nov22_noGenes2x2/summary/"

# batch number - batchrun(batchno)mcrun{mcrun}
batchno = 0

# Number of monte carlo runs - used to average
mcruns = 2

# Generation to extract summary from
#gen = range(0,130,5)
gen = [100]
#gen = [0,13,99]
# Number of patches
n = 1384

# Which grid files to do ind{}.csv or indSample{}.csv (ind or indSample)
gridformat = 'indSample'

# Calculate each statistic - Error (use qnorm(0.975)=1.959964 - couldn't find equivalent qnorm) 
qnorm = 1.959964

# ------------------
# End User Input
# ------------------

# List folders in this dir
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# --------------------
# Begin gen loop
# --------------------
size_2 = []
age_2 = []
size_1 = []
age_1 = []
size_0 = []
age_0 = []
L0_2 = []
L0_1 = []
L0_0 = []
for igen in xrange(len(gen)):
	
	# ----------------------------	
	#  Preliminary vector storage
	# ----------------------------
	# All information
	size_2.append([])
	age_2.append([])
	size_1.append([])
	age_1.append([])
	size_0.append([])
	age_0.append([])
	L0_2.append([])
	L0_1.append([])
	L0_0.append([])
	
	# -----------------------------------	
	# Read in and store metrics
	# -----------------------------------
	# Loop over folders
	for imc in xrange(mcruns): 
		
		# Add a batch spot to the vectors
		size_2[igen].append([])
		age_2[igen].append([])
		size_1[igen].append([])
		age_1[igen].append([])
		size_0[igen].append([])
		age_0[igen].append([])
		L0_2[igen].append([])
		L0_1[igen].append([])
		L0_0[igen].append([])
		
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
		age = []
		size = []
		# Extract values from this grid	
		for i in xrange(len(x)-1):
			size.append(float(x[i+1][6]))
			age.append(int(x[i+1][5]))
						
		age = np.asarray(age)
		size = np.asarray(size)
				
# --------------------------------
# Summary MC
# --------------------------------

	
# -------------------------
# Write to file
# -------------------------		

# ----------------------------
# poly fits
# ----------------------------
xp = np.linspace(0,14,100)
z_all = np.polyfit(age,size,2)
p_all = np.poly1d(z_all)

# --------------------------------------------------
# Plotting size v age
# --------------------------------------------------
fontsize=18
savedpi = 1000

figure()
plot(age,size,'.b')
plot(xp,p_all(xp),'-b',label='')
xlabel('Age',fontsize=20)
ylabel('Size',fontsize=20)
axis([-0.1,10,-0.1,250])
legend(loc=0,prop={'size':14})
# Updating fontsize on axes
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
savefig(outdir+outname+'_AgeSize.png',dpi=savedpi,bbox_inches='tight')

show()