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
outname = "_Diagnostics_WCT2016Nov22_noGenes2x2_Back_"
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
		L0 = []
		age = []
		size = []
		# Extract values from this grid	
		for i in xrange(len(x)-1):
			size.append(float(x[i+1][6]))
			age.append(int(x[i+1][5]))
			L0.append(int(x[i+1][14]))			
		L0 = np.asarray(L0)
		age = np.asarray(age)
		size = np.asarray(size)
		age_0[igen][imc].append(age[np.where(L0 == 0)[0]])
		age_1[igen][imc].append(age[np.where(L0 == 1)[0]])
		age_2[igen][imc].append(age[np.where(L0 == 2)[0]])
		size_0[igen][imc].append(size[np.where(L0 == 0)[0]])
		size_1[igen][imc].append(size[np.where(L0 == 1)[0]])
		size_2[igen][imc].append(size[np.where(L0 == 2)[0]])
		
# --------------------------------
# Summary MC
# --------------------------------
'''
age_0 = np.nansum(age_0,axis = 1)/mcruns
age_1 = np.nansum(L0_1,axis = 1)/mcruns
age_2 = np.nansum(L0_2,axis = 1)/mcruns

age_0 = np.nansum(L0_0,axis = 1)/mcruns
age_1 = np.nansum(L0_1,axis = 1)/mcruns
age_2 = np.nansum(L0_2,axis = 1)/mcruns

	
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
'''

# ----------------------------
# poly fits
# ----------------------------
z_2 = np.polyfit(age_2[0][0][0],size_2[0][0][0],2)
p_2 = np.poly1d(z_2)
z_1 = np.polyfit(age_1[0][0][0],size_1[0][0][0],2)
p_1 = np.poly1d(z_1)
z_0 = np.polyfit(age_0[0][0][0],size_0[0][0][0],2)
p_0 = np.poly1d(z_0)
xp = np.linspace(0,14,100)

z_all = np.polyfit(age,size,2)
p_all = np.poly1d(z_all)

# --------------------------------------------------
# Plotting size v age
# --------------------------------------------------
fontsize=18
savedpi = 1000
figure()
plot(age_2[0][0][0],size_2[0][0][0],'.y')
plot(xp,p_2(xp),'-y',label='AABB')
xlabel('Age',fontsize=22)
ylabel('Size',fontsize=22)
axis([-0.1,10,-0.1,250])
legend(loc=0,prop={'size':14})
# Updating fontsize on axes
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
savefig(dir+outname+'_AgeSize_L0_2.png',dpi=savedpi,bbox_inches='tight')

figure()
plot(age_1[igen][0][0],size_1[0][0][0],'.g')
plot(xp,p_1(xp),'-g',label='AaBb')
xlabel('Age',fontsize=22)
ylabel('Size',fontsize=22)
axis([-0.1,10,-0.1,250])
legend(loc=0,prop={'size':14})
# Updating fontsize on axes
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
savefig(dir+outname+'_AgeSize_L0_1.png',dpi=savedpi,bbox_inches='tight')

figure()
plot(age_0[igen][0][0],size_0[0][0][0],'.r')
plot(xp,p_0(xp),'-r',label='aabb')
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
savefig(dir+outname+'_AgeSize_L0_0.png',dpi=savedpi,bbox_inches='tight')

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
savefig(dir+outname+'_AgeSize_L0_0.png',dpi=savedpi,bbox_inches='tight')

show()