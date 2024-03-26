# -----------------------------------------------------------------------------
# diagnostics.py
# 2014 01 29: Erin Landguth
# This script grabs output.csv files, splits bars, and grabs diagnostics for checking CDPOP/CDFISH/CDmetaPOP
# v0 - Initial script
# v1 - 2014 07 11: Matching v0.79 output.
# v2 - 2014 07 30: Add in histograms of RSI and Age N at given time units
# v3 - 2014 08 20: Sensitivity run comparisons 
# v4 - 2014 08 25: Extract growth rate Nt/Nt-1, YOY/matureN and write to file: YoY = N_initial_age[0] at time+1 and mature N = total breed females + total breed males
# v5 - 2014 08 31: Summary metrics by patch and averaged.
# v6 - 2014 09 05: More summary over time
# v7 - 2014 10 08: Add in 1/2N rate of decay to He/Ho plots for checking genetic validation.
# _batch_v0: Create general for standard and +- runs (5 total).
# _batch_v1: Add in growth rate for patch mean Nt / Nt-1 - read in patch All time
# _batch_v099: Update for recent version. 
# _batch_v1.08: Update for recent version. Add in more diagnostics: mortalities.
# v2: Added more metrics: packing and move mortalities.
# v3: Added more metrics: captured individual plots. Also write each batch at given time to file for patch specific numbers
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
dir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT415_Sullivan_April2018/CompareK/"
plottitle = ''
savename = "_CompareK_"
label = ['K', '1.6K', '2K']
batchno = 3 # Total number of batches in Compare folder
linemarks = ['k--','b-o','r-','g-^','ys-']

# For Emi/Immi plots
label_D1 = ['K Emi', '1.6K Emi', '2K Emi']
label_D2 = ['K Imm', '1.6K Imm', '2K Imm']
linemarks_D1 = ['k--','b--','r--','g--','y--']
linemarks_D2 = ['k-o','b-o','r-o','g-o','y-o']
# For captured plots
label_Cap1 = ['K CapBack', '1.6K CapBack', '2K CapBack']
label_Cap2 = ['K CapOut', '1.6K CapOut', '2K CapOut']
linemarks_Cap1 = ['k--','b--','r--','g--','y--']
linemarks_Cap2 = ['k-o','b-o','r-o','g-o','y-o']

outdir = dir+"summary/"
xyfile = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT415_Sullivan_April2018/patchvars/PatchVars415_K.csv"

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 91 # Number of years 
nthfile = range(0,gen,1)
#nthfile = np.asarraty([0,1,2,3,4,5,10,19])
mcno = 2 # Number of MCs
plottime = np.asarray([90])
plotagesize = 'Y' # Plot the time specific age/size information?
#label = ['0%','25%','50%','75%','100%']
# Vertical lines
popburn = 25
barr0 = 51
barr1 = 57
#barr2 = 59
#barr3 = 64
#barr4 = 65
#barr5 = 66
#barr6 = 69
#barr7 = 71


# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# ---------------------------------
# Storage variables
# ---------------------------------
# Store values - Class
N_init_age = []
N_init_class = []
size_age = []
N_afterEmi = []
N_afterImmi = []
N_growthback = []
N_growthout = []
# Store values - Patch
PackDeaths_Emi = []
PackDeaths_Immi = []
MoveDeaths_Emi = []
MoveDeaths_Immi = []
N_init_patch = []
N_cap_back_patch = [] # Patch capture when back at spawn grounds
N_cap_back_pop = [] # total population captured when back
N_cap_out_patch = [] # patch capture when out
N_cap_out_pop = [] # total population captured when out

# Loop through batches
for ibatch in xrange(batchno):

	# Add storage spot
	N_init_age.append([])
	N_init_class.append([])
	size_age.append([])
	N_afterEmi.append([])
	N_afterImmi.append([])
	N_growthback.append([])
	N_growthout.append([])
	N_init_patch.append([])
	PackDeaths_Emi.append([])
	PackDeaths_Immi.append([])
	MoveDeaths_Emi.append([])
	MoveDeaths_Immi.append([])
	N_cap_back_patch.append([])
	N_cap_back_pop.append([])
	N_cap_out_patch.append([])
	N_cap_out_pop.append([])

	# Loop through MCs
	for imc in xrange(mcno):
	
		# -------------------------------
		# Read in classAllTime values
		# -------------------------------
		# Open output.csv file in folder
		inputfile = open(dir+'batchrun'+str(ibatch)+'mcrun'+str(imc)+'/summary_classAllTime.csv')
		
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
		
		# -------------------------------
		# Read in patchAllTime values
		# -------------------------------
		# Open output.csv file in folder
		inputfile = open(dir+'batchrun'+str(ibatch)+'mcrun'+str(imc)+'/summary_popAllTime.csv')
		
		# Read lines from the file
		lines = inputfile.readlines()

		#Close the file
		inputfile.close()

		# Create an empty matrix to append to
		values_pop = []

		# Split up each line in file and append to empty matrix for generation specified
		for i in xrange(len(lines)):
			thisline = lines[i].split(',')
			values_pop.append(thisline)
		# Delete lines
		del(lines)
		
		# Extract once ages and size classes
		ageclass = np.asarray(values[1][1].split('|'))[0:-1]
		ageclass = ageclass.astype(np.int)
		sizeclass = np.asarray(values[1][5].split('|'))[0:-1]
		sizeclass = sizeclass.astype(np.float)
		
		# Store values
		N_init_age[ibatch].append([])
		N_init_class[ibatch].append([])
		size_age[ibatch].append([])
		N_afterEmi[ibatch].append([])
		N_afterImmi[ibatch].append([])
		N_growthback[ibatch].append([])
		N_growthout[ibatch].append([])
		N_init_patch[ibatch].append([])
		PackDeaths_Emi[ibatch].append([])
		PackDeaths_Immi[ibatch].append([])
		MoveDeaths_Emi[ibatch].append([])
		MoveDeaths_Immi[ibatch].append([])
		N_cap_back_patch[ibatch].append([])
		N_cap_back_pop[ibatch].append([])
		N_cap_out_patch[ibatch].append([])
		N_cap_out_pop[ibatch].append([])
		
		# Then Loop through generations/time
		for iout in xrange(gen):
		
			N_init_age[ibatch][imc].append([])
			N_init_class[ibatch][imc].append([])
			size_age[ibatch][imc].append([])
			N_afterEmi[ibatch][imc].append([])
			N_afterImmi[ibatch][imc].append([])
			N_growthback[ibatch][imc].append([])
			N_growthout[ibatch][imc].append([])
			N_init_patch[ibatch][imc].append([])
			PackDeaths_Emi[ibatch][imc].append([])
			PackDeaths_Immi[ibatch][imc].append([])
			MoveDeaths_Emi[ibatch][imc].append([])
			MoveDeaths_Immi[ibatch][imc].append([])
			N_cap_back_patch[ibatch][imc].append([])
			N_cap_back_pop[ibatch][imc].append([])
			N_cap_out_patch[ibatch][imc].append([])
			N_cap_out_pop[ibatch][imc].append([])
			
			# Age split
			for j in xrange(len(values[1+iout][2].split('|'))-1):
				N_init_age[ibatch][imc][iout].append(int(values[1+iout][2].split('|')[j])) 
				N_init_class[ibatch][imc][iout].append(int(values[1+iout][6].split('|')[j]))
				size_age[ibatch][imc][iout].append(float(values[1+iout][3].split('|')[j])) 
				if values[1+iout][12].split('|')[j] == 'NA':
					N_afterEmi[ibatch][imc][iout].append(nan)
				else:
					N_afterEmi[ibatch][imc][iout].append(int(values[1+iout][12].split('|')[j]))
				if values[1+iout][19].split('|')[j] == 'NA':
					N_afterImmi[ibatch][imc][iout].append(nan)
				else:
					N_afterImmi[ibatch][imc][iout].append(int(values[1+iout][19].split('|')[j]))
				if values[1+iout][9].split('|')[j] == 'NA':
					N_growthback[ibatch][imc][iout].append(nan)
				else:
					N_growthback[ibatch][imc][iout].append(int(values[1+iout][9].split('|')[j]))
				if values[1+iout][15].split('|')[j] == 'NA':
					N_growthout[ibatch][imc][iout].append(nan)
				else:
					N_growthout[ibatch][imc][iout].append(int(values[1+iout][15].split('|')[j]))
				
			# Grab all patch values - patch values with total
			for j in xrange(1,len(values_pop[1+iout][3].split('|'))-1):
				N_init_patch[ibatch][imc][iout].append(float(values_pop[1+iout][3].split('|')[j])) # remove todal
				PackDeaths_Emi[ibatch][imc][iout].append(float(values_pop[1+iout][18].split('|')[j])) # remove todal
				PackDeaths_Immi[ibatch][imc][iout].append(float(values_pop[1+iout][26].split('|')[j])) # remove todal
				MoveDeaths_Emi[ibatch][imc][iout].append(float(values_pop[1+iout][17].split('|')[j])) # remove todal
				MoveDeaths_Immi[ibatch][imc][iout].append(float(values_pop[1+iout][24].split('|')[j])) # remove todal
				
				# Patch values without total
				N_cap_back_patch[ibatch][imc][iout].append(float(values_pop[1+iout][15].split('|')[j-1]))
				N_cap_out_patch[ibatch][imc][iout].append(float(values_pop[1+iout][22].split('|')[j-1]))
			# Sum the capture
			N_cap_back_pop[ibatch][imc][iout] = np.nansum(N_cap_back_patch[ibatch][imc][iout])
			N_cap_out_pop[ibatch][imc][iout] = np.nansum(N_cap_out_patch[ibatch][imc][iout])
			
			
# Turn into arrays
N_init_age = np.asarray(N_init_age)
N_init_age = N_init_age.astype(np.float)
N_init_class = np.asarray(N_init_class)
N_init_class = np.asarray(N_init_class)
size_age = np.asarray(size_age)
N_afterEmi = np.asarray(N_afterEmi)
N_afterImmi = np.asarray(N_afterImmi)
N_growthback = np.asarray(N_growthback)
N_growthout = np.asarray(N_growthout)
# Patch values
N_init_patch = np.asarray(N_init_patch)
PackDeaths_Emi = np.asarray(PackDeaths_Emi)
PackDeaths_Immi = np.asarray(PackDeaths_Immi)
MoveDeaths_Emi = np.asarray(MoveDeaths_Emi)
MoveDeaths_Immi = np.asarray(MoveDeaths_Immi)
N_cap_back_patch = np.asarray(N_cap_back_patch)
N_cap_back_pop = np.asarray(N_cap_back_pop)
N_cap_out_patch = np.asarray(N_cap_out_patch)
N_cap_out_pop = np.asarray(N_cap_out_pop)
# --------------------------------------------
# Get mean over Monte Carlosfor each batch run
# --------------------------------------------
# Capture numbers - back patch
N_cap_back_patch_m = np.nansum(N_cap_back_patch,axis=1)/mcno
N_cap_back_patch_sd = np.std(N_cap_back_patch,axis=1)	
error = qnorm*N_cap_back_patch_sd/(mcno)
N_cap_back_patch_l = N_cap_back_patch_m-error
N_cap_back_patch_r = 	N_cap_back_patch_m+error
N_cap_back_patch_min = np.min(N_cap_back_patch,axis=1)
N_cap_back_patch_max = np.max(N_cap_back_patch,axis=1)
# Capture numbers - back pop
N_cap_back_pop_m = np.nansum(N_cap_back_pop,axis=1)/mcno
N_cap_back_pop_sd = np.std(N_cap_back_pop,axis=1)	
error = qnorm*N_cap_back_pop_sd/(mcno)
N_cap_back_pop_l = N_cap_back_pop_m-error
N_cap_back_pop_r = 	N_cap_back_pop_m+error
N_cap_back_pop_min = np.min(N_cap_back_pop,axis=1)
N_cap_back_pop_max = np.max(N_cap_back_pop,axis=1)

# Capture numbers - out patch
N_cap_out_patch_m = np.nansum(N_cap_out_patch,axis=1)/mcno
N_cap_out_patch_sd = np.std(N_cap_out_patch,axis=1)	
error = qnorm*N_cap_out_patch_sd/(mcno)
N_cap_out_patch_l = N_cap_out_patch_m-error
N_cap_out_patch_r = 	N_cap_out_patch_m+error
N_cap_out_patch_min = np.min(N_cap_out_patch,axis=1)
N_cap_out_patch_max = np.max(N_cap_out_patch,axis=1)
# Capture numbers - out pop
N_cap_out_pop_m = np.nansum(N_cap_out_pop,axis=1)/mcno
N_cap_out_pop_sd = np.std(N_cap_out_pop,axis=1)	
error = qnorm*N_cap_out_pop_sd/(mcno)
N_cap_out_pop_l = N_cap_out_pop_m-error
N_cap_out_pop_r = 	N_cap_out_pop_m+error
N_cap_out_pop_min = np.min(N_cap_out_pop,axis=1)
N_cap_out_pop_max = np.max(N_cap_out_pop,axis=1)

N_init_pop_m = np.nansum(N_init_age,axis=3)
#N_init_pop_m = np.nansum(N_init_pop_m,axis=1)
#N_init_pop_sd = np.std(N_init_pop_m,axis=1)
N_init_pop_m = np.nansum(N_init_pop_m,axis=1)/mcno

N_init_age_m = np.nansum(N_init_age[:][:],axis=1)/mcno
N_init_age_sd = np.std(N_init_age[:][:],axis=1)	
error = N_init_age_sd*qnorm/mcno
N_init_age_l = N_init_age_m-error
N_init_age_r = 	N_init_age_m+error
N_init_age_min = np.min(N_init_age[:][:],axis=1)
N_init_age_max = np.max(N_init_age[:][:],axis=1)
YoY_m = N_init_age_m[:,:,0]

N_init_class_m = np.nansum(N_init_class[:][:],axis=1)/mcno
N_init_class_sd = np.std(N_init_class[:][:],axis=1)	
error = N_init_class_sd*qnorm/mcno
N_init_class_l = N_init_class_m-error
N_init_class_r = 	N_init_class_m+error
N_init_class_min = np.min(N_init_class[:][:],axis=1)
N_init_class_max = np.max(N_init_class[:][:],axis=1)

size_age_m = np.nansum(size_age[:][:],axis=1)/mcno
size_age_sd = np.std(size_age[:][:],axis=1)	
error = qnorm*size_age_sd/(mcno)
size_age_l = size_age_m-error
size_age_r = 	size_age_m+error
size_age_min = np.min(size_age[:][:],axis=1)
size_age_max = np.max(size_age[:][:],axis=1)

N_afterEmi_m = np.nansum(N_afterEmi[:][:],axis=1)/mcno
N_afterEmi_sd = np.std(N_afterEmi[:][:],axis=1)	
error = N_afterEmi_sd*qnorm/mcno
N_afterEmi_l = N_afterEmi_m-error
N_afterEmi_r = 	N_afterEmi_m+error
N_afterEmi_min = np.min(N_afterEmi[:][:],axis=1)
N_afterEmi_max = np.max(N_afterEmi[:][:],axis=1)
# Get Total
N_afterEmi_m_total = np.nansum(N_afterEmi_m,axis=2)
N_afterEmi_sd_total = np.nansum(N_afterEmi_sd,axis=2)

N_afterImmi_m = np.nansum(N_afterImmi[:][:],axis=1)/mcno
N_afterImmi_sd = np.std(N_afterImmi[:][:],axis=1)	
error = N_afterImmi_sd*qnorm/mcno
N_afterImmi_l = N_afterImmi_m-error
N_afterImmi_r = 	N_afterImmi_m+error
N_afterImmi_min = np.min(N_afterImmi[:][:],axis=1)
N_afterImmi_max = np.max(N_afterImmi[:][:],axis=1)
# Get Total
N_afterImmi_m_total = np.nansum(N_afterImmi_m,axis=2)
N_afterImmi_sd_total = np.nansum(N_afterImmi_sd,axis=2)

N_growthback_m = np.nansum(N_growthback[:][:],axis=1)/mcno
N_growthback_sd = np.std(N_growthback[:][:],axis=1)	
error = N_growthback_sd*qnorm/mcno
N_growthback_l = N_growthback_m-error
N_growthback_r = 	N_growthback_m+error
N_growthback_min = np.min(N_growthback[:][:],axis=1)
N_growthback_max = np.max(N_growthback[:][:],axis=1)

N_growthout_m = np.nansum(N_growthout[:][:],axis=1)/mcno
N_growthout_sd = np.std(N_growthout[:][:],axis=1)	
error = N_growthout_sd*qnorm/mcno
N_growthout_l = N_growthout_m-error
N_growthout_r = 	N_growthout_m+error
N_growthout_min = np.min(N_growthout[:][:],axis=1)
N_growthout_max = np.max(N_growthout[:][:],axis=1)

# Patch numbers [batch][mcrun][time][patch]
N_init_patch_m = np.nansum(N_init_patch,axis=1)/mcno
N_init_patch_sd = np.std(N_init_patch,axis=1)	
error = qnorm*N_init_patch_sd/(mcno)
N_init_patch_l = N_init_patch_m-error
N_init_patch_r = 	N_init_patch_m+error
N_init_patch_min = np.min(N_init_patch,axis=1)
N_init_patch_max = np.max(N_init_patch,axis=1)

PackDeaths_Emi_m = np.nansum(PackDeaths_Emi,axis=1)/mcno
PackDeaths_Emi_sd = np.std(PackDeaths_Emi,axis=1)
# Get Total
PackDeaths_Emi_m_total = np.nansum(PackDeaths_Emi_m,axis=2)
PackDeaths_Emi_sd_total = np.nansum(PackDeaths_Emi_sd,axis=2)

PackDeaths_Immi_m = np.nansum(PackDeaths_Immi,axis=1)/mcno
PackDeaths_Immi_sd = np.std(PackDeaths_Immi,axis=1)
# Get Total
PackDeaths_Immi_m_total = np.nansum(PackDeaths_Immi_m,axis=2)
PackDeaths_Immi_sd_total = np.nansum(PackDeaths_Immi_sd,axis=2)

MoveDeaths_Emi_m = np.nansum(MoveDeaths_Emi,axis=1)/mcno
MoveDeaths_Emi_sd = np.std(MoveDeaths_Emi,axis=1)
# Get Total
MoveDeaths_Emi_m_total = np.nansum(MoveDeaths_Emi_m,axis=2)
MoveDeaths_Emi_sd_total = np.nansum(MoveDeaths_Emi_sd,axis=2)

MoveDeaths_Immi_m = np.nansum(MoveDeaths_Immi[:][:],axis=1)/mcno
MoveDeaths_Immi_sd = np.std(MoveDeaths_Immi[:][:],axis=1)
# Get Total
MoveDeaths_Immi_m_total = np.nansum(MoveDeaths_Immi_m,axis=2)
MoveDeaths_Immi_sd_total = np.nansum(MoveDeaths_Immi_sd,axis=2)

# --------------------------------
# Other summary data 
# -------------------------------- 

# Growth rate - patch population Nt / Nt-1
N_growth_patch =  N_init_patch_m[:,1:gen,:] / N_init_patch_m[:,0:gen-1,:] 
N_growth_patch[np.where(N_growth_patch==inf)]=nan # remove inf values
# Use axis 2 here for patch mean
N_growth_patch_m = np.nansum(N_growth_patch,axis=2)/len(N_init_patch_m[0][0]) 
N_growth_patch_sd = np.nanstd(N_growth_patch,axis=2)
N_growth_patch_cv = N_growth_patch_m / 	N_growth_patch_sd
error = qnorm*N_growth_patch_sd/len(N_init_patch_m[0][0])
N_growth_patch_l = N_growth_patch_m-error
N_growth_patch_r = 	N_growth_patch_m+error
N_growth_patch_min = np.nanmin(N_growth_patch,axis=2)
N_growth_patch_max = np.nanmax(N_growth_patch,axis=2)	

# Get mean of growth rate - patch population Nt/Nt-1 across time steps
N_growth_patch_m_time_m = np.nansum(N_growth_patch_m,axis=1)/len(N_growth_patch_m[0])
N_growth_patch_sd_time_m = np.nansum(N_growth_patch_sd[:][:],axis=1)/len(N_growth_patch_m[0])
N_growth_patch_time_cv = N_growth_patch_m_time_m / N_growth_patch_sd_time_m

# Normalized by baseline
#N_init_age_m_norm = N_init_age_m / N_init_age_m[2]

# --------------------------------------------
# Write to file
# --------------------------------------------
# X,Y,patch population
X = []
Y = []
SubPop = []
# Read in PatchVars file to get X,Y
inputfile = open(xyfile)
# Read lines from the file
lines = inputfile.readlines()
#Close the file
inputfile.close()
# Create an empty matrix to append to
xyvalues = []
# Split up each line in file and append to empty matrix for generation specified
for i in xrange(len(lines)):
	thisline = lines[i].split(',')
	xyvalues.append(thisline)
# Delete lines
del(lines)
for i in xrange(len(xyvalues)-1):
	SubPop.append(xyvalues[i+1][0])
	X.append(xyvalues[i+1][1])
	Y.append(xyvalues[i+1][2])
del(xyvalues)

# Create file for each batch to write info to
for ibatch in xrange(batchno):
	for itime in xrange(len(plottime)):
		outputfile = open(outdir+'Batch'+str(ibatch)+'_'+label[ibatch]+'_'+savename+'mean_patch_pops_time'+str(plottime[itime])+'.csv','w')
		# Write out the titles
		outputfile.write('Subpopulation,X,Y,N_Back,N_Back_Captured,N_Out_Captured\n')
			
		# WRite information
		for i in xrange(len(X)):
			outputfile.write(SubPop[i]+',')
			outputfile.write(X[i]+',')
			outputfile.write(Y[i]+',')
			
			outputfile.write(str(N_init_patch_m[ibatch][plottime[itime]][i])+',')
			outputfile.write(str(N_cap_back_patch_m[ibatch][plottime[itime]][i])+',')
			outputfile.write(str(N_cap_out_patch_m[ibatch][plottime[itime]][i])+'\n')
					
		outputfile.close()		

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------

# Plot total N
# ------------
figure()
for i in xrange(len(N_init_pop_m)):
	plot(nthfile,N_init_pop_m[i],linemarks[i],label=label[i],linewidth=2)
axvline(x=popburn, color='k', linestyle='--')
axvline(x=barr0, color='k', linestyle='--')	
axvline(x=barr1, color='k', linestyle='--')
#axvline(x=barr2, color='k', linestyle='--')
#axvline(x=barr3, color='k', linestyle='--')
#axvline(x=barr4, color='k', linestyle='--')
#axvline(x=barr5, color='k', linestyle='--')
#axvline(x=barr6, color='k', linestyle='--')
#axvline(x=barr7, color='k', linestyle='--')

#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Population',fontsize=18)
title(plottitle,fontsize=21)
#axis([-0.1,gen,np.min(N_init_pop_m),np.max(N_init_pop_m)])
#axis([-0.1,gen,50000,150000])
#axis([-0.01,130,40000,120000])
#axis([-0.01,130,0,100000])
legend(loc=0)
savefig(dir+savename+'NInit_pop.png',dpi=savedpi)

# Plot captured Ns
# -----------------
figure()
for i in xrange(len(N_init_pop_m)):
	plot(nthfile,N_cap_back_pop_m[i],linemarks_Cap1[i],label=label_Cap1[i],linewidth=2)
	plot(nthfile,N_cap_out_pop_m[i],linemarks_Cap2[i],label=label_Cap2[i],linewidth=2)
axvline(x=popburn, color='k', linestyle='--')
axvline(x=barr0, color='k', linestyle='--')	
axvline(x=barr1, color='k', linestyle='--')
#axvline(x=barr2, color='k', linestyle='--')
#axvline(x=barr3, color='k', linestyle='--')
#axvline(x=barr4, color='k', linestyle='--')
#axvline(x=barr5, color='k', linestyle='--')
#axvline(x=barr6, color='k', linestyle='--')
#axvline(x=barr7, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Captured Population',fontsize=18)
title(plottitle,fontsize=21)
legend(loc=0)
savefig(dir+savename+'CapturedN.png',dpi=savedpi)



# Plot patch mean
# -------------------------------
figure()
for i in xrange(len(N_growth_patch_m)):
	plot(nthfile[1:],N_growth_patch_m[i],linemarks[i],label=label[i],linewidth=2)
xlabel('Time',fontsize=18)
ylabel('Growth rate patch mean',fontsize=18)
title(plottitle,fontsize=21)
#axis([-0.1,gen,0,2.0])
legend(loc=0)	
savefig(dir+savename+'NInit_growthrate_patchmean.png',dpi=savedpi)	

if plotagesize == 'Y':
	# Plot N initial - age class
	# -----------------------------
	ind = np.arange(len(N_init_age_m[0][0])) # the x locations for the groups
	width = 0.15                      # the width of the bars

	# Loop through each year to plot
	for it in plottime:
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		
		if batchno == 5:
			
			rects0 = ax.bar(ind,N_init_age_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_init_age_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

			rects1 = ax.bar(ind+width,N_init_age_m[1][it],width,color='blue',edgecolor='black',hatch="\\",yerr=N_init_age_sd[1][it],error_kw=dict(elinewidth=2,ecolor='black'))
			
			rects2 = ax.bar(ind+2*width,N_init_age_m[2][it],width,color='grey',edgecolor='black',hatch="0",yerr=N_init_age_sd[2][it],error_kw=dict(elinewidth=2,ecolor='black'))	
			
			rects3 = ax.bar(ind+3*width,N_init_age_m[3][it],width,color='blue',edgecolor='black',hatch="+",yerr=N_init_age_sd[3][it],error_kw=dict(elinewidth=2,ecolor='black'))
			
			rects4 = ax.bar(ind+4*width,N_init_age_m[4][it],width,color='red',edgecolor='black',hatch="-",yerr=N_init_age_sd[4][it],error_kw=dict(elinewidth=2,ecolor='black'))		
				
			# axes and labels
			ax.set_ylabel('N',fontsize=18)
			ax.set_xlabel('Age',fontsize=18)
			ax.set_title(plottitle+'N Initial (Age) Year ' +str(it),fontsize=21)
			ax.set_xlim(-width,len(ind)+width)
			ax.set_ylim(0, max(np.max(N_init_age_m[0][it]),np.max(N_init_age_m[1][it]),np.max(N_init_age_m[2][it]),np.max(N_init_age_m[3][it]),np.max(N_init_age_m[4][it])))
			xTickMarks = [str(i) for i in xrange(len(N_init_age_m[0][0]))]
			ax.set_xticks(ind+width)
			xtickNames = ax.set_xticklabels(xTickMarks)
			plt.setp(xtickNames, rotation=0)
			ax.legend((rects0[0],rects1[0],rects2[0],rects3[0],rects4[0]),label,loc=0)
			
		if batchno == 4:
		
			rects0 = ax.bar(ind,N_init_age_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_init_age_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

			rects1 = ax.bar(ind+width,N_init_age_m[1][it],width,color='blue',edgecolor='black',hatch="\\",yerr=N_init_age_sd[1][it],error_kw=dict(elinewidth=2,ecolor='black'))
			
			rects2 = ax.bar(ind+2*width,N_init_age_m[2][it],width,color='grey',edgecolor='black',hatch="0",yerr=N_init_age_sd[2][it],error_kw=dict(elinewidth=2,ecolor='black'))	
			
			rects3 = ax.bar(ind+3*width,N_init_age_m[3][it],width,color='blue',edgecolor='black',hatch="+",yerr=N_init_age_sd[3][it],error_kw=dict(elinewidth=2,ecolor='black'))
				
				
			# axes and labels
			ax.set_ylabel('N',fontsize=18)
			ax.set_xlabel('Age',fontsize=18)
			ax.set_title(plottitle+'N Initial (Age) Year ' +str(it),fontsize=21)
			ax.set_xlim(-width,len(ind)+width)
			ax.set_ylim(0, max(np.max(N_init_age_m[0][it]),np.max(N_init_age_m[1][it]),np.max(N_init_age_m[2][it]),np.max(N_init_age_m[3][it])))
			xTickMarks = [str(i) for i in xrange(len(N_init_age_m[0][0]))]
			ax.set_xticks(ind+width)
			xtickNames = ax.set_xticklabels(xTickMarks)
			plt.setp(xtickNames, rotation=0)
			ax.legend((rects0[0],rects1[0],rects2[0],rects3[0]),label,loc=0)
		
		if batchno == 3:
		
			rects0 = ax.bar(ind,N_init_age_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_init_age_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

			rects1 = ax.bar(ind+width,N_init_age_m[1][it],width,color='blue',edgecolor='black',hatch="\\",yerr=N_init_age_sd[1][it],error_kw=dict(elinewidth=2,ecolor='black'))
			
			rects2 = ax.bar(ind+2*width,N_init_age_m[2][it],width,color='grey',edgecolor='black',hatch="0",yerr=N_init_age_sd[2][it],error_kw=dict(elinewidth=2,ecolor='black'))	
			
			# axes and labels
			ax.set_ylabel('N',fontsize=18)
			ax.set_xlabel('Age',fontsize=18)
			ax.set_title(plottitle+'N Initial (Age) Year ' +str(it),fontsize=21)
			ax.set_xlim(-width,len(ind)+width)
			ax.set_ylim(0, max(np.max(N_init_age_m[0][it]),np.max(N_init_age_m[1][it]),np.max(N_init_age_m[2][it])))
			xTickMarks = [str(i) for i in xrange(len(N_init_age_m[0][0]))]
			ax.set_xticks(ind+width)
			xtickNames = ax.set_xticklabels(xTickMarks)
			plt.setp(xtickNames, rotation=0)
			ax.legend((rects0[0],rects1[0],rects2[0]),label,loc=0)
		
		
		if batchno == 2:
		
			rects0 = ax.bar(ind,N_init_age_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_init_age_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

			rects1 = ax.bar(ind+width,N_init_age_m[1][it],width,color='blue',edgecolor='black',hatch="\\",yerr=N_init_age_sd[1][it],error_kw=dict(elinewidth=2,ecolor='black'))			
				
			# axes and labels
			ax.set_ylabel('N',fontsize=18)
			ax.set_xlabel('Age',fontsize=18)
			ax.set_title(plottitle+'N Initial (Age) Year ' +str(it),fontsize=21)
			ax.set_xlim(-width,len(ind)+width)
			ax.set_ylim(0, max(np.max(N_init_age_m[0][it]),np.max(N_init_age_m[1][it])))
			xTickMarks = [str(i) for i in xrange(len(N_init_age_m[0][0]))]
			ax.set_xticks(ind+width)
			xtickNames = ax.set_xticklabels(xTickMarks)
			plt.setp(xtickNames, rotation=0)
			ax.legend((rects0[0],rects1[0]),label,loc=0)
		savefig(dir+savename+'NInit_ageclass_year_'+str(it)+'.png',dpi=savedpi)

	# Plot N initial - size class
	# -----------------------------
	ind = np.arange(len(N_init_class_m[0][0])) # the x locations for the groups
	width = 0.15                      # the width of the bars

	# Loop through each year to plot
	for it in plottime:
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		
		if batchno == 5:
		
			rects0 = ax.bar(ind,N_init_class_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_init_class_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

			rects1 = ax.bar(ind+width,N_init_class_m[1][it],width,color='blue',edgecolor='black',hatch="\\",yerr=N_init_class_sd[1][it],error_kw=dict(elinewidth=2,ecolor='black'))
			
			rects2 = ax.bar(ind+2*width,N_init_class_m[2][it],width,color='grey',edgecolor='black',hatch="0",yerr=N_init_class_sd[2][it],error_kw=dict(elinewidth=2,ecolor='black'))	
			
			rects3 = ax.bar(ind+3*width,N_init_class_m[3][it],width,color='blue',edgecolor='black',hatch="+",yerr=N_init_class_sd[3][it],error_kw=dict(elinewidth=2,ecolor='black'))
			
			rects4 = ax.bar(ind+4*width,N_init_class_m[4][it],width,color='red',edgecolor='black',hatch="-",yerr=N_init_class_sd[4][it],error_kw=dict(elinewidth=2,ecolor='black'))		
				
			# axes and labels
			ax.set_ylabel('N',fontsize=18)
			ax.set_xlabel('Size',fontsize=18)
			ax.set_title(plottitle+'N Initial (Class) Year ' +str(it),fontsize=21)
			ax.set_xlim(-width,len(ind)+width)
			ax.set_ylim(0, max(np.max(N_init_class_m[0][it]),np.max(N_init_class_m[1][it]),np.max(N_init_class_m[2][it]),np.max(N_init_class_m[3][it]),np.max(N_init_class_m[4][it])))
			xTickMarks = [str(sizeclass[i]) for i in xrange(len(N_init_class_m[0][0]))]
			ax.set_xticks(ind+width)
			xtickNames = ax.set_xticklabels(xTickMarks)
			plt.setp(xtickNames, rotation=0)
			ax.legend((rects0[0],rects1[0],rects2[0],rects3[0],rects4[0]),label,loc=0)
		
		if batchno == 4:
		
			rects0 = ax.bar(ind,N_init_class_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_init_class_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

			rects1 = ax.bar(ind+width,N_init_class_m[1][it],width,color='blue',edgecolor='black',hatch="\\",yerr=N_init_class_sd[1][it],error_kw=dict(elinewidth=2,ecolor='black'))
			
			rects2 = ax.bar(ind+2*width,N_init_class_m[2][it],width,color='grey',edgecolor='black',hatch="0",yerr=N_init_class_sd[2][it],error_kw=dict(elinewidth=2,ecolor='black'))	
			
			rects3 = ax.bar(ind+3*width,N_init_class_m[3][it],width,color='blue',edgecolor='black',hatch="+",yerr=N_init_class_sd[3][it],error_kw=dict(elinewidth=2,ecolor='black'))	
				
			# axes and labels
			ax.set_ylabel('N',fontsize=18)
			ax.set_xlabel('Size',fontsize=18)
			ax.set_title(plottitle+'N Initial (Class) Year ' +str(it),fontsize=21)
			ax.set_xlim(-width,len(ind)+width)
			ax.set_ylim(0, max(np.max(N_init_class_m[0][it]),np.max(N_init_class_m[1][it]),np.max(N_init_class_m[2][it]),np.max(N_init_class_m[3][it])))
			xTickMarks = [str(sizeclass[i]) for i in xrange(len(N_init_class_m[0][0]))]
			ax.set_xticks(ind+width)
			xtickNames = ax.set_xticklabels(xTickMarks)
			plt.setp(xtickNames, rotation=0)
			ax.legend((rects0[0],rects1[0],rects2[0],rects3[0]),label,loc=0)
			
		if batchno == 3:
		
			rects0 = ax.bar(ind,N_init_class_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_init_class_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

			rects1 = ax.bar(ind+width,N_init_class_m[1][it],width,color='blue',edgecolor='black',hatch="\\",yerr=N_init_class_sd[1][it],error_kw=dict(elinewidth=2,ecolor='black'))
			
			rects2 = ax.bar(ind+2*width,N_init_class_m[2][it],width,color='grey',edgecolor='black',hatch="0",yerr=N_init_class_sd[2][it],error_kw=dict(elinewidth=2,ecolor='black'))		
				
			# axes and labels
			ax.set_ylabel('N',fontsize=18)
			ax.set_xlabel('Size',fontsize=18)
			ax.set_title(plottitle+'N Initial (Class) Year ' +str(it),fontsize=21)
			ax.set_xlim(-width,len(ind)+width)
			ax.set_ylim(0, max(np.max(N_init_class_m[0][it]),np.max(N_init_class_m[1][it]),np.max(N_init_class_m[2][it])))
			xTickMarks = [str(sizeclass[i]) for i in xrange(len(N_init_class_m[0][0]))]
			ax.set_xticks(ind+width)
			xtickNames = ax.set_xticklabels(xTickMarks)
			plt.setp(xtickNames, rotation=0)
			ax.legend((rects0[0],rects1[0],rects2[0]),label,loc=0)

		if batchno == 2:
		
			rects0 = ax.bar(ind,N_init_class_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_init_class_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

			rects1 = ax.bar(ind+width,N_init_class_m[1][it],width,color='blue',edgecolor='black',hatch="\\",yerr=N_init_class_sd[1][it],error_kw=dict(elinewidth=2,ecolor='black'))	
				
			# axes and labels
			ax.set_ylabel('N',fontsize=18)
			ax.set_xlabel('Size',fontsize=18)
			ax.set_title(plottitle+'N Initial (Class) Year ' +str(it),fontsize=21)
			ax.set_xlim(-width,len(ind)+width)
			ax.set_ylim(0, max(np.max(N_init_class_m[0][it]),np.max(N_init_class_m[1][it])))
			xTickMarks = [str(sizeclass[i]) for i in xrange(len(N_init_class_m[0][0]))]
			ax.set_xticks(ind+width)
			xtickNames = ax.set_xticklabels(xTickMarks)
			plt.setp(xtickNames, rotation=0)
			ax.legend((rects0[0],rects1[0]),label,loc=0)		
		
		savefig(dir+savename+'NInit_sizeclass_year_'+str(it)+'.png',dpi=savedpi)

# Plot Deaths: Packing
# ------------------------------------
figure()
for i in xrange(len(PackDeaths_Emi_m_total)):
	plot(nthfile,PackDeaths_Emi_m_total[i],linemarks_D1[i],label = label_D1[i],linewidth=2)
	plot(nthfile,PackDeaths_Immi_m_total[i],linemarks_D2[i],label = label_D2[i],linewidth=2)
xlabel('Time',fontsize=18)
ylabel('Packing Deaths',fontsize=18)
title(plottitle,fontsize=21)
#axis([-0.1,gen,0,30000])
legend(loc=0)	
savefig(dir+savename+'PackingDeaths.png',dpi=savedpi)
'''
# Bar plot for first few years
n_groups = 5
fig, ax = subplots()

index = np.arange(n_groups)
bar_width = 0.25

opacity = 1.0
error_config = dict(elinewidth=2,ecolor='black')
			
rects0 = bar(index, PackDeaths_Emi_m_total[0][index], bar_width,alpha=opacity,color='red',edgecolor='black',hatch="/",yerr=PackDeaths_Emi_sd_total[0][index],error_kw=error_config,label='1km Emi')

rects1 = bar(index + bar_width, PackDeaths_Emi_m_total[1][index], bar_width,alpha=opacity,color='blue',edgecolor='black',hatch="\\",yerr=PackDeaths_Emi_sd_total[1][index],error_kw=error_config,label='6km Emi')

rects2 = bar(index + 2*bar_width, PackDeaths_Emi_m_total[2][index], bar_width,alpha=opacity,color='grey',edgecolor='black',hatch="0",yerr=PackDeaths_Emi_sd_total[2][index],error_kw=error_config,label='Max Emi')

xlabel('Time')
ylabel('Deaths')
title(plottitle)
xticks(index + bar_width / 2, ('0', '1', '2', '3', '4'))
legend()

tight_layout()
savefig(dir+savename+'PackingDeathsEmi_BarGraphs.png',dpi=savedpi)

# Bar plot with same pars from above
fig, ax = subplots()
			
rects0 = bar(index, PackDeaths_Immi_m_total[0][index], bar_width,alpha=opacity,color='red',edgecolor='black',hatch="/",yerr=PackDeaths_Immi_sd_total[0][index],error_kw=error_config,label='1km Immi')

rects1 = bar(index + bar_width, PackDeaths_Immi_m_total[1][index], bar_width,alpha=opacity,color='blue',edgecolor='black',hatch="\\",yerr=PackDeaths_Immi_sd_total[1][index],error_kw=error_config,label='6km Immi')

rects2 = bar(index + 2*bar_width, PackDeaths_Immi_m_total[2][index], bar_width,alpha=opacity,color='grey',edgecolor='black',hatch="0",yerr=PackDeaths_Immi_sd_total[2][index],error_kw=error_config,label='Max Immi')

xlabel('Time')
ylabel('Deaths')
title(plottitle)
xticks(index + bar_width / 2, ('0', '1', '2', '3', '4'))
legend()

tight_layout()
savefig(dir+savename+'PackingDeathsImmi_BarGraphs.png',dpi=savedpi)
'''
# Plot Deaths: Move
# ------------------------------------
figure()
for i in xrange(len(PackDeaths_Emi_m_total)):
	plot(nthfile,MoveDeaths_Emi_m_total[i],linemarks_D1[i],label = label_D1[i],linewidth=2)
	plot(nthfile,MoveDeaths_Immi_m_total[i],linemarks_D2[i],label = label_D2[i],linewidth=2)
xlabel('Time',fontsize=18)
ylabel('Move Deaths',fontsize=18)
title(plottitle,fontsize=21)
#axis([-0.1,gen,0,5000])
legend(loc=0)	
savefig(dir+savename+'MoveDeaths.png',dpi=savedpi)

# Plot Emi/Immi Totals
# -----------------------
figure()
for i in xrange(len(PackDeaths_Emi_m_total)):
	plot(nthfile,N_afterEmi_m_total[i],linemarks_D1[i],label=label_D1[i],linewidth=2)
	plot(nthfile,N_afterImmi_m_total[i],linemarks_D2[i],label=label_D2[i],linewidth=2)
xlabel('Time',fontsize=18)
ylabel('N',fontsize=18)
title(plottitle,fontsize=21)
#axis([-0.1,gen,0,70000])
legend(loc=0)	
savefig(dir+savename+'NEmi_NImmi_Totals.png',dpi=savedpi)	
'''
# Bar plot with same pars from above
fig, ax = subplots()
			
rects0 = bar(index, N_afterEmi_m_total[0][index], bar_width,alpha=opacity,color='red',edgecolor='black',hatch="/",yerr=N_afterEmi_sd_total[0][index],error_kw=error_config,label='1km Emi')

rects1 = bar(index + bar_width, N_afterEmi_m_total[1][index], bar_width,alpha=opacity,color='blue',edgecolor='black',hatch="\\",yerr=N_afterEmi_sd_total[1][index],error_kw=error_config,label='6km Emi')

rects2 = bar(index + 2*bar_width, N_afterEmi_m_total[2][index], bar_width,alpha=opacity,color='grey',edgecolor='black',hatch="0",yerr=N_afterEmi_sd_total[2][index],error_kw=error_config,label='Max Emi')

xlabel('Time')
ylabel('N')
title(plottitle)
xticks(index + bar_width / 2, ('0', '1', '2', '3', '4'))
legend()

tight_layout()
savefig(dir+savename+'NAfterEmi_BarGraphs.png',dpi=savedpi)

fig, ax = subplots()
			
rects0 = bar(index, N_afterImmi_m_total[0][index], bar_width,alpha=opacity,color='red',edgecolor='black',hatch="/",yerr=N_afterImmi_sd_total[0][index],error_kw=error_config,label='1km Immi')

rects1 = bar(index + bar_width, N_afterImmi_m_total[1][index], bar_width,alpha=opacity,color='blue',edgecolor='black',hatch="\\",yerr=N_afterImmi_sd_total[1][index],error_kw=error_config,label='6km Immi')

rects2 = bar(index + 2*bar_width, N_afterImmi_m_total[2][index], bar_width,alpha=opacity,color='grey',edgecolor='black',hatch="0",yerr=N_afterImmi_sd_total[2][index],error_kw=error_config,label='Max Immi')

xlabel('Time')
ylabel('N')
title(plottitle)
xticks(index + bar_width / 2, ('0', '1', '2', '3', '4'))
legend()

tight_layout()
savefig(dir+savename+'NAfterImmi_BarGraphs.png',dpi=savedpi)
'''
show()




