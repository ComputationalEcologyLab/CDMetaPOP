# -----------------------------------------------------------------------------
# population.py
# 2016 11 28: Erin Landguth
# This script grabs summary_popAllTime.csv files and plots population numbers
# _batch_v1.08: Update for recent version. Add in different populaiton plots:
#	Initial N(when they are back), Capture Back, N when they are out, Capture
# 	Out
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
dir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2016Sampling/WCT2016Nov22_noGenes2x2_6kmMove_modKv3/"
plottitle = ''
savename = "_Diagnostics_WCT2016Nov22_noGenes2x2_6kmMove_modKv3"
#label = ['100years', '200years', '300years']
label = ['Back','Captured Back','Out','Captured Out']
label2 = ['Back','Out']
batchno = 1
#linemarks = ['k--','b-o','r-','g-^','ys-']
linemarks = ['k','b','r','g','y']
linemarks_lr = ['k--','b--','r--','--.','--.']

outdir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2016Sampling/WCT2016Nov22_noGenes2x2_6kmMove_modKv3/summary/"
xyfile = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2016Sampling/PatchVars_WCT1384_2016Nov_Sampling_withGenes.csv"

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 125 # Number of years 
nthfile = range(0,gen,1)
#nthfile = np.asarraty([0,1,2,3,4,5,10,19])
mcno = 3 # Number of MCs
plottime = np.asarray([124])
Naxes = [-0.01,125,0,100000]

# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# ---------------------------------
# Storage variables
# ---------------------------------
# Store values - Patch
N_init_patch = [] # Each patch value when back
N_init_pop = [] # total population in patch when back
N_cap_back_patch = [] # Patch capture when back at spawn grounds
N_cap_back_pop = [] # total population captured when back
N_out_patch = [] # Each patch value when out
N_out_pop = [] # total population when out
N_cap_out_patch = [] # patch capture when out
N_cap_out_pop = [] # total population captured when out

# Loop through batches
for ibatch in xrange(batchno):

	# Add storage spot
	N_init_patch.append([])
	N_init_pop.append([])
	N_cap_back_patch.append([])
	N_cap_back_pop.append([])
	N_out_patch.append([])
	N_out_pop.append([])
	N_cap_out_patch.append([])
	N_cap_out_pop.append([])

	# Loop through MCs
	for imc in xrange(mcno):
	
		# Store values
		N_init_patch[ibatch].append([])
		N_init_pop[ibatch].append([])
		N_cap_back_patch[ibatch].append([])
		N_cap_back_pop[ibatch].append([])
		N_out_patch[ibatch].append([])
		N_out_pop[ibatch].append([])
		N_cap_out_patch[ibatch].append([])
		N_cap_out_pop[ibatch].append([])
		
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
		
		# ----------------------------------
		# Then Loop through generations/time
		# ----------------------------------
		for iout in xrange(gen):
		
			N_init_patch[ibatch][imc].append([])
			N_init_pop[ibatch][imc].append([])
			N_cap_back_patch[ibatch][imc].append([])
			N_cap_back_pop[ibatch][imc].append([])
			N_out_patch[ibatch][imc].append([])
			N_out_pop[ibatch][imc].append([])
			N_cap_out_patch[ibatch][imc].append([])
			N_cap_out_pop[ibatch][imc].append([])			
			
			# Grab all patch values - patch values with total
			for j in xrange(1,len(values_pop[1+iout][3].split('|'))-1):
				N_init_patch[ibatch][imc][iout].append(float(values_pop[1+iout][3].split('|')[j]))
				N_out_patch[ibatch][imc][iout].append(float(values_pop[1+iout][19].split('|')[j]))				
			# Sum totals
			N_init_pop[ibatch][imc][iout] = sum(N_init_patch[ibatch][imc][iout])
			N_out_pop[ibatch][imc][iout] = sum(N_out_patch[ibatch][imc][iout])
			
			# Grab all patch values - patch values withOUT total
			for j in xrange(len(values_pop[1+iout][15].split('|'))-1):
				N_cap_back_patch[ibatch][imc][iout].append(float(values_pop[1+iout][15].split('|')[j]))
				N_cap_out_patch[ibatch][imc][iout].append(float(values_pop[1+iout][22].split('|')[j]))
			# Sum them here
			N_cap_back_pop[ibatch][imc][iout] = np.nansum(N_cap_back_patch[ibatch][imc][iout])
			N_cap_out_pop[ibatch][imc][iout] = np.nansum(N_cap_out_patch[ibatch][imc][iout])
		
# Turn into arrays
N_init_patch = np.asarray(N_init_patch)
N_init_pop = np.asarray(N_init_pop)
N_cap_back_patch = np.asarray(N_cap_back_patch)
N_cap_back_pop = np.asarray(N_cap_back_pop)
N_out_patch = np.asarray(N_out_patch)
N_out_pop = np.asarray(N_out_pop)
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

# Patch numbers [batch][mcrun][time][patch] - when back
N_init_patch_m = np.nansum(N_init_patch,axis=1)/mcno
N_init_patch_sd = np.std(N_init_patch,axis=1)	
error = qnorm*N_init_patch_sd/(mcno)
N_init_patch_l = N_init_patch_m-error
N_init_patch_r = 	N_init_patch_m+error
N_init_patch_min = np.min(N_init_patch,axis=1)
N_init_patch_max = np.max(N_init_patch,axis=1)
# Total population
N_init_pop_m = np.nansum(N_init_pop,axis=1)/mcno
N_init_pop_sd = np.std(N_init_pop,axis=1)	
error = qnorm*N_init_pop_sd/(mcno)
N_init_pop_l = N_init_pop_m-error
N_init_pop_r = 	N_init_pop_m+error
N_init_pop_min = np.min(N_init_pop,axis=1)
N_init_pop_max = np.max(N_init_pop,axis=1)

# Patch numbers [batch][mcrun][time][patch] - when out
N_out_patch_m = np.nansum(N_out_patch,axis=1)/mcno
N_out_patch_sd = np.std(N_out_patch,axis=1)	
error = qnorm*N_out_patch_sd/(mcno)
N_out_patch_l = N_out_patch_m-error
N_out_patch_r = 	N_out_patch_m+error
N_out_patch_min = np.min(N_out_patch,axis=1)
N_out_patch_max = np.max(N_out_patch,axis=1)
# Total population
N_out_pop_m = np.nansum(N_out_pop,axis=1)/mcno
N_out_pop_sd = np.std(N_out_pop,axis=1)	
error = qnorm*N_out_pop_sd/(mcno)
N_out_pop_l = N_out_pop_m-error
N_out_pop_r = 	N_out_pop_m+error
N_out_pop_min = np.min(N_out_pop,axis=1)
N_out_pop_max = np.max(N_out_pop,axis=1)

# --------------------------------
# Other summary data 
# -------------------------------- 
# Capture proportion - back population
N_cap_back_pop_prop_m = N_cap_back_pop_m/N_init_pop_m
N_cap_back_pop_prop_l = N_cap_back_pop_l/N_init_pop_l
N_cap_back_pop_prop_r = N_cap_back_pop_r/N_init_pop_r
# Capture proportion - out population
N_cap_out_pop_prop_m = N_cap_out_pop_m/N_out_pop_m
N_cap_out_pop_prop_l = N_cap_out_pop_l/N_out_pop_l
N_cap_out_pop_prop_r = N_cap_out_pop_r/N_out_pop_r

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

# Create file to write info to
for itime in xrange(len(plottime)):
	outputfile = open(outdir+savename+'mean_patch_pops_time'+str(plottime[itime])+'.csv','w')
	# Write out the titles
	outputfile.write('Subpopulation,X,Y,Back,Back_Captured,Out,Out_Captured\n')	
	# WRite information
	for i in xrange(len(X)):
		outputfile.write(SubPop[i]+',')
		outputfile.write(X[i]+',')
		outputfile.write(Y[i]+',')
		for j in xrange(batchno):
			outputfile.write(str(N_init_patch_m[j][plottime[itime]][i])+',')
			outputfile.write(str(N_cap_back_patch_m[j][plottime[itime]][i])+',')
			outputfile.write(str(N_out_patch_m[j][plottime[itime]][i])+',')
			outputfile.write(str(N_cap_out_patch_m[j][plottime[itime]][i])+',')
		# To get return character on the end
		outputfile.write('\n')		
	outputfile.close()		

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------

# Plot total N
# ------------
figure(1)
for i in xrange(len(N_init_pop_m)):
	plot(nthfile,N_init_pop_m[i],linemarks[i],linewidth=2)
	plot(nthfile,N_init_pop_l[i],linemarks_lr[i],linewidth=1)
	plot(nthfile,N_init_pop_r[i],linemarks_lr[i],linewidth=1)
	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Population',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes)
legend(loc=0)
savefig(outdir+savename+'BackPopulation_Total.png',dpi=savedpi)

# Plot all Ns
# -----------
figure(2)
plot(nthfile,N_init_pop_m[0],linemarks[0],label=label[0],linewidth=2)
plot(nthfile,N_cap_back_pop_m[0],linemarks[1],label=label[1],linewidth=2)
plot(nthfile,N_out_pop_m[0],linemarks[2],label=label[2],linewidth=2)
plot(nthfile,N_cap_out_pop_m[0],linemarks[3],label=label[3],linewidth=2)
xlabel('Time',fontsize=18)
ylabel('Population',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes)
legend(loc=0)
savefig(outdir+savename+'AllPopulations.png',dpi=savedpi)
	
# Plot capture proportions
# -------------
figure(3)
plot(nthfile,N_cap_back_pop_prop_m[0],linemarks[0],label = label2[0],linewidth=2)
plot(nthfile,N_cap_out_pop_prop_m[0],linemarks[1],label = label2[1],linewidth=2)	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Capture',fontsize=18)
title(plottitle,fontsize=21)
axis([-0.01,125,0,1])
legend(loc=0)
savefig(outdir+savename+'NCapture_Proportion.png',dpi=savedpi)


show()




