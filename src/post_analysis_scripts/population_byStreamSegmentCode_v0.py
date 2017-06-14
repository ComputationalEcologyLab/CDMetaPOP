# -----------------------------------------------------------------------------
# population_byStreamSegmentCode.py
# April 2017: Erin Landguth
# This script grabs summary_popAllTime.csv files and plots population numbers
# by defined stream segments that are coded to PatchVars file. 
# v0: Initial script
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
	
# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
#End::count_unique()

# ---------
# User info
# ---------
dir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2017March_HarvestResponse/HarvestModels/Compare/"
dir = 'D:/projects/CDmetaPOP/Seattle/Runs/KZeller_HarvestSims/postprocessing/'
plottitle = ''
savename = "_HarvestModelsCompare_"
savename = '_testforKZ_'
mcno = 3 # Number of MCs
batchno = 1 # Number of batches to match label and linemwarks next. 
# Note label and line marks can be larger than batchno. 
label = ['E/K_Return', 'NoE/K_Return', 'E/K_Stay','NoE/K_Stay']
linemarks = ['k','b','r','g','y']
linemarks_lr = ['k--','b--','r--','g--','y--']
linemarks2 = ['k-o','b-o','r-o','g-o','y-o']

# Create a folder to store results
outdir = dir+"summary/"
xyfile = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2017March_HarvestResponse/PatchVars_WCT1384_2016Nov_Sampling_withoutGenes_modKv3_withSegmentCode.csv"

savedpi = 300 # dpi to save
gen = 150 # Number of years
timeofHarvest = 21 

Naxes = [-0.01,gen,0,80000] # Total population plots
Naxes2 = [-0.01,gen,0,40000] # Subpop/Sullivan plots
Naxes3 = [-0.01,gen,0,10000] # Slate and flume plots

# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

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

# ---------------------------------
# Storage variables
# ---------------------------------
# Store values - Patch
N_init_patch = [] # Each patch value when back
N_init_pop = [] # total population in patch when back
N_out_patch = [] # Each patch value when out
N_out_pop = [] # total population when out

# Loop through batches
for ibatch in xrange(batchno):

	# Add storage spot
	N_init_patch.append([])
	N_init_pop.append([])
	N_out_patch.append([])
	N_out_pop.append([])
	
	# Loop through MCs
	for imc in xrange(mcno):
	
		# Store values
		N_init_patch[ibatch].append([])
		N_init_pop[ibatch].append([])
		N_out_patch[ibatch].append([])
		N_out_pop[ibatch].append([])
				
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
			N_out_patch[ibatch][imc].append([])
			N_out_pop[ibatch][imc].append([])
						
			# Grab all patch values - patch values with total
			for j in xrange(1,len(values_pop[1+iout][3].split('|'))-1):
				N_init_patch[ibatch][imc][iout].append(float(values_pop[1+iout][3].split('|')[j]))
				N_out_patch[ibatch][imc][iout].append(float(values_pop[1+iout][19].split('|')[j]))
							
			# Sum totals
			N_init_pop[ibatch][imc][iout] = sum(N_init_patch[ibatch][imc][iout])
			N_out_pop[ibatch][imc][iout] = sum(N_out_patch[ibatch][imc][iout])
			
# Turn into arrays
N_init_patch = np.asarray(N_init_patch)
N_init_pop = np.asarray(N_init_pop)
N_out_patch = np.asarray(N_out_patch)
N_out_pop = np.asarray(N_out_pop)

# --------------------------------------------
# Get mean over Monte Carlosfor each batch run
# --------------------------------------------
qnorm = 1.959964 # For CIs, not in function 
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

# Get total population by SubPop
SubPop = np.asarray(SubPop,dtype=int)
uni_SubPop = count_unique(SubPop)

# Storage vectors 
N_out_subpop_m = []
N_init_subpop_m = []
N_out_subpop_l = []
N_init_subpop_l = []
N_out_subpop_r = []
N_init_subpop_r = []

# For each unique reach, loop through patches and sum
for iseg in xrange(len(uni_SubPop[0])):
	# Get index locations for each segment
	this_index = np.where(uni_SubPop[0][iseg] == SubPop)[0]
	# For OUT patches
	N_out_subpop_m.append(np.nansum(N_out_patch_m[:,:,this_index],axis=2))
	N_out_subpop_l.append(np.nansum(N_out_patch_l[:,:,this_index],axis=2))
	N_out_subpop_r.append(np.nansum(N_out_patch_r[:,:,this_index],axis=2))
	# For BACK patches
	N_init_subpop_m.append(np.nansum(N_init_patch_m[:,:,this_index],axis=2))
	N_init_subpop_l.append(np.nansum(N_init_patch_l[:,:,this_index],axis=2))
	N_init_subpop_r.append(np.nansum(N_init_patch_r[:,:,this_index],axis=2))

# --------------------------------------------
# Write to file
# --------------------------------------------

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------
nthfile = range(0,gen,1)
# Plot total N - back
# ------------
figure(1)
for i in xrange(len(N_init_pop_m)):
	plot(nthfile,N_init_pop_m[i],linemarks[i],linewidth=2,label=label[i])
	plot(nthfile,N_init_pop_l[i],linemarks_lr[i],linewidth=1)
	plot(nthfile,N_init_pop_r[i],linemarks_lr[i],linewidth=1)
axvline(x=timeofHarvest, color='k', linestyle='--')
	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Total Back Population',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes)
legend(loc=0)
savefig(outdir+savename+'BackPopulation_Total.png',dpi=savedpi)

# Plot total N - out
# -----------
figure(2)
for i in xrange(len(N_init_pop_m)):
	plot(nthfile,N_out_pop_m[i],linemarks[i],linewidth=2,label=label[i])
	plot(nthfile,N_out_pop_l[i],linemarks_lr[i],linewidth=1)
	plot(nthfile,N_out_pop_r[i],linemarks_lr[i],linewidth=1)
axvline(x=timeofHarvest, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Total Out Population',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes)
legend(loc=0)
savefig(outdir+savename+'OutPopulation_Total.png',dpi=savedpi)
	
# Plot subpopulation totals - all
# -------------------------------
figure(3)
for j in xrange(len(uni_SubPop[0])):
	plot(nthfile,N_init_subpop_m[j][0],label=str(j))
axvline(x=timeofHarvest, color='k', linestyle='--')
axvline(x=timeofHarvest, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Total SubPopulations Model 1',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes2)
legend(loc=1)
savefig(outdir+savename+'SubPopulations_Total_Model1.png',dpi=savedpi)	
figure(4)
for j in xrange(len(uni_SubPop[0])):
	plot(nthfile,N_init_subpop_m[j][1],label=str(j))
axvline(x=timeofHarvest, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Total SubPopulations Model 2',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes2)
legend(loc=1)
savefig(outdir+savename+'SubPopulations_Total_Model2.png',dpi=savedpi)
figure(5)
for j in xrange(len(uni_SubPop[0])):
	plot(nthfile,N_init_subpop_m[j][2])
axvline(x=timeofHarvest, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Total SubPopulations Model 3',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes2)
legend(loc=1)
savefig(outdir+savename+'SubPopulations_Total_Model3.png',dpi=savedpi)

'''
figure(6)
for j in xrange(len(uni_SubPop[0])):
	plot(nthfile,N_init_subpop_m[j][3])
axvline(x=timeofHarvest, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Total SubPopulations Model 4',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes2)
legend(loc=1)
savefig(outdir+savename+'SubPopulations_Total_Model4.png',dpi=savedpi)
'''
# Plot Sullivan and Slate, Flume only
# -----------------------------
figure(7)
for i in xrange(batchno):
	plot(nthfile,N_init_subpop_m[4][i],linemarks[i],label=label[i])
axvline(x=timeofHarvest, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Sullivan Population',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes2)
legend(loc=0)
savefig(outdir+savename+'Sullivan_AllModels.png',dpi=savedpi)
figure(8)
for i in xrange(batchno):
	plot(nthfile,N_init_subpop_m[2][i],linemarks[i],label=label[i])
axvline(x=timeofHarvest, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Slate Population',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes2)
legend(loc=0)
savefig(outdir+savename+'Slate_AllModels.png',dpi=savedpi)
figure(9)
for i in xrange(batchno):
	plot(nthfile,N_init_subpop_m[3][i],linemarks[i],label=label[i])
axvline(x=timeofHarvest, color='k', linestyle='--')
xlabel('Time',fontsize=18)
ylabel('Flume Population',fontsize=18)
title(plottitle,fontsize=21)
axis(Naxes2)
legend(loc=0)
savefig(outdir+savename+'Flume_AllModels.png',dpi=savedpi)	
show()




