# -----------------------------------------------------------------------------
# plotRSI.py
# 2015 12 29: Erin Landguth
# This script grabs output.csv files; plots RSI
# v1.00 - For v1.00 files.
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
#dir = 'D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v1.00_20151130/All_MovementDistances/'
plottitle = ""
savename = "_RSI_movethreshold_"
label = ['25% Max', '50% Max', '100% Max']
linemarks = ['k--','b-o','r-.','g-^','ys-']
batchno = 3

dir = 'D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v1.00_20151130/All_StrayRates/'
plottitle = ""
savename = "_RSI_strayrates_"
label = ['0.001', '0.01', '0.05']
linemarks = ['k--','b-o','r-.','g-^','ys-']
batchno = 3

dir = 'D:/projects/CDmetaPOP/Seattle/Runs/dataWCT1384_v1.00_20151130/All_Landscapes/'
savename = "_RSI_landscapes_"
label = ['Riverine', 'Ex-Barrier', 'Fut-Barrier','Remove-Barrier']
batchno = 4
linemarks = ['k--','b-o','r-','g-^','ys-']

outdir = dir

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 125 # Number of years 
nthfile = range(0,gen,1)
#nthfile = np.asarraty([0,1,2,3,4,5,10,19])
mcno = 2 # Number of MCs

# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# ---------------------------------
# Storage variables
# ---------------------------------
# Store values 
R = []
S1 = []
S2 = []
I = []
N = []

# Loop through batches
for ibatch in xrange(batchno):

	# Add storage spot
	R.append([])
	S1.append([])
	S2.append([])
	I.append([])
	N.append([])
	
	# Loop through MCs
	for imc in xrange(mcno):
	
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
		
		# Store values
		R[ibatch].append([])
		S1[ibatch].append([])
		S2[ibatch].append([])
		I[ibatch].append([])
		N[ibatch].append([])
				
		# Then Loop through generations/time
		for iout in xrange(gen):
		
			R[ibatch][imc].append([])
			S1[ibatch][imc].append([])
			S2[ibatch][imc].append([])
			I[ibatch][imc].append([])
			N[ibatch][imc].append([])
						
			# Grab all patch values - patch values with total
			for j in xrange(len(values_pop[1+iout][57].split('|'))-1):
				R[ibatch][imc][iout].append(float(values_pop[1+iout][57].split('|')[j])) 
				S1[ibatch][imc][iout].append(float(values_pop[1+iout][58].split('|')[j]))
				S2[ibatch][imc][iout].append(float(values_pop[1+iout][59].split('|')[j]))
				I[ibatch][imc][iout].append(float(values_pop[1+iout][60].split('|')[j]))
				N[ibatch][imc][iout].append(float(values_pop[1+iout][3].split('|')[j+1])) # Remove total here

# Turn into arrays
R = np.asarray(R)
S1 = np.asarray(S1)
S2 = np.asarray(S2)
I = np.asarray(I)
N = np.asarray(N)

# --------------------------------------------
# Get mean over Monte Carlosfor each batch run
# --------------------------------------------
N_m = np.nansum(N[:][:],axis=1)/mcno

R_m = np.nansum(R[:][:],axis=1)/mcno
R_sd = np.std(R[:][:],axis=1)	
error = R_sd*qnorm/mcno
R_l = R_m-error
R_r = 	R_m+error
R_min = np.min(R[:][:],axis=1)
R_max = np.max(R[:][:],axis=1)

S1_m = np.nansum(S1[:][:],axis=1)/mcno
S1_sd = np.std(S1[:][:],axis=1)	
error = S1_sd*qnorm/mcno
S1_l = S1_m-error
S1_r = 	S1_m+error
S1_min = np.min(S1[:][:],axis=1)
S1_max = np.max(S1[:][:],axis=1)

S2_m = np.nansum(S2[:][:],axis=1)/mcno
S2_sd = np.std(S2[:][:],axis=1)	
error = qnorm*S2_sd/(mcno)
S2_l = S2_m-error
S2_r = 	S2_m+error
S2_min = np.min(S2[:][:],axis=1)
S2_max = np.max(S2[:][:],axis=1)

I_m = np.nansum(I[:][:],axis=1)/mcno
I_sd = np.std(I[:][:],axis=1)	
error = I_sd*qnorm/mcno
I_l = I_m-error
I_r = 	I_m+error
I_min = np.min(I[:][:],axis=1)
I_max = np.max(I[:][:],axis=1)

# --------------------------------
# Patch mean
# -------------------------------- 

# First remove N=0 locations for averaging
# ----------------------------------------

remove0 = np.where(N_m[0][0] != 0)[0]
#remove0 = zip(*np.where(N_m != 0))

R_patch_m = np.nansum(R_m,axis=2)/len(remove0)
S1_patch_m = np.nansum(S1_m,axis=2)/len(remove0)
S2_patch_m = np.nansum(S2_m,axis=2)/len(remove0)
I_patch_m = np.nansum(I_m,axis=2)/len(remove0)
N_patch_m = np.nansum(N_m,axis=2)/len(remove0)


# --------------------------------------------
# Write to file
# --------------------------------------------

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------

# Plot total RSI through time
# ---------------------------
figure()
for i in xrange(len(R_m)):
	plot(nthfile,R_patch_m[i],linemarks[i],label=label[i],linewidth=2)
	#plot(nthfile,S1_m[i],linemarks[i],label='S1 '+label[i],linewidth=2)
	#plot(nthfile,S2_m[i],linemarks[i],label='S2 '+label[i],linewidth=2)
	#plot(nthfile,I_m[i],linemarks[i],label='I '+label[i],linewidth=2)
	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Residents (Patch Mean)',fontsize=18)
title(plottitle,fontsize=21)
axis([-0.01,gen,0,np.max(R_patch_m)])
legend(loc=0)
savefig(dir+savename+'PatchMean_R.png',dpi=savedpi)


figure()
for i in xrange(len(R_m)):
	#plot(nthfile,R_patch_m[i],linemarks[i],label='R '+label[i],linewidth=2)
	plot(nthfile,S1_patch_m[i],linemarks[i],label=label[i],linewidth=2)
	#plot(nthfile,S2_m[i],linemarks[i],label='S2 '+label[i],linewidth=2)
	#plot(nthfile,I_m[i],linemarks[i],label='I '+label[i],linewidth=2)
	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Strayers (Patch Mean)',fontsize=18)
title(plottitle,fontsize=21)
axis([-0.01,gen,0,np.max(S1_patch_m)])
legend(loc=0)
savefig(dir+savename+'PatchMean_S1.png',dpi=savedpi)

figure()
for i in xrange(len(R_m)):
	#plot(nthfile,R_patch_m[i],linemarks[i],label='R '+label[i],linewidth=2)
	#plot(nthfile,S1_m[i],linemarks[i],label='S1 '+label[i],linewidth=2)
	plot(nthfile,S2_patch_m[i],linemarks[i],label=label[i],linewidth=2)
	#plot(nthfile,I_m[i],linemarks[i],label='I '+label[i],linewidth=2)
	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Forced Strayers (Patch Mean)',fontsize=18)
title(plottitle,fontsize=21)
axis([-0.01,gen,0,np.max(S2_patch_m)])
legend(loc=0)
savefig(dir+savename+'PatchMean_S2.png',dpi=savedpi)

figure()
for i in xrange(len(R_m)):
	#plot(nthfile,R_patch_m[i],linemarks[i],label='R '+label[i],linewidth=2)
	#plot(nthfile,S1_m[i],linemarks[i],label='S1 '+label[i],linewidth=2)
	#plot(nthfile,S2_m[i],linemarks[i],label='S2 '+label[i],linewidth=2)
	plot(nthfile,I_patch_m[i],linemarks[i],label=label[i],linewidth=2)
	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Immigrants (Patch Mean)',fontsize=18)
title(plottitle,fontsize=21)
axis([-0.01,gen,0,np.max(I_patch_m)])
legend(loc=0)
savefig(dir+savename+'PatchMean_I.png',dpi=savedpi)

figure()
for i in xrange(len(R_m)):
	plot(nthfile,N_patch_m[i],linemarks[i],label=label[i],linewidth=2)
	#plot(nthfile,S1_m[i],linemarks[i],label='S1 '+label[i],linewidth=2)
	#plot(nthfile,S2_m[i],linemarks[i],label='S2 '+label[i],linewidth=2)
	#plot(nthfile,I_m[i],linemarks[i],label='I '+label[i],linewidth=2)
	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Total N (Patch Mean)',fontsize=18)
title(plottitle,fontsize=21)
axis([-0.01,gen,0,np.max(N_patch_m)])
legend(loc=0)
savefig(dir+savename+'PatchMean_N.png',dpi=savedpi)

show()




