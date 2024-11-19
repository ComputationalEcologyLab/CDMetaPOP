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
# _batch_v1.08: Update for recent version.
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
#dir = "D:/projects/CDmetaPOP/Seattle/Sampling/RivExFutBarr_100max_Straypt01_randomgenes/"
dir = "D:/projects/CDmetaPOP/Seattle/Runs/data_WCT1384_2016Sampling/Compare_modKv3/"
plottitle = ''
savename = "_Diagnostics_modKv3_"
label = ['6kmMove', 'MaxMove']
batchno = 2
linemarks = ['k--','b-o','r-','g-^','ys-']

outdir = dir

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 125 # Number of years 
nthfile = range(0,gen,1)
#nthfile = np.asarraty([0,1,2,3,4,5,10,19])
mcno = 3 # Number of MCs
plottime = np.asarray([100])
plotagesize = 'Y' # Plot the time specific age/size information?
#label = ['0%','25%','50%','75%','100%']

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
N_init_patch = []

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
			
			# Age split
			for j in xrange(len(values[1+iout][2].split('|'))-1):
				N_init_age[ibatch][imc][iout].append(int(values[1+iout][2].split('|')[j])) 
				N_init_class[ibatch][imc][iout].append(int(values[1+iout][6].split('|')[j]))
				size_age[ibatch][imc][iout].append(float(values[1+iout][3].split('|')[j])) 
				if values[1+iout][12].split('|')[j] == 'NA':
					N_afterEmi[ibatch][imc][iout].append(nan)
				else:
					N_afterEmi[ibatch][imc][iout].append(int(values[1+iout][12].split('|')[j]))
				if values[1+iout][18].split('|')[j] == 'NA':
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

# --------------------------------------------
# Get mean over Monte Carlosfor each batch run
# --------------------------------------------
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

N_afterImmi_m = np.nansum(N_afterImmi[:][:],axis=1)/mcno
N_afterImmi_sd = np.std(N_afterImmi[:][:],axis=1)	
error = N_afterImmi_sd*qnorm/mcno
N_afterImmi_l = N_afterImmi_m-error
N_afterImmi_r = 	N_afterImmi_m+error
N_afterImmi_min = np.min(N_afterImmi[:][:],axis=1)
N_afterImmi_max = np.max(N_afterImmi[:][:],axis=1)

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

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------

# Plot total N
# ------------
figure()
for i in xrange(len(N_init_pop_m)):
	plot(nthfile,N_init_pop_m[i],linemarks[i],label=label[i],linewidth=2)
	
#fill_between(nthfile, N_init_pop_m[0], N_init_pop_m[4])
xlabel('Time',fontsize=18)
ylabel('Population',fontsize=18)
title(plottitle,fontsize=21)
#axis([-0.1,gen,np.min(N_init_pop_m),np.max(N_init_pop_m)])
#axis([-0.1,gen,50000,150000])
#axis([-0.01,130,40000,120000])
axis([-0.01,130,0,140000])
legend(loc=0)
savefig(dir+savename+'NInit_pop.png',dpi=savedpi)

# Plot patch mean
# -------------------------------
figure()
for i in xrange(len(N_growth_patch_m)):
	plot(nthfile[1:],N_growth_patch_m[i],linemarks[i],label=label[i],linewidth=2)
xlabel('Time',fontsize=18)
ylabel('Growth rate patch mean',fontsize=18)
title(plottitle,fontsize=21)
axis([-0.1,gen,0,2.0])
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




show()




