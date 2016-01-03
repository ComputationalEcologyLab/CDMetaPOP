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
# v8 - 2014 10 20: Make general for a run not comparing to batches.
# v095 - 2015 01 14: Update for new format summary_class{year}.csv and summary_pop{year}.csv options.
# class_v095: plots class/age information for given time. Plot N_Age and Size_age for starters
# class_v096: option plots for different N's in file - Ns after migraiton, Ns after growth.
# class_v098: Update script for newest version. 
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
# Directory locations of folders
dir = "D:/projects/CDmetaPOP/Seattle/Runs/dataWCT406_v098_Master_20150306/WCT406_IBL_1425756970/"
plottitle = 'WCT406_IBL '
savename = 'WCT406_IBL'
outdir = dir
batchno = 0

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 110 # Number of years 
nthfile = range(0,gen,1)
#nthfile = np.asarraty([0,1,2,3,4,5,10,19])
maxA = 19*36 # loci * alleles
mcno = 1 # Number of MCs
plottime = np.asarray([10,100])
# Plot other N information Y or N
plototherN = 'Y'

# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# ---------------------------------
# Storage variables
# ---------------------------------
# Store values
N_init_age = []
N_init_class = []
size_age = []
N_afterEmi = []
N_afterImmi = []
N_growthback = []
N_growthout = []


# Loop through batches
for ibatch in xrange(1):

	# Add storage spot
	N_init_age.append([])
	N_init_class.append([])
	size_age.append([])
	N_afterEmi.append([])
	N_afterImmi.append([])
	N_growthback.append([])
	N_growthout.append([])
	
	# Loop through MCs
	for imc in xrange(mcno):
	
		# Open output.csv file in folder
		inputfile = open(dir+'batchrun'+str(batchno)+'mcrun'+str(imc)+'/summary_classAllTime.csv')
		
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
		
		# Then Loop through generations/time
		for iout in xrange(gen):
		
			N_init_age[ibatch][imc].append([])
			N_init_class[ibatch][imc].append([])
			size_age[ibatch][imc].append([])
			N_afterEmi[ibatch][imc].append([])
			N_afterImmi[ibatch][imc].append([])
			N_growthback[ibatch][imc].append([])
			N_growthout[ibatch][imc].append([])
			
			# Age split
			for j in xrange(len(values[1+iout][1].split('|'))-1):
				N_init_age[ibatch][imc][iout].append(int(values[1+iout][2].split('|')[j])) 
				N_init_class[ibatch][imc][iout].append(int(values[1+iout][6].split('|')[j]))
				size_age[ibatch][imc][iout].append(float(values[1+iout][3].split('|')[j])) 	
				N_afterEmi[ibatch][imc][iout].append(int(values[1+iout][13].split('|')[j]))
				N_afterImmi[ibatch][imc][iout].append(int(values[1+iout][20].split('|')[j]))
				N_growthback[ibatch][imc][iout].append(int(values[1+iout][10].split('|')[j]))
				N_growthout[ibatch][imc][iout].append(int(values[1+iout][17].split('|')[j]))

# Turn into arrays
N_init_age = np.asarray(N_init_age)
N_init_class = np.asarray(N_init_class)
size_age = np.asarray(size_age)
N_afterEmi = np.asarray(N_afterEmi)
N_afterImmi = np.asarray(N_afterImmi)
N_growthback = np.asarray(N_growthback)
N_growthout = np.asarray(N_growthout)

# --------------------------------------------
# Get mean over Monte Carlos for each batch run
# --------------------------------------------
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

# --------------------------------
# Other summary data 
# -------------------------------- 

# --------------------------------------------
# Write to file
# --------------------------------------------

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------

# Plot N ages at given plottime
# -----------------------------
ind = np.arange(len(N_init_age_m[0][0])) # the x locations for the groups
width = 0.25                      # the width of the bars

# Loop through each year to plot
for it in plottime:
	
	TotalN = np.nansum(N_init_age_m[0][it])
	
	fig = plt.figure()
	ax = fig.add_subplot(111)

	rects = ax.bar(ind+width,N_init_age_m[0][it],width,color='grey',edgecolor='black',hatch="//",yerr=N_init_age_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

	# axes and labels
	ax.set_ylabel('N = '+str(TotalN),fontsize=18)
	ax.set_xlabel('Age',fontsize=18)
	ax.set_title(plottitle+'Year ' +str(it),fontsize=21)
	ax.set_xlim(-width,len(ind)+width)
	xTickMarks = [str(i) for i in xrange(len(N_init_age_m[0][0]))]
	ax.set_xticks(ind+width)
	xtickNames = ax.set_xticklabels(xTickMarks)
	plt.setp(xtickNames, rotation=0)
	#ax.legend(rects[0],'Total = '+str(TotalN),loc=0)
	savefig(dir+savename+'_Ninitage_year_'+str(it)+'.png',dpi=savedpi)

# Plot Size age bar - stacked side by side
# ------------------------------------------
# Loop through each year to plot
for it in plottime:

	fig = plt.figure()
	ax = fig.add_subplot(111)

	rects = ax.bar(ind+width,size_age_m[0][it],width,color='grey',edgecolor='black',hatch="/",yerr=size_age_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))				

	# axes and labels
	ax.set_ylabel('Size',fontsize=18)
	ax.set_xlabel('Age',fontsize=18)
	ax.set_title(plottitle+'Year ' +str(it),fontsize=21)
	ax.set_xlim(-width,len(ind)+width)
	xTickMarks = [str(i) for i in xrange(len(N_init_age_m[0][0]))]
	ax.set_xticks(ind+width)
	xtickNames = ax.set_xticklabels(xTickMarks)
	plt.setp(xtickNames, rotation=0)
	savefig(dir+savename+'_SizeAge_year_'+str(it)+'.png',dpi=savedpi)

# Plot other N summary measures by class
# --------------------------------------
if plototherN == 'Y':
	width = 0.15                      # the width of the bars
	
	# All N
	# ---------------
	for it in plottime:
		
		TotalNEmi = np.nansum(N_afterEmi_m[0][it])
		TotalNImmi = np.nansum(N_afterImmi_m[0][it])
		TotalNGback = np.nansum(N_growthback_m[0][it])
		TotalNGout = np.nansum(N_growthout_m[0][it])
		TotalN = np.nansum(N_init_class_m[0][it])
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		
		rects0 = ax.bar(ind,N_init_class_m[0][it],width,color='yellow',edgecolor='black',hatch="0",yerr=N_init_class_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))			

		rects1 = ax.bar(ind+width,N_growthback_m[0][it],width,color='grey',edgecolor='black',hatch="-",yerr=N_growthback_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))
		
		rects2 = ax.bar(ind+2*width,N_afterEmi_m[0][it],width,color='white',edgecolor='black',hatch="x",yerr=N_afterEmi_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))	
		
		rects3 = ax.bar(ind+3*width,N_growthout_m[0][it],width,color='blue',edgecolor='black',hatch="+",yerr=N_growthout_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))
		
		rects4 = ax.bar(ind+4*width,N_afterImmi_m[0][it],width,color='red',edgecolor='black',hatch="\\",yerr=N_afterImmi_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))		
			
		# axes and labels
		ax.set_ylabel('N',fontsize=18)
		ax.set_xlabel('Class',fontsize=18)
		ax.set_title(plottitle+'Year ' +str(it),fontsize=21)
		ax.set_xlim(-width,len(ind)+width)
		ax.set_ylim(0, 800)
		xTickMarks = [str(sizeclass[i]) for i in xrange(len(N_afterEmi_m[0][0]))]
		ax.set_xticks(ind+width)
		xtickNames = ax.set_xticklabels(xTickMarks)
		plt.setp(xtickNames, rotation=0)
		ax.legend((rects0[0],rects1[0],rects2[0],rects3[0],rects4[0]),('Init '+str(TotalN),'GrowBack '+str(TotalNGback),'EmiOut '+str(TotalNEmi),'GrowOut '+str(TotalNGout),'ImmiBack '+str(TotalNImmi)),loc=0)
		savefig(dir+savename+'_AllN_year_'+str(it)+'.png',dpi=savedpi)
	
show()

