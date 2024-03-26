# CDmetaPOP.py
# Author: Erin L Landguth / Casey Day
# Created: February 2008
# v 1.0 Release: MARCH 2014
# v 2.0: May 2020
# ----------------------------------------------------------------------------
# General CDmetaPOP in3formation
appName = "CDmetaPOP"
appVers = "version 2.67"
appRele = "2024.03.26-11:24:01"
authorNames = "Erin L Landguth, Casey Day, Andrew Bearlin, Ryan Simmons, Travis Seaborn, et al."

# ---------------
# Global symbols
#----------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False
# File absolute paths for importing functions
SRC_PATH =  "../src/"

# ------------------------------------------
# Import Modules with Except/Try statements
# ------------------------------------------
# Python specific functions
import datetime,time,pdb,os,sys,shutil,gc,warnings

# For parallel processing
from multiprocessing import Process
from multiprocessing import Queue
import multiprocessing
import numpy as np

#Import the package specific folders
CDPOP_folder = os.path.dirname(os.path.abspath(SRC_PATH+"CDmetaPOP"))

if CDPOP_folder not in sys.path:
	sys.path.insert(0, CDPOP_folder)

# CDmetaPOP functions
from CDmetaPOP_Modules import * 
from CDmetaPOP_PreProcess import *
from CDmetaPOP_mainloop import *
		
#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 
if __name__ == '__main__':
	warnings.filterwarnings("ignore")	
	
	# ------------------------------------------------------	
	# Start timer, get script arguments, create log writeout
	# ------------------------------------------------------
	# Timing events: start
	start_time = datetime.datetime.now()
	foldertime = int(time.time())
	
	# Split PopVars len() assume noproc used
	if len(sys.argv) == 4:
		datadir = sys.argv[1]+'/'
		fileans = datadir+sys.argv[2]
		outdir = datadir+sys.argv[3]+str(foldertime)+'/'
	
	# If user did not specify .rip file
	else:
		print("User must specify data directory, input file name, and output file directory, e.g., at command line type CDmetaPOP.py ../CDmetaPOP_data/ RunVars.csv exampleout_foldername.")
		sys.exit(-1)	
	
	# If .ip file does not exist
	if not os.path.exists(fileans):
		print(("Cannot find or open runtime inputs file(%s)"%(fileans)))
		sys.exit(-1)
	
	# Create output file directory - will automatically put in the data directory
	os.mkdir(outdir)
	current_system_pid = os.getpid() # Get parent ID for terminate
		
	# ------------------------------------	
	# Call DoUserInput()
	# ------------------------------------
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Call function and store inputvariables
	batchVars1,batchVarsIndex,nSimulations = loadFile(fileans,1,',',True)
	
	# ----------------------------------------	
	# Begin Run Looping - assign processors
	# ----------------------------------------
	# This loop is defined by the number of rows in RunVars.csv
	for irun in range(nSimulations):
		
		# Timing events: start
		start_timeB = datetime.datetime.now()
		print(("On run: %s"%(str(irun))+'\n')) 
		
		# Store RunVars Here
		popvarsfile = batchVars1['Popvars'][irun].split(';')
		sizeans = batchVars1['sizecontrol'][irun]
		constMortans = batchVars1['constMortans'][irun]
		mcruns = int(batchVars1['mcruns'][irun])
		looptime = int(batchVars1['runtime'][irun])
		nthfile_out = batchVars1['output_years'][irun]
		gridformat = batchVars1['gridformat'][irun]
		gridsample = batchVars1['gridsampling'][irun]
		outputans = str(batchVars1['summaryOutput'][irun])
		cdclimgentimelist = batchVars1['cdclimgentime'][irun]
		startcomp = int(batchVars1['startcomp'][irun])
		implementcomp = str(batchVars1['implementcomp'][irun])
		
		# ------------------------------------------
		# Error check: Rows in PopVars must be equal
		# ------------------------------------------
		storerowno = []
		for ispecies in range(len(popvarsfile)):
			temp1,temp2,checkrowno = loadFile(datadir+popvarsfile[ispecies],1,',',True)
			storerowno.append(checkrowno)
		storerowno = np.asarray(storerowno)
		checkunique = np.unique(storerowno)
		if len(checkunique) != 1:
			print('Length of PopVars files do no match.')
			sys.exit(-1)
		
		# ------------------------------------------------------
		# Create Qs for multiprocessing Put/Gets 
		# ------------------------------------------------------
		XQs = [] # For multiprocessing setup, create empty list to fill with # species Qs
		extinctQ = Queue() # To track extinction. If all species extinct, exit system
		global_extinctQ = Queue() # To track global extinction
		logfHndl = []
		for ispecies in range(len(popvarsfile)):
			#Ignore queues if only one species
			if len(popvarsfile) > 1:
				XQs.append([])
				for ispecies2 in range(len(popvarsfile)): 
					XQs[ispecies].append(Queue())
			
			# This properly names log file
			logSessionPath = outdir+"CDmetaPOP"+str(ispecies)+".log"
			logfHndl.append(open(logSessionPath,'a'))
			
			msgVerbose = True
			logMsg(logfHndl[ispecies],"\n%s Release %s Version %s\n"%(appName,appRele,appVers))
			logMsg(logfHndl[ispecies],"Author(s): %s"%(authorNames)+'\n')
			logMsg(logfHndl[ispecies],"Session runtime inputs from: %s"%(fileans)+'\n\n') 
			logMsg(logfHndl[ispecies],"Session popvars inputs from: %s"%(popvarsfile[ispecies])+'\n\n')
			logMsg(logfHndl[ispecies],"On run: %s"%(str(irun))+'\n\n')
			logfHndl[ispecies].close() # in order to write out above			
			msgVerbose = False
			
		# --------------------------------------
		# Split processors here len(popvarsfile)
		# -------------------------------------- 
		#pdb.set_trace()
		if len(popvarsfile) > multiprocessing.cpu_count():
			print("PopVars files given greater than number of CPUs.")
			sys.exit(-1)
		sp = [] # list of processes for appending
		__spec__ = None #This is a fix for an Ipython error that was looking for this variable when using multiprocessing
		if len(popvarsfile) > 1:
			
			for ispecies in range(len(popvarsfile)):
				# Need to create target function main_loop
				sp.append(Process(target=main_loop, name='S'+str(ispecies), args=(ispecies,datadir+popvarsfile[ispecies],irun,datadir,sizeans,constMortans,mcruns,looptime,nthfile_out,gridformat,gridsample,outputans,cdclimgentimelist,outdir,startcomp,implementcomp,outdir+"CDmetaPOP"+str(ispecies)+".log",XQs, len(popvarsfile), extinctQ, global_extinctQ,current_system_pid)))
			# Now Start Processes
			for ispecies in range(len(popvarsfile)):
				sp[ispecies].start()
			# Now Join Processes
			for ispecies in range(len(popvarsfile)):
				sp[ispecies].join()
		elif len(popvarsfile) == 1:
			main_loop(0,datadir+popvarsfile[0],irun,datadir,sizeans,constMortans,mcruns,looptime,nthfile_out,gridformat,gridsample,outputans,cdclimgentimelist,outdir,startcomp,implementcomp,outdir+"CDmetaPOP0.log",XQs, len(popvarsfile), extinctQ, global_extinctQ,current_system_pid)
			
	# Close the logfHndl file
	for ispecies in range(len(popvarsfile)):
		logfHndl[ispecies].close()	
	#End::Batch Loop
	
# End::Main Loop	

