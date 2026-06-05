# CDmetaPOP.py
# Author: Erin L Landguth / Casey Day
# Created: February 2008
# v 1.0 Release: MARCH 2014
# v 2.0: May 2020
# ----------------------------------------------------------------------------
# General CDmetaPOP information
appName = "CDmetaPOP"
appVers = "version 3.08"
appRele = "2025.07.28-09:09:09"
authorNames = "Erin L Landguth, Casey Day, et al."

# ------------------------------------------
# Import Modules with Except/Try statements
# ------------------------------------------
# Python specific functions
import time, os, sys, warnings
from dataclasses import dataclass
from typing import List
import multiprocessing as mp

# CDmetaPOP functions
import CDmetaPOP_PreProcess as preprocess
import CDmetaPOP_mainloop as mainloop

@dataclass
class inputStruct:
    '''
    Data structure for grouping the input parameters
    '''
    sizeans: str
    constMortans: str
    mcruns: int
    looptime: int
    nthfile_out: int
    gridformat: str
    gridsample: str
    outputans: str
    cdclimgentimelist: List[int]
    startcomp: int
    implementcomp: int
    ncores: int
    @classmethod
    def from_batch(cls, batch, irun):
        '''
        A function for reading in the parameters, that interfaces with the data structure
        '''
        return cls(
            sizeans = batch['sizecontrol'][irun],
            constMortans = batch['constMortans'][irun],
            mcruns = int(batch['mcruns'][irun]),
            looptime = int(batch['runtime'][irun]),
            nthfile_out = batch['output_years'][irun],
            gridformat = batch['gridformat'][irun],
            gridsample = batch['gridsampling'][irun],
            outputans = str(batch['summaryOutput'][irun]),
            cdclimgentimelist = batch['cdclimgentime'][irun],
            startcomp = int(batch['startcomp'][irun]),
            implementcomp = str(batch['implementcomp'][irun]),
            ncores = int(batchVars1['ncores'][irun])
        )

@dataclass
class ContextStruct:
    appName: str
    appVers: str
    appRele: str
    authorNames: str
    popvarsfile: str
    spcNO: int
    fileans: str
    irun: int
    datadir: str
    outdir: str
    passlogfHndl: str
    nspecies: int
    current_system_pid: int

def check_Popvars_Lenght(datadir, popvarsfile):
    '''Error check: Rows in PopVars must be equal'''
    storerowno = [preprocess.loadFile(datadir + f, 1, ',', True)[2] for f in popvarsfile]
    if len(set(storerowno)) != 1:
        print('Length of PopVars files do no match.')
        sys.exit(-1)

def create_Qs(nspecies):
    """Create Qs for multiprocessing Put/Gets"""
    extinctQ = mp.Queue() # To track extinction. If all species extinct, exit system
    global_extinctQ = mp.Queue() # To track global extinction
    XQs = [] # For multiprocessing setup, create empty list to fill with # species Qs
    if nspecies > 1: # Ignore queues if only one species
        for _ in range(nspecies):
            XQs.append([mp.Queue() for _ in range(nspecies)])
    return XQs, extinctQ, global_extinctQ

def create_workers(inputs, popvarsfile, irun, datadir, outdir, XQs, extinctQ, global_extinctQ):
    """Creates and joins processes"""
    nspecies = len(popvarsfile)
    if nspecies > mp.cpu_count():
        print("PopVars files given greater than number of CPUs.")
        sys.exit(-1)

    sp = []
    
    for ispecies in range(nspecies):
        context = ContextStruct(
            appName=appName, 
            appVers=appVers, 
            appRele=appRele, 
            authorNames=authorNames, 
            popvarsfile=popvarsfile, #[ispecies],
            spcNO=ispecies, 
            fileans=datadir + popvarsfile[ispecies],
            irun=irun, 
            datadir=datadir, 
            outdir=outdir,
            passlogfHndl=f"{outdir}CDmetaPOP{ispecies}.log", 
            nspecies=nspecies, 
            current_system_pid=os.getpid()
        )
        
        if nspecies > 1:
            p = mp.Process(
                target=mainloop.main_loop, 
                name=f'S{ispecies}', 
                args=(inputs, context, XQs, extinctQ, global_extinctQ)
            )
            sp.append(p)
            p.start()
        else:
            mainloop.main_loop(inputs, context, XQs, extinctQ, global_extinctQ)
            
    for p in sp:
        p.join()

#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 
if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    
    # ------------------------------------------------------
    # Start timer, get script arguments, create log writeout
    # ------------------------------------------------------

    if len(sys.argv) == 4: # Split PopVars len() assume noproc used
        datadir = sys.argv[1]+'/'
        fileans = datadir+sys.argv[2]
        outdir = datadir+sys.argv[3]+str(int(time.time()))+'/'
    else: # If user did not specify .rip file
        print("User must specify data directory, input file name, and output file directory, e.g., at command line type CDmetaPOP.py ../CDmetaPOP_data/ RunVars.csv exampleout_foldername.")
        sys.exit(-1)
    
    # If .ip file does not exist
    if not os.path.exists(fileans):
        print(("Cannot find or open runtime inputs file(%s)"%(fileans)))
        sys.exit(-1)
    
    # Create output file director
    os.mkdir(outdir)
    current_system_pid = os.getpid() # Get parent ID for terminate
    
    # Call function and store inputvariables
    batchVars1,batchVarsIndex,nSimulations = preprocess.loadFile(fileans,1,',',True)

    # ----------------------------------------
    # Begin Run Looping - assign processors
    # ----------------------------------------
    # This loop is defined by the number of rows in RunVars.csv
    for irun in range(nSimulations):
        print(("On run: %s"%(str(irun))+'\n')) 
        
        # Store RunVars Here
        popvarsfile = batchVars1['Popvars'][irun].split(';')
        nspecies = len(popvarsfile)
        inputs = inputStruct.from_batch(batchVars1, irun)

        check_Popvars_Lenght(datadir, popvarsfile)

        XQs, extinctQ, global_extinctQ = create_Qs(nspecies)

        create_workers(inputs, popvarsfile, irun, datadir, outdir, XQs, extinctQ, global_extinctQ)

    #End::Batch Loop

# End::Main Loop

