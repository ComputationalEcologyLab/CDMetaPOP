======
README
======

---------------------- 
CDMetaPOP 3.00 release
----------------------
  
Welcome to the CDMetaPOP v3.0 release! This release includes installation instructions, version notes, some examples, and technical documentation. 

For the user manual that contains complete documentation, please see the doc/ directory in this repository.
  
Program Contributors: Erin Landguth, Casey Day, Andrew Bearlin, Jason Dunham, Ryan Simmons, Brenna Forrester, Kaeli Davenport, and Travis Seaborn
Link: https://github.com/ComputationalEcologyLab/CDMetaPOP
Version: 3.08 
Python: 3.8
Release Date: 2025.09.18
README Update: 2025.11.7 (ccd)
  
--------
Contents
--------
  
Included in this release are the following:

src -> CDMetaPOP source files

doc -> README.txt, user manual, history, and disclaimer

example_files -> test example files
- RunVars.csv - Runs several scenarios (Variable climate, population introduction) with varying life history parameters.
- RunVars_AddMyy.csv and RunVars_AddFyy.csv - Demonstrates the Trojan Y chromosome control strategy using YY males or YY females, respectively.
- RunVars_multispecies.csv - Demonstrates 2- and 3-species applications, including competition.
- DiseaseExamples
    - Adaptive_Run06 - Demonstrates genetic adaptation to disease via resistance and tolerance.
    - OnePatch_SIDP - SIR simulations to match simple theoretical patterns. 
  
---------------------------------------
Requirements and Pre-requisite Software
---------------------------------------

1. Baseline Requirements. CDMetaPOP requires the Python3.8.x interpreter, NumPy package, and SciPy package. Remember that Python modules usually require particular Python interpreters, so be sure the version ID for any external Python module or package (e.g. NumPy or others) matches the version of your Python interpreter (normally v3.8.x). To avoid Python installation errors, we highly recommend installing Python from any number of the freely available bundlers that include the NumPy and SciPy packages, e.g., Anaconda (recommended), Canopy, ActiveState.

---------------------------
CDMetaPOP v3.0 Installation
--------------------------- 

Linux or Windows: Download the repository from Github. Navigate to the directory on your PC where you wish to install CDMetaPOP, and unpack the zipped repository using a free archive tool like 7Zip (7z.exe), Pkunzip, Unzip, or an equivalent. Seven-Zip (7Z.exe) is highly recommended since it can handle all common formats on Windows, MAC OS X and Linux. On Windows, it is best to set up a project specific modeling subdirectory to perform your modeling outside of any folder that has spaces in its name (like "My Documents").

---------------------
Example CDMetaPOP Run
---------------------

The primary example run ('RunVars.csv') is for 7 patches representing an effective distance matrix calculated using a least-cost path algorithm through riverine distance . To run the following example, follow these steps:

1. Double check that the 3 directories provided in the Git download are in the same directory. 

2. The included file ‘RunVars.csv’ in the example_files directory specifies the parameters that apply to all species that can be changed and used in a sample CDMetaPOP run. Open ‘RunVars.csv’ in your editor of choice. A spreadsheet program like Microsoft Excel allows for easy editing of the tabular values. The location(s) of the 'PopVars' files (one for each species in the simulation) are specified in the first column.

3. The various ‘PopVars’ files define the patch files in the first column. The included ‘PatchVars’ files will also be in the same folder (../example_files/patchvars). ‘ClassVars’ files are in turn specified in the ‘PatchVars’ files and example ‘ClassVars’ files will be in the ../example_files/classvars/ folder. 

4. There will be 5 lines of information in ‘RunVars.csv’: a header line and 4 lines of information corresponding to 4 separate, single-species CDMetaPOP runs. See Table 1 in user manual which contains a breakdown for each column header and the parameters that can be defined. The ‘Input’ in the table listed is for the first row in the file. Make sure you save this file in the same format – a comma delimited file – when you make changes to the parameters. Do not change the ‘Input’ (first row) labeling. Select ‘Yes’ or ‘OK’ for any Excel questions about saving in this format. 'RunVars_multispecies.csv' contains two runs, the first for a 2-species example and the second for a 3-species example.

5. Start the program: For example, if you use python from the command line, then open a terminal window and change your shell directory to the CDMetaPOP src home directory (i.e., > cd C:\"homedirectorylocation"\src). 

6. Launch CDMetaPOP: There are a number of ways to run CDMetaPOP. If you are using a command shell you can run the program by typing

    `python CDMetaPOP.py C:/”homedirectorylocation”/example_files/ RunVars.csv output_test`

   Or a short-cut if your data is located at the same folder level as the src folder:
   
   `python CDMetaPOP.py ../example_files/ RunVars.csv output_test`

   Note that there are 5 arguments here that must be included with spaces in between: 

    1. `python` starts python, for example from the command line. Note that other python environments may have different calls here. In iPython (the IDE distributed with Anaconda) the call is “run”. 
    2. `CDMetaPOP.py` runs CDMetaPOP program.
    3. `C:/”homedirectorylocation”/example_files` is the directory location of the input test files. You can point this directory to other project files, for example. We suggest not having any spaces in your directory names. So as projects accumulate you can rename input folders that contain the project specific files (e.g., dataWestslope or dataBullTrout).
    4. `RunVars.csv` is the primary parameter file (comma delimited) which can be renamed (e.g., “RunVars_WCT.csv”). Caution should be taken when going between operating systems and saving this file as a .csv.
    5. `output_test` is the name of the directory that will be created with CDMetaPOP output in the directory specified by the third argument above.

9. Check for successful model run completion: The program will provide step-by-step output in the Shell window. Each row of RunVars.csv will run an independent simulation in sequence for each line in PopVars.csv (batches). Once completed, a simulation time will be printed out and folders run0batch0mc0species0, run0batch0mc1species0,  run0batch1mc0species0, etc. will be created in your CDMetaPOP home directory to store output from the separate runs, batches Monte-Carlo replicates, and species (each line in the RunVars file corresponds to a separate 'run' and each line in the PopVars file corresponds to a separate 'batch'. Monte Carlo runs are specified by 'mc'). These folders are located in the data folder specified in above step. The output folder will have a unique date/time stamp after the name of the output folder in case you want to run multiple CDMetaPOP runs in this same directory. The program will also provide a log file with program steps in your specified output directory. If parameters are such that all species become extinct before the specified generation time, the program will end. The program will provide error and feedback for parameters that are outside of ranges or incorrectly entered.

Happy Simulations!

Computational Ecology Laboratory
The University of Montana
32 Campus Drive
Missoula MT, 59812-1002
computationalecologylab@gmail.com
