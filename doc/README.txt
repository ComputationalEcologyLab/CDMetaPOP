======
README
======

---------------------- 
CDMetaPOP 1.00 release
----------------------
  
Welcome to the CDMetaPOP v1.0 release! This release includes installation instructions, version notes, some examples, and technical documentation. 
  
Program Contributors: Erin Landguth, Andrew Bearlin, Casey Day, Jason Dunham
Link: http://cel.dbs.umt.edu/software/CDMetaPOP/
Version: 1.00 
Python: 2.7.4
Release Date: 2015.12.14
README Update: 2015.12.14 (ell)
  
--------
Contents
--------
  
Included in this release are the following:

src -> CDMetaPOP source files
doc -> README.txt, user manual, history, and disclaimer
data -> test example files
  
---------------------------------------
Requirements and Pre-requisite Software
---------------------------------------

1. Baseline Requirements. CDMetaPOP requires the Python2.7.x interpreter, NumPy package, and SciPy package. Remember that Python modules usually require particular Python interpreters, so be sure the version ID for any external Python module or package (e.g. NumPy or others) matches the version of your Python interpreter (normally v2.7.x). aTo avoid Python installation errors, we highly recommend installing Python from any number of the freely available bundlers, e.g., Canopy, ActiveState, Anaconda.

---------------------------
CDMetaPOP v1.0 Installation
--------------------------- 

Linux or Windows: Unpack the CDMetaPOP Archive. Navigate to the directory on your PC where you wish to install CDMetaPOP, and unpack the supplied zip archive file using a free archive tool like 7Zip (7z.exe), Pkunzip, Unzip, or an equivalent. Seven-Zip (7Z.exe) is highly recommended since it can handle all common formats on Windows, MAC OS X and Linux. On Windows, it is best to setup a project specific modeling subdirectory to perform your modeling outside of any folder that has spaces in its name (like "My Documents").

---------------------
Example CDMetaPOP Run
---------------------

The example run is for 7 patches representing an effective distance matrix calculated using a least-cost path algorithm through riverine distance . To run the following example, follow these steps:

1. Double check that the 3 directories provided in the archive are in the same directory. 

2. The included file ‘PopVars.csv’ in the data directory specifies the parameters that can be changed and used in a sample CDMetaPOP run. Open ‘PopVars.csv’ in your editor of choice. A spreadsheet program like Microsoft Excel, allows for easy editing of the tabular values.

3. The file ‘PopVars.csv’ defines the patch files in the first column. The included files ‘PatchVars.csv’ will also be in the same folder (../data). ‘ClassVars.csv’ files are in turn specified in the ‘PatchVars.csv’ file and example ‘ClassVars.csv’ will be in the ../data/classvars/ folder. 

4. There will be 3 lines of information in ‘PopVars.csv’: a header line and 2 lines of information corresponding to 2 separate CDMetaPOP runs (batch process). See Table 1 in user manual which contains a breakdown for each column header and the parameters that can be changed. The ‘Input’ in the table listed is for the first row in the file. Make sure you save this file in the same format – a comma delimited file – when you make changes to the parameters. Do not change the ‘Input’ (first row) labeling. Select ‘Yes’ or ‘OK’ for any Excel questions about saving in this format.

5. Start the program: For example, if you use python from the command line, then open a terminal window and change your shell directory to the CDMetaPOP src home directory (i.e., > cd C:\"homedirectorylocation"\src). 

6. Run the program: There are a number of ways to run this program. If you are using a command shell you can run the program by typing “python CDMetaPOP.py C:/”homedirectorylocation”/data PopVars.csv output_test”. Or a short-cut if your data is located at the same folder level as the src folder: “python CDMetaPOP.py ../data/ PopVars.csv output_test”. Note that there are 5 arguments here that must be included with spaces in between: 

a.	“python” starts python, for example from the command line. Note that other python environments may have different calls here. In PyLab (the IDE distributed with Canopy), the call is “run”. 
b.	“CDMetaPOP.py” runs CDMetaPOP program.
c.	“C:/”homedirectorylocation”/data” is the directory location of the input test files. You can point this directory to other project files, for example. We suggest not having any spaces in your directory names. So as projects accumulate you can rename input folders that contain the project specific files (e.g., dataWestslope or dataBullTrout).
d.	“PopVars.csv” is the parameter file (comma delimited) which can be renamed (e.g., “PopVars_WCT.csv”). Caution should be taken when going between operating systems and saving this file as a .csv.
e.	“output_test” is the name of the directory that will be created with CDMetaPOP output in the directory specified by the third argument above.

7. Check for successful model run completion: The program will provide step-by-step output in the Shell window. Once completed, a simulation time will be printed out and folders batchrun0mcrun0, batchrun0mcrun1, batchrun0mcrun2, batchrun1mcrun0, and batchrun2mcrun0 will be created in your CDMetaPOP home directory to store output from the separate batch and/or Monte-Carlo runs (each line in the PopVars file corresponds to a separate batchrun and the specified ‘mcruns’ for each batch). These folders are located in the data folder specified in above step. The output folder will have a unique date/time stamp after the name of the output folder in case you want to run multiple CDMetaPOP runs in this same directory. The program will also provide a log file with program steps in your specified output directory. If parameters are such that population becomes extinct before specified generation time, then program will end. The program will provide error and feedback for parameters that are outside of ranges or incorrectly entered.

Happy Simulations!

Erin.

Contact Information
Erin Landguth
Computational Ecology Laboratory
The University of Montana
32 Campus Drive
Missoula MT, 59812-1002
erin.landguth@mso.umt.edu
