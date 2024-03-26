### sGD, a spatially explicit estimator of genetic diversity indices
### version 1.02
### Written by Andrew Shirk, University of Washington, Climate Impacts Group

# import required modules

import sys, string, os, csv, operator
import traceback,pdb

__filename__ = "sGD_Python_v1.02.py"


# Define workspace, genefile, costdistancematrix, outputfile name, and parameters

genefile = '/home/elandguth/projects/CDmetaPOP/data_wbp20_20150929/WBP_N_Selection_1443535893/batchrun0mcrun0/genalex_ind25.csv'
#gen = [25,35,40,125,0,50,70,90]
costdistancematrix = '/home/elandguth/projects/CDmetaPOP/data_wbp20_20150929/cdmats/Edmat_PIALXY2km4zones.csv'
workspace = '/home/elandguth/projects/CDmetaPOP/data_wbp20_20150929/sGD'
outfilename = 'N_YesSelection_sGD_gen25'
neighborhoodsize = 30000   ### make sure you enter value expressed in map units
minsamplesize = 20

# create definitions

def getuniquealleles(seq): 
    set = {} 
    map(set.__setitem__, seq, []) 
    return set.keys()
 
def raise_python_error(filename): 
    """Handle python errors and provide details to user"""
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    err = traceback.format_exc().splitlines()[-1]

    print "Python error on **" + line + "** of " + filename
    print err
    exit(0)

    
def sGD():

    try:

                
        # read in genotypes
        print "Reading input file: " + genefile
        f = open(genefile, "r")
        lines = f.readlines()
        f.close()
        numloci = int(lines[0].split(',')[0])
        numindivs = int(lines[0].split(',')[1])
        navalue = 0
        genotypes = []

        for i in xrange(numindivs):
            iline = lines[i + 3]
            n_split_iline = iline.split('\n')[0]
            comma_split_iline = n_split_iline.split(',')
            genotypes.append(comma_split_iline)

            
        # read in costdistance matrix
        print "Reading cost-distance matrix file: " + costdistancematrix
        f = open(costdistancematrix, "r")
        cdlines = f.readlines()
        f.close()

        costdistances = []
        
        for i in xrange(numindivs):
            # Get patch this individual is in
            inpatch = int(genotypes[i][1])-1
            cdline = cdlines[inpatch]
            n_split_cdline = cdline.split('\n')[0]
            comma_split_cdline = n_split_cdline.split(',')
            costdistances.append(comma_split_cdline)

        # determine what individuals are within the max cost-distance of each sample location
        maxcostdistance = 0.5 * neighborhoodsize
        neighborhoods = []
        hoodx = []
        hoody = []
        hoodname = []
        NAcount = 0

        for i in xrange(numindivs):
            # Get patch this individual is in
            #inpatch_i = int(genotypes[i][1])-1
            cdline = costdistances[i]
            templist = []
            for j in xrange(numindivs):
                #inpatch_j = int(genotypes[j][1])-1   
                jindiv = cdline[j]
                if float(jindiv) < float(maxcostdistance):
                    templist.append(j)
                    
            if len(templist) >= minsamplesize:
                neighborhoods.append(templist)
                locline = lines[i+3]
                loc_x = locline.split(',')[-2]
                loc_y = (locline.split('\n')[-2]).split(',')[-1]
                indivID = locline.split(',')[0]
                hoodx.append(loc_x)
                hoody.append(loc_y)
                hoodname.append(indivID)

            else:
                locline = lines[i+3]
                hoodx.append(locline.split(',')[-2])
                hoody.append((locline.split('\n')[-2]).split(',')[-1])
                hoodname.append(locline.split(',')[0])
            
            if len(templist) < minsamplesize:
                neighborhoods.append('NA')
                NAcount = NAcount + 1

        print str(NAcount) + " neighborhoods do not meet the minimum sample size and will have NA values written in ouput file"

        # generate genotype table for neighborhood and calculate genetic diversity indices
        print "Calculating spatial indices of genetic diversity"
        outputfile = open(workspace + os.sep + outfilename + '.csv','w')
        outputfile.write('Neighborhood,X,Y,n,Ho,Hs,Fis,Ar,A')
        outputfile.write('\n')
        ObsHet = []
        ExpHet = []
        AllelicRichness = []
        MeanAlleles = []
        InbreedingCoeff = []
        nlist = []

        # calculate the smallest number of individuals in a neighborhood (used to calculate Ar later)

        glist = []

        for i in xrange(len(neighborhoods)):
            if neighborhoods[i] !='NA':
                glist.append(len(neighborhoods[i]))

        g = min(glist)

        counter = 0

        # for the ith neighborhood...

        for i in xrange(len(neighborhoods)):

            if neighborhoods[i] != 'NA':
                neighborhoodgenes = []
                counter = counter + 1
                
                # build a table of genotypes containing all individuals in the neighborhood 

                for j in neighborhoods[i]:
                    neighborhoodgenes.append(genotypes[j])
                    
                # calculate diversity indices within the ith neighborhood 
                
                # for the kth locus...:
                Hsk = []
                Hok = []
                Pijgk = []
                Ak = []
                n = len(neighborhoods[i])
                nlist.append(n)

                for k in xrange(numloci):
                    allelelist = []
                    p2 = []
                    NAcount = 0
                    hetcount = 0
                    p2sum = 0
               
                    for j in xrange(len(neighborhoodgenes)): 
                      
                        # track number of NAs per locus
                        if neighborhoodgenes[j][2*k+2] == str(navalue):
                            NAcount = NAcount + 1
                        if neighborhoodgenes[j][2*k+3] == str(navalue):
                            NAcount = NAcount + 1

                        # track heterozygous loci 
                        if neighborhoodgenes[j][2*k+2] != neighborhoodgenes[j][2*k+3]:
                            hetcount = hetcount + 1
                            
              
                    nk = len(neighborhoods[i]) - (NAcount / 2) # n adjusted for NAs
                      
                    # for the jth individual in the ith neighborhood, identify unique alleles at the kth locus and write to list 'uniquealleles' -- used later for He
                    for j in xrange(len(neighborhoodgenes)):  
                        allelelist.append(neighborhoodgenes[j][2*k+2])
                        allelelist.append(neighborhoodgenes[j][2*k+3])

                    uniquealleles = getuniquealleles(allelelist)

                    # remove navalue from unique allele list
                    for allele in uniquealleles:
                        if allele.startswith(str(navalue)):
                            uniquealleles.remove(allele)

                    Ak.append(len(uniquealleles))
                    
                    if nk >> 0:

                        Hok.append(float(hetcount) / float(nk)) #nk is already adjusted for NAs

                    if nk == 0:           # if neighborhood is all nodata, assign NA value
                        Hok.append('NA')

                    # for the lth unique allele at the kth locus, in the ith neighborhood, calculate allele frequency (p)
                    for allele in uniquealleles: 
                    
                        allelecount = 0
                             
                        # for the jth individual in neighborhood i at locus k, add up number of occurances of allele l
                        for j in xrange(len(neighborhoodgenes)): 

                            # tally lth allele across individuals
                            if neighborhoodgenes[j][2*k+2] == allele:
                                allelecount = allelecount + 1
                            if neighborhoodgenes[j][2*k+3] == allele:
                                allelecount = allelecount + 1

                        p = float(allelecount) /  float(nk * 2)
                        p2.append (p ** 2) ## lth allele frequency squared

                    for m in xrange(len(p2)):
                        p2sum = p2sum + p2[m] 
                            
                    if nk == 0:                 ## if neighborhood is all NA's, assign value of NA
                        Hsk.append('NA')
                        
                    if nk == 1:                 ## if there is only one indiv without NAs in neighborhood...
                        Hsk.append((1.0 - float(p2sum) - (Hok[k] / (2 * float(nk)))))

                    if nk >> 1:
                        Hsk.append(float(nk) / float(nk-1) * (1.0 - float(p2sum) - (Hok[k] / (2 * float(nk)))))

                    ### calculate Ar

                    Pijg = 0 # Ar per locus across all unique alleles
                    Nj = nk # renamed to keep terms consistent

                    # for each unique allele per locus, tally the number of observations (Nij)
                    for allele in uniquealleles:
                        Nij = 0
                            
                        for j in xrange(len(neighborhoodgenes)):
                            if neighborhoodgenes[j][2*k+2]== allele:
                                Nij = Nij + 1
                            if neighborhoodgenes[j][2*k+3]== allele:
                                Nij = Nij + 1
                            
                        # now use combinatorics to calculate Pijg...
                        
                        # for the numerator:
                        if Nj > Nij:
                            nfactorial = reduce(operator.mul,xrange(1,int(Nj-Nij) + 1))
                        if Nj <= Nij:
                            nfactorial = 1
                            
                        kfactorial = reduce(operator.mul,xrange(1,int(g) + 1))
                        
                        if Nj-Nij > g:
                            nminuskfactorial = reduce(operator.mul,xrange(1,(int(Nj-Nij)-g) + 1))
                            numerator = float(nfactorial / (kfactorial * nminuskfactorial))
                        if Nj-Nij <= g:
                            numerator = 1

                        # for the denominator:
                        
                        nfactorial = reduce(operator.mul,xrange(1,int(Nj) + 1))
                        
                        if Nj > g:
                            nminuskfactorial = reduce(operator.mul,xrange(1,int(Nj-g) + 1))
                        if Nj <= g:
                            nminuskfactorial = 1
                            
                        denominator = float(nfactorial / (kfactorial * nminuskfactorial))

                        if Nj-Nij > g:
                            Pijg = Pijg + (1 - (float(numerator/denominator)))

                        if Nj-Nij <= g:
                            Pijg = Pijg + 1
                                
                    Pijgk.append(Pijg)

                #calculate per locus Fis,Ho,He,Ar

                Ho = 0
                Hs = 0
                FIS = 0
                Ar = 0
                A = 0
                NAlocicount = 0
                
                for x in xrange(numloci):
                    if Hok[x] != 'NA':
                        Ho = Ho + Hok[x]
                    if Hok[x] == 'NA':
                        NAlocicount = NAlocicount + 1
                        
                ObsHet.append(Ho /(numloci - NAlocicount))

                NAlocicount = 0
                                   
                for x in xrange(numloci):
                    if Hok[x] != 'NA':
                        Hs = Hs + Hsk[x]
                    if Hok[x] == 'NA':
                        NAlocicount = NAlocicount + 1
                        
                ExpHet.append(Hs /(numloci - NAlocicount))

                NAlocicount = 0

                if ExpHet[i] != 0:
                    FIS = 1 - float(ObsHet[i] / ExpHet[i])

                if ExpHet[i] == 0:
                    FIS = 0
                                    
                if NAlocicount < numloci:
                    InbreedingCoeff.append(FIS)

                if NAlocicount == numloci:
                    InbreedingCoeff.append('NA')

                NAlocicount = 0
                
                for x in xrange(numloci):
                    if Pijgk[x] != 'NA':
                        Ar = Ar + Pijgk[x]

                    if Pijgk[x] == 'NA':
                        NAlocicount = NAlocicount + 1

                if NAlocicount < numloci:
                    AllelicRichness.append(float(Ar) /float(numloci - NAlocicount))

                if NAlocicount == numloci:
                    AllelicRichness.append('NA')

                NAlocicount = 0
                
                for x in xrange(numloci):
                    if Ak != 'NA':
                        A = A + Ak[x]

                    if Ak[x] == 'NA':
                        NAlocicount = NAlocicount + 1

                if NAlocicount < numloci:
                    MeanAlleles.append(float(A) /float(numloci - NAlocicount))

                if NAlocicount == numloci:
                    MeanAlleles.append('NA')

            else:
                nlist.append('NA')
                ObsHet.append('NA')
                ExpHet.append('NA')
                InbreedingCoeff.append('NA')
                AllelicRichness.append('NA')
                MeanAlleles.append('NA')
            
            
        for i in xrange(len(neighborhoods)):
            outputfile.write(str(hoodname[i]) + ',')
            outputfile.write(str(hoodx[i]) + ',')
            outputfile.write(str(hoody[i]) + ',')
            outputfile.write(str(nlist[i]) + ',')
            outputfile.write(str(ObsHet[i]) + ',')
            outputfile.write(str(ExpHet[i]) + ',')
            outputfile.write(str(InbreedingCoeff[i]) + ',')
            outputfile.write(str(AllelicRichness[i]) + ',')
            outputfile.write(str(MeanAlleles[i]) + ',') 
            outputfile.write('\n')

        outputfile.close()

        print "sGD complete, output file " + outfilename + ".csv written to " + workspace

    except:
        raise_python_error(__filename__)

if __name__ == "__main__":
    sGD()
