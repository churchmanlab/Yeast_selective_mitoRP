HELP_STRING = """
SAM2mitoFPBEDplus_MTC

Author: Stirling Churchman
Date: January 15, 2014
Updated 11/2022 by M. Couvillion to handle soft-clipped reads (STAR alignment)

Here we take alignments and convert them to A site positions. Then the alignments are put back into a BED format.
I used the following to interpret the SAM file. http://samtools.sourceforge.net/SAMv1.pdf and http://picard.sourceforge.net/explain-flags.html

This will set the A site at 16 nt from the 3' end (for mitoribosome footprint reads that have been 5' end trimmed (1nt) use reads 36-40ntonly)

Only footprints sizes where shifts have been specified will be processed. 


     -h     print this help message
     -i     file input (required)
     -s     size range (e.g. 37to41)


"""
import sys
from getopt import getopt
import re

def main(argv=None):
    if argv is None:
        argv = sys.argv

    inputFile = ""
    outputFilePlus = ""
    sizeRange = ""

    try:
        optlist, args = getopt(argv[1:], "hi:s:")
    except:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    if len(optlist) == 0:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    for (opt, opt_arg) in optlist:
        #print opt
        #print opt_arg
        if opt == '-h':
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == '-i':
            inputFile = opt_arg
        elif opt == '-s':
            sizeRange = opt_arg
    if inputFile == "":
        print HELP_STRING
        sys.exit(1)
    
    def getReadSize(inString):
        if inString.count('N') == 0:
            readSize = int(inString.split('M')[0])
        elif inString.count('N') == 1:
            readSize = sum([int(i) for i in re.split('M|N', inString)[0:3:2]])
        elif inString.count('N') == 2:
            readSize = sum([int(i) for i in re.split('M|N', inString)[0:5:2]])
        else:
            readSize = sum([int(i) for i in re.split('M|N', inString)[0:7:2]])
        return readSize
        
    def getRefSize(inString):
        if inString.count('N') == 0:
            refSize = int(inString.split('M')[0])
        elif inString.count('N') == 1:
            refSize = sum([int(i) for i in re.split('M|N', inString)[0:3]])
        elif inString.count('N') == 2:
            refSize = sum([int(i) for i in re.split('M|N', inString)[0:5]])
        else:
            refSize = sum([int(i) for i in re.split('M|N', inString)[0:7]])
        return refSize


    def getAsitePlus(inString,readSize,fivePrime,refSize, lower, upper):
        if read_size in range(lower, upper + 1):
            threeP = fivePrime + refSize - 1
            if inString.count('N') == 0:
                if readSize == 37:
                    Asite = threeP - 15
                elif readSize == 38:
                    Asite = threeP - 16
                elif readSize == 39 or readSize == 40 or readSize == 41:
                    Asite = threeP - 17
                else:
                    Asite = 'OutOfRange'
            elif inString.count('N') == 1:
                e1 = int(re.split('M|N', inString)[0])
                i1 = int(re.split('M|N', inString)[1])
                e2 = int(re.split('M|N', inString)[2])
                if readSize == 37:
                    offset = 15
                    if e2 > offset:
                        Asite = threeP - offset
                    else:
                        Asite = threeP - e2 - i1 - (offset-e2)
                elif readSize == 38:
                    offset = 16
                    if e2 > offset:
                        Asite = threeP - offset
                    else:
                        Asite = threeP - e2 - i1 - (offset-e2)
                elif readSize == 39 or readSize == 40 or readSize == 41:
                    offset = 17
                    if e2 > offset:
                        Asite = threeP - offset
                    else:
                        Asite = threeP - e2 - i1 - (offset-e2)
                else:
                    Asite = 'OutOfRange'
            elif inString.count('N') == 2:
                e1 = int(re.split('M|N', inString)[0])
                i1 = int(re.split('M|N', inString)[1])
                e2 = int(re.split('M|N', inString)[2])
                i2 = int(re.split('M|N', inString)[3])
                e3 = int(re.split('M|N', inString)[4])
                if readSize == 37:
                    offset = 15
                    if e3 > offset:
                        Asite = threeP - offset
                    elif e2 + e3 > offset:
                        Asite = threeP - e3 - i2 - (offset-e3)
                    else:
                        Asite = threeP - e3 - i2 - e2 - i1 -(offset - e2 - e3)
                elif readSize == 38:
                    offset = 16
                    if e3 > offset:
                        Asite = threeP - offset
                    elif e2 + e3 > offset:
                        Asite = threeP - e3 - i2 - (offset-e3)
                    else:
                        Asite = threeP - e3 - i2 - e2 - i1 -(offset - e2 - e3)
                elif readSize == 39 or readSize == 40 or readSize == 41:
                    offset = 17
                    if e3 > offset:
                        Asite = threeP - offset
                    elif e2 + e3 > offset:
                        Asite = threeP - e3 - i2 - (offset-e3)
                    else:
                        Asite = threeP - e3 - i2 - e2 - i1 -(offset - e2 - e3)
                else:
                    Asite = 'OutOfRange'
            else:
                Asite = 'OutOfRange'
        else:
            Asite = 'OutOfRange'
        return Asite


    
    InFile = file(inputFile,'r')
    outFilePlus = open(inputFile.replace('_for'+sizeRange+'.sam', '.Asite_'+sizeRange+'_P.bed'), 'w')
    flag = 0
    minsize = int(sizeRange.split('to')[0])
    maxsize = int(sizeRange.split('to')[1])
    # iterate through and process all lines in input file
    for i,line in enumerate(open(inputFile)):
        if i%100000 == 0:
            print i


        clist = line.replace('\n','').split('\t')
        # only process line if the FLAG is 0 (plus strand primary) or 256 (plus strand not primary)
        if clist[0] != '@':
        
            if clist[1] == '0' or clist[1] == '256':

                #pull out the 5' end of the sequening read
                fiveP = int(clist[3])
                # pull out the CIGAR string
                cigar = clist[5]
                # for now ignore reads with indels and split alignments
                if 'D' not in cigar and 'I' not in cigar:
                    if 'S' not in cigar:
                        # parse CIGAR string
                        noS = cigar
                    # if there is soft clipping only at 5' end
                    elif cigar[-1] != 'S' and (cigar[1] == 'S' or cigar[2] == 'S'):
                        # Split by S then split by M
                        noS = cigar.split('S')[1]
                    # if there is soft clipping at the 3' end (doesn't matter whether 5' end is soft-clipped or not because cigar.split('S')[-2] will give the M part of the string either way
                    elif cigar[-1] == 'S':
                        noS = cigar.split('S')[-2]
                    else:
                        print(cigar)
                        continue
                    read_size = getReadSize(noS)
                    ref_size = getRefSize(noS)
                    Asite = getAsitePlus(noS,read_size,fiveP,ref_size, minsize, maxsize)

                    # Now we output the adjusted read in BED format
                    if Asite != 'OutOfRange':
                        start=Asite-1
                        outFilePlus.write(clist[2]+'\t%s\t%s\n' % (start, Asite))
                else:
                    continue



##############################################
if __name__ == "__main__":
    sys.exit(main())
