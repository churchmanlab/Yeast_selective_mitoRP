HELP_STRING = """
SAM2lengthBED.py

Author: Mary Couvillion
Date: December 21, 2016
Updated May 14, 2020 to handle soft-cipped reads: ignore the soft-clipped nucleotides unless they are part of the poly(A) tail
Updated November 8, 2022 for yeast, to handle introns. Delete minus strand

Make a text file that lists 5' position, size, and strand of read.
To interpret the SAM file: http://samtools.sourceforge.net/SAMv1.pdf and http://picard.sourceforge.net/explain-flags.html

     -h     print this help message
     -f     file input (required)
     -p     plus output file
"""
import sys
from getopt import getopt
import re

def main(argv=None):
    if argv is None:
        argv = sys.argv

    inputFile = ""
    outputFilePlus = ""
#     outputFileMinus = ""
    try:
        optlist, args = getopt(argv[1:], "hf:p:")
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
        elif opt == '-f':
            inputFile = opt_arg
        elif opt == '-p':
            outputFilePlus = opt_arg
    if inputFile == "" :
        print 'Problem with input file'
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


    InFile = file(inputFile,'r')
    outFilePlus = open(outputFilePlus, 'w')
#     outFileMinus = open(outputFileMinus, 'w')


    # iterate through and process all lines in input file
    for i,line in enumerate(open(inputFile)):
        if i%100000 == 0:
            print i
        line = line.strip()
        if line[0] == '@':
            continue
        
        clist = line.replace('\n','').split('\t')
        # only process line if the FLAG is 0, otherwise it is on the other strand (I think...) FLAG = 16 for reverse complement
        if clist[1] == '0' or clist[1] == '256':
            strand = '+'

            #pull out the 5' end of the sequening read
            fiveP = int(clist[3])
            # pull out the CIGAR string
            cigar = clist[5]
            # parse CIGAR string
            if 'D' not in cigar and 'I' not in cigar:
                if 'S' not in cigar:
                    noS = cigar
                    read_size = getReadSize(noS)
                # if there is soft clipping only at 5' end
                elif cigar[-1] != 'S' and (cigar[1] == 'S' or cigar[2] == 'S'):
                    # Split by S then split by M
                    noS = cigar.split('S')[1]
                    read_size = getReadSize(noS)
                    
                # if there is soft clipping at the 3' end (doesn't matter whether 5' end is soft-clipped or not because cigar.split('S')[-2] will give the M part of the string either way
                elif cigar[-1] == 'S':
                    tailLength = int(cigar.split('M')[-1].split('S')[0])
                    noS = cigar.split('S')[-2]
                    seq = clist[9]
                    tail = seq[-tailLength:]
                    count = 0
                    for letter in tail:
                        if letter != 'A':
                            count = count + 1                            
                        elif letter == 'A':
                            count = count
                        # If the tail is all As, add that length 
                    if count == 0:
                        read_size = getReadSize(noS) + tailLength
                    # If not, ignore tail
                    elif count != 0 :
                        read_size = getReadSize(noS)

                else:
                    continue
                lst = re.split('M|N', noS)[0:-1]
                refsize = sum([int(i) for i in lst])
                threeP = fiveP + refsize - 1

            else:
                continue
            # Now we output the adjusted read in BED format
            outFilePlus.write(clist[2]+'\t%s\t%s\t%s\n' % (threeP, read_size, strand))
            
        else:
            continue

##############################################
if __name__ == "__main__":
    sys.exit(main())
