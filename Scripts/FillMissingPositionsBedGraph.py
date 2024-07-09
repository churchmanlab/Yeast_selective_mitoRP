HELP_STRING = """
FillMissingPositionsBedGraph.py


Author: Mary Couvillion
Date: 6/16/20 updated from MakeMaxSubcodonProfiles_flat.py to split into fillMissingPositions and maxFrame
Date: 11/10/22 updated to use with yeast genome (plus file only)

     -h     print this help message
     -p     plus strand bedGraph (required)

     
command line, e.g.: python ./Scripts/FillMissingPositionsBedGraph.py -p Gal_30m_1_Mito_mRNA.noDups.Asite_31to33_P.bedGraph

"""
import sys
import itertools
from getopt import getopt



def main(argv=None):
    if argv is None:
        argv = sys.argv

    plusFile = ""
    
    try:
        optlist, args = getopt(argv[1:], "hp:")
    except:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    if len(optlist) == 0:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    for (opt, opt_arg) in optlist:
        if opt == '-h':
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == '-p':
            plusFile = opt_arg

    if plusFile == '':
        print HELP_STRING
        sys.exit(1)

    PlusFile = open(plusFile,'r')
  
    tempPlusFile = open(plusFile.replace('.bedGraph','All.bedGraph'), 'w')



    # First fill in missing positions in bedGraph file
    with PlusFile as f:
        for line in f:
            first = line.split()[0]
            if first=='chrM':
                value = float(line.split()[3])
                for position in range(int(line.split()[1])+1,int(line.split()[2])+1):
                    tempPlusFile.write(line.split()[0] + '\t' + str(position - 1) + '\t' + str(position) + '\t' + str(value) + '\n')
            else:
                tempPlusFile.write(line)



##############################################
if __name__ == "__main__":
    sys.exit(main())
