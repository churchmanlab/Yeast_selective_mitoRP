#Sort a .sam file by length of reads

#On command line write python SortSAMbyLength.py infilename
import sys 
infilename=sys.argv[1]
min=int(sys.argv[2])
max=int(sys.argv[3])

#Open the fasta file to read
InputFile = open(infilename, 'r')

#Create a file to write
targets = list(range(min, max+1))
#outputfilenames = [infilename.replace('.sam', f'_{i}.sam') for i in targets]
outputfilenames = [infilename.replace(".sam", "_{}.sam".format(i)) for i in targets]
outputfiles = [open(outputfilename, 'w') for outputfilename in outputfilenames]


#Also could do:
#outfilename=infilename.split('.')[0]+'whatever'+infilename.split('.')[1]

count = 0   

for line in InputFile.readlines():
    
    
    if line[0] == '@':
        for outputfile in outputfiles:
            outputfile.write(line)

    if line[0] != '@':
        column = line.split()
        CIGAR = column[5]
        if CIGAR != '*' and 'D' not in CIGAR and 'I' not in CIGAR and 'N' not in CIGAR:
            if 'S' not in CIGAR or (CIGAR[1] != 'S' and CIGAR[2] != 'S'):
                seq_length = int(CIGAR.split('M')[0])
            elif CIGAR[1] == 'S' or CIGAR[2] == 'S':
                seq_length = int(CIGAR.split('S')[1].split('M')[0])

            if seq_length in targets:
                idx = targets.index(seq_length)
                outputfiles[idx].write(line)



            else:
                count = count + 1
            
            
print count, 'records are out of length bounds'

for outputfile in outputfiles:
    outputfile.close()
        

InputFile.close()
