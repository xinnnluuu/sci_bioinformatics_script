# Script: extract_allmut.gz.py
# Description: Extract and filter mutation data from a gzipped input file
# Usage: python extract_allmut.gz.py <input.gz> <output.txt>
# This script processes mutation data and filters based on specific criteria:
# - Mutation frequency > 10%
# - Coverage depth >= 3
# - At least 2 supporting reads for alternative bases

import sys
import re  # For regular expression operations
import gzip  # For handling gzipped files

# Get command line arguments
infilename = sys.argv[1]  # Input gzipped file
outfilename = sys.argv[2]  # Output text file

# Open input gzipped file in text mode and output file
infile = gzip.open(infilename,'rt')
outfile = open(outfilename,'wt')

# Process the file line by line
for line in infile:
    # Skip lines that don't start with "chr"
    if not re.match("^chr",line):
        continue
    
    # Split line into fields by tab
    info = re.split("\t",line)
    
    # Skip if reference base is N
    if info[2] == "N":
        continue
    
    # Filter conditions:
    # info[8]: mutation frequency > 10%
    # info[4]: coverage depth >= 3
    if float(info[8]) > 0.10 and int(info[4]) >= 3:
        # Extract numbers from the base count field (info[6])
        # Format expected: A:x,C:y,G:z,T:w where x,y,z,w are counts
        num_list = re.findall("[0-9]+",info[6])
        
        # Check mutation type and supporting reads
        # For each reference base, sum the counts of alternative bases
        # Write to output if at least 2 supporting reads for alternatives
        if info[2] == "A" and (int(num_list[1])+int(num_list[2])+int(num_list[3]))>=2:
            # Sum of C,G,T reads when ref is A
            outfile.write(info[0]+"\t"+info[1]+"\t"+info[2]+"\t"+info[6]+"\n")
        elif info[2] == "C" and (int(num_list[0])+int(num_list[2])+int(num_list[3]))>=2:
            # Sum of A,G,T reads when ref is C
            outfile.write(info[0]+"\t"+info[1]+"\t"+info[2]+"\t"+info[6]+"\n")
        elif info[2] == "G" and (int(num_list[0])+int(num_list[1])+int(num_list[3]))>=2:
            # Sum of A,C,T reads when ref is G
            outfile.write(info[0]+"\t"+info[1]+"\t"+info[2]+"\t"+info[6]+"\n")
        elif info[2] == "T" and (int(num_list[0])+int(num_list[1])+int(num_list[2]))>=2:
            # Sum of A,C,G reads when ref is T
            outfile.write(info[0]+"\t"+info[1]+"\t"+info[2]+"\t"+info[6]+"\n")

# Close files
outfile.close()
infile.close()
