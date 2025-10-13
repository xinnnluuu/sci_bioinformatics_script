#!/usr/bin/env python3

################################################################################
# Script Name: filter_variant_sites.py
# Description: Filter variant sites from gzipped input file based on quality criteria
#              Filters sites by:
#              - Alternative allele frequency (> 10%)
#              - Coverage depth (>= 3)
#              - Alternative allele support (>= 2 reads for non-reference bases)
# Usage:       python filter_variant_sites.py input.gz output.txt
# Input:       Gzipped tab-delimited file with variant information
#              Expected columns: chr, pos, ref_base, ..., depth, ..., alt_freq, ..., base_counts
# Output:      Tab-delimited file with: chr, pos, ref_base, base_counts
# Author:      Ziang Lu
# Date:        2025.10.13
################################################################################

import sys
import re
import gzip

################################################################################
# Check command line arguments
################################################################################

if len(sys.argv) < 3:
    print("Usage: python {} input.gz output.txt".format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)

infilename = sys.argv[1]
outfilename = sys.argv[2]

################################################################################
# Open input and output files
################################################################################

# Open gzipped input file in text mode
infile = gzip.open(infilename, 'rt')

# Open output file in text mode
outfile = open(outfilename, 'wt')

################################################################################
# Process each line and apply filters
################################################################################

for line in infile:
    # Skip lines that don't start with "chr" (e.g., header lines)
    if not re.match("^chr", line):
        continue
    
    # Split line into fields
    info = re.split("\t", line)
    
    # Skip sites with 'N' as reference base (ambiguous base)
    if info[2] == "N":
        continue
    
    # Apply quality filters:
    # - Alternative allele frequency > 10% (column 8)
    # - Total coverage depth >= 3 reads (column 4)
    if float(info[8]) > 0.10 and int(info[4]) >= 3:
        # Extract base counts from column 6
        # Expected format: numbers separated by non-numeric characters
        # Typically: A_count, C_count, G_count, T_count
        num_list = re.findall("[0-9]+", info[6])
        
        # Check if alternative alleles (non-reference bases) have >= 2 supporting reads
        # num_list indices: [0]=A, [1]=C, [2]=G, [3]=T
        
        if info[2] == "A" and (int(num_list[1]) + int(num_list[2]) + int(num_list[3])) >= 2:
            # Reference is A, check if C+G+T >= 2
            outfile.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[6] + "\n")
            
        elif info[2] == "C" and (int(num_list[0]) + int(num_list[2]) + int(num_list[3])) >= 2:
            # Reference is C, check if A+G+T >= 2
            outfile.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[6] + "\n")
            
        elif info[2] == "G" and (int(num_list[0]) + int(num_list[1]) + int(num_list[3])) >= 2:
            # Reference is G, check if A+C+T >= 2
            outfile.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[6] + "\n")
            
        elif info[2] == "T" and (int(num_list[0]) + int(num_list[1]) + int(num_list[2])) >= 2:
            # Reference is T, check if A+C+G >= 2
            outfile.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[6] + "\n")

################################################################################
# Close files
################################################################################

outfile.close()
infile.close()

################################################################################
# End of script
################################################################################
