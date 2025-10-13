#!/usr/bin/env python3

################################################################################
# Script Name: filter_variant_sites_strict.py
# Description: Filter variant sites with stricter quality criteria than previous version
#              Applies more stringent filters:
#              - Alternative allele frequency > 5%
#              - Alternative allele support >= 5 reads (vs. 2 in previous script)
# Usage:       python filter_variant_sites_strict.py input.gz output.txt
# Input:       Gzipped tab-delimited file with variant information
#              Expected columns: chr, pos, ref_base, ..., ..., ..., base_counts, ..., alt_freq, ...
#              base_counts format: A_count,C_count,G_count,T_count
# Output:      Tab-delimited file with: chr, pos, ref_base, base_counts
#              Only sites passing stricter quality filters
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
# Process each line and apply strict quality filters
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
    
    # Apply strict quality filter:
    # - Alternative allele frequency > 5% (column 8)
    # Note: This is less strict than the 10% threshold in the previous script
    if float(info[8]) > 0.05:
        # Extract base counts from column 6
        # num_list indices: [0]=A_count, [1]=C_count, [2]=G_count, [3]=T_count
        num_list = re.findall("[0-9]+", info[6])
        
        # Check if alternative alleles (non-reference bases) have >= 5 supporting reads
        # This is stricter than the >= 2 reads threshold in the previous script
        
        if info[2] == "A" and (int(num_list[1]) + int(num_list[2]) + int(num_list[3])) >= 5:
            # Reference is A: check if C+G+T >= 5
            outfile.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[6] + "\n")
            
        elif info[2] == "C" and (int(num_list[0]) + int(num_list[2]) + int(num_list[3])) >= 5:
            # Reference is C: check if A+G+T >= 5
            outfile.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[6] + "\n")
            
        elif info[2] == "G" and (int(num_list[0]) + int(num_list[1]) + int(num_list[3])) >= 5:
            # Reference is G: check if A+C+T >= 5
            outfile.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[6] + "\n")
            
        elif info[2] == "T" and (int(num_list[0]) + int(num_list[1]) + int(num_list[2])) >= 5:
            # Reference is T: check if A+C+G >= 5
            outfile.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[6] + "\n")

################################################################################
# Close files
################################################################################

outfile.close()
infile.close()

################################################################################
# End of script
################################################################################
