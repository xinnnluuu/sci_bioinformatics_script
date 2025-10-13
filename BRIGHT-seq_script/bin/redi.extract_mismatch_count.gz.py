#!/usr/bin/env python3

################################################################################
# Script Name: count_base_substitutions.py
# Description: Count all possible base substitution types from gzipped input file
#              Tallies the 12 possible single nucleotide substitutions:
#              A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T, T>A, T>C, T>G
# Usage:       python count_base_substitutions.py input.gz
# Input:       Gzipped tab-delimited file with variant information
#              Expected columns: chr, pos, ref_base, ..., ..., ..., base_counts, ...
#              base_counts format: A_count,C_count,G_count,T_count
# Output:      Tab-delimited summary of substitution counts:
#              AC   count
#              AG   count
#              ... (for all 12 substitution types)
# Author:      Ziang Lu
# Date:        2025.10.13
################################################################################

import re
import sys
import gzip

################################################################################
# Check command line arguments
################################################################################

if len(sys.argv) < 2:
    print("Usage: python {} input.gz".format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)

infilename = sys.argv[1]

################################################################################
# Initialize substitution counter dictionary
################################################################################

# Dictionary to store counts for all 12 possible base substitutions
# Format: "RefBase>AltBase" -> count
count_dict = {
    "AC": 0,  # A to C
    "AG": 0,  # A to G
    "AT": 0,  # A to T
    "CA": 0,  # C to A
    "CG": 0,  # C to G
    "CT": 0,  # C to T
    "GA": 0,  # G to A
    "GC": 0,  # G to C
    "GT": 0,  # G to T
    "TA": 0,  # T to A
    "TC": 0,  # T to C
    "TG": 0   # T to G
}

################################################################################
# Process input file and accumulate substitution counts
################################################################################

# Open gzipped input file in text mode
infile = gzip.open(infilename, 'rt')

for line in infile:
    # Skip lines that don't start with "chr" (e.g., header lines)
    if not re.match("^chr", line):
        continue
    
    # Split line into fields
    info = re.split("\t", line)
    
    # Skip sites with 'N' as reference base (ambiguous base)
    if info[2] == "N":
        continue
    
    # Extract base counts from column 6
    # num_list indices: [0]=A_count, [1]=C_count, [2]=G_count, [3]=T_count
    num_list = re.findall("[0-9]+", info[6])
    
    # Accumulate substitution counts based on reference base
    # For each reference base, count substitutions to the other 3 bases
    
    if info[2] == "A":
        # Reference is A: count A>C, A>G, A>T
        count_dict["AC"] += int(num_list[1])  # A to C
        count_dict["AG"] += int(num_list[2])  # A to G
        count_dict["AT"] += int(num_list[3])  # A to T
        
    elif info[2] == "C":
        # Reference is C: count C>A, C>G, C>T
        count_dict["CA"] += int(num_list[0])  # C to A
        count_dict["CG"] += int(num_list[2])  # C to G
        count_dict["CT"] += int(num_list[3])  # C to T
        
    elif info[2] == "G":
        # Reference is G: count G>A, G>C, G>T
        count_dict["GA"] += int(num_list[0])  # G to A
        count_dict["GC"] += int(num_list[1])  # G to C
        count_dict["GT"] += int(num_list[3])  # G to T
        
    elif info[2] == "T":
        # Reference is T: count T>A, T>C, T>G
        count_dict["TA"] += int(num_list[0])  # T to A
        count_dict["TC"] += int(num_list[1])  # T to C
        count_dict["TG"] += int(num_list[2])  # T to G

infile.close()

################################################################################
# Output results
################################################################################

# Print total count for each substitution type
for key, value in count_dict.items():
    print(f"{key}\t{value}")

################################################################################
# End of script
################################################################################
