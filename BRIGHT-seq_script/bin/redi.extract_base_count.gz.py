#!/usr/bin/env python3

################################################################################
# Script Name: count_base_coverage.py
# Description: Calculate total coverage for each base type from gzipped input file
#              Sums up coverage depth for A, C, G, T bases across all positions
# Usage:       python count_base_coverage.py input.gz
# Input:       Gzipped tab-delimited file with base coverage information
#              Expected columns: chr, pos, ref_base, ..., coverage_depth, ...
# Output:      Tab-delimited summary of total coverage for each base:
#              A    total_coverage
#              C    total_coverage
#              G    total_coverage
#              T    total_coverage
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
# Initialize coverage counter dictionary
################################################################################

# Dictionary to store cumulative coverage for each base type
count_dict = {
    "A": 0,
    "C": 0,
    "G": 0,
    "T": 0
}

################################################################################
# Process input file and accumulate coverage
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
    
    # Add coverage depth (column 4) to the corresponding base counter
    # info[2]: reference base (A/C/G/T)
    # info[4]: coverage depth at this position
    count_dict[info[2]] += int(info[4])

infile.close()

################################################################################
# Output results
################################################################################

# Print total coverage for each base type
for key, value in count_dict.items():
    print(f"{key}\t{value}")

################################################################################
# End of script
################################################################################
