#!/usr/bin/env python3

################################################################################
# Script Name: calculate_base_coverage.py
# Description: Calculate base-level coverage statistics from BAM file
#              Computes coverage and frequency for each base (A, C, G, T) at 
#              every position in the reference genome
# Usage:       python calculate_base_coverage.py input.bam
# Input:       BAM file (sorted and indexed)
# Output:      Tab-delimited table with:
#              - Reference name and position
#              - Coverage count and rate for each base (A, C, G, T)
#              - Total coverage at each position
# Dependencies: pysam, numpy
# Author:      Ziang Lu
# Date:        2025.10.13
################################################################################

import sys
import pysam
import numpy

################################################################################
# Check command line arguments
################################################################################

if len(sys.argv) < 2:
    print("Usage: python {} input.bam".format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)

bam_file = sys.argv[1]

################################################################################
# Open BAM file and get reference information
################################################################################

# Open BAM file in read binary mode
inputfile = pysam.AlignmentFile(bam_file, 'rb')

# Get list of reference sequences (chromosomes/contigs)
chrlist = inputfile.references

# Define base order for output
base_list = ['A', 'C', 'G', 'T']

################################################################################
# Print header
################################################################################

print("ref", "pos", "A", "A_rate", "C", "C_rate", "G", "G_rate", "T", "T_rate", 
      "sum", sep="\t")

################################################################################
# Process each chromosome/contig
################################################################################

for mychr in chrlist:
    # Get coverage array for current chromosome
    # Returns tuple of 4 arrays (A, C, G, T coverage at each position)
    tmp_cov_array = inputfile.count_coverage(contig=mychr)
    
    # Process each position in the chromosome
    for site in range(len(tmp_cov_array[0])):
        # Print chromosome name and position (1-based)
        print(mychr, site + 1, sep='\t', end='\t')
        
        # Calculate total coverage at this position
        site_sum = 0
        for base in range(len(base_list)):
            site_sum += tmp_cov_array[base][site]
        
        # Print coverage and rate for each base
        for base in range(len(base_list)):
            # Calculate base frequency (avoid division by zero)
            if site_sum == 0:
                site_rate = 0
            else:
                site_rate = round(tmp_cov_array[base][site] / site_sum, 4)
            
            # Print base count and rate
            print(tmp_cov_array[base][site], site_rate, sep='\t', end='\t')
        
        # Print total coverage for this position
        print(site_sum)

# Close BAM file
inputfile.close()

################################################################################
# End of script
################################################################################
