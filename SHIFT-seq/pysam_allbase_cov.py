#!/home/luziang/project/my_miniconda/envs/normal_01/bin/python

# Script: pysam_allbase_cov.py
# Description: Calculate per-base coverage statistics from a BAM file
# Usage: python pysam_allbase_cov.py <input.bam>
# Output: Tab-separated values with reference, position, and coverage statistics for each base

import sys
import pysam  # For handling BAM files
import numpy  # For numerical operations

# Get BAM file path from command line argument
bam = sys.argv[1]

# Open BAM file for reading
inputfile = pysam.AlignmentFile(bam,'rb')

# Get list of reference sequences from BAM file
chrlist = inputfile.references

# Define nucleotide bases to analyze
base_list = ['A','C','G','T']

# Print header for output
print("ref","pos","A","A_rate","C","C_rate","G","G_rate","T","T_rate","sum",sep = "\t")

# Iterate through each chromosome/reference sequence
for mychr in chrlist:
    # Get coverage counts for all bases in the current chromosome
    # Returns a tuple of 4 arrays (A,C,G,T) with coverage at each position
    tmp_cov_array = inputfile.count_coverage(contig=mychr)
    
    # Iterate through each position in the chromosome
    for site in range(len(tmp_cov_array[0])):
        # Print chromosome name and 1-based position
        print(mychr,site+1,sep = '\t',end = '\t')
        
        # Calculate total coverage at this position
        site_sum = 0
        for base in range(len(base_list)):
            site_sum += tmp_cov_array[base][site]
        
        # Calculate and print coverage and frequency for each base
        for base in range(len(base_list)):
            if site_sum == 0:
                site_rate = 0
            else:
                # Calculate frequency rounded to 4 decimal places
                site_rate = round(tmp_cov_array[base][site]/site_sum,4)
            print(tmp_cov_array[base][site],site_rate,sep = '\t',end = '\t')
        
        # Print total coverage for this position
        print(site_sum)
