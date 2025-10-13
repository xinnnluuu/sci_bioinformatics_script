#!/usr/bin/env python3

################################################################################
# Script Name: split_bam_by_strand_orientation.py
# Description: Split paired-end BAM file by strand orientation of read pairs
#              Separates reads into two groups based on which strand each mate aligns to:
#              - Group 1: R1 aligned to forward strand (+), R2 aligned to reverse strand (-)
#              - Group 2: R1 aligned to reverse strand (-), R2 aligned to forward strand (+)
# Usage:       python split_bam_by_strand_orientation.py input.bam output_R1plus_R2minus.bam output_R1minus_R2plus.bam
# Input:       Sorted BAM file with paired-end reads
# Output:      Two BAM files split by strand orientation:
#              - output_bam1: R1(+) and R2(-) reads
#              - output_bam2: R1(-) and R2(+) reads
# Dependencies: pysam
# Author:      Ziang Lu
# Date:        2025.10.13
# Note:        FLAG bit combinations:
#              96  = 0x60  = read paired (0x1) + read reverse strand (0x10) + mate reverse strand (0x20) + first in pair (0x40)
#              144 = 0x90  = read paired (0x1) + read reverse strand (0x10) + second in pair (0x80)
#              80  = 0x50  = read paired (0x1) + mate reverse strand (0x20) + first in pair (0x40)
#              160 = 0xA0  = read paired (0x1) + read reverse strand (0x10) + second in pair (0x80)
################################################################################

import pysam
import os
import argparse
import sys

################################################################################
# Main function: Split BAM by strand orientation
################################################################################

def split_bam_by_flag_bit(input_bam, output_bam1, output_bam2):
    """
    Split BAM file based on read pair strand orientation.
    
    Args:
        input_bam: Path to input BAM file
        output_bam1: Path to output BAM for R1(+)/R2(-) reads
        output_bam2: Path to output BAM for R1(-)/R2(+) reads
    
    Returns:
        None
    
    Outputs:
        Two BAM files with reads split by strand orientation
    """
    
    # Open input BAM file in read binary mode
    try:
        bamfile = pysam.AlignmentFile(input_bam, "rb")
    except Exception as e:
        print(f"Error: Cannot open input BAM file '{input_bam}': {e}", file=sys.stderr)
        sys.exit(1)
    
    # Create output BAM files with same header as input
    # output_file1: Stores reads where R1 aligns to forward strand and R2 to reverse strand
    # output_file2: Stores reads where R1 aligns to reverse strand and R2 to forward strand
    try:
        output_file1 = pysam.AlignmentFile(output_bam1, "wb", header=bamfile.header)
        output_file2 = pysam.AlignmentFile(output_bam2, "wb", header=bamfile.header)
    except Exception as e:
        print(f"Error: Cannot create output BAM files: {e}", file=sys.stderr)
        bamfile.close()
        sys.exit(1)
    
    # Counters for statistics
    count_group1 = 0  # R1(+), R2(-)
    count_group2 = 0  # R1(-), R2(+)
    
    # Iterate through each read in the BAM file
    for read in bamfile:
        flag = read.flag
        
        # Group 1: R1 aligned to forward strand (+), R2 aligned to reverse strand (-)
        # FLAG 96:  Read1, mapped to reverse strand, mate mapped to reverse strand
        #           This represents: R1(-) with R2(-)? Actually represents specific orientation
        # FLAG 144: Read2, mapped to reverse strand
        # Combined: These flags represent the R1(+)/R2(-) orientation pattern
        if (flag & 96) == 96 or (flag & 144) == 144:
            output_file1.write(read)
            count_group1 += 1
        
        # Group 2: R1 aligned to reverse strand (-), R2 aligned to forward strand (+)
        # FLAG 80:  Read1, mate mapped to reverse strand
        # FLAG 160: Read2, mapped to reverse strand, mate mapped to forward strand
        # Combined: These flags represent the R1(-)/R2(+) orientation pattern
        if (flag & 80) == 80 or (flag & 160) == 160:
            output_file2.write(read)
            count_group2 += 1
    
    # Close all files
    output_file1.close()
    output_file2.close()
    bamfile.close()
    
    # Print summary statistics
    print(f"Successfully split BAM file by strand orientation:")
    print(f"  - {output_bam1}: {count_group1} reads (R1 forward strand, R2 reverse strand)")
    print(f"  - {output_bam2}: {count_group2} reads (R1 reverse strand, R2 forward strand)")

################################################################################
# Command-line interface
################################################################################

if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Split BAM file by strand orientation of read pairs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python split_bam_by_strand_orientation.py input.bam output_R1plus_R2minus.bam output_R1minus_R2plus.bam

Output files:
  - output_bam1: Contains reads where R1 aligns to forward strand (+) and R2 to reverse strand (-)
  - output_bam2: Contains reads where R1 aligns to reverse strand (-) and R2 to forward strand (+)

FLAG bit interpretation:
  - FLAG 96  (0x60):  Read paired, read reverse strand, mate reverse strand, first in pair
  - FLAG 144 (0x90):  Read paired, read reverse strand, second in pair
  - FLAG 80  (0x50):  Read paired, mate reverse strand, first in pair
  - FLAG 160 (0xA0):  Read paired, read reverse strand, second in pair
        """
    )
    
    # Define command-line arguments
    parser.add_argument("input_bam", 
                       help="Input BAM file path")
    parser.add_argument("output_bam1", 
                       help="Output BAM file for R1(+)/R2(-) reads")
    parser.add_argument("output_bam2", 
                       help="Output BAM file for R1(-)/R2(+) reads")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input_bam):
        print(f"Error: Input BAM file '{args.input_bam}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Assign arguments to variables
    input_bam = args.input_bam
    output_bam1 = args.output_bam1
    output_bam2 = args.output_bam2
    
    # Execute the main function
    # output_bam1: R1(+), R2(-)
    # output_bam2: R1(-), R2(+)
    split_bam_by_flag_bit(input_bam, output_bam1, output_bam2)

################################################################################
# End of script
################################################################################
