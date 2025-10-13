#!/usr/bin/perl

################################################################################
# Script Name: compare_pysam_with_ref.pl
# Description: Compare pysam base coverage file with reference sequence
#              Determines the most abundant base at each position and compares
#              it with the reference base
# Usage:       perl compare_pysam_with_ref.pl ref.fasta pysam.cov
# Input:       1. ref.fasta - Reference genome in FASTA format
#              2. pysam.cov - Base coverage file from pysam
# Output:      Tab-delimited file with detected base and reference base columns
# Author:      Ziang Lu
# Date:        2025.10.13
################################################################################

use strict;
use warnings;
use List::Util qw(max);

# Usage information
my $usage = <<USAGE;
Usage:
    perl $0 ref.fasta pysam.cov
    
Description:
    Compare pysam_base_cov file with reference sequence.
    
Arguments:
    ref.fasta   - Reference sequence in FASTA format
    pysam.cov   - Base coverage file from pysam
    
Output:
    Original columns from pysam.cov plus:
    - det_base: Detected base (A/C/G/T) with highest coverage
    - ref_base: Reference base at that position
USAGE

# Check command line arguments
if (@ARGV == 0) {
    die $usage;
}

################################################################################
# PART 1: Read reference genome from FASTA file
################################################################################

open(my $IN1, '<', $ARGV[0]) or die "Cannot open reference file '$ARGV[0]': $!\n";

my $chr_name;      # Current chromosome/contig name
my @seq;           # Temporary array to store sequence characters
my %ref;           # Hash to store reference sequences: {chr_name} => [base array]

while (<$IN1>) {
    if (/^>([\w\.\-]+)\s/) {
        # FASTA header line - extract chromosome name
        $chr_name = $1;
        $ref{$chr_name} = [];  # Initialize array reference for this chromosome
    } else {
        # Sequence line
        chomp;
        s/\r//g;               # Remove carriage return (for Windows format)
        @seq = split //;       # Split sequence into individual bases
        push @{$ref{$chr_name}}, @seq;  # Append bases to chromosome array
    }
}
close $IN1;

################################################################################
# PART 2: Process pysam coverage file
################################################################################

open(my $IN2, '<', $ARGV[1]) or die "Cannot open coverage file '$ARGV[1]': $!\n";

my $det_base;  # Variable to store detected base (base with highest coverage)

while (<$IN2>) {
    chomp;
    
    # Handle header line
    if (/^ref/) {
        print "$_\tdet_base\tref_base\n";
        next;
    }
    
    # Parse coverage line
    # Expected format: ref pos A count C count G count T count
    my @fields = split;
    
    # Extract coverage counts for each base
    # Assuming columns: [0]=ref, [1]=pos, [2]=A_count, [4]=C_count, [6]=G_count, [8]=T_count
    my @counts = ($fields[2], $fields[4], $fields[6], $fields[8]);
    
    # Find the maximum coverage value
    my $max_count = max @counts;
    
    # Determine which base has the maximum coverage
    if ($max_count == $fields[2]) {
        $det_base = "A";
    } elsif ($max_count == $fields[4]) {
        $det_base = "C";
    } elsif ($max_count == $fields[6]) {
        $det_base = "G";
    } elsif ($max_count == $fields[8]) {
        $det_base = "T";
    }
    
    # Output: original line + detected base + reference base
    # Note: Array index is position-1 (0-based indexing)
    print "$_\t$det_base\t$ref{$fields[0]}[$fields[1]-1]\n";
}

close $IN2;

################################################################################
# End of script
################################################################################
