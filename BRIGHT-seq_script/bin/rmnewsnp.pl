#!/usr/bin/perl

################################################################################
# Script Name: filter_sites_by_exclusion.pl
# Description: Filter out sites that exist in an exclusion list
#              Reads a list of sites to exclude, then filters the main input
#              to output only sites NOT in the exclusion list
# Usage:       perl filter_sites_by_exclusion.pl main_input.txt exclusion_list.txt
# Input:       
#              - ARGV[0]: Main input file (tab-delimited, chr and pos in columns 0-1)
#              - ARGV[1]: Exclusion list file (one site per line: "chr\tpos")
# Output:      Lines from main input that are NOT in the exclusion list
#              (printed to STDOUT)
# Author:      Ziang Lu
# Date:        2025.10.13
################################################################################

use strict;
use warnings;

################################################################################
# Check command line arguments
################################################################################

if (@ARGV < 2) {
    die "Usage: perl $0 main_input.txt exclusion_list.txt\n";
}

################################################################################
# Load exclusion list into hash
################################################################################

# Open exclusion list file (second argument)
open(IN1, '<', $ARGV[1]) or die "Cannot open exclusion file '$ARGV[1]': $!\n";

# Hash to store sites that should be excluded
# Key format: "chr\tpos"
my %newsnp;

while (<IN1>) {
    chomp;  # Remove newline character
    # Store each site (chr\tpos) as a hash key with value 1
    $newsnp{$_} = 1;
}

close IN1;

################################################################################
# Process main input file and filter out excluded sites
################################################################################

# Open main input file (first argument)
open(IN2, '<', $ARGV[0]) or die "Cannot open main input file '$ARGV[0]': $!\n";

while (<IN2>) {
    # Split line into fields
    @_ = split;
    
    # Construct site identifier from chromosome (column 0) and position (column 1)
    my $site = $_[0] . "\t" . $_[1];
    
    # Skip this line if the site exists in the exclusion list
    next if (exists $newsnp{$site});
    
    # Print the entire line if site is NOT in exclusion list
    print $_;
}

close IN2;

################################################################################
# End of script
################################################################################
