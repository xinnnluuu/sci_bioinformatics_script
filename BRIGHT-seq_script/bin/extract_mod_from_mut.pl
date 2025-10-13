#!/usr/bin/perl

################################################################################
# Script Name: filter_conversion_sites.pl
# Description: Filter conversion sites by type from methylation data
#              Extracts specific base conversion patterns based on strand orientation
# Usage:       perl filter_conversion_sites.pl input_file conversion_type
# Input:       1. input_file - Tab-delimited file with conversion information
#                   Expected to have strand (+/-), reference base, and converted base
#              2. conversion_type - Filter type: "CtoA" or "CtoT"
# Output:      Filtered lines matching the specified conversion pattern
# Author:      Ziang Lu
# Date:        2025.10.13
################################################################################

use strict;
use warnings;

################################################################################
# Process input file and filter by conversion type
################################################################################

open(my $IN, '<', $ARGV[0]) or die "Cannot open input file '$ARGV[0]': $!\n";

while (<$IN>) {
    # Skip header line
    next if (/^chr\t/);
    
    # Filter based on conversion type
    if ($ARGV[1] =~ /CtoA/) {
        # C-to-A conversion filter:
        # - Forward strand (+): C > A
        # - Reverse strand (-): G > T (complementary to C > A)
        print if (/(\+\tC\tA|\-\tG\tT)/);
        
    } elsif ($ARGV[1] =~ /CtoT/) {
        # C-to-T conversion filter:
        # - Forward strand (+): C > T
        # - Reverse strand (-): G > A (complementary to C > T)
        print if (/(\+\tC\tT|\-\tG\tA)/);
    }
}

close $IN;

################################################################################
# End of script
################################################################################
