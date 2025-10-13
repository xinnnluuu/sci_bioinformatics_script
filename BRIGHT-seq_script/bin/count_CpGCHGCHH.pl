#!/usr/bin/perl

################################################################################
# Script Name: count_methylation_context.pl
# Description: Count methylation sites in different sequence contexts (CpG, CHG, CHH)
#              Supports filtering by conversion type (C-to-A, C-to-T) or counting all
# Usage:       perl count_methylation_context.pl input_file conversion_type
# Input:       1. input_file - Tab-delimited file with methylation information
#              2. conversion_type - Filter type: "CtoA", "CtoT", or "ALL"
# Output:      Count summary of CpG, CHG, CHH sites and total
# Author:      Ziang Lu
# Date:        2025.10.13
################################################################################

use strict;
use warnings;

# Initialize counters for different methylation contexts
my $cpg_c = 0;  # CpG context counter
my $chg_c = 0;  # CHG context counter
my $chh_c = 0;  # CHH context counter

################################################################################
# Process input file
################################################################################

open(my $IN2, '<', $ARGV[0]) or die "Cannot open input file '$ARGV[0]': $!\n";

while (<$IN2>) {
    # Filter by conversion type specified in second argument
    if ($ARGV[1] =~ /CtoA/) {
        # C-to-A conversion: count C>A on forward strand or G>T on reverse strand
        if ($_ =~ /\tC\tA\t/ or $_ =~ /\tG\tT\t/) {
            $cpg_c++ if ($_ =~ /CpG/);
            $chg_c++ if ($_ =~ /CHG/);
            $chh_c++ if ($_ =~ /CHH/);
        }
    } elsif ($ARGV[1] =~ /CtoT/) {
        # C-to-T conversion: count C>T on forward strand or G>A on reverse strand
        if ($_ =~ /\tC\tT\t/ or $_ =~ /\tG\tA\t/) {
            $cpg_c++ if ($_ =~ /CpG/);
            $chg_c++ if ($_ =~ /CHG/);
            $chh_c++ if ($_ =~ /CHH/);
        }
    } elsif ($ARGV[1] =~ /ALL/) {
        # Count all sites regardless of conversion type
        $cpg_c++ if ($_ =~ /CpG/);
        $chg_c++ if ($_ =~ /CHG/);
        $chh_c++ if ($_ =~ /CHH/);
    }
}

close $IN2;

################################################################################
# Calculate total and output results
################################################################################

# Calculate total count across all contexts
my $total = $cpg_c + $chg_c + $chh_c;

# Print summary statistics
print "CpG\t$cpg_c\n";
print "CHG\t$chg_c\n";
print "CHH\t$chh_c\n";
print "Total\t$total\n";

################################################################################
# End of script
################################################################################
