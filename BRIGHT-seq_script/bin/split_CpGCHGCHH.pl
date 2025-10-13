#!/usr/bin/perl

################################################################################
# Script Name: annotate_cytosine_mutations_with_context.pl
# Description: Annotate cytosine mutation sites with methylation context information
#              Identifies the most frequent non-reference base and classifies
#              cytosine sites into CpG, CHG, or CHH contexts based on reference files
# Usage:       perl annotate_cytosine_mutations_with_context.pl cpg_ref.txt chg_ref.txt plus_strand.txt minus_strand.txt
# Input:       
#              - ARGV[0]: CpG context reference file (chr, pos, strand, context, sequence)
#              - ARGV[1]: CHG context reference file (chr, pos, strand, context, sequence)
#              - ARGV[2]: Plus strand mutation sites (chr, pos, ref_base, base_counts)
#              - ARGV[3]: Minus strand mutation sites (chr, pos, ref_base, base_counts)
#              - base_counts format: [A_count, C_count, G_count, T_count]
# Output:      Tab-delimited file with columns:
#              chr, pos, strand, ref, mut, mut_rate, ACGT_count, C_type
#              (printed to STDOUT)
# Author:      Ziang Lu
# Date:        2025.10.13
# Note:        For plus strand (+): C mutations, most frequent among A/G/T
#              For minus strand (-): G mutations, most frequent among T/C/A
################################################################################

use strict;
use warnings;

################################################################################
# Check command line arguments
################################################################################

if (@ARGV < 4) {
    die "Usage: perl $0 cpg_ref.txt chg_ref.txt plus_strand.txt minus_strand.txt\n";
}

################################################################################
# Load CpG context reference sites
################################################################################

# Open CpG reference file
open(IN1, '<', $ARGV[0]) or die "Cannot open CpG reference file '$ARGV[0]': $!\n";

# Hash to store CpG sites
# Structure: $cpg{chromosome}{position} = 0
my %cpg;

while (<IN1>) {
    chomp;
    # Parse reference file format: chr, pos, strand, context, sequence
    my ($chr, $pos, $strand, $context, $sequence) = split(/\t/);
    
    # If context contains "CG", mark this position as CpG site
    if ($context =~ /CG/) {
        $cpg{$chr}{$pos} = 0;
    }
}

close IN1;

################################################################################
# Load CHG context reference sites
################################################################################

# Open CHG reference file
open(IN2, '<', $ARGV[1]) or die "Cannot open CHG reference file '$ARGV[1]': $!\n";

# Hash to store CHG sites
# Structure: $chg{chromosome}{position} = 0
my %chg;

while (<IN2>) {
    chomp;
    # Parse reference file format: chr, pos, strand, context, sequence
    my ($chr, $pos, $strand, $context, $sequence) = split(/\t/);
    
    # If context contains "CHG", mark this position as CHG site
    if ($context =~ /CHG/) {
        $chg{$chr}{$pos} = 0;
    }
}

close IN2;

################################################################################
# Print output header
################################################################################

print "chr\tpos\tstrand\tref\tmut\tmut_rate\tACGT_count\tC_type\n";

################################################################################
# Declare variables for mutation analysis
################################################################################

my ($a_count, $c_count, $g_count, $t_count, $mut_base, $mut_rate, $c_type);

################################################################################
# Process plus strand (+) mutation sites (C mutations)
################################################################################

# Open plus strand mutation file
open(IN4, '<', $ARGV[2]) or die "Cannot open plus strand file '$ARGV[2]': $!\n";

while (<IN4>) {
    chomp;
    @_ = split(/\t/);
    my $chr = $_[0];
    my $pos = $_[1];
    
    # Extract base counts from format: [A_count, C_count, G_count, T_count]
    if (/\[(\d*)\,\s(\d*)\,\s(\d*)\,\s(\d*)\]/) {
        $a_count = $1;
        $c_count = $2;
        $g_count = $3;
        $t_count = $4;
    }
    
    # Determine the most frequent non-reference base (excluding C)
    # For plus strand C sites, check A, G, T frequencies
    if ($a_count >= $g_count and $a_count >= $t_count) {
        $mut_base = "A";
        $mut_rate = sprintf("%.4f", $a_count / ($a_count + $c_count + $g_count + $t_count));
    } elsif ($g_count >= $a_count and $g_count >= $t_count) {
        $mut_base = "G";
        $mut_rate = sprintf("%.4f", $g_count / ($a_count + $c_count + $g_count + $t_count));
    } elsif ($t_count >= $a_count and $t_count >= $g_count) {
        $mut_base = "T";
        $mut_rate = sprintf("%.4f", $t_count / ($a_count + $c_count + $g_count + $t_count));
    }
    
    # Classify cytosine context and output
    # Only process if reference base is C
    if (exists $cpg{$chr}{$pos} and $_[2] =~ /C/) {
        $c_type = "CpG";
        print "$_[0]\t$_[1]\t+\t$_[2]\t$mut_base\t$mut_rate\t$_[3]\t$c_type\n";
    } elsif (exists $chg{$chr}{$pos} and $_[2] =~ /C/) {
        $c_type = "CHG";
        print "$_[0]\t$_[1]\t+\t$_[2]\t$mut_base\t$mut_rate\t$_[3]\t$c_type\n";
    } elsif ($_[2] =~ /C/) {
        $c_type = "CHH";
        print "$_[0]\t$_[1]\t+\t$_[2]\t$mut_base\t$mut_rate\t$_[3]\t$c_type\n";
    }
}

close IN4;

################################################################################
# Process minus strand (-) mutation sites (G mutations)
################################################################################

# Open minus strand mutation file
open(IN5, '<', $ARGV[3]) or die "Cannot open minus strand file '$ARGV[3]': $!\n";

while (<IN5>) {
    chomp;
    @_ = split(/\t/);
    my $chr = $_[0];
    my $pos = $_[1];
    
    # Extract base counts from format: [A_count, C_count, G_count, T_count]
    if (/\[(\d*)\,\s(\d*)\,\s(\d*)\,\s(\d*)\]/) {
        $a_count = $1;
        $c_count = $2;
        $g_count = $3;
        $t_count = $4;
    }
    
    # Determine the most frequent non-reference base (excluding G)
    # For minus strand G sites, check T, C, A frequencies
    # Note: Priority order is T, C, A (different from plus strand)
    if ($t_count >= $c_count and $t_count >= $a_count) {
        $mut_base = "T";
        $mut_rate = sprintf("%.4f", $t_count / ($a_count + $c_count + $g_count + $t_count));
    } elsif ($c_count >= $t_count and $c_count >= $a_count) {
        $mut_base = "C";
        $mut_rate = sprintf("%.4f", $c_count / ($a_count + $c_count + $g_count + $t_count));
    } elsif ($a_count >= $t_count and $a_count >= $c_count) {
        $mut_base = "A";
        $mut_rate = sprintf("%.4f", $a_count / ($a_count + $c_count + $g_count + $t_count));
    }
    
    # Classify cytosine context and output
    # Only process if reference base is G (complement of C on minus strand)
    if (exists $cpg{$chr}{$pos} and $_[2] =~ /G/) {
        $c_type = "CpG";
        print "$_[0]\t$_[1]\t-\t$_[2]\t$mut_base\t$mut_rate\t$_[3]\t$c_type\n";
    } elsif (exists $chg{$chr}{$pos} and $_[2] =~ /G/) {
        $c_type = "CHG";
        print "$_[0]\t$_[1]\t-\t$_[2]\t$mut_base\t$mut_rate\t$_[3]\t$c_type\n";
    } elsif ($_[2] =~ /G/) {
        $c_type = "CHH";
        print "$_[0]\t$_[1]\t-\t$_[2]\t$mut_base\t$mut_rate\t$_[3]\t$c_type\n";
    }
}

close IN5;

################################################################################
# End of script
################################################################################
