#!/usr/bin/perl
################################################################################
# Script Name: annotate_cytosine_mutations_with_context.pl
# Description: Annotate cytosine mutation sites with methylation context information.
#              This script performs the following tasks:
#              1. Parses base counts from input, supporting '*' (indicating 0).
#              2. Identifies the most frequent non-reference mutation.
#              3. Recalculates mutation rate using: Max_Alt_Count / (Max_Alt_Count + Ref_Count).
#              4. Classifies cytosine sites into CpG, CHG, or CHH contexts.
# Usage:       perl annotate_cytosine_mutations_with_context.pl <cpg_ref> <chg_ref> <plus_strand_input> <minus_strand_input>
# Input:       
#              - ARGV[0]: CpG context reference file (chr, pos, strand, context, sequence)
#              - ARGV[1]: CHG context reference file (chr, pos, strand, context, sequence)
#              - ARGV[2]: Plus strand mutation sites (chr, pos, ref_base, base_counts...)
#              - ARGV[3]: Minus strand mutation sites (chr, pos, ref_base, base_counts...)
#              - base_counts format: [A_count, C_count, G_count, T_count] e.g. [*, 5, *, 1]
# Output:      Tab-delimited file to STDOUT with columns:
#              chr, pos, strand, ref, mut, mut_rate, ACGT_count, C_type
# Author:      Luziang
# Date:        2025.12.30
# Note:        For plus strand (+): Ref=C, Mutation calculated against A/G/T
#              For minus strand (-): Ref=G, Mutation calculated against T/C/A
################################################################################

use strict;
use warnings;

# Check if correct number of arguments are provided
if (@ARGV < 4) {
    die "Usage: perl annotate_cytosine_mutations_with_context.pl <CpG_Ref> <CHG_Ref> <Plus_Input> <Minus_Input>\n";
}

# ------------------------------------------------------------------
# Step 1: Load Context Reference Data (CpG)
# ------------------------------------------------------------------
open IN1, $ARGV[0] or die "Error opening CpG reference file $ARGV[0]: $!";
my %cpg;
while(<IN1>){
    chomp;
    my ($chr, $pos, $strand, $context, $sequence) = split(/\t/);
    if($context =~ /CG/){
        $cpg{$chr}{$pos} = 0; # Store position as key
    }
}
close IN1;

# ------------------------------------------------------------------
# Step 2: Load Context Reference Data (CHG)
# ------------------------------------------------------------------
open IN2, $ARGV[1] or die "Error opening CHG reference file $ARGV[1]: $!";
my %chg;
while(<IN2>){
    chomp;
    my ($chr, $pos, $strand, $context, $sequence) = split(/\t/);
    if($context =~/CHG/){
        $chg{$chr}{$pos} = 0;
    }
}
close IN2;

# Print Output Header
print "chr\tpos\tstrand\tref\tmut\tmut_rate\tACGT_count\tC_type\n";

my ($a_count, $c_count, $g_count, $t_count, $mut_base, $mut_rate, $c_type);

# ------------------------------------------------------------------
# Step 3: Process Plus Strand Input (Ref = C)
# ------------------------------------------------------------------
open IN4, $ARGV[2] or die "Error opening Plus Strand file $ARGV[2]: $!";
while(<IN4>){
    chomp;
    # Skip empty lines or comments
    next if /^\s*$/ || /^#/;

    my @cols = split(/\t/);
    my $chr = $cols[0];
    my $pos = $cols[1];
    my $ref_base_in_file = $cols[2];
    
    # Parse base counts: [A, C, G, T]
    # Supports format like [*, 5, *, 1]
    if(/\[(.*?)\]/){
        my $content = $1;
        $content =~ s/\s//g; # Remove whitespace
        my @counts = split(/,/, $content);
        
        # Replace '*' with 0
        for (@counts) { $_ = 0 if $_ eq '*'; }
        
        $a_count = $counts[0];
        $c_count = $counts[1];
        $g_count = $counts[2];
        $t_count = $counts[3];
    } else {
        next; # Skip lines without valid count format
    }

    # Determine dominant mutation base (Ref is C, so check A, G, T)
    my $max_alt_count = 0;
    
    if($a_count >= $g_count and $a_count >= $t_count){
        $mut_base = "A";
        $max_alt_count = $a_count;
    } elsif($g_count >= $a_count and $g_count >= $t_count){
        $mut_base = "G";
        $max_alt_count = $g_count;
    } else { # T is max
        $mut_base = "T";
        $max_alt_count = $t_count;
    }
    
    # Calculate Mutation Rate: Max_Alt / (Max_Alt + Ref_Count)
    # Ref count for plus strand is C_count
    my $denominator = $max_alt_count + $c_count;
    if ($denominator > 0) {
        $mut_rate = sprintf("%.4f", $max_alt_count / $denominator);
    } else {
        $mut_rate = 0.0000;
    }

    # Determine Context Type (CpG > CHG > CHH)
    # Keep output column 3 (original counts) as is
    my $col3 = defined $cols[3] ? $cols[3] : "NA";
    
    if(exists $cpg{$chr}{$pos} and $ref_base_in_file =~ /C/){
        $c_type="CpG";
        print "$chr\t$pos\t+\t$ref_base_in_file\t$mut_base\t$mut_rate\t$col3\t$c_type\n";
    }elsif(exists $chg{$chr}{$pos} and $ref_base_in_file =~ /C/){
        $c_type="CHG";
        print "$chr\t$pos\t+\t$ref_base_in_file\t$mut_base\t$mut_rate\t$col3\t$c_type\n";
    }elsif($ref_base_in_file =~ /C/){
        $c_type="CHH";
        print "$chr\t$pos\t+\t$ref_base_in_file\t$mut_base\t$mut_rate\t$col3\t$c_type\n";
    }
}
close IN4;

# ------------------------------------------------------------------
# Step 4: Process Minus Strand Input (Ref = G)
# ------------------------------------------------------------------
open IN5, $ARGV[3] or die "Error opening Minus Strand file $ARGV[3]: $!";
while(<IN5>){
    chomp;
    next if /^\s*$/ || /^#/;

    my @cols = split(/\t/);
    my $chr = $cols[0];
    my $pos = $cols[1];
    my $ref_base_in_file = $cols[2];

    # Parse base counts
    if(/\[(.*?)\]/){
        my $content = $1;
        $content =~ s/\s//g; 
        my @counts = split(/,/, $content);
        for (@counts) { $_ = 0 if $_ eq '*'; }
        
        $a_count = $counts[0];
        $c_count = $counts[1];
        $g_count = $counts[2];
        $t_count = $counts[3];
    } else {
        next;
    }

    # Determine dominant mutation base (Ref is G, so check T, C, A)
    my $max_alt_count = 0;
    
    if($t_count >= $c_count and $t_count >= $a_count){
        $mut_base = "T";
        $max_alt_count = $t_count;
    } elsif($c_count >= $t_count and $c_count >= $a_count){
        $mut_base = "C";
        $max_alt_count = $c_count;
    } else { # A is max
        $mut_base = "A";
        $max_alt_count = $a_count;
    }

    # Calculate Mutation Rate: Max_Alt / (Max_Alt + Ref_Count)
    # Ref count for minus strand is G_count
    my $denominator = $max_alt_count + $g_count;
    if ($denominator > 0) {
        $mut_rate = sprintf("%.4f", $max_alt_count / $denominator);
    } else {
        $mut_rate = 0.0000;
    }

    # Determine Context Type
    my $col3 = defined $cols[3] ? $cols[3] : "NA";
    
    if(exists $cpg{$chr}{$pos} and $ref_base_in_file =~ /G/){
        $c_type="CpG";
        print "$chr\t$pos\t-\t$ref_base_in_file\t$mut_base\t$mut_rate\t$col3\t$c_type\n";
    }elsif(exists $chg{$chr}{$pos} and $ref_base_in_file =~ /G/){
        $c_type="CHG";
        print "$chr\t$pos\t-\t$ref_base_in_file\t$mut_base\t$mut_rate\t$col3\t$c_type\n";
    }elsif($ref_base_in_file =~ /G/){
        $c_type="CHH";
        print "$chr\t$pos\t-\t$ref_base_in_file\t$mut_base\t$mut_rate\t$col3\t$c_type\n";
    }
}
close IN5;
