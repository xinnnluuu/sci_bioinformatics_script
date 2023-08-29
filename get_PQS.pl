#!/usr/bin/perl
#
#It's a script to analysis PQS in fasta sequencce. It's the origin version. PQS type output following prority rules.
#perl *.pl  *.fasta > result.txt 
#output text: sequence_name\tpqs_type\tpqs_sequence\tpqs_pos

use strict;

# make loop parts
my $loop3 = ".{1,3}?";
my $loop5 = ".{4,5}?";
my $loop7 = ".{6,7}?";
my $loop_normal = ".{1,7}?";
my $loopL_any = ".{8,12}?";
my $loopL_middle = ".{8,21}?";
# make G-repeat pars
my $G3 = "G{3,}";
my $G2 = "G{2}";

# build 4 expression with different loop length
my $PQS_loop3 = "$G3($loop3$G3)($loop_normal$G3){2}|$G3($loop_normal$G3)($loop3$G3)($loop_normal$G3)|$G3($loop_normal$G3){2}($loop3$G3)";
my $PQS_loop5 = "$G3($loop5$G3)($loop_normal$G3){2}|$G3($loop_normal$G3)($loop5$G3)($loop_normal$G3)|$G3($loop_normal$G3){2}($loop5$G3)";
my $PQS_loop7 = "$G3($loop7$G3)($loop_normal$G3){2}|$G3($loop_normal$G3)($loop7$G3)($loop_normal$G3)|$G3($loop_normal$G3){2}($loop7$G3)";
my $PQS_loop_long = "$G3($loopL_any$G3)($loop_normal$G3){2}|$G3($loop_normal$G3)($loopL_middle$G3)($loop_normal$G3)|$G3($loop_normal$G3){2}($loopL_any$G3)";

# make bulge parts
my $bulge1G = "(GG[^G]GG?|GG?[^G]GG)";
my $bulge5G = "(GG[^G]{1,5}GG?|GG?[^G]{1,5}GG)";
my $bulge7G = "(GG[^G]{1,7}GG?|GG?[^G]{1,7}GG)";

# build simple bulge expression
my $PQS_bulge_single_b7 = "$bulge7G($loop_normal$G3){3}|$G3$bulge7G($loop_normal$G3){2}|$G3($loop_normal$G3)$bulge7G($loop_normal$G3)|$G3($loop_normal$G3){2}$bulge7G";

my $PQS_bulge_mutiple_b1_n2 = "$bulge1G($loop_normal$bulge1G)($loop_normal$G3){2}|$bulge1G($loop_normal$G3)($loop_normal$bulge1G)($loop_normal$G3)|$bulge1G($loop_normal$G3){2}($loop_normal$bulge1G)|$G3($loop_normal$bulge1G){2}($loop_normal$G3)|$G3($loop_normal$bulge1G)($loop_normal$G3)($loop_normal$bulge1G)|$G3($loop_normal$G3)($loop_normal$bulge1G){2}";
my $PQS_bulge_mutiple_b1_n3 = "$bulge1G($loop_normal$bulge1G){2}($loop_normal$G3)|$bulge1G($loop_normal$bulge1G)($loop_normal$G3)($loop_normal$bulge1G)|$bulge1G($loop_normal$G3)($loop_normal$bulge1G){2}|$G3($loop_normal$bulge1G){3}";
my $PQS_bulge_mutiple_b1_n4 = "$bulge1G($loop_normal$bulge1G){3}";
my $PQS_bulge_mutiple_b1 = "$PQS_bulge_mutiple_b1_n2|$PQS_bulge_mutiple_b1_n3|$PQS_bulge_mutiple_b1_n4";

# build complex bulge expression
my $PQS_bulge_mutiple_b5_n2 = "$bulge5G($loop_normal$bulge5G)($loop_normal$G3){2}|$bulge5G($loop_normal$G3)($loop_normal$bulge5G)($loop_normal$G3)|$bulge5G($loop_normal$G3){2}($loop_normal$bulge5G)|$G3($loop_normal$bulge5G){2}($loop_normal$G3)|$G3($loop_normal$bulge5G)($loop_normal$G3)($loop_normal$bulge5G)|$G3($loop_normal$G3)($loop_normal$bulge5G){2}";
my $PQS_bulge_mutiple_b5_n3 = "$bulge5G($loop_normal$bulge5G){2}($loop_normal$G3)|$bulge5G($loop_normal$bulge5G)($loop_normal$G3)($loop_normal$bulge5G)|$bulge5G($loop_normal$G3)($loop_normal$bulge5G){2}|$G3($loop_normal$bulge5G){3}";
my $PQS_bulge_mutiple_b5_n4 = "$bulge5G($loop_normal$bulge5G){3}";
my $PQS_bulge_mutiple_b5 = "$PQS_bulge_mutiple_b5_n2|$PQS_bulge_mutiple_b5_n3|$PQS_bulge_mutiple_b5_n4";

# build two-tetrad expression
my $PQS_tetrad_G2 = "$G2($loop_normal$G2){3}";

#match and output

open IN,@ARGV[0] or die$!;
my $pos;
while (<IN>) {
	# get name and sequence
	if (/^>(.*)/){
	$pos = $1;
	next;
	};
	my @loop3 = ();
	my @loop5 = ();
	my @loop7 = ();
	my @loop_long = ();
	my @bulge_single_b7 = ();
	my @bulge_mutiple_b1 = ();
	my @bulge_mutiple_b5 = ();
	my @tetrad_G2 = ();
	my @site_pos = ();
	#loop3 match
	while ($_ =~ m/($PQS_loop3)/gi) {
		push @loop3, $1;
		push @site_pos,pos($_);
	}
	#loop5 match and loop3 output
	if (not @loop3) {
		while ($_ =~ m/($PQS_loop5)/gi) {
			push @loop5, $1;
			push @site_pos,pos($_);
		}
	}else{
		my $outbase = join ",",@loop3;
		my $outpos = join ",",@site_pos;
		print "$pos\tloop3\t$outbase\t$outpos\n";
		next
	}
	#loop7 match and loop5 output
	if (not @loop5) {
		while ($_ =~ m/($PQS_loop7)/gi) {
			push @loop7, $1;
			push @site_pos,pos($_);
		}
	}else{
		my $outbase = join ",",@loop5;
		my $outpos = join ",",@site_pos;
		print "$pos\tloop5\t$outbase\t$outpos\n";
		next
	}
	#long loop match and loop7 output
	if (not @loop7) {
		while ($_ =~ m/($PQS_loop_long)/gi) {
			push @loop_long, $1;
			push @site_pos,pos($_);
		}
	}else{
		my $outbase = join ",",@loop7;
		my $outpos = join ",",@site_pos;
		print "$pos\tloop7\t$outbase\t$outpos\n";
		next
	}
	#simple bulge match and long loop output

	if (not @loop_long) {
		while ($_ =~ m/($PQS_bulge_single_b7)/gi) {
			push @bulge_single_b7, $1;
			push @site_pos,pos($_);
		}
		while ($_ =~ m/($PQS_bulge_mutiple_b1)/gi) {
			push @bulge_mutiple_b1, $1;
			push @site_pos,pos($_);
		}
	}else{
		my $outbase = join ",",@loop_long;
		my $outpos = join ",",@site_pos;
		print "$pos\tlong_loop\t$outbase\t$outpos\n";
		next
	}
	#complex bulge/two-tetrad match and simple bulge output
	if (not((@bulge_single_b7) or (@bulge_mutiple_b1))) {
		while ($_ =~ m/($PQS_bulge_mutiple_b5)/gi) {
			push @bulge_mutiple_b5, $1;
			push @site_pos,pos($_);
		}
		while ($_ =~ m/($PQS_tetrad_G2)/gi) {
			push @tetrad_G2, $1;
			push @site_pos,pos($_);
		}
	}elsif (@bulge_single_b7){
		my $outbase = join ",",@bulge_single_b7;
		my $outpos = join ",",@site_pos;
		print "$pos\tbulge_single_b7\t$outbase\t$outpos\n";
		next
	}elsif (@bulge_mutiple_b1){
		my $outbase = join ",",@bulge_mutiple_b1;
		my $outpos = join ",",@site_pos;
		print "$pos\tbulge_mutiple_b1\t$outbase\t$outpos\n";
		next
	}
	#complex bulge and two-tetrad output
	if (not((@bulge_mutiple_b5) or (@tetrad_G2))) {
		print "$pos\tother\n";
	}elsif (@bulge_mutiple_b5) {
		my $outbase = join ",",@bulge_mutiple_b5;
		my $outpos = join ",",@site_pos;
		print "$pos\tbulge_mutiple_b5\t$outbase\t$outpos\n";
		next
	}elsif (@tetrad_G2) {
		my $outbase = join ",",@tetrad_G2;
		my $outpos = join ",",@site_pos;
		print "$pos\ttetrad_G2\t$outbase\t$outpos\n";
		next
	}
}


