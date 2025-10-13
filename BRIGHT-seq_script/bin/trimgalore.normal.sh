#!/usr/bin/bash
#
# Program: use trimgalore to trim reads
#
# History: 2022.01.14 created
#
# Usage: bash *.sh $fq.r1 $fq.r2 $outputdir  $threads
#	output $prefix_val_1.fq.gz $prefix_val_2.fq.gz

if test $4
then
	trim_galore -j $4 -q 20 --stringency 3 --length 20  -e 0.1 --gzip --paired $1 $2 -o $3 
else
	echo "parameter needed is incomplete"
fi
