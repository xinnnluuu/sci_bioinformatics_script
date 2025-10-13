#!/usr/bin/bash
#
# Program: use hisat3 to align reads
#
# History: 2024.10.1  created
#
# Usage:  bash *.sh $genomeindex $fq.r1 $fq.r2 $outputsam $muttype $threads

if test $6
then
	TYPE=$5 
	hisat-3n -x $1 -1 $2 -2 $3 -S $4 --base-change ${TYPE/to/,} --repeat -p $6 --no-spliced-alignment
else
	echo "parameter needed is incomplete"
	echo "bash *.sh $genomeindex $fq.r1 $fq.r2 $outputsam $muttype $threads"
fi
