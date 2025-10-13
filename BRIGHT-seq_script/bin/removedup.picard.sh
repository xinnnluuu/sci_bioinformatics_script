#!/usr/bin/bash
#
# Program: use picard to remove PCR duplicates 
#
# History: 2022.01.14  created
# 2024.4.10  fix  path of picard
# 2025.5.08  fix  path of picard
#
# Usage: bash *.sh $sorted.inputbam $outputbam $prefix

if test $3
then
	java -jar /home/luziang/project/my_miniconda/pkgs/picard-2.21.1-0/share/picard-2.21.1-0/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$1 O=$2 M=$3
	#picard MarkDuplicates REMOVE_DUPLICATES=true I=$1 O=$2 M=$3
else
	echo "parameter needed is incomplete"
fi
