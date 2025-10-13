#!/usr/bin/bash
#
# Program: use bowtie2 to align read with nomal setting 
#
# History: 2022.03.13  created
#
# Usage:  bash *.sh ${threads} ${genome_index} ${infq1} ${infq2} ${outbam}

if test $5
then
	bowtie2 --no-unal -N 1 -p $1 -x $2 -1 $3 -2 $4 | samtools view -b -o $5 -@ $1 -
else
	echo "parameter needed is incomplete"
fi
