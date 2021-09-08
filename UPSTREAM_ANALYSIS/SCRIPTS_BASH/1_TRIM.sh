#!/bin/bash

for i in {101..139};
do
	for k in {1..3};
	do
		x=100
		s=$(($i-$x))
		STF="P9503_${i}_S${s}_L00${k}_R1_001.fastq.gz"
		STR="P9503_${i}_S${s}_L00${k}_R2_001.fastq.gz"
		STF_P="P9503_${i}_S${s}_L00${k}_R1_001_FP.fastq.gz"
		STF_UP="P9503_${i}_S${s}_L00${k}_R1_001_FUP.fastq.gz"
		STR_P="P9503_${i}_S${s}_L00${k}_R2_001_RP.fastq.gz"
		STR_UP="P9503_${i}_S${s}_L00${k}_R2_001_RUP.fastq.gz"
		echo Running
		echo $STF
		echo $STR
		echo $STF_P
		echo $STF_UP
		echo $STR_P
		echo $STR_UP
		java -jar /Users/vaskar/Documents/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 \
/Volumes/5TB_Vaskar/Fastq_files/$STF \
/Volumes/5TB_Vaskar/Fastq_files/$STR \
/Volumes/5TB_Vaskar/FORWARD_READ_PAIRED/$STF_P \
/Volumes/5TB_Vaskar/FORWARD_READ_UNPAIRED/$STF_UP \
/Volumes/5TB_Vaskar/REVERSE_READ_PAIRED/$STR_P \
/Volumes/5TB_Vaskar/REVERSE_READ_UNPAIRED/$STR_UP \
ILLUMINACLIP:/Users/vaskar/Documents/Trimmomatic-0.38/TruSeq3-PE.fa:2:30:10 HEADCROP:12 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		echo Done
	done
done
