#!/bin/bash

for i in {101..118};
do
        for k in {1..3};
        do
                x=100
                s=$(($i-$x))
		STF="P9503_${i}_S${s}_L00${k}_R1_001_FP.fastq.gz"
		STR="P9503_${i}_S${s}_L00${k}_R2_001_RP.fastq.gz"
		STO="P9503_${i}_S${s}_L00${k}.sam"
		SUM="${i}_L00${k}"
		echo Running
		echo $STF
		echo $STR
		echo $STO
		hisat2 --rna-strandness RF --new-summary --summary-file /Volumes/5TB_Vaskar/Summary_HISAT2_SC/$SUM -x genome_tran \
-1 /Volumes/5TB_Vaskar/Trimmed_files_fq/$STF \
-2 /Volumes/5TB_Vaskar/Trimmed_files_fq/$STR \
-S /Volumes/5TB_Vaskar/ALIGN_SAM_SC/$STO
		echo Done
        done
done
