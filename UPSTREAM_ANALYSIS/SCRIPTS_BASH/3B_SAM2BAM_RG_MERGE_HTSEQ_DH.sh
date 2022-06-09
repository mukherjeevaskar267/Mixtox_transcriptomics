#!/bin/bash

r=1
for i in {119..139}        
do
	for k in {1..3};
        do
                x=100
                s=$(($i-$x))
		r=$(($r+1))
		STI="P9503_${i}_S${s}_L00${k}.sam"
		STO="P9503_${i}_S${s}_L00${k}.bam"
		STAT="P9503_${i}_S${s}_L00${k}.txt"
		echo Running
                echo $STI
                echo $STO
		echo $STAT
		samtools view -Su /Volumes/5TB_Vaskar/ALIGN_SAM_DH/$STI | samtools sort -o /Volumes/5TB_Vaskar/ALIGN_BAM_DH/$STO
		samtools index /Volumes/5TB_Vaskar/ALIGN_BAM_DH/$STO
		samtools stats /Volumes/5TB_Vaskar/ALIGN_BAM_DH/$STO > /Volumes/5TB_Vaskar/BAM_STAT_DH/$STAT
                echo "Done sam to bam"
                STI="P9503_${i}_S${s}_L00${k}.bam"
                STO="P9503_${i}_S${s}_L00${k}_RG.bam"
                ID=$(awk -v n=$r 'NR==n {print $2}' /Volumes/5TB_Vaskar/RG_assignment_DH.txt)
                PU=$(awk -v n=$r 'NR==n {print $5}' /Volumes/5TB_Vaskar/RG_assignment_DH.txt)
                SM=$(awk -v n=$r 'NR==n {print $6}' /Volumes/5TB_Vaskar/RG_assignment_DH.txt)
                RG="${ID}_${PU}_${SM}"
                echo Running
                echo $STI
                echo $STO
                echo $RG
                samtools addreplacerg -r "ID:${ID}" -r "LB:LIB1" -r "SM:${SM}" -r "PL:ILLUMINA" -r "PU:${PU}" -r "CN:NGI.SE" -r "PM:HISEQ2500" -o /Volumes/5TB_Vaskar/RG_BAM_DH/$STO /Volumes/5TB_Vaskar/ALIGN_BAM_DH/$STI
	done
	echo "RG files generated"
	j1=1
        j2=2
        j3=3
        STI1="P9503_${i}_S${s}_L00${j1}_RG.bam"
        STI2="P9503_${i}_S${s}_L00${j2}_RG.bam"
        STI3="P9503_${i}_S${s}_L00${j3}_RG.bam"
        STOM="P9503_${i}_MRG.bam"
	MSTAT="P9503_${i}_MRG_STAT.txt"
        echo $STI1
        echo $STI2
        echo $STI3
        echo $STOM
	echo $MSTAT
        echo "merging BAM files"
	samtools merge /Volumes/5TB_Vaskar/MRG_BAM_DH/$STOM /Volumes/5TB_Vaskar/RG_BAM_DH/$STI1 /Volumes/5TB_Vaskar/RG_BAM_DH/$STI2 /Volumes/5TB_Vaskar/RG_BAM_DH/$STI3
	OF="P9503_${i}.txt"
	echo "Running HTSeq_count for $STOM"
	echo "Output file = $OF"
	samtools index /Volumes/5TB_Vaskar/MRG_BAM_DH/$STOM
	samtools stats /Volumes/5TB_Vaskar/MRG_BAM_DH/$STOM > /Volumes/5TB_Vaskar/BAM_MER_STAT_DH/$MSTAT
	htseq-count -f bam -r pos -s reverse --idattr=gene_id \
/Volumes/5TB_Vaskar/MRG_BAM_DH/$STOM \
/Volumes/5TB_Vaskar/Debaryomyces_hansenii_cbs767_gca_000006445.ASM644v2.42_string_replaced.gtf > /Volumes/5TB_Vaskar/HTSEQ_DH/$OF
done
