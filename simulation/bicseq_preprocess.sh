#!/bin/bash
reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa

mkdir bicseq_samtools
modifiedSamtools view -U BWA,bicseq_samtools/,N,N simulated_sorted.bam
#mkdir bicseq_samtools_q10
#modifiedSamtools view -U BWA,bicseq_samtools_q10/,N,N -q 10 simulated_sorted.bam

mkdir bicseq_norm
echo -e "chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm" > norm_configFile;
for i in {1..23}
do
        if [ $i == 23 ]
        then
                chr="X";
        else
                chr=$i;
        fi
        echo -e "$chr\t$reference.$chr\t/home/qlalf1457/NBICseq-norm_v0.2.4/hg19.CRG.50bp/hg19.50mer.CRC.chr$chr.txt\t$PWD/bicseq_samtools/$chr.seq\t$PWD/bicseq_norm/$chr.norm.bin" >> norm_configFile;

done

mkdir tmp;
perl ~/NBICseq-norm_v0.2.4/NBICseq-norm.pl -l 100 -s 350 norm_configFile ./NB_parameters --tmp tmp

