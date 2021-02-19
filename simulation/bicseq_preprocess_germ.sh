#!/bin/bash

reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa

mkdir bicseq_samtools_germ
modifiedSamtools view -U BWA,bicseq_samtools_germ/,N,N germ_simulated_sorted.bam
#mkdir bicseq_samtools_q10
#modifiedSamtools view -U BWA,bicseq_samtools_q10/,N,N -q 10 simulated_sorted.bam


mkdir bicseq_norm_germ
echo -e "chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm" > norm_configFile_germ;
for i in {1..23}
do
        if [ $i == 23 ]
        then
                chr="X";
        else
                chr=$i;
        fi
        echo -e "$chr\t$reference.$chr\t/home/qlalf1457/NBICseq-norm_v0.2.4/hg19.CRG.50bp/hg19.50mer.CRC.chr$chr.txt\t$PWD/bicseq_samtools_germ/$chr.seq\t$PWD/bicseq_norm_germ/$chr.norm.bin" >> norm_configFile_germ;

done

mkdir tmp
perl ~/NBICseq-norm_v0.2.4/NBICseq-norm.pl -l 100 -s 350 norm_configFile_germ ./NB_parameters_germ --tmp tmp

