#!/bin/bash

delly_v0.7.6_CentOS5.4_x86_64bit call -t DEL -g /DAS_Storage5/leeyh/hg19/ucsc.hg19.fasta -o DEL.bcf simulated_sorted.bam &
pids[1]=$!;
delly_v0.7.6_CentOS5.4_x86_64bit call -t DUP -g /DAS_Storage5/leeyh/hg19/ucsc.hg19.fasta -o DUP.bcf simulated_sorted.bam &
pids[2]=$!;
delly_v0.7.6_CentOS5.4_x86_64bit call -t INV -g /DAS_Storage5/leeyh/hg19/ucsc.hg19.fasta -o INV.bcf simulated_sorted.bam &
pids[3]=$!;
delly_v0.7.6_CentOS5.4_x86_64bit call -t TRA -g /DAS_Storage5/leeyh/hg19/ucsc.hg19.fasta -o TRA.bcf simulated_sorted.bam &
pids[4]=$!;


for pid in ${pids[*]}; do
        wait $pid;
done


bcftools view DEL.bcf > DEL.vcf

bcftools view DUP.bcf > DUP.vcf
 bcftools view INV.bcf > INV.vcf
 bcftools view TRA.bcf > TRA.vcf


mkdir samtools_bicseq

modifiedSamtools view -U BWA,samtools_bicseq/,N,N -q 10 simulated_sorted.bam

for i in {1..23}
do
	if($i == 23)
		chr = "X";
	else
		chr = $i;

	echo $chr;

done


mkdir simulated_q10
echo -e "chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm" > norm_configFile;
for i in {1..23}
do
        if [ $i == 23 ]
        then
                chr="X";
        else
                chr=$i;
        fi
        echo -e "chr$chr\t/DATA1/AITL/hg19/chr$chr.fa\t/home/qlalf1457/NBICseq-norm_v0.2.4/hg19.CRG.50bp/hg19.50mer.CRC.chr$chr.txt\t$PWD/samtools_bicseq/chr$chr.seq\t$PWD/simulated_q10/$chr.norm.bin" >> norm_configFile;

done

perl ~/NBICseq-norm_v0.2.4/NBICseq-norm.pl -l 101 -s 500 norm_configFile ./NB_parameters

