#!/bin/bash
mkdir mpileup
for i in `seq 1 22`; do
	samtools mpileup -B -f /DATA1/AITL/hg19/ucsc.hg19.fasta simulated_sorted.bam -r chr$i > mpileup/$i.mpileup &
        pids[${i}]=$!;
done
samtools mpileup -B -f /DATA1/AITL/hg19/ucsc.hg19.fasta simulated_sorted.bam -r chrX > mpileup/X.mpileup &
pids[23]=$!;

for pid in ${pids[*]}; do
        wait $pid;
done

