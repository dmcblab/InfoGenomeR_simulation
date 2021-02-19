#!/bin/bash
for i in `seq 1 5`; do
        bwa-0.7.15 mem  -t 4 /DATA1/AITL/hg19/ucsc.hg19.fasta $1$i\1.fq $1$i\2.fq | samtools view -bS >  $1$i.bam &
        pids[${i}]=$!;
done
        bwa-0.7.15 mem  -t 8 /DATA1/AITL/hg19/ucsc.hg19.fasta germ_simulated1.fq germ_simulated2.fq | samtools view -bS >  germ_simulated.bam &
        pids[6]=$!;

for pid in ${pids[*]}; do
        wait $pid;
done

for i in `seq 1 5`; do
        samtools sort --threads 4 simulated$i.bam > simulated$i\_sorted.bam &
        pids[${i}]=$!;
done
        samtools sort --threads 8 germ_simulated.bam > germ_simulated_sorted.bam &
        pids[6]=$!;

for pid in ${pids[*]}; do
        wait $pid;
done

for i in `seq 1 5`; do
        mv simulated$i\_sorted.bam simulated$i.bam
done
	mv simulated1.bam simulated_f5.bam

	i=1;
        samtools index germ_simulated_sorted.bam &
        pids[${i}]=$!;
        i=2;
        samtools index $1\_f5.bam &
        pids[${i}]=$!;
	samtools merge $1\_f10.bam $1\_f5.bam $1\2.bam --threads 30
        i=3;
        samtools index $1\_f10.bam &
        pids[${i}]=$!;
        samtools merge $1\_f15.bam $1\_f10.bam $1\3.bam --threads 30
        i=4;
        samtools index $1\_f15.bam &
        pids[${i}]=$!;
        samtools merge $1\_f20.bam $1\_f15.bam $1\4.bam --threads 30
        i=5;
        samtools index $1\_f20.bam &
        pids[${i}]=$!;
        samtools merge $1\_f25.bam $1\_f20.bam $1\5.bam --threads 30
        i=6;
        samtools index $1\_f25.bam &
        pids[${i}]=$!;


for pid in ${pids[*]}; do
        wait $pid;
done

