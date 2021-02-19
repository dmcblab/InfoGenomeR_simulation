#!/bin/bash

for i in {1..10}
do
	echo $i
	mkdir simul$i
	cp *.R *.sh *.pl *.awk simul$i/
	cd simul$i
	perl SV_chr_rand_sub_function.pl > results
	cat results | awk -f SV.results_to_SVs.awk > true_SV_sets
	wgsim -N 455445570 -1 101 -2 101 simulated.fa simulated.fq1 simulated.fq2
	bwa-0.7.15 mem  -t 40 /DATA1/AITL/hg19/ucsc.hg19.fasta simulated.fq1 simulated.fq2 | samtools view -bS >  simulated.bam
	samtools sort simulated.bam > simulated_sorted.bam
	samtools index simulated_sorted.bam
	cd ../
done
