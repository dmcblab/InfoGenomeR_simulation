#!/bin/bash
reference="/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa"

j=1;
	cd set5
	bwa-0.7.15 mem -t 6 $reference normal21.fq normal22.fq | samtools view -bS > normal2.bam && samtools sort --threads 6 normal2.bam > normal2_sorted.bam && samtools index normal2_sorted.bam && rm normal2.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        bwa-0.7.15 mem -t 15 $reference normal5.11.fq normal5.12.fq | samtools view -bS > normal5.1.bam && samtools sort --threads 15 normal5.1.bam > normal5.1_sorted.bam && samtools index normal5.1_sorted.bam && rm normal5.1.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        bwa-0.7.15 mem -t 15 $reference normal5.21.fq normal5.22.fq | samtools view -bS > normal5.2.bam && samtools sort --threads 15 normal5.2.bam > normal5.2_sorted.bam && samtools index normal5.2_sorted.bam && rm normal5.2.bam &
        pids[${j}]=$!;
        j=$(($j+1));
	cd../

for pid in ${pids[*]}; do
        wait $pid;
done
