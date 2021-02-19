#!/bin/bash
for k in `seq 0 5`; do
	rm set$k\/*.fq
	rm set$k\/*.aln
done

j=1;
i=0;

	cd set$i;
        samtools merge tumor3_0.6.bam tumor1.8_sorted.bam normal0.3_sorted.bam normal0.45.1_sorted.bam normal0.45.2_sorted.bam --threads 1 && samtools index tumor3_0.6.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        samtools merge tumor3_0.75.bam tumor1.8_sorted.bam tumor0.45.1_sorted.bam normal0.3_sorted.bam  normal0.45.1_sorted.bam  --threads 1 && samtools index tumor3_0.75.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        samtools merge tumor3_0.9.bam tumor1.8_sorted.bam tumor0.45.1_sorted.bam tumor0.45.2_sorted.bam normal0.3_sorted.bam  --threads 1 && samtools index tumor3_0.9.bam &
        pids[${j}]=$!;
        j=$(($j+1));

	cd ../

for i in `seq 1 4`; do
	cd set$i
	samtools merge tumor5_0.6.bam tumor3_sorted.bam normal0.5_sorted.bam normal0.75.1_sorted.bam normal0.75.2_sorted.bam --threads 2 && samtools index tumor5_0.6.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        samtools merge tumor5_0.75.bam tumor3_sorted.bam tumor0.75.1_sorted.bam normal0.5_sorted.bam normal0.75.1_sorted.bam  --threads 2 && samtools index tumor5_0.75.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        samtools merge tumor5_0.9.bam tumor3_sorted.bam tumor0.75.1_sorted.bam tumor0.75.2_sorted.bam normal0.5_sorted.bam --threads 2 && samtools index tumor5_0.9.bam &
        pids[${j}]=$!;
        j=$(($j+1));
	cd ../
done

	mkdir f3;
	mkdir f5;
	mkdir f10;
	mkdir f15;
	mkdir f20;
        samtools merge f3/germ_simulated_sorted.bam set1/normal0.75.1_sorted.bam set2/normal0.75.1_sorted.bam set3/normal0.75.1_sorted.bam set4/normal0.75.1_sorted.bam --threads 1 && samtools index f3/germ_simulated_sorted.bam &
        pids[${j}]=$!;
        j=$(($j+1));
	samtools merge f5/germ_simulated_sorted.bam set1/normal0.5_sorted.bam set2/normal0.5_sorted.bam set3/normal0.5_sorted.bam set4/normal0.5_sorted.bam set1/normal0.75.1_sorted.bam set2/normal0.75.1_sorted.bam set3/normal0.75.1_sorted.bam set4/normal0.75.1_sorted.bam --threads 4 && samtools index f5/germ_simulated_sorted.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        samtools merge f10/germ_simulated_sorted.bam set1/normal0.5_sorted.bam set2/normal0.5_sorted.bam set3/normal0.5_sorted.bam set4/normal0.5_sorted.bam set1/normal0.75.1_sorted.bam set1/normal0.75.2_sorted.bam set2/normal0.75.1_sorted.bam set2/normal0.75.2_sorted.bam set3/normal0.75.1_sorted.bam set3/normal0.75.2_sorted.bam set4/normal0.75.1_sorted.bam set4/normal0.75.2_sorted.bam set5/normal2_sorted.bam --threads 12 && samtools index f10/germ_simulated_sorted.bam &
        pids[${j}]=$!;
        j=$(($j+1));

for pid in ${pids[*]}; do
        wait $pid;
done
