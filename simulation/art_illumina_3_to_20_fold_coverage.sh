#!/bin/bash

i=1;

	samtools merge f15/germ_simulated_sorted.bam f10/germ_simulated_sorted.bam set5/normal5.1_sorted.bam --threads 4 && samtools index f15/germ_simulated_sorted.bam &
        pids[${i}]=$!;
        i=$(($i+1));
        samtools merge f20/germ_simulated_sorted.bam f10/germ_simulated_sorted.bam set5/normal5.1_sorted.bam set5/normal5.2_sorted.bam --threads 4 && samtools index f20/germ_simulated_sorted.bam &
       pids[${i}]=$!;
        i=$(($i+1));


for j in 0.6 0.75 0.9; do
	cd f3;
	mkdir p$j;
	cp ../set0/tumor3_$j.bam p$j/simulated_sorted.bam && cp ../set0/tumor3_$j.bam.bai p$j/simulated_sorted.bam.bai &
        pids[${i}]=$!;
        i=$(($i+1));
	cd ../

	cd f5;
	mkdir p$j;
	cp ../set1/tumor5_$j.bam p$j/simulated_sorted.bam && cp ../set1/tumor5_$j.bam.bai p$j/simulated_sorted.bam.bai &
        pids[${i}]=$!;
        i=$(($i+1));

	cd ../

	cd f10;
	mkdir p$j;
	samtools merge p$j/simulated_sorted.bam ../set1/tumor5_$j.bam ../set2/tumor5_$j.bam --threads 4 && samtools index p$j/simulated_sorted.bam &
        pids[${i}]=$!;
        i=$(($i+1));
	cd ../

        cd f15;
        mkdir p$j;
        samtools merge p$j/simulated_sorted.bam ../set1/tumor5_$j.bam ../set2/tumor5_$j.bam ../set3/tumor5_$j.bam --threads 4 && samtools index p$j/simulated_sorted.bam &
        pids[${i}]=$!;
        i=$(($i+1));
        cd ../

        cd f20;
        mkdir p$j;
        samtools merge p$j/simulated_sorted.bam ../set1/tumor5_$j.bam ../set2/tumor5_$j.bam ../set3/tumor5_$j.bam ../set4/tumor5_$j.bam --threads 4 && samtools index p$j/simulated_sorted.bam &
        pids[${i}]=$!;
        i=$(($i+1));
        cd ../


done


for pid in ${pids[*]}; do
        wait $pid;
done

