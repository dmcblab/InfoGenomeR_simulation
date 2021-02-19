#!/bin/bash
reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa
ref_bam=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19_sorted.bam

export PATH=/DASstorage6/leeyh/novoBreak_distribution_v1.1.3rc:$PATH

k=1
        for j in 3 5 10 15 20; do
                cd f$j
                cp ../*.sh ./
		k=1
                for i in 0.6 0.75 0.9; do
                        cd p$i;
 	                cp ../*.sh ./
			mkdir novobreak;
			/DASstorage6/leeyh/novoBreak_distribution_v1.1.3rc/run_novoBreak_orient_partial_corrected.sh /DASstorage6/leeyh/novoBreak_distribution_v1.1.3rc $reference simulated_sorted.bam $ref_bam 8 novobreak & 
                        pids[${k}]=$!;
                        k=$(($k+1));

                        cd ../
                done
                cd ../

		for pid in ${pids[*]}; do
			wait $pid;
		done


        done
