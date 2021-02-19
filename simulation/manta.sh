#!/bin/bash
reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa
k=1
        for j in 3 5 10 15 20; do
                cd f$j
                for i in 0.6 0.75 0.9; do
                        cd p$i;
		        ~/manta-1.1.0.centos5_x86_64/bin/configManta.py --tumorBam simulated_sorted.bam --referenceFasta $reference --runDir manta
			~/manta-1.1.0.centos5_x86_64/bin/configManta.py --normalBam ../germ_simulated_sorted.bam --tumorBam simulated_sorted.bam --referenceFasta $reference --runDir manta_somatic

			cd manta
			./runWorkflow.py -m local -j 1 &
		        pids[${k}]=$!;	
		        k=$(($k+1));
			cd ../

                        cd manta_somatic
                        ./runWorkflow.py -m local -j 1 &
                        pids[${k}]=$!;
                        k=$(($k+1));
                        cd ../


                        cd ../
                done
                cd ../
        done



        for pid in ${pids[*]}; do
                wait $pid;
        done

