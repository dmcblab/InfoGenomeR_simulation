#!/bin/bash
reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa
k=1
        for j in 3 5 10 15 20; do
		k=1;
                cd f$j
                cp ../*.sh ./
                for i in 0.6 0.75 0.9; do
                        cd p$i;
	                cp ../*.sh ./

			mkdir vcf
			for l in `seq 1 22`; do
				samtools mpileup -t DP,AD,ADF,ADR,SP -q 10 -d 5000 -I -uf $reference simulated_sorted.bam -r $l | bcftools call -c -v > vcf/$l.vcf &
				pids[${k}]=$!;
				k=$(($k+1));
			done
			l="X";
			      samtools mpileup -t DP,AD,ADF,ADR,SP -q 10 -d 5000 -I -uf $reference simulated_sorted.bam -r $l | bcftools call -c -v > vcf/$l.vcf &
                                pids[${k}]=$!;
                                k=$(($k+1));

			for pid in ${pids[*]}; do
				wait $pid;
			done



                        cd ../
                done
                cd ../
        done
