#!/bin/bash
reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa
declare -a arr;
k=0;

	for j in 3 5 10 15 20; do
                for i in 0.6 0.75 0.9; do
			mkdir f$j/p$i/vcf_somatic;
			for l in `seq 1 22`; do
				arr[$k]="samtools mpileup -t DP,AD,ADF,ADR,SP -q 10 -d 5000 -I -uf $reference f$j/germ_simulated_sorted.bam f$j/p$i/simulated_sorted.bam -r $l | bcftools call -c -v > f$j/p$i/vcf_somatic/$l.vcf";
				k=$(($k+1));
			done
			l="X";
                               arr[$k]="samtools mpileup -t DP,AD,ADF,ADR,SP -q 10 -d 5000 -I -uf $reference f$j/germ_simulated_sorted.bam f$j/p$i/simulated_sorted.bam -r $l | bcftools call -c -v > f$j/p$i/vcf_somatic/$l.vcf";
                                k=$(($k+1));

		done
        done

	k=$(($k-1));
	
        for z in `seq 0 23`; do
                eval ${arr[$k]} &
                pids[$z]=$!;
                k=$(($k -1));
        done

        while [ $k -ge 0 ];do

                for l in `seq 0 23`;do
                        n=`ps | grep ${pids[$l]} | wc -l`;
                        if [ $n -eq 0 ] && [ $k -ge 0 ];then
                                eval ${arr[$k]} &
                                pids[$l]=$!;
                                k=$(($k -1));
                        fi
                done
                sleep 600;

        done

        for z in `seq 0 23`;do
                wait ${pids[$z]};
        done
