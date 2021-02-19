#!/bin/bash

for i in {1..1}
do
	echo $i
	mkdir simul$i
	cp *.R *.sh *.pl *.awk simul$i/
	cd simul$i
	perl SV_chr_rand_sub_function_N_excluded.pl > results
	cat results | awk -f SV.results_to_SVs.awk > true_SV_sets
	./art_illumina.sh simulated
	./bwa_mem.sh simulated

	./SV_germ.sh &
        ./bicseq_preprocess_germ.sh &

	k=1;
	for j in 5 10 15 20 25
	do
		mkdir f$j
		mv simulated_f$j.bam  f$j/simulated_sorted.bam
		mv simulated_f$j.bam.bai f$j/simulated_sorted.bam.bai
		cp *.R *.sh *.pl *.awk f$j/
		cd f$j
		./SV.sh &
	        pids[${k}]=$!;
		k=$(($k+1));
	        ./bicseq_preprocess.sh &
                pids[${k}]=$!;
                k=$(($k+1));
 	       ./bam_to_NPE.sh &
                pids[${k}]=$!;
                k=$(($k+1));
		cd ../
	done

	for pid in ${pids[*]}; do
	        wait $pid;
	done
#        art_illumina -ss HS20 -i germ_simulated.fa -p -l 100 -f 15 -m 350 -s 20 -o germ_simulated
#	bwa-0.7.15 mem  -t 24 /DATA1/AITL/hg19/ucsc.hg19.fasta germ_simulated1.fq germ_simulated2.fq | samtools view -bS >  germ_simulated.bam
#	samtools sort germ_simulated.bam > germ_simulated_sorted.bam
#	samtools index germ_simulated_sorted.bam
#        ./SV_germ.sh &
#       ./bicseq_preprocess_germ.sh &
	cd ../

done
