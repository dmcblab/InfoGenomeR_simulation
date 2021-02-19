#!/bin/bash
#IND=1762
k=1
	./art_illumina_3_to_20.sh simulated
	./art_illumina_3_to_20_purity.sh
	./art_illumina_3_to_20_fold_coverage.sh

        for j in 3 5 10 15 20; do
                cp *.R *.sh *.pl *.awk f$j/
                cd f$j
                for i in 0.6 0.75 0.9; do
                        cp *.R *.sh *.pl *.awk p$i/
                        cd p$i;
                        cd ../
                done
                cd ../
        done

	./weaver.sh &
	weaver_pid=$!;

	
	for j in 3 5 10 15 20; do
                cp *.R *.sh *.pl *.awk f$j/
		cd f$j
		./SV_germ.sh &
	        pids[${k}]=$!;
                k=$(($k+1));
               ./bicseq_preprocess_germ.sh &
                pids[${k}]=$!;
                k=$(($k+1));

		for i in 0.6 0.75 0.9; do
 	        	cp *.R *.sh *.pl *.awk p$i/
			cd p$i;
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
		cd ../
	done



	for pid in ${pids[*]}; do
	        wait $pid;
	done

	./manta_somatic.sh
	./manta_germline.sh
	./novobreak_somatic.sh

	k=0;
	gfServer start localhost 6667 /DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.2bit & 
	
	./crest_extract.sh && ./crest_job.sh &
        pids[${k}]=$!;
	consertc=${pids[$k]};
        k=$(($k+1));
##	./weaver.sh &
##       pids[${k}]=$!;
##        k=$(($k+1));
	./conserting_pre_job1.sh && ./conserting_pre_job3.sh && ./conserting_pre_job4.sh &
	pids[${k}]=$!;
	consertp=${pids[$k]};
        k=$(($k+1));
 
	wait $consertc;
	wait $consertp;
        ./conserting_job.sh && ./conserting_after_job_for_no_ai.sh &
        pids[${k}]=$!;
        k=$(($k+1));

        for pid in ${pids[*]}; do
                wait $pid;
        done

######  ./conserting_pre_job2.sh############################# 필요 없ì
	./het_SNP_detection_somatic.sh  ################################# 다 해결됨
	#./het_SNP_detection.sh

	./main_job1.sh
	./main_job2.sh
	./main_job3.sh
	./main_job4.sh

	wait $weaver_pid
