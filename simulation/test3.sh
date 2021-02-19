#!/bin/bash
reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa

export PATH=/DASstorage6/leeyh/novoBreak_distribution_v1.1.3rc:$PATH
declare -a arr;

k=0;

	for j in 3 5 10 15 20; do
		for i in 0.6 0.75 0.9; do
			arr[$k]="echo $j && echo $i && sleep 2";
			k=$(($k+1));
                done
	done
	
	k=$(($k-1));

	for z in 0 1 2 3 4 5; do
		eval ${arr[$k]} &
                pids[$z]=$!;
                k=$(($k -1));
	done
 
        while [ $k -ge 0 ];do

                for l in 0 1 2 3 4 5;do
                        n=`ps | grep ${pids[$l]} | wc -l`;
                        if [ $n -eq 0 ] && [ $k -ge 0 ];then
                                eval ${arr[$k]} &
                                pids[$l]=$!;
                                k=$(($k -1));
                        fi
                done
                sleep 1;

        done

	for z in 0 1 2 3 4 5;do
        	wait ${pids[$z]};
	done
