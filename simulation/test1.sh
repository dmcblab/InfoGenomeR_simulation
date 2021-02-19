#!/bin/bash
declare -a arr;
k=0;
for i in 15 14 13 12 11
do
	arr[$k]="echo $i &&sleep $i";
	k=$(($k+1));
done
	k=$(($k-1));
		eval ${arr[$k]} &
		pids[0]=$!;
		k=$(($k -1));

                eval ${arr[$k]} &
                pids[1]=$!;
                k=$(($k -1));

        while [ $k -ge 0 ];do

		for l in 0 1;do
			n=`ps | grep ${pids[$l]} | wc -l`;
			if [ $n -eq 0 ] && [ $k -ge 0 ];then
				eval ${arr[$k]} &
				pids[$l]=$!;
				k=$(($k -1));
			fi
		done
		sleep 1;

	done
	
	wait ${pids[0]};
	wait ${pids[1]};
