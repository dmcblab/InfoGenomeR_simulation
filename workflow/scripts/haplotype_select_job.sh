#!/bin/bash

IND=2342;


for i in {1..23}
do
	perl individual_select.pl 1 $i $IND &
        pids[${i}]=$!;
done



for pid in ${pids[*]}; do
        wait $pid;
done


for i in {1..23}
do
	perl individual_select.pl 2 $i $IND &
	pids[${i}]=$!;
done



for pid in ${pids[*]}; do
	wait $pid;
done

cat $IND.g*.fa.*.SVs > germline_initial_SVs
