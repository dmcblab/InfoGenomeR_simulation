#!/bin/bash

IND=1762;


for i in {1..23}
do

	perl individual_select.pl 2 $i $IND &
        pids[${i}]=$!;

done



for pid in ${pids[*]}; do
        wait $pid;
done

cat $IND.g*.fa.*.SVs > germline_initial_SVs
