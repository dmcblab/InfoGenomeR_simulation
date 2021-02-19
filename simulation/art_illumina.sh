#!/bin/bash

for i in `seq 1 5`; do
	art_illumina -ss HS20 -i $1.fa -p -l 100 -f 5 -m 350 -s 20 -o $1$i -rs $RANDOM &
        pids[${i}]=$!;
done
	art_illumina -ss HS20 -i germ_simulated.fa -p -l 100 -f 10 -m 350 -s 20 -o germ_simulated -rs $RANDOM &
        pids[6]=$!;

for pid in ${pids[*]}; do
        wait $pid;
done
