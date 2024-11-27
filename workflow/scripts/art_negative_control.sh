#!/bin/bash
declare -a pid;

hap_coverage=20


out_f=p0.0
mkdir -p $out_f;


art_illumina -1 PrecisionFDAv2R1.txt -2 PrecisionFDAv2R2.txt  -i simulated.fa -p -l 150 -f $hap_coverage -m 400 -s 150 -o $out_f\/clonal_simulated  &
pids[0]=$!;

wait ${pid[0]}

rm $out_f\/*.aln


unset pid;
declare -a pid;

for read_idx in `seq 1 2`;do
        fq_idx=1;
        cat $out_f\/clonal_simulated$read_idx.fq | awk -v fq_idx_initial=$fq_idx 'BEGIN{fq_idx=fq_idx_initial;}{
                if(NR % 4 == 1){
                        print "@"fq_idx;
                        fq_idx=fq_idx+1;
                }else{
                        print $0;
                }
        }' |  pigz -p 20 > $out_f\/sample_simulated$read_idx.fq.gz &
        pid[$(($read_idx-1))]=$!;
done


wait ${pid[0]}
wait ${pid[1]}


rm $out_f\/*.fq
