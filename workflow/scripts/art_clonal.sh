#!/bin/bash
declare -a pid;

simul_idx=$1;
mosaic_batch=$2;

hap_coverage=20

simul_name=simul$simul_idx

out_f=p1.0/batch$mosaic_batch
mkdir -p $out_f;


art_illumina -na -1 PrecisionFDAv2R1.txt -2 PrecisionFDAv2R2.txt  -i InfoGenomeR_simulation/simulation/${simul_name}/mosaic_batch${mosaic_batch}/sub_clonal/simulated.fa  -p -l 150 -f $hap_coverage -m 400 -s 150 -o $out_f\/sub_clonal_simulated  &
pids[0]=$!;

wait ${pid[0]}


unset pid;
declare -a pid;

for read_idx in `seq 1 2`;do
        fq_idx=1;
        cat $out_f\/sub_clonal_simulated$read_idx.fq | awk -v fq_idx_initial=$fq_idx 'BEGIN{fq_idx=fq_idx_initial;}{
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


#rm $out_f\/*.fq
