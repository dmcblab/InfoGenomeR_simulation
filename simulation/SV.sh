#!/bin/bash
reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa
delly_v0.7.6_CentOS5.4_x86_64bit call -t DEL -g $reference -o DEL.bcf simulated_sorted.bam && delly_v0.7.6_CentOS5.4_x86_64bit call -t DUP -g $reference  -o DUP.bcf simulated_sorted.bam && delly_v0.7.6_CentOS5.4_x86_64bit call -t INV -g $reference -o INV.bcf simulated_sorted.bam && delly_v0.7.6_CentOS5.4_x86_64bit call -t TRA -g $reference -o TRA.bcf simulated_sorted.bam &
pids[1]=$!;

for pid in ${pids[*]}; do
        wait $pid;
done


bcftools view DEL.bcf > DEL.vcf

bcftools view DUP.bcf > DUP.vcf
 bcftools view INV.bcf > INV.vcf
 bcftools view TRA.bcf > TRA.vcf


