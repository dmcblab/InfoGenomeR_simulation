#!/bin/bash
reference=/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa
delly_v0.7.6_CentOS5.4_x86_64bit call -t DEL -g $reference -o germ_DEL.bcf germ_simulated_sorted.bam && delly_v0.7.6_CentOS5.4_x86_64bit call -t DUP -g $reference  -o germ_DUP.bcf germ_simulated_sorted.bam && delly_v0.7.6_CentOS5.4_x86_64bit call -t INV -g $reference -o germ_INV.bcf germ_simulated_sorted.bam && delly_v0.7.6_CentOS5.4_x86_64bit call -t TRA -g $reference -o germ_TRA.bcf germ_simulated_sorted.bam &
pids[1]=$!;


for pid in ${pids[*]}; do
        wait $pid;
done


bcftools view germ_DEL.bcf > germ_DEL.vcf
bcftools view germ_DUP.bcf > germ_DUP.vcf
bcftools view germ_INV.bcf > germ_INV.vcf
bcftools view germ_TRA.bcf > germ_TRA.vcf

