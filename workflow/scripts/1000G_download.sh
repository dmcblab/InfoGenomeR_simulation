#!/bin/bash
declare -a pids;
let pid_idx=0;

for i in `seq 1 22`;do
	wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz && 
	wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi &
	pids[$pid_idx]=$!;
	 ((++pid_idx))
done

#for X
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz &&
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi &
pids[$pid_idx]=$!;


for pid_idx in `0 22`;do
	wait ${pids[$pid_idx]};
done


ln -s 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz  1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz 
ln -s 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi  1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
