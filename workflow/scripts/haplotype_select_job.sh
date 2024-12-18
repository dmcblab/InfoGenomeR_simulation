#!/bin/bash
SCRIPT_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")

humandb=$(readlink -f $SCRIPT_DIR/../../humandb)

haplotype_selected_output_dir=$(readlink -f $1);
haplotype_name=$2;
output_dir=$(readlink -f $3);

ref_version=$4;
echo "ref version: $ref_version";

mkdir -p $output_dir
cd "$output_dir" || { echo "Failed to change directory to $output_dir"; exit 1; }


IND=$(bcftools view -h  ${humandb}/${ref_version}/1000G/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
grep '#CHROM' | \
awk -v haplotype="$haplotype_name" '{for(i=1;i<=NF;i++){if($i==haplotype) print i}}')


for i in {1..23}
do
	perl ${SCRIPT_DIR}/haplotype_select.pl 1 $i $IND $humandb $haplotype_selected_output_dir $output_dir $ref_version
done


for i in {1..23}
do
	perl ${SCRIPT_DIR}/haplotype_select.pl 2 $i $IND $humandb $haplotype_selected_output_dir $output_dir $ref_version 
done


cat $IND.g*.fa.*.SVs > germline_initial_SVs
