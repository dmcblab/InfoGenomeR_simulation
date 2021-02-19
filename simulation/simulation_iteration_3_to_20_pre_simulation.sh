#!/bin/bash
sample_name="HG00732"
IND=`cat /DASstorage6/leeyh/1000G_SNP_indel_SV_integrated/samples.index | awk '{if($1=="'$sample_name'") print NR}'`
##cat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf  |  awk '{if($1762 != "0|0" && $1762 !=".") print $0}' | grep 'SVTYPE'  | grep -v 'CNV' | wc -l
#SVs=2953;  for 1762
SVs=3074;
new_TEI=380;
new_TEI_per_chr_hap=8;

temp=$PWD
cd /DASstorage6/leeyh/1000G_SNP_indel_SV_integrated/
./individual_select.sh $IND $SVs $new_TEI $new_TEI_per_chr_hap

./individual_select_loop_hap_1.sh $IND
./individual_select_loop_hap_2.sh $IND

cd $temp

for i in {1..32}
do
	echo $i
	mkdir simul$i
	cp *.R *.sh *.pl *.awk simul$i/


	cd simul$i
	cp /DASstorage6/leeyh/1000G_SNP_indel_SV_integrated/$IND* ./


	perl SV_chr_rand_sub_function_N_excluded.pl $IND > results && cat results | awk -f SV.results_to_SVs.awk > true_SV_sets && cat germ_results | awk -f SV.results_to_SVs.awk > true_SV_sets_germline &
        pids[${i}]=$!;
	cd ../



done

for pid in ${pids[*]}; do
        wait $pid;
done

