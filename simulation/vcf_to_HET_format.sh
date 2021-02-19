#!/bin/bash
rm het_snps.format
for i in {1..22}
do
	cat ./mpileup/$i.vcf | awk -f vcf_to_HET_format.awk >> het_snps.format
done
        cat ./mpileup/X.vcf | awk -f vcf_to_HET_format.awk >> het_snps.format

