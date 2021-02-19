#!/bin/bash
for i in `seq 1 23`; do
	if [ $i -eq 23 ];then
		i="X";
	fi
	cat vcf_somatic/$i.vcf | tail -n +130 | awk '{

		if(length($4)==1 && length($5) == 1){
			split($10,f1,":");
			split($11,f2,":");
			split(f1[7], normal, ",");
			split(f2[7], tumor, ",");
			print "'$i'""\t"$2"\t"$3"\t"$4"\t"$5"\t"normal[1]"\t"normal[2]"\t"tumor[1]"\t"tumor[2];
		}
	}'	> vcf_somatic/$i.vcf.format




	cat vcf_somatic/$i.vcf.format | awk '{

		if($6+$7 > 10 && $8+$9 > 10 && $6/($6+$7)> 0.2 && $6/($6+$7) < 0.8){
			print $0;
		}

	}' > vcf_somatic/$i.vcf.format.normal_het


        cat vcf_somatic/$i.vcf.format | awk '{

                if($8+$9 > 10 && $8/($8+$9)> 0.2 && $8/($8+$9) < 0.8){
                        print $0;
                }

        }' > vcf_somatic/$i.vcf.format.tumor_het

	cat vcf_somatic/*.vcf.format.tumor_het > vcf_somatic/het_snps.format
	cat  vcf_somatic/*.vcf.format.normal_het > vcf_somatic/het_snps.format.somatic


done
