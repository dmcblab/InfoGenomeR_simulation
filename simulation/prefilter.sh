#!/bin/bash

#echo "Type SRR ID"
#read tumorIDini

thres=$1;
rm filtered.format

for i in {1..4}
do
	if [ $i == 1 ]
	then
		tumorID=DEL.vcf
	elif [ $i == 2 ]
	then
                tumorID=INV.vcf
        elif [ $i == 3 ]
	then
                tumorID=DUP.vcf
        else
                tumorID=TRA.vcf
	fi      

		cat $tumorID | awk -F "\t" '{
					        n=split($8,f,";");
						if(n>10){
							split(f[13],SR, "=");
							SR_number=SR[2];
						}else{
							SR_number=0;
						}
						split(f[2],SVTYPE,"=");
					        split(f[4],chr2,"=");
					        split(f[5],pos2,"=");
					        split(f[6],PE,"=");
						split(f[7], MAPQ, "=");
						split(f[8], ori, "=");
					        if($7=="PASS" && (PE[2] > '$thres' && SR_number > -1)&& MAPQ[2] > 20){
					                print "<"SVTYPE[2]">""\t"$1"\t"$2"\t"chr2[2]"\t"pos2[2]"\t"ori[2]"\t"PE[2]"\t"SR_number"\t"MAPQ[2]"\t"PE[2]+SR_number"\t0\t0\t1\t1\t"
					        }
					}' >> filtered.format

done





