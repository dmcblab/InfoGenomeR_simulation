#!/bin/bash

IND=2342;
SVs=0;
new_TEI=0;
new_TEI_per_chr_hap=0;

numt_p=0.01;
alu_p=0.76;
l1_p=0.18;
sva_p=0.05;

numt_cul=0.01;
alu_cul=0.77;
l1_cul=0.95;
sva_cul=1;

TEI_p=0;
SVs_p=1;
for i in `seq 1 23`
do
	if [ $i != 23 ]; then
		j=$i;
	else
		j="X";
	fi
	chr_end=300000000
	zcat 1000G_phased_high_coverage/1kGP_high_coverage_Illumina.chr$j.filtered.SNV_INDEL_SV_phased_panel.vcf.gz |  grep -v '#' |  awk -F "\t" 'BEGIN{srand();}{
		if($'$IND' != "."){
                split($'$IND',g, "|");

		if(length($4)==1 && length($5)==1){
			if(g[1]==0){
				g1=$4;
			}else{
				split($5,a,",");
				g1=a[g[1]];
			}
	                if(g[2]==0){
	                        g2=$4;
	                }else{
	                        split($5,a,",");
	                        g2=a[g[2]];
	                }
	
	                if(g[1]!=0) print $1"\t"$2"\t"$4"\t"$5"\t"$"'$IND'""\t"g1"\t"g[1]"\tm" > "'$IND'.g1.vcf.'$j'"
			if(g[2]!=0) print $1"\t"$2"\t"$4"\t"$5"\t"$"'$IND'""\t"g2"\t"g[2]"\tm" > "'$IND'.g2.vcf.'$j'"
		}else{
			split($8, s, "SVTYPE=");
			split(s[2], SVTYPE, ";");
			if(rand() < '$TEI_p'){
				if(SVTYPE[1]=="ALU" || SVTYPE[1]=="SVA" || SVTYPE[1]=="LINE1"){
					if(SVTYPE[1]=="ALU"){svname="alu";}else if(SVTYPE[1]=="SVA"){svname="sva";}else{svname="L1";};
					split($8, m, "MEINFO=");
					split(m[2],MEINFO, ";");
					if(g[1]==0){
						g1=$4;
					}else{
						g1=MEINFO[1];
					}
					if(g[2]==0){
						g2=$4;
					}else{
						g2=MEINFO[1];
					}
					if(g[1]!=0) print $1"\t"$2"\t"$4"\t"$5"\t"$"'$IND'""\t"g1"\t"g[1]"\t"svname > "'$IND'.g1.vcf.'$j'"
					if(g[2]!=0) print $1"\t"$2"\t"$4"\t"$5"\t"$"'$IND'""\t"g2"\t"g[2]"\t"svname > "'$IND'.g2.vcf.'$j'"

				}else if(SVTYPE[1]=="INS"){
					split($8, ms, "MSTART=");
					split(ms[2],MSTART, ";");
					split($8, me, "MEND=");
					split(me[2],MEND, ";");
					split($8,e,";END=");
					split(e[2],end,";");
					if(g[1]==0){
						g1=$4;
					}else{
						g1="numt,"MSTART[1]","MEND[1]",+";
					}
					if(g[2]==0){
						g2=$4;
					}else{
						g2="numt,"MSTART[1]","MEND[1]",+";
					}
					if(g[1]!=0) print $1"\t"$2"\t"$4"\t"end[1]"\t"$"'$IND'""\t"g1"\t"g[1]"\tnumt" > "'$IND'.g1.vcf.'$j'"
					if(g[2]!=0) print $1"\t"$2"\t"$4"\t"end[1]"\t"$"'$IND'""\t"g2"\t"g[2]"\tnumt" > "'$IND'.g2.vcf.'$j'"
				}
			}
			if(rand() < '$SVs_p'){
				if(SVTYPE[1]=="DUP" || SVTYPE[1]=="DEL" || SVTYPE[1]=="INV"){
					if(SVTYPE[1]=="DUP"){svname="dup";}else if(SVTYPE[1]=="DEL"){svname="del";}else{svname="inv";};
					split($8,e,";END=");
					split(e[2],end,";");
				       if(g[1]==0){
						g1=$4;
					}else{
						split($5,a,",");
						g1=a[g[1]];
					}
					if(g[2]==0){
						g2=$4;
					}else{
						split($5,a,",");
						g2=a[g[2]];
					}
					if(g[1]!=0) print $1"\t"$2"\t"$4"\t"end[1]"\t"$"'$IND'""\t"g1"\t"g[1]"\t"svname > "'$IND'.g1.vcf.'$j'"
					if(g[2]!=0) print $1"\t"$2"\t"$4"\t"end[1]"\t"$"'$IND'""\t"g2"\t"g[2]"\t"svname > "'$IND'.g2.vcf.'$j'"
				}
			}

		}
		}


	}END{
		for(SVs_i=1;SVs_i<='$new_TEI_per_chr_hap';SVs_i++){
			p=rand();
			coor=1;
			while(coor<=1 || coor>='$chr_end'){
				coor=int(rand()*'$chr_end')+1;
			}
				coor_end=coor+int(rand()*15000);
			if(p<'$numt_cul'){
				print '$j'"\t"coor"\t""N""\t"coor_end"\t""0|0""\t""numt,1,"coor_end-coor",+""\t""1""\t""numt" > "'$IND'.g1.vcf.'$j'"
			}else if(p<'$alu_cul'){
				print '$j'"\t"coor"\t""N""\t""<INS:ME:ALU>""\t""0|0""\t""AluUndef,1,281,+""\t""1""\t""alu" > "'$IND'.g1.vcf.'$j'"
			}else if(p<'$l1_cul'){
				print '$j'"\t"coor"\t""N""\t""<INS:ME:LINE1>""\t""0|0""\t""LINE1,26,6019,+""\t""1""\t""L1" > "'$IND'.g1.vcf.'$j'"
			}else if(p<'$sva_cul'){
				print '$j'"\t"coor"\t""N""\t""<INS:ME:SVA>""\t""0|0""\t""SVA,304,1627,-""\t""1""\t""sva" > "'$IND'.g1.vcf.'$j'"
			}
		}

                for(SVs_i=1;SVs_i<='$new_TEI_per_chr_hap';SVs_i++){
                        p=rand();
                        coor=1;
                        while(coor<=1 || coor>='$chr_end'){
                                coor=int(rand()*'$chr_end')+1;
                        }
                                coor_end=coor+int(rand()*15000);
                        if(p<'$numt_cul'){
                                print '$j'"\t"coor"\t""N""\t"coor_end"\t""0|0""\t""numt,1,"coor_end-coor",+""\t""1""\t""numt" > "'$IND'.g2.vcf.'$j'"
                        }else if(p<'$alu_cul'){
                                print '$j'"\t"coor"\t""N""\t""<INS:ME:ALU>""\t""0|0""\t""AluUndef,1,281,+""\t""1""\t""alu" > "'$IND'.g2.vcf.'$j'"
                        }else if(p<'$l1_cul'){
                                print '$j'"\t"coor"\t""N""\t""<INS:ME:LINE1>""\t""0|0""\t""LINE1,26,6019,+""\t""1""\t""L1" > "'$IND'.g2.vcf.'$j'"
                        }else if(p<'$sva_cul'){
                                print '$j'"\t"coor"\t""N""\t""<INS:ME:SVA>""\t""0|0""\t""SVA,304,1627,-""\t""1""\t""sva" > "'$IND'.g2.vcf.'$j'"
                        }
                }


	}' &
        pids[${i}]=$!;

done


for pid in ${pids[*]}; do
        wait $pid;
done


for i in {1..23}
do
        if [ $i != 23 ]; then
                j=$i;
        else
                j="X";
        fi
	sort -nk 2,2 $IND.g1.vcf.$j > $IND.g1.vcf.$j.swp && mv $IND.g1.vcf.$j.swp $IND.g1.vcf.$j
        sort -nk 2,2 $IND.g2.vcf.$j > $IND.g2.vcf.$j.swp && mv $IND.g2.vcf.$j.swp $IND.g2.vcf.$j
done
