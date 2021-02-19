#!/bin/bash


		cat novoBreak.pass.flt.vcf | awk -F "\t" '{
					        n=split($8,f,";");
					        split(f[8],chr2,"=");
					        split(f[9],pos2,"=");
						split(f[10],len,"=");
						split(f[2], ori, "=");


						if($5=="<TRA>" || $5=="<INV>"){
						        print $5"\t"$1"\t"$2"\t"chr2[2]"\t"pos2[2]"\t"ori[2]"\t"($36+$38)/2"\t"$19"\t"$20"\t"($36+$38)/2+$19"\t0\t0\t1\t1";
						}else{
							if($5=="<DUP>"){
								if(len[2]>0){
		                                                        print $5"\t"$1"\t"$2"\t"chr2[2]"\t"pos2[2]"\t"ori[2]"\t"($36+$38)/2"\t"$19"\t"$20"\t"($36+$38)/2+$19"\t0\t0\t1\t1"
								}else{
									split(ori[2], cori, "to");
	                                                                print $5"\t"$1"\t"$2"\t"chr2[2]"\t"pos2[2]"\t"cori[2]"to"cori[1]"\t"($36+$38)/2"\t"$19"\t"$20"\t"($36+$38)/2+$19"\t0\t0\t1\t1";
								}
							}
							if($5=="<DEL>"){
                                                               if(len[2]>0){
                                                                        split(ori[2], cori, "to");
                                                                        print $5"\t"$1"\t"$2"\t"chr2[2]"\t"pos2[2]"\t"cori[2]"to"cori[1]"\t"($36+$38)/2"\t"$19"\t"$20"\t"($36+$38)/2+$19"\t0\t0\t1\t1";
                                                                }else{
                                                                        print $5"\t"$1"\t"$2"\t"chr2[2]"\t"pos2[2]"\t"ori[2]"\t"($36+$38)/2"\t"$19"\t"$20"\t"($36+$38)/2+$19"\t0\t0\t1\t1";
                                                                }

							}

						
						}
					}' > novoBreak.pass.flt.vcf.format



		#cat novoBreak.pass.flt.vcf.format | awk '{if($7+$8 > 3 && $9 !=0 ) print $0}' > novoBreak.pass.flt.vcf.format.swp 
		#mv novoBreak.pass.flt.vcf.format.swp novoBreak.pass.flt.vcf.format
