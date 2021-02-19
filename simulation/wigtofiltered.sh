#!/bin/bash
cat wgEncodeCrgMapabilityAlign100mer.wig  | awk '{split($1,chr,"chr");if($3-$2>20 && $4<0.2){
 if(chr[2]==1 || chr[2]==2 || chr[2]==3 || chr[2]==4 || chr[2]==5 || chr[2]==6 ||chr[2]==7 || chr[2]==8 || chr[2]==9 || chr[2]==10 || chr[2]==11 || chr[2]==12||chr[2]==13 || chr[2]==14 || chr[2]==15 || chr[2]==16 || chr[2]==17 || chr[2]==18 ||chr[2]==19 || chr[2]==20 || chr[2]==21 || chr[2]==22 || chr[2]=="X"){
 print chr[2]"\t"$2"\t"$3"\t"$4}}}' > wgEncodeCrgMapabilityAlign100mer.wig.filtered;
