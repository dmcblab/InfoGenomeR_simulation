{if(NR>24){

	split($1,chr,"chr");
	split($10,f,":");
	split($8,t,";");
	split(t[3],het,"HET=");
	if(het[2]==1)
		print chr[2]"\t"$2"\t"$3"\t"$4"\t"$5"\t""NA""\t""NA""\t"f[5]"\t"f[6];



}}

