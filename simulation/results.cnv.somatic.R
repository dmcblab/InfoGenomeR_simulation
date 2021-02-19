OUT=data.frame(matrix(ncol=6, nrow=0));

for(i in 1:23){
	t=read.table("results.cnv.edited");
	n=read.table("germline_results.cnv.edited");
	t=t[t$V1==i,];
	n=n[n$V1==i,];


	n=n[n$V4!=2,];
 
        t_remove=c();
        for(t_i in 1:nrow(t)){
                if(length(which((t$V2[t_i]>=n$V2 & t$V3[t_i]<=n$V3)))!=0){
                        t_remove[length(t_remove)+1]=t_i;
                }
        }
        t=t[-t_remove,];

	OUT=rbind(OUT,t);
#	breaks=c(t[,2],t[,3],n[,2],n[,3]);
#	breaks=sort(unique(breaks));
#
#	somatic=data.frame(matrix(ncol=6, nrow=0));
#
#	b_i=1;
#	while(b_i <=length(breaks)-2){
#	        somatic[nrow(somatic)+1,1]=i;
#		somatic[nrow(somatic),2]=breaks[b_i];
#                somatic[nrow(somatic),3]=breaks[b_i+1];
#		if(breaks[b_i+2]==breaks[b_i+1]+1){
#			b_i=b_i+2;
#		}else{
#			b_i=b_i+1;
#		}
#	}
#	
#	for(s_i in 1:nrow(somatic)){
#		for(t_i in 1:nrow(t)){
#			if(t[t_i,2]<=somatic[s_i,2] && t[t_i,3]>=somatic[s_i,3]){
#				tumor_coverage=t[t_i,4];
#			}
#		}
#		for(n_i in 1:nrow(n)){
#                        if(n[n_i,2]<=somatic[s_i,2] && n[n_i,3]>=somatic[s_i,3]){
#				normal_coverage=n[n_i,4];
#                        }
#		}
#
#		somatic[s_i,4]=tumor_coverage/normal_coverage;
#	}
#	OUT=rbind(OUT,somatic);
#}
}
for(i in 2:nrow(OUT)){
	if(OUT[i-1,1] == OUT[i,1] &&  (OUT[i-1,3] != OUT[i,2] || OUT[i-1,3] +1 != OUT[i,2] )){
		OUT[i,2] = OUT[i-1,3]
	}

}

i=1;
while(i <nrow(OUT)){
	if(OUT[i,1] == OUT[i+1,1] && OUT[i,4] == OUT[i+1,4] && OUT[i,5] == OUT[i+1,5] && OUT[i,6] == OUT[i+1,6]){
		OUT[i,3]=OUT[i+1,3];
		OUT=OUT[-(i+1),];
	}else{
		i=i+1;
	}
}
write.table(OUT,"somatic_results.cnv", quote=F, sep="\t", row.names=F, col.names=F);



