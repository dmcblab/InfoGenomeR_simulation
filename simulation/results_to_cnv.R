args=commandArgs(TRUE);

results=list();
for(i in 0:args[1]){
	if(args[2]=="total"){
		results[[i+1]]=read.table(paste("SVs.results.",i,".ref_coor",sep=""),stringsAsFactors=F);	
	}else{
                results[[i+1]]=read.table(paste("germline_SVs.results.",i,".ref_coor",sep=""),stringsAsFactors=F);
	}
	 results[[i+1]][results[[i+1]][,5]=="X",5]=23;
	results[[i+1]][,5]=as.integer(results[[i+1]][,5]);
        results[[i+1]][,4]=as.integer(results[[i+1]][,4]);
}

chrs=list();
for( i in 1:23){
	chrs[[i]]=data.frame(matrix(,nrow=0,ncol=5));
}

for(i in 0:args[1]){
	for(j in 1:nrow(results[[i+1]])){
		if(results[[i+1]][j,5]=="X"){
			chr_index=23;
		}else{
			chr_index=results[[i+1]][j,5];
		}
		if(results[[i+1]][j,1]>0){
			chrs[[chr_index]]=rbind(chrs[[chr_index]], results[[i+1]][j,]);
		}else{
			swp=abs(results[[i+1]][j,1]);
			results[[i+1]][j,1]=abs(results[[i+1]][j,2])
			results[[i+1]][j,2]=swp;
                        chrs[[chr_index]]=rbind(chrs[[chr_index]], results[[i+1]][j,]);
		}

	}
}


chr_cnvs=list();
for( i in 1:23){
        chr_cnvs[[i]]=data.frame(matrix(,nrow=0,ncol=5));
}


for( i in 1:23){

	t=chrs[[i]];
	t_break=sort(unique(c(t[,1],t[,2])));

	for(j in 1:(length(t_break)-1)){
		k=j+1;
		if(t_break[k]-t_break[j]>1){
			chr_cnvs[[i]]=rbind(chr_cnvs[[i]],c(t_break[j],t_break[k],0,0,0));
		}

	}
}


for( i in 1:23){
	for(j in 1:nrow(chr_cnvs[[i]])){
		for(k in 1:nrow(chrs[[i]])){
			if(chr_cnvs[[i]][j,1]>=chrs[[i]][k,1] &&  chr_cnvs[[i]][j,2]<=chrs[[i]][k,2]){
				chr_cnvs[[i]][j,3]=chr_cnvs[[i]][j,3]+1;
				if(chrs[[i]][k,4]==1){
	                                chr_cnvs[[i]][j,4]=chr_cnvs[[i]][j,4]+1;
				}else{
                                        chr_cnvs[[i]][j,5]=chr_cnvs[[i]][j,5]+1;
				}
			}
		}
	}
}


out=data.frame(matrix(ncol=6,nrow=0));
names(out)=c("V1","V2","V3","V4","V5","V6");

for( i in 1:23){
	for( j in 1:nrow(chr_cnvs[[i]])){
		out[nrow(out)+1,]=c(i,chr_cnvs[[i]][j,1], chr_cnvs[[i]][j,2],chr_cnvs[[i]][j,3],chr_cnvs[[i]][j,4],chr_cnvs[[i]][j,5]);
	}
}
#out[out$V1==23,1]="X";
if(args[2]=="total"){
	write.table(out,"results.cnv", col.names=F, row.names=F, sep="\t", quote=F);
}else{
        write.table(out,"germline_results.cnv", col.names=F, row.names=F, sep="\t", quote=F);
}
