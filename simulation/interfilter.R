args=commandArgs(TRUE)
#args=c();
#args[1]="synthetic.challenge.set1.tumor.v2.bam";
#args[2]="synthetic.challenge.set1.normal.v2.bam";

for(name in c("DEL", "DUP", "TRA", "INV")){
        tumor=read.table(paste(name,".vcf",sep=""),stringsAsFactors=F);
	normal=read.table(paste("../germ_",name,".vcf",sep=""),stringsAsFactors=F);
	tumor$count=0;
	for(i in 1:nrow(tumor)){
		tumor$pos2[i]=strsplit(tumor$V8[i],split=";END=")[[1]][2];
		tumor$pos2[i]=strsplit(tumor$pos2[i],split=";")[[1]][1];

		tumor$chr2[i]=strsplit(tumor$V8[i],split=";CHR2=")[[1]][2];
		tumor$chr2[i]=strsplit(tumor$chr2[i],split=";")[[1]][1];

	}

	for(i in 1:nrow(normal)){
		normal$pos2[i]=strsplit(normal$V8[i],split=";END=")[[1]][2];
		normal$pos2[i]=strsplit(normal$pos2[i],split=";")[[1]][1];

		normal$chr2[i]=strsplit(normal$V8[i],split=";CHR2=")[[1]][2];
		normal$chr2[i]=strsplit(normal$chr2[i],split=";")[[1]][1];

	}
	tumor$pos2=as.integer(tumor$pos2);
	normal$pos2=as.integer(normal$pos2);


	for(i in 1:nrow(tumor)){
		if(length(which(normal$V1==tumor$V1[i] & abs(normal$V2-tumor$V2[i])<500 & normal$chr2==tumor$chr2[i] & abs(normal$pos2-tumor$pos2[i])<500))!=0)
			tumor$count[i]=1;
	}
	tumor=tumor[tumor$count==0,];
	write.table(tumor[,1:10], paste("somatic.",name,".vcf",sep=""), sep="\t", quote=F, col.names=F, row.names=F);
}
