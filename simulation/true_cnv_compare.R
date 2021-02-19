true=read.table("results.cnv",stringsAsFactors=F)
test=read.table("job/iter1/iter2/iter3/iter1/iter2/iter3/iter4/iter5/iter6/iter7/SNP6_level2.txt.copynumber.filtered.segmean.CN_opt.ACN",stringsAsFactors=F, header=T);
test=test[,c("Chromosome","Start.bp","End.bp","modal_cn","allele_1","allele_2")];
test[test$Chromosome==23,1]="X";
names(test)=c("V1","V2","V3","V4","allele_1","allele_2");
test=test[complete.cases(test),];
ploidy=3;
thres=0.9;
test$count=0;
true$count=0;
test=test[test$V3-test$V2>5,];

#swp=test;
#test=true;
#true=swp;



overlap<-function(a,b){
        if(a[1]<=b[1] && a[2]>=b[2]){
                return(c(b[1],b[2]));
        }else if(b[1]<=a[1] && a[2]<=b[2]){
                return(c(a[1],a[2]));
        }else if(a[1]<=b[1] && a[2]>=b[1]){
                return(c(b[1],a[2]));
        }else if(a[1]<=b[2] && a[2]>=b[2]){
                return(c(a[1],b[2]));
        }else{
                return(0);
        }

}

for(true_i in 1:nrow(true)){
	overlap_sum=0;
	for(test_i in 1:nrow(test)){
		if(test[test_i,1]==true[true_i,1]){
			if((test[test_i,4]== ploidy && true[true_i,4] == ploidy) || (test[test_i,4]-ploidy)*(true[true_i,4]-ploidy)>0 && test[test_i,4]-true[true_i,4]==0){
				r=overlap(c(true[true_i,2],true[true_i,3]),c(test[test_i,2],test[test_i,3]));
				if(length(r)==2){
					overlap_sum=overlap_sum+r[2]-r[1]+1;
				}
			}
		}
	}
	if(overlap_sum/(true[true_i,3]-true[true_i,2]+1)>0.9)
		true[true_i,"count"]=1;


}


#precision=length(which(test$count!=0))/nrow(test)
recall=length(which(true$count!=0))/nrow(true)

#print(precision)
print(recall)
#Fscore=2*(precision*recall/(precision+recall));
#print(Fscore)
