break_thres=10000;
true=read.table("../simul5/results.cnv");
#true=read.table("somatic_results.cnv");
#true=true[true$V3-true$V2>1e6,]
true_breaks=data.frame(matrix(ncol=3, nrow=0));
for(i in 1:(nrow(true)-1)){
	if(true[i,1]==true[i+1,1]){
		if(!is.nan(true[i,4]) && !is.nan(true[i+1,4]) && !is.na(true[i,4]) && !is.na(true[i+1,4]) ){
			if(true[i,4]<true[i+1,4]){
				true_breaks[nrow(true_breaks)+1,]=c(true[i,1], true[i,3], "up");
			}else if(true[i,4]>true[i+1,4]){
	                        true_breaks[nrow(true_breaks)+1,]=c(true[i,1], true[i,3], "down");
			}
		}
	}
}


#test=read.table("job/iter1/SNP6_level2.txt.copynumber.filtered.segmean.CN_opt", header=T, stringsAsFactors=F);
test=read.table("../simul5/Weaver.cnv",stringsAsFactors=F);
#test=test[test[,3]-test[,2]>1e6,];
names(test)[6]="modal_cn";
#test=read.table("job/CNV.output_8", header=T, stringsAsFactors=F);
#names(test)[7]="modal_cn";
######################################################
#test=test[abs(test$End.bp-test$Start.bp)>1e4,]
#######################################################


test_breaks=data.frame(matrix(ncol=3, nrow=0));
for(i in 1:(nrow(test)-1)){
        if(test[i,1]==test[i+1,1]){
                if(test[i,"modal_cn"]<test[i+1,"modal_cn"]){
                        test_breaks[nrow(test_breaks)+1,]=c(test[i,1], test[i,3], "up");
                }else if(test[i,"modal_cn"]>test[i+1,"modal_cn"]){
                        test_breaks[nrow(test_breaks)+1,]=c(test[i,1], test[i,3], "down");
                }

        }
}


true_breaks$count=0;
test_breaks$count=0;
true=true_breaks;
test=test_breaks;
true$X2=as.integer(true$X2);
test$X2=as.integer(test$X2);

for(i in 1:nrow(test)){
	t=which(true$X1==test[i,1] & true$X3==test[i,3] & true$count==0 & abs(true$X2-test[i,2]) < break_thres);
	if(length(t)!=0){
		test[i,"count"]=length(t);
		true[t[1],"count"]=1;
	}
}





precision=length(which(test$count!=0))/nrow(test)
recall=length(which(true$count!=0))/nrow(true)

print(precision)
print(recall)
Fscore=2*(precision*recall/(precision+recall));
print(Fscore)

