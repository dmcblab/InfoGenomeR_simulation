true=read.table("true_SV_sets_dup_checked",stringsAsFactors=F)
true$count=0;
test=read.table("weaver/CROSS.SV_format", stringsAsFactors=F)
test=read.table("job_bicseq_raw/iter1/iter2/iter3/iter4/iter5/iter6/filtered.format.truncated.break_adjusted",stringsAsFactors=F)
test$count=0;
break_thres=2e5;
#true=true[-which(true$V3==true$V6 & abs(true$V4-true$V7)<1e6),]
#test=test[-which(test$V2==test$V4 & abs(test$V3-test$V5)<1e6),]
#test=test[test$V7>14,]
for(true_i in 1:nrow(true)){
        reverse=strsplit(true[true_i,1],split="to");
        reverse=paste(reverse[[1]][2],"to", reverse[[1]][1],sep="");

        m=which((test$V2==true[true_i,3] & abs(test$V3-true[true_i,4]) < break_thres & test$V4==true[true_i,6] & abs(test$V5-true[true_i,7]) < break_thres & test$V6== true[true_i,1] )|
                (test$V4==true[true_i,3] & abs(test$V5-true[true_i,4]) < break_thres & test$V2==true[true_i,6] & abs(test$V3-true[true_i,7]) < break_thres & test$V6 == reverse)        )

        if(length(m)!=0){
                true[true_i,"count"]=true[true_i,"count"]+1;
                test[m,"count"]=test[m,"count"]+1;
        }

}


precision=length(which(test$count!=0))/nrow(test)
recall=length(which(true$count!=0))/nrow(true)

print(precision)
print(recall)
Fscore=2*(precision*recall/(precision+recall));
print(Fscore)
