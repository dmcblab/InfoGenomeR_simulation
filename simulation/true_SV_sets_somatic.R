true=read.table("true_SV_sets_dup_checked",stringsAsFactors=F)
true$count=0;
test=read.table("true_SV_sets_germline_dup_checked", stringsAsFactors=F);
test$count=0;

break_thres=1;

for(true_i in 1:nrow(true)){
        reverse=strsplit(true[true_i,1],split="to");
        reverse=paste(reverse[[1]][2],"to", reverse[[1]][1],sep="");

        m=which((test$V2==true[true_i,2] & test$V3==true[true_i,3] & abs(test$V4-true[true_i,4]) < break_thres & test$V5==true[true_i,5] & test$V6==true[true_i,6] & abs(test$V7-true[true_i,7]) < break_thres & test$V1== true[true_i,1] )|
                (test$V5==true[true_i,2] & test$V6==true[true_i,3] & abs(test$V7-true[true_i,4]) < break_thres & test$V2==true[true_i,5] & test$V3==true[true_i,6] & abs(test$V4-true[true_i,7]) < break_thres & test$V1 == reverse)        )

        if(length(m)!=0){
                true[true_i,"count"]=true[true_i,"count"]+1;
                test[m,"count"]=test[m,"count"]+1;
        }

}
write.table(true[true$count==0,], "true_SV_sets_somatic_dup_checked",quote=F, sep="\t", row.names=F, col.names=F);
