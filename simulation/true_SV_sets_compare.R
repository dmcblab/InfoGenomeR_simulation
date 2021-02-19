true=read.table("true_SV_sets_dup_checked.refcoor.germline_initial_SVs_changed_size1000",stringsAsFactors=F)
true$count=0;
#true=true[true$V3==true$V6 & (true$V1=="3to3" |true$V1=="5to5"),]
#test=read.table("../simul5/Weaver", stringsAsFactors=F)
test=read.table("/DATA2/leeyh/BG_job/2017.11.02.simulation_job/f15/germ/iter4/filtered.format.truncated.break_adjusted.CN_opt.filtered");
#test=read.table("/DAS_Storage5/leeyh/novoBreak/simul2/novoBreak.pass.flt.vcf.format",stringsAsFactors=F);
#test=read.table("CREST");
#test=read.table("/DASstorage6/leeyh/BG_job/2017.10.10.simulations/simul1/job/iter1/iter2/t", stringsAsFactors=F);

#test=read.table("/DASstorage6/leeyh/BG_job/2017.08.29.simulations_NGS/simul1/job_bicseq_raw/iter1/iter2/iter3/iter4/iter5/iter6/filtered.format.truncated.break_adjusted",stringsAsFactors=F)
#test=read.table("simulated_sorted.bam.Weaver.FINAL_SV.format",stringsAsFactors=F)
#test=test[test$V1=="<INV>",];
test$count=0;
break_thres=2e4;
#test=test[test$V7>7,]
for(true_i in 1:nrow(true)){
        reverse=strsplit(true[true_i,1],split="to");
        reverse=paste(reverse[[1]][2],"to", reverse[[1]][1],sep="");

        m=which((test$V2==true[true_i,3] & abs(test$V3-true[true_i,4]) < break_thres & test$V4==true[true_i,6] & abs(test$V5-true[true_i,7]) < break_thres & test$V6== true[true_i,1] )|
                (test$V4==true[true_i,3] & abs(test$V5-true[true_i,4]) < break_thres & test$V2==true[true_i,6] & abs(test$V3-true[true_i,7]) < break_thres & test$V6 == reverse)        )

  #      m=which((test$V2==true[true_i,3] & abs(test$V3-true[true_i,4]) < break_thres & test$V4==true[true_i,6] & abs(test$V5-true[true_i,7]) < break_thres)|
  #              (test$V4==true[true_i,3] & abs(test$V5-true[true_i,4]) < break_thres & test$V2==true[true_i,6] & abs(test$V3-true[true_i,7]) < break_thres)        )

        if(length(m)!=0){
                true[true_i,"count"]=true[true_i,"count"]+1;
                test[m,"count"]=test[m,"count"]+1;
        }

}
true=true[-which(true$V3==true$V6 & abs(true$V4-true$V7)<1e5),]
test=test[-which(test$V2==test$V4 & abs(test$V3-test$V5)<1e5),]



precision=length(which(test$count!=0))/nrow(test)
recall=length(which(true$count!=0))/nrow(true)

print(precision)
print(recall)
Fscore=2*(precision*recall/(precision+recall));
print(Fscore)
