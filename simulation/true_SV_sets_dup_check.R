t=read.table("true_SV_sets",stringsAsFactors=F)
t$V4=abs(t$V4);
t$V7=abs(t$V7);

t$V8=0;
new_t=t[0,]
i=1;
while(nrow(t)!=0){


        reverse=strsplit(t[i,1],split="to");
        reverse=paste(reverse[[1]][2],"to", reverse[[1]][1],sep="");
        m=which((t$V1==t[i,1] & t$V2==t[i,2] &t$V3==t[i,3] &t$V4==t[i,4] &t$V5==t[i,5] &t$V6==t[i,6] &t$V7==t[i,7])     |
                (t$V1==reverse & t$V5==t[i,2] &t$V6==t[i,3] &t$V7==t[i,4] &t$V2==t[i,5] &t$V3==t[i,6] &t$V4==t[i,7])            )

	t[1,8]=length(m);
	new_t=rbind(new_t,t[1,]);
        t=t[-m,];

}

new_t$V4=abs(new_t$V4)
new_t$V7=abs(new_t$V7)


write.table(new_t, "true_SV_sets_dup_checked", col.names=F, row.names=F, quote=F, sep="\t")


