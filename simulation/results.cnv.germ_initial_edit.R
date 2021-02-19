args=commandArgs(T);
if(args[1]=="germline"){
	t=read.table("germline_results.cnv",stringsAsFactors=F)
}else{
        t=read.table("results.cnv",stringsAsFactors=F)

}
s=read.table("germline_initial_SVs_changed_size1000",stringsAsFactors=F);
s[s[,3]=="X",3]=23;
s[s[,6]=="X",6]=23;

for(i in 1:nrow(s)){
	cor=which(t[,1]==s[i,3] & t[,2]<=s[i,4] & t[,3] >=s[i,7]);
	if(length(cor)==1){

		t[seq(cor+2,nrow(t)+2),]=t[seq(cor,nrow(t)),];
		t[cor+1,] = t[cor,];
		if(s[i,1]=="3to5"){
			t[cor,3] = s[i,4]-1;
			t[cor+1,2]=s[i,4];
			t[cor+1,3]=s[i,7];
			t[cor+2,2]=s[i,7]+1;
			if(s[i,2]==1){
				t[cor+1,5]=0;
			}else{
				t[cor+1,6]=0;
			}
		}
                if(s[i,1]=="5to3"){
                        t[cor,3] = s[i,4]-1;
                        t[cor+1,2]=s[i,4];
                        t[cor+1,3]=s[i,7];
                        t[cor+2,2]=s[i,7]+1;
                        if(s[i,2]==1){
                                t[cor+1,5]=t[cor,5]*(s[i,8]+1);
                        }else{
                                t[cor+1,6]=t[cor,6]*(s[i,8]+1);
                        }
                }
                t[cor+1,4]=t[cor+1,5]+t[cor+1,6];


	}
}

if(args[1]=="germline"){
	write.table(t,"germline_results.cnv.edited", quote=F, sep="\t", row.names=F, col.names=F)
}else{
        write.table(t,"results.cnv.edited", quote=F, sep="\t", row.names=F, col.names=F)
}
