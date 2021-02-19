BEGIN{
	index1="initial";
}
{
        if(split($0,f,"\t")>2){
	#if(length($0)>2){
		if(index1 == "initial"){
		        index1=$2;
			hap1=$4;
			chr1=$5;
		}else{
			if($5!=chr1){
				if(index1 >0){
					ori1=3;
				}else{
					ori1=5;
				}
				if($1 > 0){
					ori2=5;
				}else{
					ori2=3;
				}
				type=ori1"to"ori2;
				print type"\t"hap1"\t"chr1"\t"index1"\t"$4"\t"$5"\t"$1;

			}
			else if(index1+1 != $1){
                                if(index1 > 0 &&  $1 > 0)
					type="3to5"
				if(index1 < 0 && $1 < 0)
					type="5to3"
				#if(index1 * $1 > 0 && index1 >= $1)
				#	type="5to3"
				#if(index1 * $1 > 0 && index1 < $1)
				#	type="3to5"
				if(index1 * $1 < 0 && index1 < 0)
					type="5to5"
				if(index1 * $1 < 0 && index1 > 0)
					type="3to3"
                                print type"\t"hap1"\t"chr1"\t"index1"\t"$4"\t"$5"\t"$1;

			}
			hap1=$4;
			chr1=$5;
			index1=$2;
		}
	}else{

		index1="initial";
	}

}

