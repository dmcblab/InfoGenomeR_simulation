args=commandArgs(TRUE);
somatic_min=2e4;
somatic_max=50000000;
insertion_max=10000000;
#inv_max=21332326;
inv_a=0.428702;
inv_b=1.262533;
#del_max=212123331;
del_a=0.265586;
del_b=1.905053;
#dup_max=223968675;
dup_a=0.2629577;
dup_b=2.686045;

ins_a=0.4284392;
ins_b=1.979344;
ins_max=1e8;

if(args[1]=="dup"){
	size=min(somatic_min+round(rbeta(1,dup_a, dup_b)*somatic_max), somatic_max)
}else if(args[1]=="del"){
        size=min(somatic_min+round(rbeta(1,del_a, del_b)*somatic_max), somatic_max)
}else if(args[1]=="simple_inv"){
        size=min(somatic_min+round(rbeta(1,inv_a, inv_b)*somatic_max), somatic_max)

}else if(args[1]=="insertion"){
	size=min(somatic_min+round(rbeta(1,ins_a, ins_b)*10000000), 10000000)
}

print(size);
