args=commandArgs(TRUE);
germline_min=2e3;
inv_max=48367;
inv_a=0.5229973;
inv_b=2.667276;
del_max=2258237;
del_a=0.4599259;
del_b=62.90983;
dup_max=801035;
dup_a=0.8964938;
dup_b=11.24391;

if(args[1]=="dup"){
	size=min(germline_min+round(rbeta(1,dup_a, dup_b)*dup_max), dup_max)
}else if(args[1]=="del"){
        size=min(germline_min+round(rbeta(1,del_a, del_b)*del_max), del_max)
}else if(args[1]=="simple_inv"){
        size=min(germline_min+round(rbeta(1,inv_a, inv_b)*inv_max), inv_max)
}
print(size);
